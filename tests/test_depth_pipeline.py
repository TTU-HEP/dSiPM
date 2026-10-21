from pathlib import Path
import json
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from dsipm.depth_reconstruction import DepthContract, PhotonObservations, fit_calibration, reconstruct
from dsipm.depth_geometry import FIBER_OFFSETS_MM, rotation, dreamsim_front_contract
from dsipm.depth_io import save_cache, load_cache, corrected_event_from_arrays
from dsipm.depth_observations import make_observations

REPO = Path(__file__).resolve().parents[1]


def fixture(prefix, events=4, n=120):
    """Synthetic unit-test arrays, NOT generated Geant4/shower events."""
    rng = np.random.default_rng(8)
    c = DepthContract()
    keys = np.repeat([f"{prefix}:{i}" for i in range(events)], n)
    ids = np.tile(np.arange(n, dtype=np.int64), events)
    r = rng.uniform(.05, .8, len(keys))
    k = 1.01+.015*r+.005*rng.random(len(keys))
    x = rng.integers(-2, 3, len(keys)); y = rng.integers(-2, 3, len(keys))
    # Two separated populations test that the method returns a candidate, not
    # the unoccupied midpoint between depth modes.
    z = np.where(np.arange(len(keys)) % 2, 450., 1450.)+rng.normal(0, 15, len(keys))
    obs = PhotonObservations(keys, ids, x, y, c.forward_time(z, k), r, 4*x.astype(float), 4*y.astype(float))
    return obs, c, z, k


class DepthPipelineTests(unittest.TestCase):
    def test_cli_calibration_and_all_resolution_methods(self):
        with TemporaryDirectory() as d:
            d = Path(d)
            a, c, z, k = fixture("train")
            b, _, zb, _ = fixture("test")
            save_cache(d/"train", a, c, provenance={"truth_validation": "synthetic fixture exact kappa", "timestamp_sigma_ns": 0},
                       truth_z_mm=z, kappa_truth=k, calibration_eligible=True)
            save_cache(d/"test", b, c, provenance={"timestamp_sigma_ns": 0}, truth_z_mm=zb)
            subprocess.run([sys.executable, str(REPO/"studyDepthReconstruction.py"), "calibrate", "--input", str(d/"train"),
                            "--output", str(d/"cal.json")], check=True, capture_output=True, text=True)
            cmd = [sys.executable, str(REPO/"studyDepthReconstruction.py"), "scan", "--input", str(d/"test"),
                   "--calibration", str(d/"cal.json"), "--output", str(d/"scan"), "--bootstrap-draws", "8"]
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            report = json.loads((d/"scan/summary.json").read_text())
            self.assertEqual(len(report["results"]), 12)
            for r in report["results"]:
                self.assertTrue(r["evaluation"]["closure"])
                self.assertEqual(r["evaluation"]["input"], len(zb))
                with np.load(d/"scan"/(r["name"]+".npz"), allow_pickle=False) as f:
                    np.testing.assert_equal(f["photon_id"], b.photon_id)
                    np.testing.assert_equal(f["event_key"], b.event_key)
            self.assertNotEqual(subprocess.run(cmd, capture_output=True).returncode, 0)

    def test_unvalidated_truth_cannot_be_trained_by_cli(self):
        with TemporaryDirectory() as d:
            d = Path(d); a, c, z, k = fixture("x")
            save_cache(d/"input", a, c, provenance={"timestamp_sigma_ns": 0}, truth_z_mm=z, kappa_truth=k)
            p = subprocess.run([sys.executable, str(REPO/"studyDepthReconstruction.py"), "calibrate", "--input", str(d/"input"),
                                "--output", str(d/"cal.json")], capture_output=True, text=True)
            self.assertNotEqual(p.returncode, 0); self.assertIn("not validated", p.stderr)
            self.assertFalse((d/"cal.json").exists())

    def test_cache_corruption_rejected(self):
        with TemporaryDirectory() as d:
            d = Path(d); a, c, z, k = fixture("a")
            save_cache(d/"input", a, c, provenance={})
            with (d/"input/photons.npz").open("ab") as f:
                f.write(b"changed")
            with self.assertRaisesRegex(ValueError, "hash mismatch"): load_cache(d/"input")

    def test_two_depth_populations_do_not_force_a_mean_between_modes(self):
        a, c, z, k = fixture("cal")
        b, _, zb, _ = fixture("eval")
        cal = fit_calibration(a, k, c, source_sha256="b"*64)
        r = reconstruct(b, cal, c)
        # Loose structural regression, NOT a detector-performance threshold.
        self.assertLess(np.mean(np.isfinite(r.z_mm) & (r.z_mm > 800) & (r.z_mm < 1100)), .05)

    def test_raw_corrected_adapter_geometry_and_observable_separation(self):
        c, tface, _ = dreamsim_front_contract(pde=1.)
        f = np.arange(4, dtype=np.int32)
        calo = np.column_stack([2+FIBER_OFFSETS_MM[:, 0], 2+FIBER_OFFSETS_MM[:, 1], np.full(4, 1000.)])
        world = calo @ rotation()/10
        birth = calo.copy(); birth[:, 2] = -500.
        birth_world = birth @ rotation()/10
        mom = np.tile([0., 0., 1.], (4, 1)) @ rotation()
        data = dict(OP_time_final=c.forward_time(np.full(4, 500.), 1., calo[:, 0], calo[:, 1])+tface,
            OP_time_produced=np.zeros(4), OP_productionFiber=f, OP_finalFiber=f.copy(), OP_isCoreC=np.ones(4, bool))
        for j, axis in enumerate("xyz"):
            data["OP_pos_final_"+axis] = world[:, j]
            data["OP_pos_produced_"+axis] = birth_world[:, j]
            data["OP_mom_final_"+axis] = mom[:, j]
        e, ids, diagnostic, counts = corrected_event_from_arrays(data, c)
        np.testing.assert_allclose(e.local_x_mm, 0, atol=1e-10)
        np.testing.assert_allclose(diagnostic["truth_z_mm"], 500, atol=1e-10)
        np.testing.assert_allclose(diagnostic["exit_direction_kappa_proxy"], 1, atol=1e-10)
        batch = make_observations(e, "hash:0", c, entrance_time_ns=tface, photon_ids=ids)
        self.assertEqual(counts["selected"], 4); self.assertEqual(len(batch.photons.time_ns), 4)
        data["OP_finalFiber"][:] = -99
        _, _, _, unknown = corrected_event_from_arrays(data, c)
        self.assertEqual(unknown["unavailable_final_fiber"], 4)
        data["OP_finalFiber"][0] = 3
        with self.assertRaisesRegex(ValueError, "fiber mismatch"): corrected_event_from_arrays(data, c)


if __name__ == "__main__":
    unittest.main()
