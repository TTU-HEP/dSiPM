"""Strict compact-cache and existing-ROOT adapters for the depth study.

ROOT work on TTU HPCC is Slurm-only. No global monkeypatching of historical
readers and no changes to production/locked schema-v3 outputs. Truth fields
are returned separately, never included in PhotonObservations.
"""
from dataclasses import asdict
import hashlib
import json
import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from .geometry import localize_endpoints, world_to_calorimeter, rotation
from .depth import DepthContract, PhotonObservations


def file_sha256(path):
    with Path(path).open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def save_cache(directory, obs, contract, *, provenance, truth_z_mm=None, kappa_truth=None,
               calibration_eligible=False):
    """New directory only. An eligible kappa requires external physics validation.

    Do NOT mark a saved-direction proxy or (arrival-emission)/(L-z) from the
    known-bug transport as validated true path factor.
    """
    target = Path(directory)
    target.mkdir(parents=True, exist_ok=False)
    columns = {k: getattr(obs, k) for k in obs.__dataclass_fields__}
    for key, value in [("truth_z_mm", truth_z_mm), ("kappa_truth", kappa_truth)]:
        if value is not None:
            a = np.asarray(value, float)
            if a.shape != obs.time_ns.shape or not np.isfinite(a).all():
                raise ValueError("Misaligned/nonfinite cache truth")
            columns[key] = a
    if calibration_eligible and (kappa_truth is None or not provenance.get("truth_validation")):
        raise ValueError("Eligible calibration needs a documented truth validation")
    np.savez_compressed(target/"photons.npz", **columns)
    meta = dict(schema="depth-cache-v1", contract=asdict(contract), contract_sha256=contract.sha256,
                input_sha256=file_sha256(target/"photons.npz"), provenance=provenance,
                calibration_eligible=calibration_eligible, events=sorted(set(obs.event_key.tolist())))
    with (target/"metadata.json").open("x") as f:
        json.dump(meta, f, indent=2, allow_nan=False)
        f.write("\n")


def load_cache(directory):
    directory = Path(directory)
    meta = json.loads((directory/"metadata.json").read_text())
    if meta["schema"] != "depth-cache-v1" or file_sha256(directory/"photons.npz") != meta["input_sha256"]:
        raise ValueError("Cache schema/hash mismatch")
    contract = DepthContract(**meta["contract"])
    if contract.sha256 != meta["contract_sha256"]:
        raise ValueError("Cache contract mismatch")
    with np.load(directory/"photons.npz", allow_pickle=False) as f:
        obs = PhotonObservations(**{k: f[k] for k in PhotonObservations.__dataclass_fields__})
        truth = {k: f[k] for k in ("truth_z_mm", "kappa_truth") if k in f}
    if meta["events"] != sorted(set(obs.event_key.tolist())):
        raise ValueError("Cache event manifest mismatch")
    return obs, contract, truth, meta


def corrected_event_from_arrays(arrays, contract):
    """Convert one tree entry; keep ORIGINAL photon indices and raw truth apart.

    Production fiber and corrected endpoints identify the physical fiber.
    The pinned producer leaves exitFiber=-99; that sentinel is counted as
    unavailable, NOT treated as either a mismatch or a measured identity.
    Stored production time is a first-step global time, NOT exact photon birth.
    Final direction yields an ideal-cylinder kappa PROXY, not measured truth.
    """
    if not contract.geometry_id.startswith("DREAMSim-54b2447-RyRx-2-2-"):
        raise ValueError("ROOT adapter supports only the pinned 2-degree production geometry")
    names = ["OP_time_final", "OP_time_produced", "OP_productionFiber", "OP_finalFiber", "OP_isCoreC"]
    names += ["OP_pos_"+where+"_"+a for where in ("final", "produced") for a in "xyz"]
    names += ["OP_mom_final_"+a for a in "xyz"]
    data = {k: np.asarray(arrays[k]) for k in names}
    n = len(data["OP_time_final"])
    if any(a.shape != (n,) or not np.isfinite(a).all() for a in data.values()):
        raise ValueError("Malformed/nonfinite ROOT photon branches")
    fiber = data["OP_productionFiber"]
    t = data["OP_time_final"]
    keep = ((fiber >= 0) & (fiber < 4) & data["OP_isCoreC"].astype(bool) &
            (data["OP_pos_final_z"] > 0) & (t > contract.gate_ns[0]) & (t < contract.gate_ns[1]))
    ids = np.flatnonzero(keep)
    final_fiber = data["OP_finalFiber"][ids]
    unknown_final = final_fiber == -99
    if np.any((~unknown_final) & (final_fiber != fiber[ids])):
        raise ValueError("Production/final fiber mismatch: investigate transport, do not remap")
    points = np.column_stack([data["OP_pos_final_"+a][ids] for a in "xyz"])*10
    local, sensor, rods, sensor_ids, _ = localize_endpoints(points, fiber[ids])
    birth = world_to_calorimeter(np.column_stack([data["OP_pos_produced_"+a][ids] for a in "xyz"])*10)
    if np.any(np.linalg.norm(birth[:, :2]-sensor, axis=1) > contract.core_radius_mm+1e-5):
        raise ValueError("Production point not in localized physical core")
    depth = birth[:, 2]+contract.length_mm/2
    if np.any((depth < -1e-5) | (depth > contract.length_mm+1e-5)):
        raise ValueError("Production depth outside detector")
    # Only geometry roundoff is clamped here, NOT any reconstructed estimate.
    depth = np.clip(depth, 0, contract.length_mm)
    momentum = np.column_stack([data["OP_mom_final_"+a][ids] for a in "xyz"]) @ rotation().T
    norm = np.linalg.norm(momentum, axis=1)
    cosine = np.divide(momentum[:, 2], norm, out=np.full(len(ids), np.nan), where=norm > 0)
    kproxy = np.divide(1., cosine, out=np.full(len(ids), np.nan), where=cosine > 0)
    event = SimpleNamespace(time_ns=t[ids], local_x_mm=local[:, 0], local_y_mm=local[:, 1],
        sensor_x_mm=sensor[:, 0], sensor_y_mm=sensor[:, 1], rod_ix=rods[:, 0], rod_iy=rods[:, 1], sensor_id=sensor_ids)
    diagnostics = dict(truth_z_mm=depth, production_xyz_mm=birth,
        recorded_first_step_time_ns=data["OP_time_produced"][ids], exit_direction_kappa_proxy=kproxy)
    return event, ids, diagnostics, dict(raw=n, selected=len(ids), rejected=n-len(ids),
                                       unavailable_final_fiber=int(unknown_final.sum()))


def iter_corrected_root(path, expected_sha256, contract, *, entry_start=0, entry_stop=None):
    """Read existing ROOT under Slurm and verify full input hash, before AND after.

    Deliberately no calibration auto-fit: the transport/birth-time truth gate
    needs a separate decision. uproot/awkward are optional read dependencies.
    """
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Raw HPCC photon work must run under Slurm")
    if entry_start < 0 or (entry_stop is not None and entry_stop <= entry_start):
        raise ValueError("Invalid entry range")
    if file_sha256(path) != expected_sha256:
        raise ValueError("ROOT input checksum mismatch")
    import awkward as ak
    import uproot
    with uproot.open(path) as f:
        tree = f["tree"]
        names = ["OP_time_final", "OP_time_produced", "OP_productionFiber", "OP_finalFiber", "OP_isCoreC"]
        names += ["OP_pos_"+where+"_"+a for where in ("final", "produced") for a in "xyz"]
        names += ["OP_mom_final_"+a for a in "xyz"]
        names += ["beamX", "beamY", "beamZ"]
        entry = entry_start
        for batch in tree.iterate(names, entry_start=entry_start, entry_stop=entry_stop, step_size="32 MB", library="ak"):
            for j in range(len(batch)):
                if not np.allclose([float(batch[b][j]) for b in ("beamX", "beamY", "beamZ")], [25., -25., -1050.], atol=1e-8, rtol=0):
                    raise ValueError("Beam-front mismatch")
                data = {k: ak.to_numpy(batch[k][j]) for k in names if k.startswith("OP_")}
                yield f"{expected_sha256}:{entry}", corrected_event_from_arrays(data, contract)
                entry += 1
    if file_sha256(path) != expected_sha256:
        raise ValueError("ROOT changed during read")
