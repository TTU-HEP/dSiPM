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


UPSTREAM_PHOTON_FIELDS = ['OP_time_final', 'OP_time_produced', 'OP_productionFiber', 'OP_finalFiber', 'OP_isCoreC']
UPSTREAM_PHOTON_FIELDS += ['OP_pos_'+p+'_'+a for p in ('final', 'produced') for a in 'xyz']
UPSTREAM_PHOTON_FIELDS += ['OP_mom_final_'+a for a in 'xyz']
UPSTREAM_PHOTON_FIELDS += ['OP_isCoreQ', 'OP_isCerenkov', 'OP_productionRod', 'OP_productionLayer', 'OP_trackid']


def upstream_event_from_arrays(arrays, contract, detector, material):
    """Actual transported downstream exits only; NEVER substitute analytic times.

    Materials are separate reconstruction/calibration populations. Raw counts
    include all records so single-rod sampling cannot masquerade as full optical
    coverage. Geant4 copy numbers independently cross-check localized endpoints.
    """
    detector.check_contract(contract, material)
    data = {k: np.asarray(arrays[k]) for k in UPSTREAM_PHOTON_FIELDS}
    n = len(data['OP_time_final'])
    if any(a.shape != (n,) or not np.isfinite(a).all() for a in data.values()):
        raise ValueError('Malformed/nonfinite upstream photon branches')
    for name in ('OP_productionFiber', 'OP_finalFiber', 'OP_productionRod', 'OP_productionLayer', 'OP_trackid'):
        if data[name].dtype.kind not in 'iu':
            raise ValueError('Nonintegral upstream photon identity')
    for name in ('OP_isCoreC', 'OP_isCoreQ', 'OP_isCerenkov'):
        if not np.isin(data[name], [0, 1]).all():
            raise ValueError('Malformed upstream photon flags')
    if len(np.unique(data['OP_trackid'])) != n or np.any(data['OP_trackid'] <= 0):
        raise ValueError('Duplicate/invalid upstream track identity')
    if np.any(data['OP_isCoreC'].astype(bool) & data['OP_isCoreQ'].astype(bool)):
        raise ValueError('Ambiguous plastic/quartz identity')
    core = data['OP_isCoreC' if material == 'plastic' else 'OP_isCoreQ'].astype(bool)
    wanted = core & data['OP_isCerenkov'].astype(bool)
    sample = (data['OP_productionRod'] == detector.sample_rod_layer[0]) & (data['OP_productionLayer'] == detector.sample_rod_layer[1])
    recorded = data['OP_time_final'] > 0
    if np.any(recorded & ~sample):
        raise ValueError('Recorded optical exit outside source-defined sample rod/layer')
    copy = data['OP_productionFiber']
    allowed = [f[1] for f in detector.fibers if f[0] == material]
    if not np.isin(copy[core], allowed).all():
        raise ValueError('Core flag/fiber material mismatch')
    endpoints = np.column_stack([data['OP_pos_final_'+a] for a in 'xyz'])*10
    calo = world_to_calorimeter(endpoints, x_deg=detector.rotation_deg[0], y_deg=detector.rotation_deg[1])
    downstream = np.abs(calo[:, 2]-detector.length_mm/2) <= 1e-4
    gate = (data['OP_time_final'] > contract.gate_ns[0]) & (data['OP_time_final'] < contract.gate_ns[1])
    incident = wanted & recorded & downstream
    ids = np.flatnonzero(incident & gate)
    final_copy = data['OP_finalFiber'][ids]
    if np.any((final_copy != -99) & (final_copy != copy[ids])):
        raise ValueError('Production/final fiber mismatch')
    local, sensor, rods, sensor_ids, _ = localize_endpoints(endpoints[ids], copy[ids], detector=detector, material=material)
    if (not np.array_equal(rods[:, 0]+detector.n_rods//2, data['OP_productionRod'][ids]) or
            not np.array_equal(rods[:, 1]+detector.n_layers//2, data['OP_productionLayer'][ids])):
        raise ValueError('Localized endpoint and source replica identity disagree')
    birth = world_to_calorimeter(np.column_stack([data['OP_pos_produced_'+a][ids] for a in 'xyz'])*10,
        x_deg=detector.rotation_deg[0], y_deg=detector.rotation_deg[1])
    z = birth[:, 2]+detector.length_mm/2
    if (np.any(np.linalg.norm(birth[:, :2]-sensor, axis=1) > contract.core_radius_mm+1e-5) or
            np.any((z < -1e-5) | (z > contract.length_mm+1e-5))):
        raise ValueError('Source production point outside localized fiber core/length')
    momentum = np.column_stack([data['OP_mom_final_'+a][ids] for a in 'xyz']) @ rotation(*detector.rotation_deg).T
    norm = np.linalg.norm(momentum, axis=1)
    cosine = np.divide(momentum[:, 2], norm, out=np.full(len(ids), np.nan), where=norm > 0)
    event = SimpleNamespace(time_ns=data['OP_time_final'][ids], local_x_mm=local[:, 0], local_y_mm=local[:, 1],
        sensor_x_mm=sensor[:, 0], sensor_y_mm=sensor[:, 1], rod_ix=rods[:, 0], rod_iy=rods[:, 1], sensor_id=sensor_ids)
    truth = dict(truth_z_mm=np.clip(z, 0, detector.length_mm), production_xyz_mm=birth,
        recorded_first_step_time_ns=data['OP_time_produced'][ids],
        exit_direction_kappa_proxy=np.divide(1., cosine, out=np.full(len(ids), np.nan), where=cosine > 0))
    counts = dict(raw=n, selected=len(ids), rejected=n-len(ids), material_cherenkov_core=int(wanted.sum()),
        unrecorded_core=int((wanted & ~recorded).sum()), other_records=int((~wanted).sum()),
        recorded_core_exits=int((wanted & recorded).sum()), non_downstream_exits=int((wanted & recorded & ~downstream).sum()),
        downstream_before_gate=int(incident.sum()), outside_gate=int((incident & ~gate).sum()),
        non_sample_core=int((wanted & ~sample).sum()),
        unavailable_final_fiber=int((final_copy == -99).sum()))
    if (counts['raw'] != counts['other_records']+counts['unrecorded_core']+counts['recorded_core_exits'] or
            counts['recorded_core_exits'] != counts['non_downstream_exits']+counts['downstream_before_gate'] or
            counts['downstream_before_gate'] != counts['outside_gate']+counts['selected']):
        raise ValueError('Upstream raw/transport/window accounting closure failed')
    return event, ids, truth, counts


def corrected_event_from_arrays(arrays, contract, *, detector=None, material='plastic'):
    """Convert one tree entry; keep ORIGINAL photon indices and raw truth apart.

    Production fiber and corrected endpoints identify the physical fiber.
    The pinned producer leaves exitFiber=-99; that sentinel is counted as
    unavailable, NOT treated as either a mismatch or a measured identity.
    Stored production time is a first-step global time, NOT exact photon birth.
    Final direction yields an ideal-cylinder kappa PROXY, not measured truth.
    """
    if detector is not None:
        return upstream_event_from_arrays(arrays, contract, detector, material)
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


def iter_corrected_root(path, expected_sha256, contract, *, entry_start=0, entry_stop=None,
                        detector=None, material='plastic', provenance=None):
    """Read existing ROOT under Slurm and verify full input hash, before AND after.

    Deliberately no calibration auto-fit: the transport/birth-time truth gate
    needs a separate decision. uproot/awkward are optional read dependencies.
    """
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Raw HPCC photon work must run under Slurm")
    if entry_start < 0 or (entry_stop is not None and entry_stop <= entry_start):
        raise ValueError("Invalid entry range")
    if detector is not None:
        detector.check_contract(contract, material)
        if (not provenance or provenance.get('root_sha256') != expected_sha256 or
                provenance.get('detector_sha256') != detector.sha256 or
                provenance.get('producer_revision') != detector.revision):
            raise ValueError('Need ROOT hash tied to the exact DREAMSim source/run detector profile')
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
        if detector is not None:
            names = list(dict.fromkeys(names+UPSTREAM_PHOTON_FIELDS+['beamPX', 'beamPY', 'beamPZ']))
            # Current upstream assigns these members but does not export them.
            # Rotations therefore come from the pinned producer run profile;
            # if a later schema exports them, also cross-check every event.
            rotations = ['caloRotationX', 'caloRotationY']
            if any(k in tree.keys() for k in rotations):
                names += rotations
        if not set(names).issubset(tree.keys()):
            raise ValueError('Required ROOT branches missing')
        entry = entry_start
        for batch in tree.iterate(names, entry_start=entry_start, entry_stop=entry_stop, step_size="32 MB", library="ak"):
            for j in range(len(batch)):
                beam = detector.gun_world_mm if detector is not None else [25., -25., -1050.]
                if not np.allclose([float(batch[b][j]) for b in ("beamX", "beamY", "beamZ")], beam, atol=1e-8, rtol=0):
                    raise ValueError("Beam-front mismatch")
                if detector is not None:
                    if 'caloRotationX' in names and not np.allclose([float(batch[b][j]) for b in ('caloRotationX', 'caloRotationY')], detector.rotation_deg, atol=1e-8, rtol=0):
                        raise ValueError('Calorimeter rotation mismatch')
                    direction = np.array([float(batch[b][j]) for b in ('beamPX', 'beamPY', 'beamPZ')])
                    if not np.isfinite(direction).all() or direction[2] <= 0 or not np.allclose(direction/np.linalg.norm(direction), [0, 0, 1], atol=1e-8, rtol=0):
                        raise ValueError('Beam-direction mismatch')
                data = {k: ak.to_numpy(batch[k][j]) for k in names if k.startswith("OP_")}
                yield f"{expected_sha256}:{entry}", corrected_event_from_arrays(data, contract, detector=detector, material=material)
                entry += 1
    if file_sha256(path) != expected_sha256:
        raise ValueError("ROOT changed during read")
