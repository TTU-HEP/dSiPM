"""Matched sensor/timestamp observations for corrected FullDetector events.

This adapter does not read production_z_mm (truth). Supply corrected physical
fiber-local and sensor coordinates: legacy display shifts are NOT a geometry.
All sensor masks are formed BEFORE timestamp smearing. Entrance-clock error
is one common draw per event; photon errors are independent conditional draws.
"""
from dataclasses import dataclass
import hashlib

import numpy as np

from .depth_reconstruction import PhotonObservations


def _uniform(event_key, ids, seed, stream):
    """Counter-keyed variates, invariant to order, subsets and batching.

    SplitMix64 finalizer, high 53 bits; distinct SHA-256 stream keys. These are
    Monte Carlo variates, not cryptographic random numbers. Unlike one shared
    RNG, skipping zero-jitter draws cannot change another photon's PDE.
    """
    if seed < 0 or not isinstance(seed, int) or not event_key:
        raise ValueError("Need nonnegative integer seed and full event key")
    key = int.from_bytes(hashlib.sha256(f"{seed}:{stream}:{event_key}".encode()).digest()[:8], "little")
    x = np.asarray(ids, dtype=np.uint64) ^ np.uint64(key)
    with np.errstate(over="ignore"):
        x = x + np.uint64(0x9E3779B97F4A7C15)
        x = (x ^ (x >> 30)) * np.uint64(0xBF58476D1CE4E5B9)
        x = (x ^ (x >> 27)) * np.uint64(0x94D049BB133111EB)
        x ^= x >> 31
    return ((x >> 11).astype(float)+.5) / 2**53


def _normal(event_key, ids, seed, stream):
    u = _uniform(event_key, ids, seed, stream+"/u")
    v = _uniform(event_key, ids, seed, stream+"/v")
    return np.sqrt(-2*np.log(u))*np.cos(2*np.pi*v)


def smear_observations(obs, *, sigma_t_ns=0., entrance_sigma_ns=0., seed=20261021, draw=0):
    """Add jitter to a NOISELESS, already digitized cache, with fixed identities.

    This does not rerun PDE/holdoff or re-gate after smearing. The caller must
    not label an already smeared input as noiseless. The same standard-normal
    draws are scaled at every resolution (including zero) for a matched scan.
    """
    if (not np.isfinite([sigma_t_ns, entrance_sigma_ns]).all() or min(sigma_t_ns, entrance_sigma_ns) < 0 or
            not isinstance(draw, int) or draw < 0):
        raise ValueError("Invalid jitter")
    time = obs.time_ns.copy()
    for event in np.unique(obs.event_key):
        m = obs.event_key == event
        clock = entrance_sigma_ns*_normal(event, [0], seed, f"clock/{draw}")[0]
        time[m] += sigma_t_ns*_normal(event, obs.photon_id[m], seed, f"photon/{draw}")-clock
    columns = {k: getattr(obs, k) for k in obs.__dataclass_fields__}
    columns["time_ns"] = time
    return PhotonObservations(**columns)


@dataclass(frozen=True)
class ObservationBatch:
    photons: PhotonObservations
    accepted_indices: np.ndarray
    detected_mask: np.ndarray | None
    accepted_mask: np.ndarray
    in_gate_mask: np.ndarray
    entrance_error_ns: float

    @property
    def counts(self):
        detected = None if self.detected_mask is None else int(self.detected_mask.sum())
        return dict(input=len(self.in_gate_mask), in_gate=int(self.in_gate_mask.sum()),
                    detected=detected, accepted=int(self.accepted_mask.sum()),
                    lost_holdoff=None if detected is None else detected-int(self.accepted_mask.sum()))


def make_observations(event, event_key, contract, *, entrance_time_ns=0.,
                      sigma_t_ns=0., entrance_sigma_ns=0., seed=20261021, draw=0,
                      photon_ids=None, accepted_mask=None, external_detected_mask=None):
    """Digitize one full event, or reuse a previously validated accepted mask.

    event is a corrected FullDetectorPhotonEvent-compatible object. Only time,
    local x/y, sensor ID/centre, rod indices are read, never truth. Unique
    event_key convention: '<source ROOT SHA-256>:<original tree entry>'.
    draw affects jitter only; PDE and holdoff are identical for all resolutions
    AND jitter realizations at a fixed pitch. Input gate is strict (lo, hi).
    Existing digitization may be passed via accepted_mask. Unless its detected
    mask is also supplied, detected/lost counts are UNKNOWN, never fabricated.
    """
    if (not np.isfinite([entrance_time_ns, sigma_t_ns, entrance_sigma_ns]).all() or
            sigma_t_ns < 0 or entrance_sigma_ns < 0 or not isinstance(draw, int) or draw < 0):
        raise ValueError("Invalid timing options")
    t = np.asarray(event.time_ns, float)
    n = len(t)
    arrays = {name: np.asarray(getattr(event, name)) for name in
              ("local_x_mm", "local_y_mm", "sensor_id", "sensor_x_mm", "sensor_y_mm", "rod_ix", "rod_iy")}
    if t.ndim != 1 or not np.isfinite(t).all() or any(a.shape != (n,) for a in arrays.values()):
        raise ValueError("Malformed event columns")
    if not all(np.isfinite(a).all() for a in arrays.values()):
        raise ValueError("Nonfinite event")
    for name in ("sensor_id", "rod_ix", "rod_iy"):
        if arrays[name].dtype.kind not in "iu":
            raise ValueError("Sensor and rod identities must be integers")
    ids = np.arange(n, dtype=np.int64) if photon_ids is None else np.asarray(photon_ids)
    if ids.shape != (n,) or ids.dtype.kind not in "iu" or np.any(ids < 0) or len(np.unique(ids)) != n:
        raise ValueError("Need unique original integer photon IDs")
    x, y = arrays["local_x_mm"], arrays["local_y_mm"]
    if np.any(np.hypot(x, y) > contract.core_radius_mm+1e-5):
        raise ValueError("Not physical fiber-local endpoints; check geometry")
    gate = (t > contract.gate_ns[0]) & (t < contract.gate_ns[1])
    pitch = contract.pitch_um/1000
    px = np.floor((x-contract.pixel_origin_mm)/pitch).astype(np.int64)
    py = np.floor((y-contract.pixel_origin_mm)/pitch).astype(np.int64)
    if accepted_mask is None:
        if external_detected_mask is not None:
            raise ValueError("External detected mask requires an accepted mask")
        detected = gate & (_uniform(event_key, ids, seed, "pde-v1") < contract.pde)
        accepted = np.zeros(n, bool)
        di = np.flatnonzero(detected)
        order = di[np.lexsort((ids[di], t[di], py[di], px[di], arrays["sensor_id"][di]))]
        last_key, last_time = None, None
        for j in order:
            key = (int(arrays["sensor_id"][j]), int(px[j]), int(py[j]))
            if key != last_key or t[j]-last_time >= contract.recovery_ns:
                accepted[j] = True
                last_key, last_time = key, t[j]
    else:
        accepted = np.asarray(accepted_mask)
        if accepted.dtype != bool or accepted.shape != (n,) or np.any(accepted & ~gate):
            raise ValueError("Invalid external accepted mask/gate")
        detected = None if external_detected_mask is None else np.asarray(external_detected_mask)
        if detected is not None and (detected.dtype != bool or detected.shape != (n,) or
                                     np.any(accepted & ~detected) or np.any(detected & ~gate)):
            raise ValueError("External detected/accepted closure failure")
    ai = np.flatnonzero(accepted)
    exq = contract.pixel_origin_mm+(px[ai]+.5)*pitch
    eyq = contract.pixel_origin_mm+(py[ai]+.5)*pitch
    clock = float(entrance_sigma_ns*_normal(event_key, [0], seed, f"clock/{draw}")[0])
    time = t[ai]-entrance_time_ns-clock+sigma_t_ns*_normal(event_key, ids[ai], seed, f"photon/{draw}")
    obs = PhotonObservations(np.full(len(ai), event_key), ids[ai], arrays["rod_ix"][ai], arrays["rod_iy"][ai],
            time, np.hypot(exq, eyq)/contract.core_radius_mm,
            arrays["sensor_x_mm"][ai], arrays["sensor_y_mm"][ai])
    return ObservationBatch(obs, ai, detected, accepted, gate, clock)
