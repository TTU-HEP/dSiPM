"""Sensor response: independent-fiber FLASH, analog recharge and readout proxies.

Historical seeded reference APIs are preserved for numerical closure; matched
depth observations use event/photon-keyed streams. Arrival times are not voltage
waveforms; scope_waveform is an explicitly convolved, uncalibrated proxy.
"""
from __future__ import annotations
from collections.abc import Sequence
from dataclasses import dataclass
import hashlib
import numpy as np
from .depth import PhotonObservations

def _apply_nonparalyzable_dead_time(
    cell_ids: np.ndarray,
    time_ns: np.ndarray,
    dead_time_ns: float | None,
) -> np.ndarray:
    """Return a mask of avalanches accepted by a non-paralyzable SPAD."""

    if dead_time_ns is None or dead_time_ns <= 0 or time_ns.size == 0:
        return np.ones(time_ns.size, dtype=bool)

    accepted = np.zeros(time_ns.size, dtype=bool)
    last_accepted: dict[int, float] = {}
    for index, (cell_id, arrival_time) in enumerate(zip(cell_ids, time_ns, strict=True)):
        previous = last_accepted.get(int(cell_id))
        if previous is None or float(arrival_time) - previous >= dead_time_ns:
            accepted[index] = True
            last_accepted[int(cell_id)] = float(arrival_time)
    return accepted


SENSOR_HALF_WIDTH_MM = 1.5


@dataclass(frozen=True)
class FullDetectorPhotonEvent:
    """Selected Cherenkov photons and event-level calorimeter observables."""

    event_index: int
    beam_energy_gev: float
    c_signal: float
    s_signal: float
    ecalo_truth_gev: float
    eleak_truth_gev: float
    raw_photon_count: int
    local_x_mm: np.ndarray
    local_y_mm: np.ndarray
    sensor_x_mm: np.ndarray
    sensor_y_mm: np.ndarray
    rod_ix: np.ndarray
    rod_iy: np.ndarray
    sensor_id: np.ndarray
    time_ns: np.ndarray
    wavelength_nm: np.ndarray
    production_z_mm: np.ndarray

    @property
    def selected_photon_count(self) -> int:
        return int(self.time_ns.size)

    @property
    def active_sensor_count(self) -> int:
        return int(np.unique(self.sensor_id).size)

    @property
    def max_sensor_photons(self) -> int:
        if self.sensor_id.size == 0:
            return 0
        _, counts = np.unique(self.sensor_id, return_counts=True)
        return int(np.max(counts))


@dataclass(frozen=True)
class FullDetectorScenario:
    """Constant-PDE response point used in the full-detector scan."""

    pde: float
    pitch_um: float
    recovery_ns: float

    @property
    def name(self) -> str:
        return (
            f"pde{int(round(100 * self.pde)):02d}_"
            f"pitch{self.pitch_um:g}um_recovery{self.recovery_ns:g}ns"
        )


@dataclass(frozen=True)
class AnalogSiPMScenario:
    """Representative conventional analog-SiPM response parameters.

    ``microcell_size_um`` is the manufacturer's nominal microcell size.  The
    spatial model uses the nearest square grid to ``n_microcells`` across the
    active area, which preserves the measured cell count and the local photon
    density more faithfully than treating the nominal active-cell size as the
    lattice pitch.
    """

    name: str
    pde: float
    microcell_size_um: float
    n_microcells: int
    recovery_tau_ns: float
    sensor_width_mm: float = 3.0
    fast_output_fwhm_ns: float | None = None
    pde_wavelength_nm: float | None = None
    overvoltage_v: float | None = None
    source_url: str | None = None

    @property
    def grid_side(self) -> int:
        return int(round(np.sqrt(self.n_microcells)))

    @property
    def modeled_microcells(self) -> int:
        return self.grid_side**2

    @property
    def effective_pitch_um(self) -> float:
        return 1000.0 * self.sensor_width_mm / self.grid_side


@dataclass(frozen=True)
class FullDetectorDigitization:
    """Counts and masks for one full-detector response point."""

    detected_mask: np.ndarray
    accepted_mask: np.ndarray
    n_detected: int
    n_accepted: int

    @property
    def recovery_survival(self) -> float:
        return self.n_accepted / self.n_detected if self.n_detected else 0.0


@dataclass(frozen=True)
class AnalogSiPMDigitization:
    """Charge response for a conventional summed analog SiPM."""

    detected_mask: np.ndarray
    charge_per_detected_photon: np.ndarray
    n_detected: int
    charge_pe: float

    @property
    def charge_survival(self) -> float:
        return self.charge_pe / self.n_detected if self.n_detected else 0.0


def global_cell_ids(
    event: FullDetectorPhotonEvent,
    pitch_um: float,
    *,
    sensor_half_width_mm: float = SENSOR_HALF_WIDTH_MM,
    mask: np.ndarray | None = None,
) -> np.ndarray:
    """Return a unique SPAD-cell id across every fiber in the detector."""

    if pitch_um <= 0:
        raise ValueError("pitch_um must be positive")
    if sensor_half_width_mm <= 0:
        raise ValueError("sensor_half_width_mm must be positive")
    select = np.ones(event.time_ns.size, dtype=bool) if mask is None else np.asarray(mask, dtype=bool)
    if select.shape != event.time_ns.shape:
        raise ValueError("mask must match the photon arrays")
    x = event.local_x_mm[select]
    y = event.local_y_mm[select]
    sensor = event.sensor_id[select]
    pitch_mm = pitch_um / 1000.0
    width_mm = 2.0 * sensor_half_width_mm
    inside = (
        (x >= -sensor_half_width_mm)
        & (x < sensor_half_width_mm)
        & (y >= -sensor_half_width_mm)
        & (y < sensor_half_width_mm)
    )
    if not np.all(inside):
        raise ValueError("localized photon falls outside the modeled sensor")
    n_axis = int(np.ceil(width_mm / pitch_mm))
    ix = np.floor((x + sensor_half_width_mm) / pitch_mm).astype(np.int64)
    iy = np.floor((y + sensor_half_width_mm) / pitch_mm).astype(np.int64)
    local_cell = iy * n_axis + ix
    return sensor * np.int64(n_axis * n_axis) + local_cell


def analog_microcell_ids(
    event: FullDetectorPhotonEvent,
    scenario: AnalogSiPMScenario,
    *,
    mask: np.ndarray | None = None,
) -> np.ndarray:
    """Return spatially local microcell ids for a conventional analog SiPM."""

    if scenario.n_microcells <= 0:
        raise ValueError("n_microcells must be positive")
    if scenario.sensor_width_mm <= 0:
        raise ValueError("sensor_width_mm must be positive")
    select = np.ones(event.time_ns.size, dtype=bool) if mask is None else np.asarray(mask, dtype=bool)
    if select.shape != event.time_ns.shape:
        raise ValueError("mask must match the photon arrays")
    half_width = 0.5 * scenario.sensor_width_mm
    x = event.local_x_mm[select]
    y = event.local_y_mm[select]
    sensor = event.sensor_id[select]
    inside = (
        (x >= -half_width)
        & (x < half_width)
        & (y >= -half_width)
        & (y < half_width)
    )
    if not np.all(inside):
        raise ValueError("localized photon falls outside the modeled analog SiPM")
    side = scenario.grid_side
    ix = np.floor((x + half_width) * side / scenario.sensor_width_mm).astype(np.int64)
    iy = np.floor((y + half_width) * side / scenario.sensor_width_mm).astype(np.int64)
    local_cell = iy * side + ix
    return sensor * np.int64(side * side) + local_cell


def digitize_analog_sipm_event(
    event: FullDetectorPhotonEvent,
    scenario: AnalogSiPMScenario,
    *,
    seed: int = 4307,
) -> AnalogSiPMDigitization:
    """Apply PDE and exponential partial recharge to a summed analog SiPM.

    The first avalanche in a microcell contributes one photoelectron of
    charge.  A later avalanche after ``dt`` contributes
    ``1 - exp(-dt / tau)`` and resets that microcell.  PDE selection is made
    before recharge, so this is an optimistic charge-only model: the reduced
    triggering probability of a partially recharged cell, optical crosstalk,
    afterpulsing, electronics noise, and waveform integration are not yet
    included.
    """

    if not 0.0 <= scenario.pde <= 1.0:
        raise ValueError("pde must be between zero and one")
    if scenario.recovery_tau_ns < 0:
        raise ValueError("recovery_tau_ns cannot be negative")
    random = np.random.default_rng(
        np.random.SeedSequence([seed, event.event_index, 1729])
    )
    detected_mask = random.random(event.time_ns.size) < scenario.pde
    detected_time = event.time_ns[detected_mask]
    n_detected = int(detected_time.size)
    if n_detected == 0:
        return AnalogSiPMDigitization(
            detected_mask=detected_mask,
            charge_per_detected_photon=np.zeros(0, dtype=np.float64),
            n_detected=0,
            charge_pe=0.0,
        )

    cell_ids = analog_microcell_ids(event, scenario, mask=detected_mask)
    charge = np.ones(n_detected, dtype=np.float64)
    if scenario.recovery_tau_ns > 0:
        order = np.lexsort((detected_time, cell_ids))
        ordered_cells = cell_ids[order]
        ordered_time = detected_time[order]
        same_cell = ordered_cells[1:] == ordered_cells[:-1]
        delta_time = ordered_time[1:] - ordered_time[:-1]
        ordered_charge = np.ones(n_detected, dtype=np.float64)
        repeated_indices = np.flatnonzero(same_cell) + 1
        ordered_charge[repeated_indices] = -np.expm1(
            -delta_time[same_cell] / scenario.recovery_tau_ns
        )
        charge[order] = ordered_charge
    return AnalogSiPMDigitization(
        detected_mask=detected_mask,
        charge_per_detected_photon=charge,
        n_detected=n_detected,
        charge_pe=float(np.sum(charge)),
    )


def digitize_full_detector_event(
    event: FullDetectorPhotonEvent,
    scenario: FullDetectorScenario,
    *,
    seed: int = 4307,
) -> FullDetectorDigitization:
    """Apply constant PDE and non-paralyzable recovery to one event."""

    if not 0.0 <= scenario.pde <= 1.0:
        raise ValueError("pde must be between zero and one")
    if scenario.recovery_ns < 0:
        raise ValueError("recovery_ns cannot be negative")
    random = np.random.default_rng(
        np.random.SeedSequence([seed, event.event_index, 1729])
    )
    detected_mask = random.random(event.time_ns.size) < scenario.pde
    detected_indices = np.flatnonzero(detected_mask)
    accepted_mask = np.zeros(event.time_ns.size, dtype=bool)
    if detected_indices.size:
        cell_ids = global_cell_ids(event, scenario.pitch_um, mask=detected_mask)
        accepted_detected = _apply_nonparalyzable_dead_time(
            cell_ids,
            event.time_ns[detected_mask],
            scenario.recovery_ns,
        )
        accepted_mask[detected_indices[accepted_detected]] = True
    return FullDetectorDigitization(
        detected_mask=detected_mask,
        accepted_mask=accepted_mask,
        n_detected=int(detected_indices.size),
        n_accepted=int(np.count_nonzero(accepted_mask)),
    )


def scan_full_detector_event(
    event: FullDetectorPhotonEvent,
    *,
    pde: float,
    pitches_um: Sequence[float],
    recovery_times_ns: Sequence[float],
    seed: int = 4307,
) -> tuple[int, np.ndarray]:
    """Return detected count and accepted counts for a pitch/recovery grid."""

    pitches = np.asarray(pitches_um, dtype=np.float64)
    recovery = np.asarray(recovery_times_ns, dtype=np.float64)
    if not 0.0 <= pde <= 1.0:
        raise ValueError("pde must be between zero and one")
    if np.any(pitches <= 0):
        raise ValueError("pitches must be positive")
    if np.any(recovery < 0):
        raise ValueError("recovery times cannot be negative")

    random = np.random.default_rng(
        np.random.SeedSequence([seed, event.event_index, 1729])
    )
    detected_mask = random.random(event.time_ns.size) < pde
    detected_count = int(np.count_nonzero(detected_mask))
    accepted = np.zeros((pitches.size, recovery.size), dtype=np.int64)
    if detected_count == 0:
        return detected_count, accepted
    detected_time = event.time_ns[detected_mask]
    for pitch_index, pitch_um in enumerate(pitches):
        cell_ids = global_cell_ids(event, float(pitch_um), mask=detected_mask)
        accepted[pitch_index] = accepted_counts_by_recovery(
            cell_ids,
            detected_time,
            recovery,
        )
    return detected_count, accepted


def accepted_counts_by_recovery(
    cell_ids: np.ndarray,
    time_ns: np.ndarray,
    recovery_times_ns: Sequence[float],
) -> np.ndarray:
    """Return accepted photon counts for several non-paralyzable recoveries.

    Photons are grouped by SPAD.  Sorting is stable, so a globally
    time-ordered event remains time ordered within every SPAD.  Only cells
    with multiple photons need a recovery calculation; singleton cells are
    accepted for every recovery time.
    """

    cell_ids = np.asarray(cell_ids, dtype=np.int64)
    time_ns = np.asarray(time_ns, dtype=np.float64)
    recovery = np.asarray(recovery_times_ns, dtype=np.float64)
    if cell_ids.shape != time_ns.shape:
        raise ValueError("cell_ids and time_ns must have the same shape")
    if np.any(recovery < 0):
        raise ValueError("recovery times cannot be negative")
    if cell_ids.size == 0:
        return np.zeros(recovery.size, dtype=np.int64)

    order = np.argsort(cell_ids, kind="stable")
    sorted_cells = cell_ids[order]
    sorted_times = time_ns[order]
    starts = np.r_[0, np.flatnonzero(sorted_cells[1:] != sorted_cells[:-1]) + 1]
    ends = np.r_[starts[1:], sorted_cells.size]

    accepted = np.full(recovery.size, starts.size, dtype=np.int64)
    accepted[recovery <= 0] = cell_ids.size
    finite_indices = np.flatnonzero(recovery > 0)

    for start, end in zip(starts, ends, strict=True):
        if end - start < 2:
            continue
        times = sorted_times[start:end]
        if not np.all(times[:-1] <= times[1:]):
            times = np.sort(times, kind="stable")
        for recovery_index in finite_indices:
            holdoff = recovery[recovery_index]
            n_cell = 1
            last_accepted = times[0]
            for arrival in times[1:]:
                if arrival - last_accepted >= holdoff:
                    n_cell += 1
                    last_accepted = arrival
            accepted[recovery_index] += n_cell - 1
    return accepted


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


def partial_recharge_weights(
    cell_ids: np.ndarray,
    time_ns: np.ndarray,
    recovery_tau_ns: float,
) -> np.ndarray:
    """Return analog charge weights for avalanches in recovering microcells.

    The first avalanche in a cell has unit charge.  Each later avalanche
    contributes ``1 - exp(-dt / tau)`` where ``dt`` is measured from the
    preceding avalanche in that cell.  Inputs need not be ordered.
    """

    cell_ids = np.asarray(cell_ids, dtype=np.int64)
    time_ns = np.asarray(time_ns, dtype=np.float64)
    if cell_ids.shape != time_ns.shape:
        raise ValueError("cell_ids and time_ns must have the same shape")
    if recovery_tau_ns < 0:
        raise ValueError("recovery_tau_ns cannot be negative")
    if time_ns.size == 0 or recovery_tau_ns == 0:
        return np.ones(time_ns.size, dtype=np.float64)

    order = np.lexsort((time_ns, cell_ids))
    ordered_cells = cell_ids[order]
    ordered_time = time_ns[order]
    ordered_weights = np.ones(time_ns.size, dtype=np.float64)
    same_cell = ordered_cells[1:] == ordered_cells[:-1]
    repeated = np.flatnonzero(same_cell) + 1
    delta_time = ordered_time[1:] - ordered_time[:-1]
    ordered_weights[repeated] = -np.expm1(
        -delta_time[same_cell] / recovery_tau_ns
    )
    weights = np.empty(time_ns.size, dtype=np.float64)
    weights[order] = ordered_weights
    return weights


def scope_kernel(
    bin_width_ns: float,
    *,
    rise_10_90_ns: float = 0.6,
    fwhm_ns: float = 1.5,
    duration_ns: float = 12.0,
) -> np.ndarray:
    """Return a causal, unit-area fast-output display kernel.

    The leading side is Gaussian and the falling side exponential.  Its
    parameters are chosen analytically so that the sampled shape approximates
    the requested 10--90% rise time and FWHM.  The default values are the
    MicroFC-30035-SMT fast-output specifications.  This is not an electronics
    simulation and therefore carries no voltage calibration.
    """

    if bin_width_ns <= 0:
        raise ValueError("bin_width_ns must be positive")
    if rise_10_90_ns <= 0 or fwhm_ns <= 0 or duration_ns <= 0:
        raise ValueError("kernel timing parameters must be positive")

    gaussian_10 = np.sqrt(2.0 * np.log(10.0))
    gaussian_90 = np.sqrt(2.0 * np.log(1.0 / 0.9))
    sigma = rise_10_90_ns / (gaussian_10 - gaussian_90)
    gaussian_half_width = sigma * np.sqrt(2.0 * np.log(2.0))
    tail_half_width = fwhm_ns - gaussian_half_width
    if tail_half_width <= 0:
        raise ValueError("FWHM is too short for the requested rise time")
    decay_tau = tail_half_width / np.log(2.0)

    # Put the 10% rising crossing after t=0 so the kernel is causal without
    # clipping the specified leading edge.
    peak_time = gaussian_10 * sigma + 0.5 * bin_width_ns
    time = np.arange(0.0, duration_ns + 0.5 * bin_width_ns, bin_width_ns)
    before_peak = time <= peak_time
    values = np.empty_like(time)
    values[before_peak] = np.exp(
        -0.5 * ((time[before_peak] - peak_time) / sigma) ** 2
    )
    values[~before_peak] = np.exp(
        -(time[~before_peak] - peak_time) / decay_tau
    )
    values /= np.sum(values)
    return values


def scope_waveform(
    impulse_histogram: np.ndarray,
    *,
    bin_width_ns: float,
    rise_10_90_ns: float = 0.6,
    fwhm_ns: float = 1.5,
    duration_ns: float = 12.0,
) -> np.ndarray:
    """Convolve a binned avalanche/charge stream with the display kernel."""

    impulse_histogram = np.asarray(impulse_histogram, dtype=np.float64)
    if impulse_histogram.ndim != 1:
        raise ValueError("impulse_histogram must be one-dimensional")
    kernel = scope_kernel(
        bin_width_ns,
        rise_10_90_ns=rise_10_90_ns,
        fwhm_ns=fwhm_ns,
        duration_ns=duration_ns,
    )
    return np.convolve(impulse_histogram, kernel, mode="full")[: impulse_histogram.size]


def histogram_observables(
    histogram: np.ndarray,
    edges_ns: np.ndarray,
    *,
    prompt_end_ns: float,
) -> dict[str, float]:
    """Return integral, prompt fraction, mean time, and central-80% width."""

    histogram = np.asarray(histogram, dtype=np.float64)
    edges_ns = np.asarray(edges_ns, dtype=np.float64)
    if edges_ns.size != histogram.size + 1:
        raise ValueError("edges_ns must contain one more value than histogram")
    if np.any(histogram < 0):
        raise ValueError("histogram cannot contain negative weights")
    centers = 0.5 * (edges_ns[:-1] + edges_ns[1:])
    total = float(np.sum(histogram))
    if total == 0:
        return {
            "integral": 0.0,
            "prompt_fraction": 0.0,
            "mean_time_ns": float("nan"),
            "central80_width_ns": float("nan"),
        }
    cumulative = np.cumsum(histogram)
    q10_index = int(np.searchsorted(cumulative, 0.10 * total, side="left"))
    q90_index = int(np.searchsorted(cumulative, 0.90 * total, side="left"))
    return {
        "integral": total,
        "prompt_fraction": float(np.sum(histogram[centers < prompt_end_ns]) / total),
        "mean_time_ns": float(np.sum(histogram * centers) / total),
        "central80_width_ns": float(centers[q90_index] - centers[q10_index]),
    }


@dataclass(frozen=True)
class ChannelReadoutConfig:
    """Resources shared by a square group of SPADs."""

    spads_per_channel_side: int
    timestamp_buffer_depth: int = 8
    channel_resolving_time_ns: float = 0.05
    count_bits: int = 8
    timestamp_bits: int = 12
    address_bits: int = 12
    event_overhead_bits: int = 64

    def validate(self) -> None:
        if self.spads_per_channel_side < 1:
            raise ValueError("spads_per_channel_side must be positive")
        if self.timestamp_buffer_depth < 0:
            raise ValueError("timestamp_buffer_depth cannot be negative")
        if self.channel_resolving_time_ns < 0:
            raise ValueError("channel_resolving_time_ns cannot be negative")
        for field_name in (
            "count_bits",
            "timestamp_bits",
            "address_bits",
            "event_overhead_bits",
        ):
            if getattr(self, field_name) < 0:
                raise ValueError(f"{field_name} cannot be negative")


@dataclass(frozen=True)
class ChannelReadoutResult:
    """Count and timing output from one event after SPAD recovery."""

    accepted_hits: int
    active_channels: int
    counted_hits: int
    count_saturated_channels: int
    timestamp_candidates: int
    timestamp_records: int
    resolving_lost_hits: int
    buffer_lost_records: int
    overflowed_channels: int
    max_hits_in_channel: int
    payload_bits: int

    @property
    def count_survival(self) -> float:
        return self.counted_hits / self.accepted_hits if self.accepted_hits else 1.0

    @property
    def timestamp_survival(self) -> float:
        return (
            self.timestamp_records / self.accepted_hits
            if self.accepted_hits
            else 1.0
        )

    @property
    def any_buffer_overflow(self) -> bool:
        return self.overflowed_channels > 0


def _keep_first_per_group(
    group_ids: np.ndarray,
    maximum_per_group: int,
) -> np.ndarray:
    """Return a mask retaining the first ``maximum_per_group`` group entries."""

    if group_ids.size == 0 or maximum_per_group <= 0:
        return np.zeros(group_ids.size, dtype=bool)
    order = np.argsort(group_ids, kind="stable")
    sorted_groups = group_ids[order]
    starts = np.r_[0, np.flatnonzero(sorted_groups[1:] != sorted_groups[:-1]) + 1]
    ends = np.r_[starts[1:], sorted_groups.size]
    keep = np.zeros(group_ids.size, dtype=bool)
    for start, end in zip(starts, ends, strict=True):
        keep[order[start : min(end, start + maximum_per_group)]] = True
    return keep


def channel_readout(
    spad_ix: np.ndarray,
    spad_iy: np.ndarray,
    time_ns: np.ndarray,
    *,
    n_spads_axis: int,
    config: ChannelReadoutConfig,
) -> ChannelReadoutResult:
    """Apply shared-channel resolving time, counters, buffers, and payload.

    Inputs are avalanches already accepted by the individual-SPAD recovery
    model.  Each active channel exports one address and count word plus up to
    ``timestamp_buffer_depth`` timestamp words.  The counter integrates every
    accepted avalanche until its finite bit depth saturates.
    """

    config.validate()
    spad_ix = np.asarray(spad_ix, dtype=np.int64)
    spad_iy = np.asarray(spad_iy, dtype=np.int64)
    time_ns = np.asarray(time_ns, dtype=np.float64)
    if spad_ix.shape != spad_iy.shape or spad_ix.shape != time_ns.shape:
        raise ValueError("SPAD indices and time_ns must have the same shape")
    if n_spads_axis < 1:
        raise ValueError("n_spads_axis must be positive")
    if spad_ix.size and (
        np.any(spad_ix < 0)
        or np.any(spad_iy < 0)
        or np.any(spad_ix >= n_spads_axis)
        or np.any(spad_iy >= n_spads_axis)
    ):
        raise ValueError("SPAD indices fall outside the sensor grid")

    accepted_hits = int(time_ns.size)
    if not accepted_hits:
        return ChannelReadoutResult(
            accepted_hits=0,
            active_channels=0,
            counted_hits=0,
            count_saturated_channels=0,
            timestamp_candidates=0,
            timestamp_records=0,
            resolving_lost_hits=0,
            buffer_lost_records=0,
            overflowed_channels=0,
            max_hits_in_channel=0,
            payload_bits=config.event_overhead_bits,
        )

    order = np.argsort(time_ns, kind="stable")
    time_ns = time_ns[order]
    spad_ix = spad_ix[order]
    spad_iy = spad_iy[order]

    group_side = config.spads_per_channel_side
    n_channels_axis = int(np.ceil(n_spads_axis / group_side))
    channel_ids = (
        (spad_iy // group_side) * n_channels_axis + (spad_ix // group_side)
    )
    _, channel_hit_counts = np.unique(channel_ids, return_counts=True)
    active_channels = int(channel_hit_counts.size)
    max_hits = int(np.max(channel_hit_counts))

    maximum_count = (1 << config.count_bits) - 1 if config.count_bits else 0
    counted_hits = int(np.minimum(channel_hit_counts, maximum_count).sum())
    count_saturated = int(np.count_nonzero(channel_hit_counts > maximum_count))

    resolving_mask = _apply_nonparalyzable_dead_time(
        channel_ids,
        time_ns,
        config.channel_resolving_time_ns,
    )
    candidate_channels = channel_ids[resolving_mask]
    timestamp_candidates = int(candidate_channels.size)
    _, candidate_counts = np.unique(candidate_channels, return_counts=True)
    overflowed_channels = int(
        np.count_nonzero(candidate_counts > config.timestamp_buffer_depth)
    )
    buffer_mask = _keep_first_per_group(
        candidate_channels,
        config.timestamp_buffer_depth,
    )
    timestamp_records = int(np.count_nonzero(buffer_mask))

    payload_bits = (
        config.event_overhead_bits
        + active_channels * (config.address_bits + config.count_bits)
        + timestamp_records * config.timestamp_bits
    )
    return ChannelReadoutResult(
        accepted_hits=accepted_hits,
        active_channels=active_channels,
        counted_hits=counted_hits,
        count_saturated_channels=count_saturated,
        timestamp_candidates=timestamp_candidates,
        timestamp_records=timestamp_records,
        resolving_lost_hits=accepted_hits - timestamp_candidates,
        buffer_lost_records=timestamp_candidates - timestamp_records,
        overflowed_channels=overflowed_channels,
        max_hits_in_channel=max_hits,
        payload_bits=int(payload_bits),
    )
