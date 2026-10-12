"""Deterministic photon-to-digital response model.

This module contains no ROOT code, which keeps the detector behavior testable.
Every device sees the same per-photon random variates for a given seed and
event number.  That makes pairwise device comparisons deterministic and
removes the extra Monte Carlo noise present in the original student workflow.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .config import DeviceConfig


HC_GEV_NM = 1.2398e-6
LEGACY_QE_COEFFICIENTS = np.asarray(
    [
        8.3089805235e-14,
        -1.7718571478e-10,
        1.4724264010e-07,
        -6.0692493107e-05,
        1.3083151701e-02,
        -1.0675291442,
    ],
    dtype=np.float64,
)


@dataclass(frozen=True)
class PhotonEvent:
    """Selected photons incident on the simulated detector plane."""

    event_index: int
    raw_photon_count: int
    x_mm: np.ndarray
    y_mm: np.ndarray
    time_ns: np.ndarray
    wavelength_nm: np.ndarray
    production_z_mm: np.ndarray

    @property
    def selected_photon_count(self) -> int:
        return int(self.time_ns.size)

    def sorted_by_time(self) -> "PhotonEvent":
        if self.time_ns.size < 2 or np.all(self.time_ns[:-1] <= self.time_ns[1:]):
            return self
        order = np.argsort(self.time_ns, kind="stable")
        return PhotonEvent(
            event_index=self.event_index,
            raw_photon_count=self.raw_photon_count,
            x_mm=self.x_mm[order],
            y_mm=self.y_mm[order],
            time_ns=self.time_ns[order],
            wavelength_nm=self.wavelength_nm[order],
            production_z_mm=self.production_z_mm[order],
        )


@dataclass(frozen=True)
class DigitizedEvent:
    """Response and compact observables for one device and one event."""

    device_name: str
    event_index: int
    n_raw: int
    n_selected: int
    n_detected: int
    n_lost_pde: int
    n_accepted: int
    n_lost_deadtime: int
    n_active_spads: int
    n_active_timing_groups: int
    n_timing_records: int
    retained_fraction: float
    deadtime_loss_fraction: float
    mean_pde: float
    mean_wavelength_nm: float
    first_timestamp_ns: float | None
    last_timestamp_ns: float | None
    payload_bits_proxy: int
    detected_x_mm: np.ndarray
    detected_y_mm: np.ndarray
    detected_production_z_mm: np.ndarray
    detected_time_ns: np.ndarray
    accepted_x_mm: np.ndarray
    accepted_y_mm: np.ndarray
    accepted_production_z_mm: np.ndarray
    accepted_time_true_ns: np.ndarray
    accepted_time_ns: np.ndarray
    timestamp_ns: np.ndarray
    cell_multiplicities: np.ndarray

    def as_row(self) -> dict[str, Any]:
        return {
            "device": self.device_name,
            "event": self.event_index,
            "n_raw": self.n_raw,
            "n_selected": self.n_selected,
            "n_detected": self.n_detected,
            "n_lost_pde": self.n_lost_pde,
            "n_accepted": self.n_accepted,
            "n_lost_deadtime": self.n_lost_deadtime,
            "n_active_spads": self.n_active_spads,
            "n_active_timing_groups": self.n_active_timing_groups,
            "n_timing_records": self.n_timing_records,
            "retained_fraction": self.retained_fraction,
            "deadtime_loss_fraction": self.deadtime_loss_fraction,
            "mean_pde": self.mean_pde,
            "mean_wavelength_nm": self.mean_wavelength_nm,
            "first_timestamp_ns": ""
            if self.first_timestamp_ns is None
            else self.first_timestamp_ns,
            "last_timestamp_ns": ""
            if self.last_timestamp_ns is None
            else self.last_timestamp_ns,
            "payload_bits_proxy": self.payload_bits_proxy,
        }


def legacy_pde(wavelength_nm: np.ndarray) -> np.ndarray:
    """Student-study wavelength response, valid only from 300 to 620 nm."""

    wavelength_nm = np.asarray(wavelength_nm, dtype=np.float64)
    efficiency = np.polyval(LEGACY_QE_COEFFICIENTS, wavelength_nm)
    valid = (wavelength_nm >= 300.0) & (wavelength_nm <= 620.0)
    efficiency = np.where(valid, efficiency, 0.0)
    return np.clip(efficiency, 0.0, 1.0)


def pde_probability(wavelength_nm: np.ndarray, config: DeviceConfig) -> np.ndarray:
    if config.pde_model == "unity":
        probability = np.ones_like(wavelength_nm, dtype=np.float64)
    elif config.pde_model == "constant":
        probability = np.full_like(
            wavelength_nm,
            float(config.pde_value),
            dtype=np.float64,
        )
    elif config.pde_model == "legacy_polynomial":
        probability = legacy_pde(wavelength_nm)
    else:  # protected by configuration validation
        raise ValueError(f"Unsupported PDE model {config.pde_model!r}")
    return np.clip(probability * config.pde_scale, 0.0, 1.0)


def _cell_indices(
    x_mm: np.ndarray,
    y_mm: np.ndarray,
    pitch_mm: float,
    detector_min_mm: float,
    detector_max_mm: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    if detector_max_mm <= detector_min_mm:
        raise ValueError("detector_max_mm must be greater than detector_min_mm")
    if pitch_mm <= 0:
        raise ValueError("pitch_mm must be positive")
    if x_mm.shape != y_mm.shape:
        raise ValueError("x_mm and y_mm must have the same shape")
    in_sensor = (
        (x_mm >= detector_min_mm)
        & (x_mm < detector_max_mm)
        & (y_mm >= detector_min_mm)
        & (y_mm < detector_max_mm)
    )
    if not np.all(in_sensor):
        raise ValueError(
            "Photon coordinates fall outside the sensor bounds; filter them "
            "before assigning SPAD cells"
        )
    n_axis = int(np.ceil((detector_max_mm - detector_min_mm) / pitch_mm))
    ix = np.floor((x_mm - detector_min_mm) / pitch_mm).astype(np.int64)
    iy = np.floor((y_mm - detector_min_mm) / pitch_mm).astype(np.int64)
    cell_ids = iy * n_axis + ix
    return ix, iy, cell_ids, n_axis


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


def _quantize_time(time_ns: np.ndarray, tdc_bin_ps: float | None) -> np.ndarray:
    if tdc_bin_ps is None or tdc_bin_ps <= 0:
        return time_ns
    step_ns = tdc_bin_ps / 1000.0
    return np.rint(time_ns / step_ns) * step_ns


def _timestamp_records(
    accepted_time_ns: np.ndarray,
    group_ids: np.ndarray,
    config: DeviceConfig,
    gate_ns: tuple[float, float],
) -> tuple[np.ndarray, int]:
    if accepted_time_ns.size == 0:
        return np.empty(0, dtype=np.float64), 0

    policy = config.timestamp_policy
    active_groups = int(np.unique(group_ids).size)
    if policy == "all":
        records = np.sort(accepted_time_ns, kind="stable")
    elif policy == "first_global":
        records = np.asarray([np.min(accepted_time_ns)])
    elif policy == "first_last_global":
        records = (
            np.asarray([accepted_time_ns[0]])
            if accepted_time_ns.size == 1
            else np.asarray([np.min(accepted_time_ns), np.max(accepted_time_ns)])
        )
    elif policy == "first_per_group":
        unique_groups = np.unique(group_ids)
        records = np.asarray(
            [np.min(accepted_time_ns[group_ids == group]) for group in unique_groups],
            dtype=np.float64,
        )
        records.sort(kind="stable")
    elif policy == "histogram_per_group":
        gate_low, gate_high = gate_ns
        width = (gate_high - gate_low) / config.histogram_bins
        time_bin = np.floor((accepted_time_ns - gate_low) / width).astype(np.int64)
        time_bin = np.clip(time_bin, 0, config.histogram_bins - 1)
        packed = group_ids * config.histogram_bins + time_bin
        _, first_indices = np.unique(packed, return_index=True)
        records = accepted_time_ns[first_indices]
        records.sort(kind="stable")
    else:  # protected by configuration validation
        raise ValueError(f"Unsupported timestamp policy {policy!r}")
    return records, active_groups


def _payload_bits(config: DeviceConfig, n_records: int, n_active_groups: int) -> int:
    """Calculate a transparent output-size proxy, not a link-level prediction."""

    bits = config.event_overhead_bits + config.count_bits
    if config.timestamp_policy == "histogram_per_group":
        bits += n_active_groups * (
            config.address_bits + config.histogram_bins * config.count_bits
        )
    else:
        bits += n_records * (config.address_bits + config.timestamp_bits)
    return int(bits)


def digitize_event(
    event: PhotonEvent,
    config: DeviceConfig,
    *,
    seed: int,
    detector_min_mm: float = -100.0,
    detector_max_mm: float = 100.0,
    gate_ns: tuple[float, float] = (0.0, 40.0),
) -> DigitizedEvent:
    """Pass one selected photon event through one device configuration."""

    event = event.sorted_by_time()
    n_selected = event.selected_photon_count
    if n_selected == 0:
        return DigitizedEvent(
            device_name=config.name,
            event_index=event.event_index,
            n_raw=event.raw_photon_count,
            n_selected=0,
            n_detected=0,
            n_lost_pde=0,
            n_accepted=0,
            n_lost_deadtime=0,
            n_active_spads=0,
            n_active_timing_groups=0,
            n_timing_records=0,
            retained_fraction=0.0,
            deadtime_loss_fraction=0.0,
            mean_pde=0.0,
            mean_wavelength_nm=float("nan"),
            first_timestamp_ns=None,
            last_timestamp_ns=None,
            payload_bits_proxy=_payload_bits(config, 0, 0),
            detected_x_mm=np.empty(0, dtype=np.float64),
            detected_y_mm=np.empty(0, dtype=np.float64),
            detected_production_z_mm=np.empty(0, dtype=np.float64),
            detected_time_ns=np.empty(0, dtype=np.float64),
            accepted_x_mm=np.empty(0, dtype=np.float64),
            accepted_y_mm=np.empty(0, dtype=np.float64),
            accepted_production_z_mm=np.empty(0, dtype=np.float64),
            accepted_time_true_ns=np.empty(0, dtype=np.float64),
            accepted_time_ns=np.empty(0, dtype=np.float64),
            timestamp_ns=np.empty(0, dtype=np.float64),
            cell_multiplicities=np.empty(0, dtype=np.int64),
        )

    # Recreating the generator for every device is intentional.  It gives all
    # devices identical per-photon uniforms and Gaussian deviates.
    random = np.random.default_rng(np.random.SeedSequence([seed, event.event_index, 3409]))
    detection_uniform = random.random(n_selected)
    timing_standard_normal = random.standard_normal(n_selected)

    probabilities = pde_probability(event.wavelength_nm, config)
    detected_mask = detection_uniform < probabilities
    detected_indices = np.flatnonzero(detected_mask)

    detected_x = event.x_mm[detected_mask]
    detected_y = event.y_mm[detected_mask]
    detected_time = event.time_ns[detected_mask]
    detected_production_z = event.production_z_mm[detected_mask]
    detected_jitter = timing_standard_normal[detected_mask]

    pitch_mm = config.spad_pitch_um / 1000.0
    cell_x = detected_x * config.coordinate_scale
    cell_y = detected_y * config.coordinate_scale
    ix, iy, cell_ids, n_axis = _cell_indices(
        cell_x,
        cell_y,
        pitch_mm,
        detector_min_mm * config.coordinate_scale,
        detector_max_mm * config.coordinate_scale,
    )
    if cell_ids.size:
        _, cell_multiplicities = np.unique(cell_ids, return_counts=True)
    else:
        cell_multiplicities = np.empty(0, dtype=np.int64)

    accepted_mask = _apply_nonparalyzable_dead_time(
        cell_ids,
        detected_time,
        config.dead_time_ns,
    )
    accepted_time_true = detected_time[accepted_mask]
    accepted_x = detected_x[accepted_mask]
    accepted_y = detected_y[accepted_mask]
    accepted_production_z = detected_production_z[accepted_mask]
    accepted_ix = ix[accepted_mask]
    accepted_iy = iy[accepted_mask]
    accepted_jitter = detected_jitter[accepted_mask]

    if config.timing_fwhm_ps is None or config.timing_fwhm_ps <= 0:
        observed_time = accepted_time_true.astype(np.float64, copy=True)
    else:
        sigma_ns = config.timing_fwhm_ps / 2354.820045
        observed_time = accepted_time_true + sigma_ns * accepted_jitter
    observed_time = _quantize_time(observed_time, config.tdc_bin_ps)

    group_nx = int(np.ceil(n_axis / config.timing_group_spads_x))
    group_ix = accepted_ix // config.timing_group_spads_x
    group_iy = accepted_iy // config.timing_group_spads_y
    group_ids = group_iy * group_nx + group_ix
    timestamp_ns, n_active_timing_groups = _timestamp_records(
        observed_time,
        group_ids,
        config,
        gate_ns,
    )

    n_detected = int(detected_indices.size)
    n_accepted = int(accepted_time_true.size)
    n_lost_deadtime = n_detected - n_accepted
    finite_wavelength = np.isfinite(event.wavelength_nm)
    mean_wavelength = (
        float(np.mean(event.wavelength_nm[finite_wavelength]))
        if finite_wavelength.any()
        else float("nan")
    )

    return DigitizedEvent(
        device_name=config.name,
        event_index=event.event_index,
        n_raw=event.raw_photon_count,
        n_selected=n_selected,
        n_detected=n_detected,
        n_lost_pde=n_selected - n_detected,
        n_accepted=n_accepted,
        n_lost_deadtime=n_lost_deadtime,
        n_active_spads=int(np.unique(cell_ids[accepted_mask]).size),
        n_active_timing_groups=n_active_timing_groups,
        n_timing_records=int(timestamp_ns.size),
        retained_fraction=n_accepted / n_selected,
        deadtime_loss_fraction=(n_lost_deadtime / n_detected) if n_detected else 0.0,
        mean_pde=float(np.mean(probabilities)),
        mean_wavelength_nm=mean_wavelength,
        first_timestamp_ns=float(np.min(timestamp_ns)) if timestamp_ns.size else None,
        last_timestamp_ns=float(np.max(timestamp_ns)) if timestamp_ns.size else None,
        payload_bits_proxy=_payload_bits(
            config,
            int(timestamp_ns.size),
            n_active_timing_groups,
        ),
        detected_x_mm=detected_x,
        detected_y_mm=detected_y,
        detected_production_z_mm=detected_production_z,
        detected_time_ns=detected_time,
        accepted_x_mm=accepted_x,
        accepted_y_mm=accepted_y,
        accepted_production_z_mm=accepted_production_z,
        accepted_time_true_ns=accepted_time_true,
        accepted_time_ns=observed_time,
        timestamp_ns=timestamp_ns,
        cell_multiplicities=cell_multiplicities,
    )
