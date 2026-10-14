"""Mergeable DREAMSim reductions for HPCC production files.

The raw optical-photon vectors stay on HPCC.  This module converts each event
to the response quantities needed by the CPAD/JINST plots and accumulates
histograms whose bin contents can be summed exactly across independent jobs.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np

from .full_detector import (
    AnalogSiPMScenario,
    FullDetectorPhotonEvent,
    FullDetectorScenario,
    digitize_analog_sipm_event,
    digitize_full_detector_event,
    global_cell_ids,
    scan_full_detector_event,
)


ANALOG_SIPM = AnalogSiPMScenario(
    name="onsemi_microfc_30035",
    pde=0.31,
    microcell_size_um=35.0,
    n_microcells=4774,
    recovery_tau_ns=82.0,
    fast_output_fwhm_ns=1.5,
    pde_wavelength_nm=420.0,
    overvoltage_v=2.5,
    source_url="https://www.onsemi.com/pdf/datasheet/microc-series-d.pdf",
)

TIME_EDGES_NS = np.linspace(5.0, 40.0, 701, dtype=np.float64)
LOCAL_EDGES_MM = np.linspace(-1.5, 1.5, 121, dtype=np.float64)
GLOBAL_EDGES_MM = np.linspace(-100.0, 100.0, 201, dtype=np.float64)
DEPTH_EDGES_MM = np.linspace(0.0, 2000.0, 201, dtype=np.float64)
DEPTH_TIME_EDGES_NS = np.linspace(5.0, 40.0, 351, dtype=np.float64)
FIBER_OCCUPANCY_MAX = 8192
CELL_OCCUPANCY_MAX = 64
FIBER_OCCUPANCY_EDGES = np.concatenate(
    (
        np.arange(0.5, FIBER_OCCUPANCY_MAX + 1.5, dtype=np.float64),
        np.asarray([np.inf]),
    )
)
CELL_OCCUPANCY_EDGES = np.concatenate(
    (
        np.arange(0.5, CELL_OCCUPANCY_MAX + 1.5, dtype=np.float64),
        np.asarray([np.inf]),
    )
)


def scenario_column(pitch_um: float, recovery_ns: float) -> str:
    return f"accepted_pitch{pitch_um:g}um_recovery{recovery_ns:g}ns"


def _quantile(values: np.ndarray, probability: float) -> float | str:
    return float(np.quantile(values, probability)) if values.size else ""


def _mean(values: np.ndarray) -> float | str:
    return float(np.mean(values)) if values.size else ""


def _rms(values: np.ndarray) -> float | str:
    return float(np.std(values)) if values.size else ""


def _occupancies(ids: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return sorted active ids and their positive integer occupancies."""

    values = np.asarray(ids, dtype=np.int64)
    if values.size == 0:
        return np.asarray([], dtype=np.int64), np.asarray([], dtype=np.int64)
    return np.unique(values, return_counts=True)


def _fraction_in_top(occupancies: np.ndarray, count: int) -> float:
    if occupancies.size == 0:
        return 0.0
    return float(np.sum(np.sort(occupancies)[-count:]) / np.sum(occupancies))


def reduce_event(
    event: FullDetectorPhotonEvent,
    *,
    sample: str,
    source_root: str,
    source_entry: int,
    pde: float,
    pitches_um: Sequence[float],
    recovery_times_ns: Sequence[float],
    reference_pitch_um: float,
    reference_recovery_ns: float,
    seed: int,
    analog_sipm: AnalogSiPMScenario = ANALOG_SIPM,
) -> tuple[dict[str, object], np.ndarray, np.ndarray, np.ndarray]:
    """Return one compact row plus reference detected/accepted masks."""

    pitches = np.asarray(pitches_um, dtype=np.float64)
    recovery = np.asarray(recovery_times_ns, dtype=np.float64)
    detected, accepted = scan_full_detector_event(
        event,
        pde=pde,
        pitches_um=pitches,
        recovery_times_ns=recovery,
        seed=seed,
    )
    reference = digitize_full_detector_event(
        event,
        FullDetectorScenario(
            pde=pde,
            pitch_um=reference_pitch_um,
            recovery_ns=reference_recovery_ns,
        ),
        seed=seed,
    )
    pitch_matches = np.flatnonzero(np.isclose(pitches, reference_pitch_um))
    recovery_matches = np.flatnonzero(
        np.isclose(recovery, reference_recovery_ns)
    )
    if pitch_matches.size != 1 or recovery_matches.size != 1:
        raise ValueError("Reference pitch and recovery must appear once in the scan")
    reference_scan_count = int(
        accepted[int(pitch_matches[0]), int(recovery_matches[0])]
    )
    if reference.n_detected != detected or reference.n_accepted != reference_scan_count:
        raise AssertionError("Reference digitization does not close to the scan result")

    analog = digitize_analog_sipm_event(event, analog_sipm, seed=seed)
    accepted_times = event.time_ns[reference.accepted_mask]
    detected_times = event.time_ns[reference.detected_mask]
    lost_mask = reference.detected_mask & ~reference.accepted_mask
    _, fiber_occupancies = _occupancies(event.sensor_id)
    if fiber_occupancies.size != event.active_sensor_count:
        raise AssertionError("Active-fiber count does not close to sensor ids")
    if int(np.sum(fiber_occupancies)) != event.selected_photon_count:
        raise AssertionError("Per-fiber occupancy does not close to selected photons")

    row: dict[str, object] = {
        "sample": sample,
        "source_root": source_root,
        "source_entry": source_entry,
        "event": event.event_index,
        "beam_energy_gev": event.beam_energy_gev,
        "sum3dCC": event.c_signal,
        "sum3dSS": event.s_signal,
        "ecalo_truth_gev": event.ecalo_truth_gev,
        "eleak_truth_gev": event.eleak_truth_gev,
        "raw_optical_photons": event.raw_photon_count,
        "selected_cherenkov_photons": event.selected_photon_count,
        "active_cherenkov_fibers": event.active_sensor_count,
        "max_photons_in_one_fiber": event.max_sensor_photons,
        "active_fiber_occupancy_q50": _quantile(fiber_occupancies, 0.50),
        "active_fiber_occupancy_q90": _quantile(fiber_occupancies, 0.90),
        "active_fiber_occupancy_q99": _quantile(fiber_occupancies, 0.99),
        "top_fiber_photon_fraction": _fraction_in_top(fiber_occupancies, 1),
        "top5_fiber_photon_fraction": _fraction_in_top(fiber_occupancies, 5),
        "active_fibers_ge2_photons": int(np.count_nonzero(fiber_occupancies >= 2)),
        "active_fibers_ge5_photons": int(np.count_nonzero(fiber_occupancies >= 5)),
        "active_fibers_ge10_photons": int(np.count_nonzero(fiber_occupancies >= 10)),
        "active_fibers_ge20_photons": int(np.count_nonzero(fiber_occupancies >= 20)),
        "active_fibers_ge50_photons": int(np.count_nonzero(fiber_occupancies >= 50)),
        "active_fibers_ge100_photons": int(np.count_nonzero(fiber_occupancies >= 100)),
        "detected_reference": reference.n_detected,
        "accepted_reference": reference.n_accepted,
        "lost_reference_recovery": int(np.count_nonzero(lost_mask)),
        "reference_recovery_survival": (
            reference.n_accepted / reference.n_detected
            if reference.n_detected
            else 1.0
        ),
        "analog_detected_at_31pct_pde": analog.n_detected,
        "analog_charge_pe": analog.charge_pe,
        "analog_charge_survival": analog.charge_survival,
        "arrival_q10_ns": _quantile(event.time_ns, 0.10),
        "arrival_median_ns": _quantile(event.time_ns, 0.50),
        "arrival_q90_ns": _quantile(event.time_ns, 0.90),
        "accepted_arrival_q10_ns": _quantile(accepted_times, 0.10),
        "accepted_arrival_median_ns": _quantile(accepted_times, 0.50),
        "accepted_arrival_q90_ns": _quantile(accepted_times, 0.90),
        "detected_arrival_mean_ns": _mean(detected_times),
        "sensor_x_mean_mm": _mean(event.sensor_x_mm),
        "sensor_y_mean_mm": _mean(event.sensor_y_mm),
        "sensor_x_rms_mm": _rms(event.sensor_x_mm),
        "sensor_y_rms_mm": _rms(event.sensor_y_mm),
        "production_depth_mean_mm": _mean(event.production_z_mm),
        "production_depth_rms_mm": _rms(event.production_z_mm),
    }
    for pitch_index, pitch_um in enumerate(pitches):
        for recovery_index, recovery_ns in enumerate(recovery):
            row[scenario_column(float(pitch_um), float(recovery_ns))] = int(
                accepted[pitch_index, recovery_index]
            )
    return row, reference.detected_mask, reference.accepted_mask, lost_mask


@dataclass
class HistogramAccumulator:
    """Additive histograms saved once per ROOT file and merged later."""

    arrays: dict[str, np.ndarray] = field(default_factory=dict)
    pitches_um: tuple[float, ...] = ()

    def __post_init__(self) -> None:
        if self.arrays:
            return
        self.arrays = {
            "time_incident": np.zeros(TIME_EDGES_NS.size - 1, dtype=np.int64),
            "time_detected": np.zeros(TIME_EDGES_NS.size - 1, dtype=np.int64),
            "time_accepted": np.zeros(TIME_EDGES_NS.size - 1, dtype=np.int64),
            "time_lost_recovery": np.zeros(TIME_EDGES_NS.size - 1, dtype=np.int64),
            "local_xy_incident": np.zeros(
                (LOCAL_EDGES_MM.size - 1, LOCAL_EDGES_MM.size - 1),
                dtype=np.int64,
            ),
            "local_xy_detected": np.zeros(
                (LOCAL_EDGES_MM.size - 1, LOCAL_EDGES_MM.size - 1),
                dtype=np.int64,
            ),
            "local_xy_accepted": np.zeros(
                (LOCAL_EDGES_MM.size - 1, LOCAL_EDGES_MM.size - 1),
                dtype=np.int64,
            ),
            "local_xy_lost_recovery": np.zeros(
                (LOCAL_EDGES_MM.size - 1, LOCAL_EDGES_MM.size - 1),
                dtype=np.int64,
            ),
            "global_xy_incident": np.zeros(
                (GLOBAL_EDGES_MM.size - 1, GLOBAL_EDGES_MM.size - 1),
                dtype=np.int64,
            ),
            "global_xy_accepted": np.zeros(
                (GLOBAL_EDGES_MM.size - 1, GLOBAL_EDGES_MM.size - 1),
                dtype=np.int64,
            ),
            "depth_time_incident": np.zeros(
                (DEPTH_EDGES_MM.size - 1, DEPTH_TIME_EDGES_NS.size - 1),
                dtype=np.int64,
            ),
            "fiber_active_count_by_incident_occupancy": np.zeros(
                FIBER_OCCUPANCY_EDGES.size - 1,
                dtype=np.int64,
            ),
            "fiber_incident_photons_by_incident_occupancy": np.zeros(
                FIBER_OCCUPANCY_EDGES.size - 1,
                dtype=np.int64,
            ),
            "fiber_detected_photons_by_incident_occupancy": np.zeros(
                FIBER_OCCUPANCY_EDGES.size - 1,
                dtype=np.int64,
            ),
            "fiber_accepted_photons_by_incident_occupancy": np.zeros(
                FIBER_OCCUPANCY_EDGES.size - 1,
                dtype=np.int64,
            ),
            "fiber_lost_photons_by_incident_occupancy": np.zeros(
                FIBER_OCCUPANCY_EDGES.size - 1,
                dtype=np.int64,
            ),
        }
        if self.pitches_um:
            shape = (len(self.pitches_um), CELL_OCCUPANCY_EDGES.size - 1)
            self.arrays["cell_active_count_by_detected_occupancy"] = np.zeros(
                shape, dtype=np.int64
            )
            self.arrays["cell_detected_photons_by_detected_occupancy"] = np.zeros(
                shape, dtype=np.int64
            )
            self.arrays["cell_extra_photons_by_detected_occupancy"] = np.zeros(
                shape, dtype=np.int64
            )

    def add(
        self,
        event: FullDetectorPhotonEvent,
        detected_mask: np.ndarray,
        accepted_mask: np.ndarray,
        lost_mask: np.ndarray,
    ) -> None:
        masks = {
            "incident": np.ones(event.time_ns.size, dtype=bool),
            "detected": detected_mask,
            "accepted": accepted_mask,
            "lost_recovery": lost_mask,
        }
        for name, mask in masks.items():
            self.arrays[f"time_{name}"] += np.histogram(
                event.time_ns[mask], bins=TIME_EDGES_NS
            )[0]
            self.arrays[f"local_xy_{name}"] += np.histogram2d(
                event.local_x_mm[mask],
                event.local_y_mm[mask],
                bins=(LOCAL_EDGES_MM, LOCAL_EDGES_MM),
            )[0].astype(np.int64)
        self.arrays["global_xy_incident"] += np.histogram2d(
            event.sensor_x_mm,
            event.sensor_y_mm,
            bins=(GLOBAL_EDGES_MM, GLOBAL_EDGES_MM),
        )[0].astype(np.int64)
        self.arrays["global_xy_accepted"] += np.histogram2d(
            event.sensor_x_mm[accepted_mask],
            event.sensor_y_mm[accepted_mask],
            bins=(GLOBAL_EDGES_MM, GLOBAL_EDGES_MM),
        )[0].astype(np.int64)
        self.arrays["depth_time_incident"] += np.histogram2d(
            event.production_z_mm,
            event.time_ns,
            bins=(DEPTH_EDGES_MM, DEPTH_TIME_EDGES_NS),
        )[0].astype(np.int64)

        sensor_ids, incident_occupancy = _occupancies(event.sensor_id)
        detected_occupancy = np.zeros(sensor_ids.size, dtype=np.int64)
        accepted_occupancy = np.zeros(sensor_ids.size, dtype=np.int64)
        lost_occupancy = np.zeros(sensor_ids.size, dtype=np.int64)
        for target, mask in (
            (detected_occupancy, detected_mask),
            (accepted_occupancy, accepted_mask),
            (lost_occupancy, lost_mask),
        ):
            selected_ids, selected_counts = _occupancies(event.sensor_id[mask])
            if selected_ids.size:
                positions = np.searchsorted(sensor_ids, selected_ids)
                if not np.array_equal(sensor_ids[positions], selected_ids):
                    raise AssertionError("Digitized fiber id is absent from incident photons")
                target[positions] = selected_counts
        fiber_histograms = {
            "fiber_active_count_by_incident_occupancy": None,
            "fiber_incident_photons_by_incident_occupancy": incident_occupancy,
            "fiber_detected_photons_by_incident_occupancy": detected_occupancy,
            "fiber_accepted_photons_by_incident_occupancy": accepted_occupancy,
            "fiber_lost_photons_by_incident_occupancy": lost_occupancy,
        }
        for name, weights in fiber_histograms.items():
            self.arrays[name] += np.histogram(
                incident_occupancy,
                bins=FIBER_OCCUPANCY_EDGES,
                weights=weights,
            )[0].astype(np.int64)

        for pitch_index, pitch_um in enumerate(self.pitches_um):
            _, cell_occupancy = _occupancies(
                global_cell_ids(event, float(pitch_um), mask=detected_mask)
            )
            cell_histograms = {
                "cell_active_count_by_detected_occupancy": None,
                "cell_detected_photons_by_detected_occupancy": cell_occupancy,
                "cell_extra_photons_by_detected_occupancy": np.maximum(
                    cell_occupancy - 1, 0
                ),
            }
            for name, weights in cell_histograms.items():
                self.arrays[name][pitch_index] += np.histogram(
                    cell_occupancy,
                    bins=CELL_OCCUPANCY_EDGES,
                    weights=weights,
                )[0].astype(np.int64)

    def payload(self) -> dict[str, np.ndarray]:
        return {
            **self.arrays,
            "time_edges_ns": TIME_EDGES_NS,
            "local_x_edges_mm": LOCAL_EDGES_MM,
            "local_y_edges_mm": LOCAL_EDGES_MM,
            "global_x_edges_mm": GLOBAL_EDGES_MM,
            "global_y_edges_mm": GLOBAL_EDGES_MM,
            "depth_edges_mm": DEPTH_EDGES_MM,
            "depth_time_edges_ns": DEPTH_TIME_EDGES_NS,
            "fiber_occupancy_edges": FIBER_OCCUPANCY_EDGES,
            "cell_occupancy_edges": CELL_OCCUPANCY_EDGES,
            "cell_occupancy_pitches_um": np.asarray(
                self.pitches_um, dtype=np.float64
            ),
        }

    def closure(self) -> dict[str, int]:
        incident = int(np.sum(self.arrays["time_incident"]))
        detected = int(np.sum(self.arrays["time_detected"]))
        accepted = int(np.sum(self.arrays["time_accepted"]))
        lost = int(np.sum(self.arrays["time_lost_recovery"]))
        if detected != accepted + lost:
            raise AssertionError("Detected time histogram does not close")
        fiber_incident = int(
            np.sum(self.arrays["fiber_incident_photons_by_incident_occupancy"])
        )
        fiber_detected = int(
            np.sum(self.arrays["fiber_detected_photons_by_incident_occupancy"])
        )
        fiber_accepted = int(
            np.sum(self.arrays["fiber_accepted_photons_by_incident_occupancy"])
        )
        fiber_lost = int(
            np.sum(self.arrays["fiber_lost_photons_by_incident_occupancy"])
        )
        if fiber_incident != incident:
            raise AssertionError("Fiber occupancy does not close to incident photons")
        if fiber_detected != detected:
            raise AssertionError("Fiber occupancy does not close to detected photons")
        if fiber_accepted != accepted or fiber_lost != lost:
            raise AssertionError("Fiber occupancy does not close to recovery response")
        if not np.array_equal(
            self.arrays["fiber_detected_photons_by_incident_occupancy"],
            self.arrays["fiber_accepted_photons_by_incident_occupancy"]
            + self.arrays["fiber_lost_photons_by_incident_occupancy"],
        ):
            raise AssertionError("Per-occupancy detected photons do not close")
        if self.pitches_um:
            cell_photons = np.sum(
                self.arrays["cell_detected_photons_by_detected_occupancy"],
                axis=1,
            )
            if np.any(cell_photons != detected):
                raise AssertionError("SPAD-cell occupancy does not close to detected photons")
        return {
            "incident_photons_in_time_histogram": incident,
            "detected_photons_in_time_histogram": detected,
            "accepted_photons_in_time_histogram": accepted,
            "lost_recovery_photons_in_time_histogram": lost,
            "active_fibers_in_occupancy_histogram": int(
                np.sum(self.arrays["fiber_active_count_by_incident_occupancy"])
            ),
            "incident_photons_in_fiber_occupancy_histogram": fiber_incident,
            "detected_photons_in_fiber_occupancy_histogram": fiber_detected,
            "accepted_photons_in_fiber_occupancy_histogram": fiber_accepted,
            "lost_recovery_photons_in_fiber_occupancy_histogram": fiber_lost,
        }
