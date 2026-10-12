"""FLASH design-requirement scans on the historical HG-DREAM samples.

The first study separates physical SPAD pitch from recovery time.  Chris
Madrid's original occupancy plots used one accepted hit per SPAD for the
entire 40 ns analysis gate.  That remains a useful worst case, but it is not
the behavior of devices with nanosecond-scale recharge or hold-off.

This module keeps the historical photon selection and spatial mapping intact
while scanning the two parameters independently.  PDE is deliberately unity
so the result answers an occupancy/recovery question, not a full sensor-yield
question.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass

import numpy as np

from .digitizer import PhotonEvent, _cell_indices


@dataclass(frozen=True)
class RequirementPoint:
    """Aggregate photon retention at one pitch and recovery-time point."""

    pitch_um: float
    recovery_ns: float
    events: int
    incident_photons: int
    accepted_photons: int

    @property
    def retained_fraction(self) -> float:
        if self.incident_photons == 0:
            return 0.0
        return self.accepted_photons / self.incident_photons


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


def scan_pitch_and_recovery(
    events: Iterable[PhotonEvent],
    *,
    pitches_um: Sequence[float],
    recovery_times_ns: Sequence[float],
    coordinate_scale: float = 3.0,
    detector_min_mm: float = -100.0,
    detector_max_mm: float = 100.0,
) -> list[RequirementPoint]:
    """Aggregate occupancy retention for a pitch-by-recovery grid."""

    pitches = np.asarray(pitches_um, dtype=np.float64)
    recovery = np.asarray(recovery_times_ns, dtype=np.float64)
    if np.any(pitches <= 0):
        raise ValueError("pitches must be positive")
    if coordinate_scale <= 0:
        raise ValueError("coordinate_scale must be positive")

    accepted = np.zeros((pitches.size, recovery.size), dtype=np.int64)
    incident = 0
    n_events = 0
    for event in events:
        event = event.sorted_by_time()
        incident += event.selected_photon_count
        n_events += 1
        scaled_x = event.x_mm * coordinate_scale
        scaled_y = event.y_mm * coordinate_scale
        for pitch_index, pitch_um in enumerate(pitches):
            _, _, cell_ids, _ = _cell_indices(
                scaled_x,
                scaled_y,
                pitch_um / 1000.0,
                detector_min_mm * coordinate_scale,
                detector_max_mm * coordinate_scale,
            )
            accepted[pitch_index] += accepted_counts_by_recovery(
                cell_ids,
                event.time_ns,
                recovery,
            )

    return [
        RequirementPoint(
            pitch_um=float(pitch),
            recovery_ns=float(holdoff),
            events=n_events,
            incident_photons=incident,
            accepted_photons=int(accepted[pitch_index, recovery_index]),
        )
        for pitch_index, pitch in enumerate(pitches)
        for recovery_index, holdoff in enumerate(recovery)
    ]


def interpolated_pitch_limit(
    points: Sequence[RequirementPoint],
    *,
    recovery_ns: float,
    target_retention: float,
) -> float | None:
    """Estimate the largest pitch meeting a target using linear interpolation."""

    selected = sorted(
        (
            point
            for point in points
            if np.isclose(point.recovery_ns, recovery_ns, rtol=0.0, atol=1.0e-12)
        ),
        key=lambda point: point.pitch_um,
    )
    if not selected:
        raise ValueError(f"no points found for recovery_ns={recovery_ns:g}")
    retention = np.asarray([point.retained_fraction for point in selected])
    pitch = np.asarray([point.pitch_um for point in selected])
    passing = np.flatnonzero(retention >= target_retention)
    if passing.size == 0:
        return None
    last = int(passing[-1])
    if last == len(selected) - 1:
        return float(pitch[last])
    x0, x1 = pitch[last], pitch[last + 1]
    y0, y1 = retention[last], retention[last + 1]
    if np.isclose(y0, y1):
        return float(x0)
    fraction = (target_retention - y0) / (y1 - y0)
    return float(x0 + fraction * (x1 - x0))
