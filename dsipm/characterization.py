"""Standard dSiPM characterization studies for the FLASH simulation.

The historical HG-DREAM study isolates pixel occupancy at unity PDE.  This
module adds two pieces needed to interpret that result as a detector study:

* a paired PDE-by-pitch-by-recovery scan on the same incident photons;
* the exact response of a finite array to a simultaneous uniform light pulse.

The PDE scan uses common random numbers.  A photon accepted at a lower PDE is
therefore also accepted at every higher PDE in the same event.  This removes
unnecessary Monte Carlo noise from comparisons and makes the PDE trend easy
to test.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass

import numpy as np

from .digitizer import PhotonEvent, _cell_indices
from .requirements import accepted_counts_by_recovery


@dataclass(frozen=True)
class PDERequirementPoint:
    """Aggregate response and event-bootstrap intervals at one scan point."""

    pde: float
    pitch_um: float
    recovery_ns: float
    events: int
    incident_photons: int
    detected_photons: int
    accepted_photons: int
    usable_fraction_low68: float
    usable_fraction_high68: float
    recovery_retention_low68: float
    recovery_retention_high68: float

    @property
    def usable_fraction(self) -> float:
        """Accepted photons divided by photons incident on the sensor."""

        if self.incident_photons == 0:
            return 0.0
        return self.accepted_photons / self.incident_photons

    @property
    def recovery_retention(self) -> float:
        """Accepted avalanches divided by PDE-selected photons."""

        if self.detected_photons == 0:
            return 0.0
        return self.accepted_photons / self.detected_photons


def expected_unique_cells(
    incident_photons: np.ndarray | Sequence[float] | float,
    n_cells: int,
) -> np.ndarray:
    """Return the exact mean fired-cell count for a uniform prompt pulse.

    Every photon independently chooses one of ``n_cells`` pixels and every
    pixel records at most one hit.  The result is the standard finite-pixel
    SiPM saturation curve, with no recovery during the pulse.
    """

    if n_cells <= 0:
        raise ValueError("n_cells must be positive")
    photons = np.asarray(incident_photons, dtype=np.float64)
    if np.any(photons < 0):
        raise ValueError("incident_photons cannot be negative")
    if n_cells == 1:
        return np.where(photons > 0, 1.0, 0.0)
    return -float(n_cells) * np.expm1(
        photons * np.log1p(-1.0 / float(n_cells))
    )


def _ratio_interval(
    numerator: np.ndarray,
    denominator: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    ratio = np.divide(
        numerator,
        denominator,
        out=np.zeros_like(numerator, dtype=np.float64),
        where=denominator > 0,
    )
    return (
        np.quantile(ratio, 0.16, axis=0),
        np.quantile(ratio, 0.84, axis=0),
    )


def scan_pde_pitch_recovery(
    events: Iterable[PhotonEvent],
    *,
    pdes: Sequence[float],
    pitches_um: Sequence[float],
    recovery_times_ns: Sequence[float],
    seed: int = 4307,
    coordinate_scale: float = 3.0,
    detector_min_mm: float = -100.0,
    detector_max_mm: float = 100.0,
    bootstrap_samples: int = 400,
    bootstrap_seed: int = 20261010,
) -> list[PDERequirementPoint]:
    """Scan usable signal and recovery survival on selected photon events.

    Confidence intervals resample whole GEANT4 events, rather than individual
    photons, so the quoted spread retains the event-to-event shower
    fluctuations.  The 16th and 84th percentiles are reported.
    """

    pde_values = np.asarray(pdes, dtype=np.float64)
    pitches = np.asarray(pitches_um, dtype=np.float64)
    recovery = np.asarray(recovery_times_ns, dtype=np.float64)
    if pde_values.size == 0 or np.any((pde_values < 0) | (pde_values > 1)):
        raise ValueError("pdes must be a non-empty sequence between zero and one")
    if pitches.size == 0 or np.any(pitches <= 0):
        raise ValueError("pitches must be a non-empty positive sequence")
    if recovery.size == 0 or np.any(recovery < 0):
        raise ValueError("recovery times must be a non-empty non-negative sequence")
    if coordinate_scale <= 0:
        raise ValueError("coordinate_scale must be positive")
    if bootstrap_samples < 0:
        raise ValueError("bootstrap_samples cannot be negative")

    incident_by_event: list[int] = []
    detected_by_event: list[np.ndarray] = []
    accepted_by_event: list[np.ndarray] = []

    for event in events:
        event = event.sorted_by_time()
        n_incident = event.selected_photon_count
        incident_by_event.append(n_incident)
        event_detected = np.zeros(pde_values.size, dtype=np.int64)
        event_accepted = np.zeros(
            (pde_values.size, pitches.size, recovery.size),
            dtype=np.int64,
        )
        if n_incident:
            rng = np.random.default_rng(
                np.random.SeedSequence([seed, event.event_index, 3409])
            )
            detection_uniform = rng.random(n_incident)
            scaled_x = event.x_mm * coordinate_scale
            scaled_y = event.y_mm * coordinate_scale

            for pde_index, pde in enumerate(pde_values):
                detected_mask = detection_uniform < pde
                event_detected[pde_index] = int(np.count_nonzero(detected_mask))
                if not np.any(detected_mask):
                    continue
                detected_x = scaled_x[detected_mask]
                detected_y = scaled_y[detected_mask]
                detected_time = event.time_ns[detected_mask]
                for pitch_index, pitch_um in enumerate(pitches):
                    _, _, cell_ids, _ = _cell_indices(
                        detected_x,
                        detected_y,
                        pitch_um / 1000.0,
                        detector_min_mm * coordinate_scale,
                        detector_max_mm * coordinate_scale,
                    )
                    event_accepted[pde_index, pitch_index] = (
                        accepted_counts_by_recovery(
                            cell_ids,
                            detected_time,
                            recovery,
                        )
                    )
        detected_by_event.append(event_detected)
        accepted_by_event.append(event_accepted)

    if not incident_by_event:
        return []

    incident_array = np.asarray(incident_by_event, dtype=np.int64)
    detected_array = np.stack(detected_by_event)
    accepted_array = np.stack(accepted_by_event)
    incident_total = int(incident_array.sum())
    detected_total = detected_array.sum(axis=0)
    accepted_total = accepted_array.sum(axis=0)

    usable_low = np.divide(
        accepted_total,
        incident_total,
        out=np.zeros_like(accepted_total, dtype=np.float64),
        where=incident_total > 0,
    )
    usable_high = usable_low.copy()
    detected_denominator = detected_total[:, None, None]
    retention_low = np.divide(
        accepted_total,
        detected_denominator,
        out=np.zeros_like(accepted_total, dtype=np.float64),
        where=detected_denominator > 0,
    )
    retention_high = retention_low.copy()

    if bootstrap_samples and incident_array.size > 1:
        bootstrap_rng = np.random.default_rng(bootstrap_seed)
        weights = bootstrap_rng.multinomial(
            incident_array.size,
            np.full(incident_array.size, 1.0 / incident_array.size),
            size=bootstrap_samples,
        ).astype(np.float64)
        bootstrap_incident = weights @ incident_array.astype(np.float64)
        bootstrap_detected = weights @ detected_array.astype(np.float64)
        bootstrap_accepted = (
            weights @ accepted_array.reshape(incident_array.size, -1).astype(np.float64)
        ).reshape(
            bootstrap_samples,
            pde_values.size,
            pitches.size,
            recovery.size,
        )
        usable_low, usable_high = _ratio_interval(
            bootstrap_accepted,
            bootstrap_incident[:, None, None, None],
        )
        retention_low, retention_high = _ratio_interval(
            bootstrap_accepted,
            bootstrap_detected[:, :, None, None],
        )

    return [
        PDERequirementPoint(
            pde=float(pde),
            pitch_um=float(pitch),
            recovery_ns=float(holdoff),
            events=int(incident_array.size),
            incident_photons=incident_total,
            detected_photons=int(detected_total[pde_index]),
            accepted_photons=int(
                accepted_total[pde_index, pitch_index, recovery_index]
            ),
            usable_fraction_low68=float(
                usable_low[pde_index, pitch_index, recovery_index]
            ),
            usable_fraction_high68=float(
                usable_high[pde_index, pitch_index, recovery_index]
            ),
            recovery_retention_low68=float(
                retention_low[pde_index, pitch_index, recovery_index]
            ),
            recovery_retention_high68=float(
                retention_high[pde_index, pitch_index, recovery_index]
            ),
        )
        for pde_index, pde in enumerate(pde_values)
        for pitch_index, pitch in enumerate(pitches)
        for recovery_index, holdoff in enumerate(recovery)
    ]
