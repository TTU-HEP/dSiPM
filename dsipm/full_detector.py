"""Full-detector DREAMSim photon extraction and dSiPM response utilities.

The older studies used optical photons from one reference rod.  The
``FullDetector`` samples retain the sensor-plane coordinates for all
rods while reusing the four local Cherenkov-fiber indices in each rod.  This
module unwraps that repeated geometry so dead time is applied independently
to every physical SPAD in every physical fiber.

Only Cherenkov optical photons are explicit in the available files.  The
scintillation observable stored in ``sum3dSS`` is retained as the companion
dual-readout signal, but it is not passed through a photon-level dSiPM model.
"""

from __future__ import annotations

from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from .digitizer import HC_GEV_NM, _apply_nonparalyzable_dead_time
from .requirements import accepted_counts_by_recovery
from .root_io import OLD_2025_X_SHIFT, OLD_2025_Y_SHIFT


ROD_PITCH_CM = 0.4
SENSOR_HALF_WIDTH_MM = 1.5
SENSOR_INDEX_OFFSET = 2048
ORIGINAL_DREAM_FIBER_DIAMETER_MM = 0.8

# Current HG-DREAM Cherenkov-channel placement used for the 2025 CERN beam
# test.  Coordinates follow the board map in CaloXDataAnalysis: the outer
# readout has 12 x 16 mm2 calorimeter channels coupled to 6 x 6 mm2 SiPMs,
# while the core has 12 x 4 mm2 channels coupled to 3 x 3 mm2 SiPMs.
HGDREAM_OUTER_BOARD_ORIGINS = (
    (-13.5, 3.5),
    (-9.5, 7.5),
    (-5.5, 7.5),
    (-1.5, 9.5),
    (2.5, 7.5),
    (6.5, 7.5),
    (10.5, 3.5),
    (-9.5, -0.5),
    (-5.5, -0.5),
    (-1.5, -2.5),
    (2.5, -0.5),
    (6.5, -0.5),
)
HGDREAM_CORE_BOARD_ORIGINS = ((-1.5, 1.875), (0.5, 1.875))

# The original DREAM towers were approximately regular hexagons, contained
# 270 rods, and measured about 80 mm from apex to apex.  The 37.1 mm value
# quoted in DREAM publications is the equal-area effective radius, not the
# center-to-center pitch.  A regular 80 mm apex-to-apex hexagon has a 40 mm
# circumradius and a sqrt(3) * 40 mm pitch between adjacent tower centers.
ORIGINAL_DREAM_TOWER_APEX_TO_APEX_MM = 80.0
ORIGINAL_DREAM_TOWER_RADIUS_MM = 0.5 * ORIGINAL_DREAM_TOWER_APEX_TO_APEX_MM
ORIGINAL_DREAM_TOWER_SPACING_MM = (
    np.sqrt(3.0) * ORIGINAL_DREAM_TOWER_RADIUS_MM
)

REQUIRED_BRANCHES = {
    "beamE",
    "sum3dSS",
    "sum3dCC",
    "eCalotruth",
    "eLeaktruth",
    "OP_time_final",
    "OP_productionFiber",
    "OP_isCoreC",
    "OP_pos_final_x",
    "OP_pos_final_y",
    "OP_pos_final_z",
    "OP_pos_produced_z",
    "OP_mom_produced_x",
    "OP_mom_produced_y",
    "OP_mom_produced_z",
}


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


@dataclass(frozen=True)
class TransverseChannelLayout:
    """Rectangular calorimeter channels and their coupled sensor sizes."""

    center_x_mm: np.ndarray
    center_y_mm: np.ndarray
    width_x_mm: np.ndarray
    width_y_mm: np.ndarray
    sensor_width_mm: np.ndarray
    rods_per_channel: np.ndarray
    is_core: np.ndarray

    @property
    def n_channels(self) -> int:
        return int(self.center_x_mm.size)


def hgdream_cherenkov_channel_layout() -> TransverseChannelLayout:
    """Return the 448-channel Cherenkov-side HG-DREAM segmentation.

    The layout contains 384 outer channels that each collect 12 rods onto a
    6 x 6 mm2 SiPM and 64 core channels that each collect three rods onto a
    3 x 3 mm2 SiPM.  Board origins are copied from the current test-beam
    channel map; the positions are converted to millimetres at the absorber
    face.
    """

    rows: list[tuple[float, float, float, float, float, int, bool]] = []
    for origin_x, origin_y in HGDREAM_OUTER_BOARD_ORIGINS:
        for local_x in range(4):
            for local_y in range(8):
                rows.append(
                    (
                        12.0 * (origin_x + local_x),
                        16.0 * (origin_y - local_y),
                        12.0,
                        16.0,
                        6.0,
                        12,
                        False,
                    )
                )
    for origin_x, origin_y in HGDREAM_CORE_BOARD_ORIGINS:
        for local_x in range(2):
            for local_y in range(16):
                rows.append(
                    (
                        12.0 * (origin_x + local_x),
                        16.0 * (origin_y - local_y / 4.0),
                        12.0,
                        4.0,
                        3.0,
                        3,
                        True,
                    )
                )
    values = np.asarray(rows, dtype=np.float64)
    return TransverseChannelLayout(
        center_x_mm=values[:, 0],
        center_y_mm=values[:, 1],
        width_x_mm=values[:, 2],
        width_y_mm=values[:, 3],
        sensor_width_mm=values[:, 4],
        rods_per_channel=values[:, 5].astype(np.int64),
        is_core=values[:, 6].astype(bool),
    )


def hgdream_cherenkov_channel_ids(
    event: FullDetectorPhotonEvent,
    layout: TransverseChannelLayout | None = None,
) -> np.ndarray:
    """Map each selected photon to the current HG-DREAM C channel."""

    layout = hgdream_cherenkov_channel_layout() if layout is None else layout
    photon_x_mm = 10.0 * ROD_PITCH_CM * event.rod_ix
    photon_y_mm = 10.0 * ROD_PITCH_CM * event.rod_iy
    channel_ids = np.full(event.time_ns.size, -1, dtype=np.int64)
    for channel_id in range(layout.n_channels):
        half_x = 0.5 * layout.width_x_mm[channel_id]
        half_y = 0.5 * layout.width_y_mm[channel_id]
        inside = (
            (photon_x_mm >= layout.center_x_mm[channel_id] - half_x)
            & (photon_x_mm < layout.center_x_mm[channel_id] + half_x)
            & (photon_y_mm >= layout.center_y_mm[channel_id] - half_y)
            & (photon_y_mm < layout.center_y_mm[channel_id] + half_y)
        )
        if np.any(inside & (channel_ids >= 0)):
            raise ValueError("HG-DREAM channel rectangles overlap")
        channel_ids[inside] = channel_id
    return channel_ids


def original_dream_pmt_centers(
    *, spacing_mm: float = ORIGINAL_DREAM_TOWER_SPACING_MM,
) -> np.ndarray:
    """Return the 19 tower centers of the original DREAM PMT readout."""

    if spacing_mm <= 0:
        raise ValueError("PMT tower spacing must be positive")
    axial_coordinates = [
        (q, r)
        for q in range(-2, 3)
        for r in range(-2, 3)
        if max(abs(q), abs(r), abs(q + r)) <= 2
    ]
    return np.asarray(
        [
            (
                spacing_mm * (q + 0.5 * r),
                spacing_mm * np.sqrt(3.0) * 0.5 * r,
            )
            for q, r in axial_coordinates
        ],
        dtype=np.float64,
    )


def original_dream_pmt_channel_ids(
    event: FullDetectorPhotonEvent,
    centers_mm: np.ndarray | None = None,
    *,
    tower_radius_mm: float = ORIGINAL_DREAM_TOWER_RADIUS_MM,
) -> np.ndarray:
    """Project the full-detector photons onto the original 19 C PMTs.

    The available event was generated with the HG-DREAM geometry, so this is
    a transverse nearest-tower projection rather than a separate simulation
    of the original DREAM module.
    """

    centers = original_dream_pmt_centers() if centers_mm is None else np.asarray(
        centers_mm, dtype=np.float64
    )
    if centers.shape != (19, 2):
        raise ValueError("original DREAM layout must contain 19 x-y centers")
    if tower_radius_mm <= 0:
        raise ValueError("original DREAM tower radius must be positive")
    photon_xy = np.column_stack(
        (
            10.0 * ROD_PITCH_CM * event.rod_ix,
            10.0 * ROD_PITCH_CM * event.rod_iy,
        )
    )
    distance_squared = np.sum(
        (photon_xy[:, np.newaxis, :] - centers[np.newaxis, :, :]) ** 2,
        axis=2,
    )
    channel_ids = np.argmin(distance_squared, axis=1).astype(np.int64)

    # Keep only photons inside the nearest pointy-top regular hexagon.  This
    # prevents an edge tower from absorbing photons that fall outside the
    # physical 19-tower DREAM footprint.
    local_xy = photon_xy - centers[channel_ids]
    abs_x = np.abs(local_xy[:, 0])
    abs_y = np.abs(local_xy[:, 1])
    inside = (
        (abs_x <= 0.5 * np.sqrt(3.0) * tower_radius_mm)
        & (abs_y <= tower_radius_mm - abs_x / np.sqrt(3.0))
    )
    channel_ids[~inside] = -1
    return channel_ids


def localize_full_detector_coordinates(
    x_cm: np.ndarray,
    y_cm: np.ndarray,
    production_fiber: np.ndarray,
    *,
    rod_pitch_cm: float = ROD_PITCH_CM,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Map absolute DREAMSim endpoints to a fiber-local sensor coordinate.

    The four offsets are the measured centers used in Chris Madrid's 2025
    one-rod study.  In the full geometry those centers repeat at the 4 mm rod
    pitch.  Rounding the shifted coordinate to the nearest pitch therefore
    identifies the parent rod and leaves the local coordinate inside the
    fiber-end sensor.
    """

    x_cm = np.asarray(x_cm, dtype=np.float64)
    y_cm = np.asarray(y_cm, dtype=np.float64)
    fiber = np.asarray(production_fiber, dtype=np.int64)
    if x_cm.shape != y_cm.shape or x_cm.shape != fiber.shape:
        raise ValueError("x, y, and production_fiber must have the same shape")
    if rod_pitch_cm <= 0:
        raise ValueError("rod_pitch_cm must be positive")
    if fiber.size and (fiber.min() < 0 or fiber.max() >= 4):
        raise ValueError("production-fiber indices must be between zero and three")

    shifted_x_cm = x_cm + np.take(OLD_2025_X_SHIFT, fiber)
    shifted_y_cm = y_cm + np.take(OLD_2025_Y_SHIFT, fiber)
    rod_ix = np.rint(shifted_x_cm / rod_pitch_cm).astype(np.int64)
    rod_iy = np.rint(shifted_y_cm / rod_pitch_cm).astype(np.int64)
    local_x_mm = 10.0 * (shifted_x_cm - rod_pitch_cm * rod_ix)
    local_y_mm = 10.0 * (shifted_y_cm - rod_pitch_cm * rod_iy)

    # OLD_2025_*_SHIFT converts the raw DREAMSim frame into a fiber-local
    # frame.  Its common component is therefore a global detector translation,
    # while the differences between the four values encode the fiber centers
    # within a rod.  Remove the common component so sensor centers use the same
    # rod-centered detector frame as rod_ix/rod_iy, HG-DREAM, and DREAM towers.
    sensor_x_mm = 10.0 * (
        np.mean(OLD_2025_X_SHIFT)
        - np.take(OLD_2025_X_SHIFT, fiber)
        + rod_pitch_cm * rod_ix
    )
    sensor_y_mm = 10.0 * (
        np.mean(OLD_2025_Y_SHIFT)
        - np.take(OLD_2025_Y_SHIFT, fiber)
        + rod_pitch_cm * rod_iy
    )
    sensor_id = (
        ((rod_ix + SENSOR_INDEX_OFFSET) << 16)
        | ((rod_iy + SENSOR_INDEX_OFFSET) << 2)
        | fiber
    ).astype(np.int64)
    return (
        local_x_mm,
        local_y_mm,
        sensor_x_mm,
        sensor_y_mm,
        rod_ix,
        rod_iy,
        sensor_id,
    )


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


def _event_from_tree_entry(
    event: object,
    event_index: int,
    *,
    gate_ns: tuple[float, float],
) -> FullDetectorPhotonEvent:
    raw_time = np.asarray(event.OP_time_final)
    raw_count = int(raw_time.size)
    order = np.argsort(raw_time, kind="stable")
    time_ns = raw_time[order]
    fiber = np.asarray(event.OP_productionFiber)[order].astype(np.int64)
    core = np.asarray(event.OP_isCoreC)[order].astype(bool)
    final_x = np.asarray(event.OP_pos_final_x)[order]
    final_y = np.asarray(event.OP_pos_final_y)[order]
    final_z = np.asarray(event.OP_pos_final_z)[order]
    produced_z = np.asarray(event.OP_pos_produced_z)[order]
    momentum_x = np.asarray(event.OP_mom_produced_x)[order]
    momentum_y = np.asarray(event.OP_mom_produced_y)[order]
    momentum_z = np.asarray(event.OP_mom_produced_z)[order]

    gate_low, gate_high = gate_ns
    valid = (
        (fiber >= 0)
        & (fiber < 4)
        & core
        & (final_z > 0.0)
        & (time_ns > gate_low)
        & (time_ns < gate_high)
        & np.isfinite(final_x)
        & np.isfinite(final_y)
        & np.isfinite(time_ns)
    )
    selected_fiber = fiber[valid]
    (
        local_x_mm,
        local_y_mm,
        sensor_x_mm,
        sensor_y_mm,
        rod_ix,
        rod_iy,
        sensor_id,
    ) = localize_full_detector_coordinates(
        final_x[valid],
        final_y[valid],
        selected_fiber,
    )
    momentum = np.sqrt(
        momentum_x[valid] ** 2 + momentum_y[valid] ** 2 + momentum_z[valid] ** 2
    )
    wavelength_nm = np.divide(
        HC_GEV_NM,
        momentum,
        out=np.full_like(momentum, np.nan, dtype=np.float64),
        where=momentum > 0,
    )
    return FullDetectorPhotonEvent(
        event_index=event_index,
        beam_energy_gev=float(event.beamE) / 1000.0,
        c_signal=float(event.sum3dCC),
        s_signal=float(event.sum3dSS),
        ecalo_truth_gev=float(event.eCalotruth),
        eleak_truth_gev=float(event.eLeaktruth),
        raw_photon_count=raw_count,
        local_x_mm=local_x_mm,
        local_y_mm=local_y_mm,
        sensor_x_mm=sensor_x_mm,
        sensor_y_mm=sensor_y_mm,
        rod_ix=rod_ix,
        rod_iy=rod_iy,
        sensor_id=sensor_id,
        time_ns=time_ns[valid].astype(np.float64, copy=False),
        wavelength_nm=wavelength_nm,
        production_z_mm=(20.0 * produced_z[valid] + 2000.0).astype(
            np.float64,
            copy=False,
        ),
    )


def iter_full_detector_events(
    path: str | Path,
    *,
    max_events: int | None = None,
    gate_ns: tuple[float, float] = (0.0, 40.0),
) -> Iterator[FullDetectorPhotonEvent]:
    """Yield full-detector Cherenkov photon events from a DREAMSim ROOT file."""

    import ROOT

    source = Path(path).expanduser().resolve()
    root_file = ROOT.TFile.Open(str(source), "READ")
    if not root_file or root_file.IsZombie():
        raise OSError(f"Could not open ROOT file: {source}")
    try:
        tree = root_file.Get("tree")
        if not tree:
            raise ValueError(f"ROOT file has no 'tree': {source}")
        branches = {branch.GetName() for branch in tree.GetListOfBranches()}
        missing = sorted(REQUIRED_BRANCHES - branches)
        if missing:
            raise ValueError(
                f"ROOT file is missing required branches: {', '.join(missing)}"
            )
        tree.SetBranchStatus("*", 0)
        for branch in REQUIRED_BRANCHES:
            tree.SetBranchStatus(branch, 1)
        for event_index, event in enumerate(tree):
            if max_events is not None and event_index >= max_events:
                break
            yield _event_from_tree_entry(
                event,
                event_index,
                gate_ns=gate_ns,
            )
    finally:
        root_file.Close()


def iter_full_detector_events_uproot(
    path: str | Path,
    *,
    max_events: int | None = None,
    gate_ns: tuple[float, float] = (0.0, 40.0),
    event_index_offset: int = 0,
    step_size: str = "64 MB",
) -> Iterator[FullDetectorPhotonEvent]:
    """Yield full-detector events with uproot instead of PyROOT.

    This reader is intended for HPCC reduction jobs, where uproot and Awkward
    are available but the ROOT Python bindings are not.  Both readers call the
    same ``_event_from_tree_entry`` conversion so coordinate localization,
    photon selection, wavelength reconstruction, and units cannot drift.

    ``event_index_offset`` gives every event in a split production sample a
    stable sample-global index.  That index also enters the seeded sensor
    response, preventing repeated random streams when many ROOT files are
    reduced independently.
    """

    import awkward as ak
    import uproot

    source = Path(path).expanduser().resolve()
    if event_index_offset < 0:
        raise ValueError("event_index_offset cannot be negative")
    with uproot.open(source) as root_file:
        if "tree" not in root_file:
            raise ValueError(f"ROOT file has no 'tree': {source}")
        tree = root_file["tree"]
        branches = set(tree.keys())
        missing = sorted(REQUIRED_BRANCHES - branches)
        if missing:
            raise ValueError(
                f"ROOT file is missing required branches: {', '.join(missing)}"
            )

        emitted = 0
        for arrays in tree.iterate(
            sorted(REQUIRED_BRANCHES),
            step_size=step_size,
            library="ak",
        ):
            batch_entries = len(arrays["beamE"])
            for batch_index in range(batch_entries):
                if max_events is not None and emitted >= max_events:
                    return
                values: dict[str, object] = {}
                for branch in REQUIRED_BRANCHES:
                    item = arrays[branch][batch_index]
                    values[branch] = (
                        ak.to_numpy(item)
                        if branch.startswith("OP_")
                        else ak.to_list(item)
                    )
                yield _event_from_tree_entry(
                    SimpleNamespace(**values),
                    event_index_offset + emitted,
                    gate_ns=gate_ns,
                )
                emitted += 1
