"""Read DREAMSim optical photons and apply documented legacy geometry maps."""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

import numpy as np

from .digitizer import HC_GEV_NM, PhotonEvent


GEOMETRY_MODES = {"old_2025", "student_2026", "original_mickey", "raw"}

OLD_2025_X_SHIFT = np.asarray([3.2908, 3.2205, 3.3597, 3.2898])
OLD_2025_Y_SHIFT = np.asarray([-3.6878, -3.7287, -3.7284, -3.6075])
STUDENT_X_SHIFT = np.asarray([3.7308, 3.5768, 3.8008, 3.6468])
STUDENT_Y_SHIFT = np.asarray([-3.7293, -3.6878, -3.6878, -3.6488])
MICKEY_X_SHIFT = np.asarray([3.8523, 3.8523, 3.8523, 3.8523])
MICKEY_Y_SHIFT = np.asarray([-3.8588, -3.8588, -3.8588, -3.8588])

SHRINK_LIMITS = np.asarray([0.1 + 0.4 * index for index in range(40)])
SHRINK_AMOUNTS = np.asarray([round(0.23 * index, 2) for index in range(40)])

REQUIRED_BRANCHES = {
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


def shrink_toward_center(values: np.ndarray) -> np.ndarray:
    """Vectorized copy of the current student-repository geometry correction."""

    indices = np.searchsorted(SHRINK_LIMITS, np.abs(values), side="right")
    indices = np.clip(indices, 0, len(SHRINK_AMOUNTS) - 1)
    return values - SHRINK_AMOUNTS[indices] * np.sign(values)


def transform_coordinates(
    x: np.ndarray,
    y: np.ndarray,
    production_fiber: np.ndarray,
    mode: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Return detector-plane coordinates in millimeters."""

    if mode not in GEOMETRY_MODES:
        available = ", ".join(sorted(GEOMETRY_MODES))
        raise ValueError(f"Unknown geometry mode {mode!r}. Choose from {available}")
    if mode == "raw":
        return 10.0 * x, 10.0 * y

    if production_fiber.size and (
        production_fiber.min() < 0 or production_fiber.max() >= 4
    ):
        raise ValueError("Production-fiber indices must be between zero and three")

    if mode == "old_2025":
        shifted_x = x + np.take(OLD_2025_X_SHIFT, production_fiber)
        shifted_y = y + np.take(OLD_2025_Y_SHIFT, production_fiber)
    elif mode == "original_mickey":
        shifted_x = x + np.take(MICKEY_X_SHIFT, production_fiber)
        shifted_y = y + np.take(MICKEY_Y_SHIFT, production_fiber)
    else:
        shifted_x = x + np.take(STUDENT_X_SHIFT, production_fiber)
        shifted_y = y + np.take(STUDENT_Y_SHIFT, production_fiber)
        shifted_x = shrink_toward_center(shifted_x)
        shifted_y = shrink_toward_center(shifted_y)
    return 10.0 * shifted_x, 10.0 * shifted_y


def _array(event: object, branch: str, order: np.ndarray) -> np.ndarray:
    return np.asarray(getattr(event, branch))[order]


def _photon_event_from_tree_entry(
    event: object,
    event_index: int,
    *,
    geometry_mode: str,
    gate_ns: tuple[float, float],
    detector_min_mm: float,
    detector_max_mm: float,
    fiber_index: int | None,
) -> PhotonEvent:
    raw_time = np.asarray(event.OP_time_final)
    raw_count = int(raw_time.size)
    order = np.argsort(raw_time, kind="stable")

    time_ns = raw_time[order]
    production_fiber = _array(event, "OP_productionFiber", order).astype(np.int64)
    is_core_cherenkov = _array(event, "OP_isCoreC", order).astype(bool)
    final_x = _array(event, "OP_pos_final_x", order)
    final_y = _array(event, "OP_pos_final_y", order)
    final_z = _array(event, "OP_pos_final_z", order)
    produced_z = _array(event, "OP_pos_produced_z", order)

    momentum_x = _array(event, "OP_mom_produced_x", order)
    momentum_y = _array(event, "OP_mom_produced_y", order)
    momentum_z = _array(event, "OP_mom_produced_z", order)
    momentum = np.sqrt(momentum_x**2 + momentum_y**2 + momentum_z**2)
    wavelength_nm = np.divide(
        HC_GEV_NM,
        momentum,
        out=np.full_like(momentum, np.nan, dtype=np.float64),
        where=momentum > 0,
    )

    valid_fiber = (production_fiber >= 0) & (production_fiber < 4)
    gate_low, gate_high = gate_ns
    fiber_mask = (
        np.ones(production_fiber.size, dtype=bool)
        if fiber_index is None
        else production_fiber == fiber_index
    )
    physics_mask = (
        valid_fiber
        & fiber_mask
        & is_core_cherenkov
        & (final_z > 0)
        & (time_ns > gate_low)
        & (time_ns < gate_high)
    )

    selected_x, selected_y = transform_coordinates(
        final_x[physics_mask],
        final_y[physics_mask],
        production_fiber[physics_mask],
        geometry_mode,
    )
    in_detector = (
        (selected_x >= detector_min_mm)
        & (selected_x < detector_max_mm)
        & (selected_y >= detector_min_mm)
        & (selected_y < detector_max_mm)
    )

    return PhotonEvent(
        event_index=event_index,
        raw_photon_count=raw_count,
        x_mm=selected_x[in_detector].astype(np.float64, copy=False),
        y_mm=selected_y[in_detector].astype(np.float64, copy=False),
        time_ns=time_ns[physics_mask][in_detector].astype(np.float64, copy=False),
        wavelength_nm=wavelength_nm[physics_mask][in_detector].astype(
            np.float64,
            copy=False,
        ),
        production_z_mm=(20.0 * produced_z[physics_mask][in_detector] + 2000.0).astype(
            np.float64,
            copy=False,
        ),
    )


def inspect_root_file(path: str | Path) -> dict[str, object]:
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
        return {
            "path": str(source),
            "size_bytes": source.stat().st_size,
            "entries": int(tree.GetEntries()),
            "missing_required_branches": missing,
            "branch_count": len(branches),
        }
    finally:
        root_file.Close()


def iter_root_events(
    path: str | Path,
    *,
    max_events: int | None,
    geometry_mode: str,
    gate_ns: tuple[float, float],
    detector_min_mm: float,
    detector_max_mm: float,
    fiber_index: int | None = None,
) -> Iterator[PhotonEvent]:
    import ROOT

    source = Path(path).expanduser().resolve()
    if fiber_index is not None and fiber_index not in range(4):
        raise ValueError("fiber_index must be between zero and three")
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

        for event_index, event in enumerate(tree):
            if max_events is not None and event_index >= max_events:
                break
            yield _photon_event_from_tree_entry(
                event,
                event_index,
                geometry_mode=geometry_mode,
                gate_ns=gate_ns,
                detector_min_mm=detector_min_mm,
                detector_max_mm=detector_max_mm,
                fiber_index=fiber_index,
            )
    finally:
        root_file.Close()
