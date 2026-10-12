"""Active-area and spectral-response helpers for FLASH requirement studies.

Pitch, photosensitive junction size, and photon-detection probability are
different design variables.  Keeping them separate prevents geometric loss
from being hidden inside a single effective PDE and makes the proposal's fill-
factor target directly testable.
"""

from __future__ import annotations

import numpy as np

from .digitizer import _cell_indices, legacy_pde


ACTIVE_AREA_SHAPES = {"square", "circle"}


def active_width_um(
    pitch_um: float,
    fill_factor: float,
    *,
    shape: str = "square",
) -> float:
    """Return square width or circular diameter for a nominal fill factor."""

    if pitch_um <= 0:
        raise ValueError("pitch_um must be positive")
    if not 0.0 <= fill_factor <= 1.0:
        raise ValueError("fill_factor must be between zero and one")
    if shape not in ACTIVE_AREA_SHAPES:
        raise ValueError(f"Unsupported active-area shape {shape!r}")
    if shape == "square":
        return float(pitch_um * np.sqrt(fill_factor))
    maximum_circle_fill = np.pi / 4.0
    if fill_factor > maximum_circle_fill + 1e-12:
        raise ValueError(
            "A circle contained in a square cell cannot exceed pi/4 fill factor"
        )
    return float(pitch_um * np.sqrt(4.0 * fill_factor / np.pi))


def active_area_mask(
    x_mm: np.ndarray,
    y_mm: np.ndarray,
    *,
    pitch_um: float,
    fill_factor: float,
    detector_min_mm: float,
    detector_max_mm: float,
    shape: str = "square",
) -> np.ndarray:
    """Return photons falling inside the centered active region of each cell.

    The square option is the default requirement proxy because it can span the
    complete zero-to-unity fill-factor range.  A fabricated layout, guard ring,
    microlens, and edge response should replace this idealized centered region
    once prototype geometry is available.
    """

    x_mm = np.asarray(x_mm, dtype=np.float64)
    y_mm = np.asarray(y_mm, dtype=np.float64)
    pitch_mm = pitch_um / 1000.0
    ix, iy, _, _ = _cell_indices(
        x_mm,
        y_mm,
        pitch_mm,
        detector_min_mm,
        detector_max_mm,
    )
    width_mm = active_width_um(pitch_um, fill_factor, shape=shape) / 1000.0
    center_x = detector_min_mm + (ix.astype(np.float64) + 0.5) * pitch_mm
    center_y = detector_min_mm + (iy.astype(np.float64) + 0.5) * pitch_mm
    dx = x_mm - center_x
    dy = y_mm - center_y
    if shape == "square":
        half_width = 0.5 * width_mm
        return (np.abs(dx) <= half_width) & (np.abs(dy) <= half_width)
    radius = 0.5 * width_mm
    return dx * dx + dy * dy <= radius * radius


def normalized_legacy_spectral_shape(wavelength_nm: np.ndarray) -> np.ndarray:
    """Return the capstone PDE curve normalized to unit peak.

    This is a traceable spectral-shape proxy, not a FLASH measurement.  It is
    zero outside the polynomial's documented 300--620 nm range.
    """

    wavelength_nm = np.asarray(wavelength_nm, dtype=np.float64)
    reference_grid = np.linspace(300.0, 620.0, 3201)
    peak = float(np.max(legacy_pde(reference_grid)))
    if peak <= 0:
        raise RuntimeError("Legacy spectral proxy has no positive support")
    return np.clip(legacy_pde(wavelength_nm) / peak, 0.0, 1.0)


def junction_pde_probability(
    wavelength_nm: np.ndarray,
    *,
    peak_probability: float,
    spectral_model: str = "legacy_shape",
) -> np.ndarray:
    """Detection probability conditional on a photon reaching active area."""

    if not 0.0 <= peak_probability <= 1.0:
        raise ValueError("peak_probability must be between zero and one")
    wavelength_nm = np.asarray(wavelength_nm, dtype=np.float64)
    if spectral_model == "flat":
        probability = np.full(wavelength_nm.shape, peak_probability)
    elif spectral_model == "legacy_shape":
        probability = peak_probability * normalized_legacy_spectral_shape(
            wavelength_nm
        )
    else:
        raise ValueError(f"Unsupported spectral model {spectral_model!r}")
    return np.clip(probability, 0.0, 1.0)
