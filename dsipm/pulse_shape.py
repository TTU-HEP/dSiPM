"""Pulse-shape utilities shared by the one-fiber sensor studies.

The functions in this module deliberately distinguish three quantities:

* photon/avalanche arrival histograms;
* the charge retained by a conventional analog SiPM during recharge; and
* a waveform-equivalent rendering obtained by convolving accepted avalanches
  with one common display kernel.

A dSiPM does not produce the analog waveform drawn by ``scope_waveform``.
That curve is only a convenient way to compare its time-tagged hits with the
mental picture supplied by an oscilloscope.
"""

from __future__ import annotations

import numpy as np


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
