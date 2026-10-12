#!/usr/bin/env python3
"""Study the event-level effect of replacing HG-DREAM SiPMs with dSiPMs.

This is intentionally a detector-response study rather than another
single-fiber occupancy plot.  It uses the absolute fiber endpoints in the
full-detector DREAMSim samples, applies PDE and SPAD recovery independently
inside every physical Cherenkov fiber, calibrates each response on electrons,
and propagates it through a fixed dual-readout estimator for pions.

The available ROOT files contain explicit Cherenkov optical photons but not
explicit scintillation photons.  The scintillation branch is therefore kept
at the DREAMSim signal level.  All figures and tables label that limitation.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from dsipm.cms_style import (
    draw_flash_label,
    draw_panel_title,
    make_canvas,
    make_legend,
    save_canvas,
    set_flash_style,
    style_axes,
    style_pad,
)
from dsipm.full_detector import (
    AnalogSiPMScenario,
    FullDetectorPhotonEvent,
    FullDetectorScenario,
    ORIGINAL_DREAM_FIBER_DIAMETER_MM,
    ORIGINAL_DREAM_TOWER_RADIUS_MM,
    digitize_analog_sipm_event,
    digitize_full_detector_event,
    hgdream_cherenkov_channel_ids,
    hgdream_cherenkov_channel_layout,
    iter_full_detector_events,
    original_dream_pmt_centers,
    original_dream_pmt_channel_ids,
    scan_full_detector_event,
)
from dsipm.workflow import load_sensor_workflow


REPO_DIR = Path(__file__).resolve().parent
DEFAULT_ELECTRON = REPO_DIR / "mc_dreamsim_e+_0_run100_0_FullDetector_10evt_e+_100_101.root"
DEFAULT_PION = REPO_DIR / "sim_200evt_pi+_90+rod.root"
WORKFLOW = load_sensor_workflow()
PITCHES_UM = np.asarray(WORKFLOW.scans.full_detector_spad_pitch_um)
RECOVERY_NS = np.asarray(WORKFLOW.scans.full_detector_recovery_ns)
PDE = WORKFLOW.reference.effective_pde
TARGET_PITCH_UM = WORKFLOW.reference.spad_pitch_um
TARGET_RECOVERY_NS = WORKFLOW.reference.spad_recovery_ns
COARSE_PITCH_UM = 50.0
# The physical 0.8 mm fibers are sub-pixel in a four-panel, full-detector
# figure.  Use a larger marker only for the event-display rendering; the
# positions, photon counts, and 4 mm2 assigned-area normalization are unchanged.
IDEAL_FIBER_DISPLAY_DIAMETER_MM = 2.0
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


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Cannot write empty table: {path}")
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _read_csv(path: Path) -> list[dict[str, object]]:
    with path.open(newline="") as handle:
        rows = [dict(row) for row in csv.DictReader(handle)]
    if not rows:
        raise ValueError(f"Cached event table is empty: {path}")
    return rows


def _validate_cached_workflow(output_dir: Path) -> None:
    """Reject event tables produced with different sensor assumptions."""

    summary_path = output_dir / "full_detector_summary.json"
    if not summary_path.exists():
        raise FileNotFoundError(
            "Cannot reuse event tables without full_detector_summary.json; "
            "run once without --reuse-event-tables"
        )
    with summary_path.open() as handle:
        summary = json.load(handle)
    cached_hash = summary.get("sensor_workflow", {}).get("sha256")
    if cached_hash != WORKFLOW.sha256:
        raise ValueError(
            "Cached event tables were produced with a different or unrecorded "
            "sensor workflow; run without --reuse-event-tables"
        )


def _accepted_cube(rows: list[dict[str, object]]) -> np.ndarray:
    """Rebuild the pitch-by-recovery count cube from a cached event table."""

    return np.asarray(
        [
            [
                [
                    float(row[_scenario_column(float(pitch), float(recovery))])
                    for recovery in RECOVERY_NS
                ]
                for pitch in PITCHES_UM
            ]
            for row in rows
        ],
        dtype=np.float64,
    )


def _scenario_column(pitch_um: float, recovery_ns: float) -> str:
    return f"accepted_pitch{pitch_um:g}um_recovery{recovery_ns:g}ns"


def _read_and_scan(
    path: Path,
    *,
    max_events: int | None,
    seed: int,
) -> tuple[list[dict[str, object]], np.ndarray]:
    rows: list[dict[str, object]] = []
    accepted_rows: list[np.ndarray] = []
    for event in iter_full_detector_events(path, max_events=max_events):
        detected, accepted = scan_full_detector_event(
            event,
            pde=PDE,
            pitches_um=PITCHES_UM,
            recovery_times_ns=RECOVERY_NS,
            seed=seed,
        )
        analog = digitize_analog_sipm_event(event, ANALOG_SIPM, seed=seed)
        row: dict[str, object] = {
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
            "detected_reference": detected,
            "analog_detected_at_31pct_pde": analog.n_detected,
            "analog_charge_pe": analog.charge_pe,
            "analog_charge_survival": analog.charge_survival,
        }
        for pitch_index, pitch_um in enumerate(PITCHES_UM):
            for recovery_index, recovery_ns in enumerate(RECOVERY_NS):
                row[_scenario_column(float(pitch_um), float(recovery_ns))] = int(
                    accepted[pitch_index, recovery_index]
                )
        rows.append(row)
        accepted_rows.append(accepted)
    if not rows:
        raise ValueError(f"No events were read from {path}")
    return rows, np.stack(accepted_rows)


def _column(rows: list[dict[str, object]], name: str) -> np.ndarray:
    return np.asarray([float(row[name]) for row in rows], dtype=np.float64)


def _response_summary(values: np.ndarray) -> dict[str, float]:
    values = np.asarray(values, dtype=np.float64)
    q16, median, q84 = np.quantile(values, [0.16, 0.50, 0.84])
    mean = float(np.mean(values))
    rms = float(np.std(values, ddof=1)) if values.size > 1 else 0.0
    return {
        "mean": mean,
        "rms": rms,
        "rms_over_mean": rms / mean if mean else 0.0,
        "median": float(median),
        "central68_half_width": float(0.5 * (q84 - q16)),
        "central68_over_median": float(0.5 * (q84 - q16) / median) if median else 0.0,
    }


def _bootstrap_metrics(
    values: np.ndarray,
    beam_energy: float,
    *,
    samples: int,
    seed: int,
) -> dict[str, float]:
    values = np.asarray(values, dtype=np.float64)
    if samples <= 0 or values.size < 2:
        return {"mean_bias_error": 0.0, "resolution_error": 0.0}
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, values.size, size=(samples, values.size))
    draws = values[indices]
    means = np.mean(draws, axis=1)
    widths = np.std(draws, axis=1, ddof=1) / means
    return {
        "mean_bias_error": float(np.std(means / beam_energy - 1.0, ddof=1)),
        "resolution_error": float(np.std(widths, ddof=1)),
    }


def _bootstrap_paired_c_metrics(
    values: np.ndarray,
    baseline: np.ndarray,
    *,
    samples: int,
    seed: int,
) -> dict[str, float]:
    """Bootstrap a response shift and width using paired pion events."""

    values = np.asarray(values, dtype=np.float64)
    baseline = np.asarray(baseline, dtype=np.float64)
    if values.shape != baseline.shape:
        raise ValueError("values and baseline must have the same shape")
    if samples <= 0 or values.size < 2:
        return {"mean_shift_error": 0.0, "resolution_error": 0.0}
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, values.size, size=(samples, values.size))
    draws = values[indices]
    baseline_draws = baseline[indices]
    means = np.mean(draws, axis=1)
    baseline_means = np.mean(baseline_draws, axis=1)
    shifts = means / baseline_means - 1.0
    widths = np.std(draws, axis=1, ddof=1) / means
    return {
        "mean_shift_error": float(np.std(shifts, ddof=1)),
        "resolution_error": float(np.std(widths, ddof=1)),
    }


def _calibrated_responses(
    electron_rows: list[dict[str, object]],
    pion_rows: list[dict[str, object]],
    electron_accepted: np.ndarray,
    pion_accepted: np.ndarray,
    *,
    bootstrap_samples: int,
) -> tuple[list[dict[str, object]], dict[str, np.ndarray], dict[str, object]]:
    electron_beam = float(np.mean(_column(electron_rows, "beam_energy_gev")))
    pion_beam = float(np.mean(_column(pion_rows, "beam_energy_gev")))
    electron_s = _column(electron_rows, "sum3dSS")
    electron_c = _column(electron_rows, "sum3dCC")
    electron_selected = _column(electron_rows, "selected_cherenkov_photons")
    pion_s = _column(pion_rows, "sum3dSS")
    pion_c = _column(pion_rows, "sum3dCC")
    pion_selected = _column(pion_rows, "selected_cherenkov_photons")

    s_scale = electron_beam / float(np.mean(electron_s))
    current_c_scale = electron_beam / float(np.mean(electron_c))
    pion_s_energy = pion_s * s_scale
    pion_current_c_energy = pion_c * current_c_scale
    current_denominator = float(np.mean(pion_current_c_energy) - pion_beam)
    if np.isclose(current_denominator, 0.0):
        raise ValueError("Cannot derive dual-readout chi from the current response")
    current_chi = float(
        (np.mean(pion_s_energy) - pion_beam) / current_denominator
    )
    if np.isclose(current_chi, 1.0):
        raise ValueError("Derived dual-readout chi is singular")
    pion_current_dr = (
        pion_s_energy - current_chi * pion_current_c_energy
    ) / (1.0 - current_chi)

    # The explicit photon branch and sum3dCC are different DREAMSim
    # observables.  Derive a second, internally consistent baseline from the
    # explicit photons and keep its chi fixed for every dSiPM scenario.  This
    # prevents a difference between ROOT branches from being misidentified as
    # a sensor nonlinearity.
    ideal_c_scale = electron_beam / float(np.mean(electron_selected))
    pion_ideal_c_energy = pion_selected * ideal_c_scale
    ideal_denominator = float(np.mean(pion_ideal_c_energy) - pion_beam)
    if np.isclose(ideal_denominator, 0.0):
        raise ValueError("Cannot derive dual-readout chi from the photon baseline")
    photon_chi = float((np.mean(pion_s_energy) - pion_beam) / ideal_denominator)
    if np.isclose(photon_chi, 1.0):
        raise ValueError("Photon-baseline dual-readout chi is singular")
    pion_ideal_dr = (
        pion_s_energy - photon_chi * pion_ideal_c_energy
    ) / (1.0 - photon_chi)

    electron_detected = _column(electron_rows, "detected_reference")
    pion_detected = _column(pion_rows, "detected_reference")
    pde_only_scale = electron_beam / float(np.mean(electron_detected))
    pion_pde_only_c_energy = pion_detected * pde_only_scale
    pion_pde_only_dr = (
        pion_s_energy - photon_chi * pion_pde_only_c_energy
    ) / (1.0 - photon_chi)

    electron_analog_detected = _column(
        electron_rows, "analog_detected_at_31pct_pde"
    )
    pion_analog_detected = _column(pion_rows, "analog_detected_at_31pct_pde")
    electron_analog_charge = _column(electron_rows, "analog_charge_pe")
    pion_analog_charge = _column(pion_rows, "analog_charge_pe")
    analog_pde_only_scale = electron_beam / float(np.mean(electron_analog_detected))
    pion_analog_pde_only_c_energy = pion_analog_detected * analog_pde_only_scale
    pion_analog_pde_only_dr = (
        pion_s_energy - photon_chi * pion_analog_pde_only_c_energy
    ) / (1.0 - photon_chi)
    analog_scale = electron_beam / float(np.mean(electron_analog_charge))
    pion_analog_c_energy = pion_analog_charge * analog_scale
    pion_analog_dr = (
        pion_s_energy - photon_chi * pion_analog_c_energy
    ) / (1.0 - photon_chi)

    metrics: list[dict[str, object]] = []
    response_arrays: dict[str, np.ndarray] = {
        "pion_s_energy": pion_s_energy,
        "pion_current_c_energy": pion_current_c_energy,
        "pion_current_dr_energy": pion_current_dr,
        "pion_ideal_c_energy": pion_ideal_c_energy,
        "pion_ideal_dr_energy": pion_ideal_dr,
        "pion_pde_only_c_energy": pion_pde_only_c_energy,
        "pion_pde_only_dr_energy": pion_pde_only_dr,
        "pion_analog_pde_only_c_energy": pion_analog_pde_only_c_energy,
        "pion_analog_pde_only_dr_energy": pion_analog_pde_only_dr,
        "pion_analog_c_energy": pion_analog_c_energy,
        "pion_analog_dr_energy": pion_analog_dr,
    }

    current_summary = _response_summary(pion_current_dr)
    current_bootstrap = _bootstrap_metrics(
        pion_current_dr,
        pion_beam,
        samples=bootstrap_samples,
        seed=20261011,
    )
    metrics.append(
        {
            "scenario": "current_analog_C_plus_S",
            "pde": "",
            "pitch_um": "",
            "recovery_ns": "",
            "electron_C_calibration_per_signal": current_c_scale,
            "pion_C_mean_gev": float(np.mean(pion_current_c_energy)),
            "pion_DR_mean_gev": current_summary["mean"],
            "pion_DR_mean_bias": current_summary["mean"] / pion_beam - 1.0,
            "pion_DR_mean_bias_error": current_bootstrap["mean_bias_error"],
            "pion_DR_rms_over_mean": current_summary["rms_over_mean"],
            "pion_DR_rms_over_mean_error": current_bootstrap["resolution_error"],
            "pion_DR_central68_over_median": current_summary[
                "central68_over_median"
            ],
        }
    )

    analog_summary = _response_summary(pion_analog_dr)
    analog_c_summary = _response_summary(pion_analog_c_energy)
    analog_bootstrap = _bootstrap_metrics(
        pion_analog_dr,
        pion_beam,
        samples=bootstrap_samples,
        seed=20261014,
    )
    analog_c_bootstrap = _bootstrap_paired_c_metrics(
        pion_analog_c_energy,
        pion_analog_pde_only_c_energy,
        samples=bootstrap_samples,
        seed=20261015,
    )
    metrics.append(
        {
            "scenario": ANALOG_SIPM.name,
            "source_url": ANALOG_SIPM.source_url,
            "pde": ANALOG_SIPM.pde,
            "pde_wavelength_nm": ANALOG_SIPM.pde_wavelength_nm,
            "overvoltage_v": ANALOG_SIPM.overvoltage_v,
            "pitch_um": "",
            "microcell_size_um": ANALOG_SIPM.microcell_size_um,
            "effective_pitch_um": ANALOG_SIPM.effective_pitch_um,
            "n_microcells_datasheet": ANALOG_SIPM.n_microcells,
            "n_microcells_modeled": ANALOG_SIPM.modeled_microcells,
            "recovery_ns": "",
            "recovery_tau_ns": ANALOG_SIPM.recovery_tau_ns,
            "fast_output_fwhm_ns": ANALOG_SIPM.fast_output_fwhm_ns,
            "electron_C_calibration_per_signal": analog_scale,
            "electron_recovery_survival": float(
                np.sum(electron_analog_charge) / np.sum(electron_analog_detected)
            ),
            "pion_recovery_survival": float(
                np.sum(pion_analog_charge) / np.sum(pion_analog_detected)
            ),
            "pion_C_mean_gev": analog_c_summary["mean"],
            "pion_C_mean_shift_vs_pde_only": analog_c_summary["mean"]
            / float(np.mean(pion_analog_pde_only_c_energy))
            - 1.0,
            "pion_C_mean_shift_error": analog_c_bootstrap["mean_shift_error"],
            "pion_C_rms_over_mean": analog_c_summary["rms_over_mean"],
            "pion_C_rms_over_mean_error": analog_c_bootstrap["resolution_error"],
            "pion_DR_mean_gev": analog_summary["mean"],
            "pion_DR_mean_bias": analog_summary["mean"] / pion_beam - 1.0,
            "pion_DR_mean_bias_error": analog_bootstrap["mean_bias_error"],
            "pion_DR_rms_over_mean": analog_summary["rms_over_mean"],
            "pion_DR_rms_over_mean_error": analog_bootstrap["resolution_error"],
            "pion_DR_central68_over_median": analog_summary[
                "central68_over_median"
            ],
        }
    )

    pde_only_summary = _response_summary(pion_pde_only_dr)
    pde_only_c_summary = _response_summary(pion_pde_only_c_energy)
    pde_only_bootstrap = _bootstrap_metrics(
        pion_pde_only_dr,
        pion_beam,
        samples=bootstrap_samples,
        seed=20261013,
    )
    metrics.append(
        {
            "scenario": "reference_pde_only_no_recovery",
            "pde": PDE,
            "pitch_um": "",
            "recovery_ns": 0.0,
            "electron_C_calibration_per_signal": pde_only_scale,
            "electron_recovery_survival": 1.0,
            "pion_recovery_survival": 1.0,
            "pion_C_mean_gev": pde_only_c_summary["mean"],
            "pion_C_mean_shift_vs_pde_only": 0.0,
            "pion_C_mean_shift_error": 0.0,
            "pion_C_rms_over_mean": pde_only_c_summary["rms_over_mean"],
            "pion_C_rms_over_mean_error": 0.0,
            "pion_DR_mean_gev": pde_only_summary["mean"],
            "pion_DR_mean_bias": pde_only_summary["mean"] / pion_beam - 1.0,
            "pion_DR_mean_bias_error": pde_only_bootstrap["mean_bias_error"],
            "pion_DR_rms_over_mean": pde_only_summary["rms_over_mean"],
            "pion_DR_rms_over_mean_error": pde_only_bootstrap["resolution_error"],
            "pion_DR_central68_over_median": pde_only_summary[
                "central68_over_median"
            ],
        }
    )

    ideal_summary = _response_summary(pion_ideal_dr)
    ideal_bootstrap = _bootstrap_metrics(
        pion_ideal_dr,
        pion_beam,
        samples=bootstrap_samples,
        seed=20261012,
    )
    metrics.append(
        {
            "scenario": "ideal_explicit_photon_C_plus_S",
            "pde": 1.0,
            "pitch_um": "",
            "recovery_ns": "",
            "electron_C_calibration_per_signal": ideal_c_scale,
            "pion_C_mean_gev": float(np.mean(pion_ideal_c_energy)),
            "pion_DR_mean_gev": ideal_summary["mean"],
            "pion_DR_mean_bias": ideal_summary["mean"] / pion_beam - 1.0,
            "pion_DR_mean_bias_error": ideal_bootstrap["mean_bias_error"],
            "pion_DR_rms_over_mean": ideal_summary["rms_over_mean"],
            "pion_DR_rms_over_mean_error": ideal_bootstrap["resolution_error"],
            "pion_DR_central68_over_median": ideal_summary[
                "central68_over_median"
            ],
        }
    )

    for pitch_index, pitch_um in enumerate(PITCHES_UM):
        for recovery_index, recovery_ns in enumerate(RECOVERY_NS):
            electron_counts = electron_accepted[:, pitch_index, recovery_index].astype(
                np.float64
            )
            pion_counts = pion_accepted[:, pitch_index, recovery_index].astype(np.float64)
            digital_scale = electron_beam / float(np.mean(electron_counts))
            pion_digital_c_energy = pion_counts * digital_scale
            pion_digital_dr = (
                pion_s_energy - photon_chi * pion_digital_c_energy
            ) / (1.0 - photon_chi)
            summary = _response_summary(pion_digital_dr)
            c_summary = _response_summary(pion_digital_c_energy)
            bootstrap = _bootstrap_metrics(
                pion_digital_dr,
                pion_beam,
                samples=bootstrap_samples,
                seed=20261100 + 10 * pitch_index + recovery_index,
            )
            c_bootstrap = _bootstrap_paired_c_metrics(
                pion_digital_c_energy,
                pion_pde_only_c_energy,
                samples=bootstrap_samples,
                seed=20261200 + 10 * pitch_index + recovery_index,
            )
            scenario = FullDetectorScenario(
                pde=PDE,
                pitch_um=float(pitch_um),
                recovery_ns=float(recovery_ns),
            )
            metrics.append(
                {
                    "scenario": scenario.name,
                    "pde": PDE,
                    "pitch_um": float(pitch_um),
                    "recovery_ns": float(recovery_ns),
                    "electron_C_calibration_per_signal": digital_scale,
                    "electron_recovery_survival": float(
                        np.sum(electron_counts) / np.sum(electron_detected)
                    ),
                    "pion_recovery_survival": float(
                        np.sum(pion_counts) / np.sum(pion_detected)
                    ),
                    "pion_C_mean_gev": c_summary["mean"],
                    "pion_C_mean_shift_vs_pde_only": c_summary["mean"]
                    / pde_only_c_summary["mean"]
                    - 1.0,
                    "pion_C_mean_shift_error": c_bootstrap["mean_shift_error"],
                    "pion_C_rms_over_mean": c_summary["rms_over_mean"],
                    "pion_C_rms_over_mean_error": c_bootstrap[
                        "resolution_error"
                    ],
                    "pion_DR_mean_gev": summary["mean"],
                    "pion_DR_mean_bias": summary["mean"] / pion_beam - 1.0,
                    "pion_DR_mean_bias_error": bootstrap["mean_bias_error"],
                    "pion_DR_rms_over_mean": summary["rms_over_mean"],
                    "pion_DR_rms_over_mean_error": bootstrap["resolution_error"],
                    "pion_DR_central68_over_median": summary[
                        "central68_over_median"
                    ],
                }
            )
            response_arrays[f"pion_c_{scenario.name}"] = pion_digital_c_energy
            response_arrays[f"pion_dr_{scenario.name}"] = pion_digital_dr

    summary = {
        "electron_beam_energy_gev_mean": electron_beam,
        "pion_beam_energy_gev_mean": pion_beam,
        "scintillation_electron_calibration_gev_per_signal": s_scale,
        "current_cherenkov_electron_calibration_gev_per_signal": current_c_scale,
        "ideal_photon_electron_calibration_gev_per_photon": ideal_c_scale,
        "reference_pde_only_electron_calibration_gev_per_detected_photon": pde_only_scale,
        "analog_sipm_electron_calibration_gev_per_charge_pe": analog_scale,
        "analog_sipm_pde_only_electron_calibration_gev_per_detected_photon": analog_pde_only_scale,
        "dual_readout_chi_current_stored_response": current_chi,
        "fixed_dual_readout_chi_explicit_photon_baseline": photon_chi,
        "selected_photon_vs_sum3dCC_correlation": {
            "electron": float(np.corrcoef(electron_selected, electron_c)[0, 1]),
            "pion": float(np.corrcoef(pion_selected, pion_c)[0, 1]),
        },
    }
    return metrics, response_arrays, summary


def _histogram(
    ROOT,
    name: str,
    values: np.ndarray,
    bins: int,
    low: float,
    high: float,
    *,
    color: int,
    style: int = 1,
):
    hist = ROOT.TH1D(name, ";;Fraction of events / bin", bins, low, high)
    hist.SetDirectory(0)
    hist.Sumw2()
    for value in values:
        hist.Fill(float(value))
    if hist.Integral() > 0:
        hist.Scale(1.0 / hist.Integral())
    hist.SetLineColor(color)
    hist.SetLineWidth(3)
    hist.SetLineStyle(style)
    style_axes(hist)
    return hist


def _plot_validation(ROOT, palette, electron_rows, pion_rows, output_dir: Path):
    canvas = make_canvas(ROOT, "c_validation", width=1400, height=650)
    canvas.Divide(2, 1, 0.015, 0.0)
    keep = []
    for pad_index, (label, rows, color, right) in enumerate(
        [
            ("electron", electron_rows, palette["blue"], "100 GeV e^{+}"),
            ("pion", pion_rows, palette["red"], "90 GeV #pi^{+}"),
        ],
        start=1,
    ):
        pad = canvas.cd(pad_index)
        style_pad(pad, left=0.155, right=0.045, top=0.135, bottom=0.16)
        x = _column(rows, "selected_cherenkov_photons")
        y = _column(rows, "sum3dCC")
        graph = ROOT.TGraph(len(x), np.ascontiguousarray(x), np.ascontiguousarray(y))
        graph.SetMarkerStyle(20)
        graph.SetMarkerSize(0.55 if label == "electron" else 0.75)
        graph.SetMarkerColorAlpha(color, 0.55)
        graph.SetLineColor(color)
        graph.SetTitle(";Selected optical photons;Stored Cherenkov signal [a.u.]")
        graph.Draw("AP")
        style_axes(graph)
        fit_values = np.polyfit(x, y, 1)
        line = ROOT.TLine(float(np.min(x)), float(np.polyval(fit_values, np.min(x))), float(np.max(x)), float(np.polyval(fit_values, np.max(x))))
        line.SetLineColor(palette["black"])
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.Draw()
        keep.extend([graph, line])
        keep.extend(draw_flash_label(ROOT, pad, extra="Simulation", right=right))
        correlation = float(np.corrcoef(x, y)[0, 1])
        text = ROOT.TLatex()
        text.SetNDC(True)
        text.SetTextFont(42)
        text.SetTextSize(0.043)
        text.DrawLatex(0.19, 0.82, f"r = {correlation:.2f},  N = {len(rows)}")
        keep.append(text)
    canvas.cd()
    save_canvas(canvas, output_dir / "full_detector_validation")
    return keep


def _plot_occupancy(ROOT, palette, electron_rows, pion_rows, output_dir: Path):
    canvas = make_canvas(ROOT, "c_occupancy", width=1400, height=650)
    canvas.Divide(2, 1, 0.015, 0.0)
    keep = []
    pad = canvas.cd(1)
    style_pad(pad, left=0.16, right=0.045, top=0.23, bottom=0.16)
    maximum = max(
        float(np.max(_column(electron_rows, "max_photons_in_one_fiber"))),
        float(np.max(_column(pion_rows, "max_photons_in_one_fiber"))),
    )
    high = 50.0 * np.ceil(maximum / 50.0)
    bins = int(high / 50.0)
    h_e = _histogram(
        ROOT,
        "h_max_e",
        _column(electron_rows, "max_photons_in_one_fiber"),
        bins,
        0.0,
        high,
        color=palette["blue"],
    )
    h_pi = _histogram(
        ROOT,
        "h_max_pi",
        _column(pion_rows, "max_photons_in_one_fiber"),
        bins,
        0.0,
        high,
        color=palette["red"],
    )
    h_e.GetXaxis().SetTitle("Hottest-fiber occupancy [photons / event]")
    h_e.SetMaximum(1.45 * max(h_e.GetMaximum(), h_pi.GetMaximum()))
    h_e.Draw("HIST")
    h_pi.Draw("HIST SAME")
    legend = make_legend(ROOT, 0.55, 0.58, 0.93, 0.73, text_size=0.037)
    legend.AddEntry(h_e, f"100 GeV e^{{+}}  (N={len(electron_rows)})", "l")
    legend.AddEntry(h_pi, f"90 GeV #pi^{{+}}  (N={len(pion_rows)})", "l")
    legend.Draw()
    keep.extend([h_e, h_pi, legend])
    keep.extend(draw_flash_label(ROOT, pad, extra="Simulation", right="HG-DREAM full detector"))
    note = ROOT.TLatex()
    note.SetNDC(True)
    note.SetTextFont(42)
    note.SetTextSize(0.034)
    note.DrawLatex(0.19, 0.86, "Brightest Cherenkov fiber in each shower")
    note.DrawLatex(
        0.19,
        0.81,
        f"Median: e^{{+}} {np.median(_column(electron_rows, 'max_photons_in_one_fiber')):.0f},  "
        f"#pi^{{+}} {np.median(_column(pion_rows, 'max_photons_in_one_fiber')):.0f} photons",
    )
    keep.append(note)

    pad = canvas.cd(2)
    style_pad(pad, left=0.16, right=0.045, top=0.135, bottom=0.16)
    scatter_graphs = []
    series = []
    for rows, color, marker, label, response_name, detected_name in [
        (
            electron_rows,
            palette["blue"],
            20,
            "FLASH, 100 GeV e^{+}",
            _scenario_column(TARGET_PITCH_UM, TARGET_RECOVERY_NS),
            "detected_reference",
        ),
        (
            pion_rows,
            palette["red"],
            21,
            "FLASH, 90 GeV #pi^{+}",
            _scenario_column(TARGET_PITCH_UM, TARGET_RECOVERY_NS),
            "detected_reference",
        ),
        (
            electron_rows,
            palette["blue"],
            24,
            "Analog SiPM, 100 GeV e^{+}",
            "analog_charge_pe",
            "analog_detected_at_31pct_pde",
        ),
        (
            pion_rows,
            palette["red"],
            25,
            "Analog SiPM, 90 GeV #pi^{+}",
            "analog_charge_pe",
            "analog_detected_at_31pct_pde",
        ),
    ]:
        x = _column(rows, "max_photons_in_one_fiber")
        detected = _column(rows, detected_name)
        response = _column(rows, response_name)
        y = np.divide(response, detected, out=np.zeros_like(response), where=detected > 0)
        series.append((x, 100.0 * y, color, marker, label))

    nonzero_survival = np.concatenate([y[y > 0] for _, y, _, _, _ in series])
    lower_survival = min(95.0, max(70.0, float(np.floor(np.min(nonzero_survival)) - 1.0)))
    for index, (x, y, color, marker, label) in enumerate(series):
        graph = ROOT.TGraph(len(x), np.ascontiguousarray(x), np.ascontiguousarray(y))
        graph.SetMarkerStyle(marker)
        graph.SetMarkerSize(0.48 if marker in (20, 21) else 0.62)
        graph.SetMarkerColorAlpha(color, 0.42 if index < 2 else 0.70)
        graph.SetLineColor(color)
        if index == 0:
            graph.SetTitle(";Hottest-fiber occupancy [incident photons / event];Event signal retained [%]")
            graph.Draw("AP")
            style_axes(graph)
            graph.GetXaxis().SetLimits(0.0, 1.03 * high)
            graph.GetYaxis().SetRangeUser(lower_survival, 100.35)
        else:
            graph.Draw("P SAME")
        keep.append(graph)
        scatter_graphs.append(graph)
    unity = ROOT.TLine(0.0, 100.0, 1.03 * high, 100.0)
    unity.SetLineColor(palette["gray"])
    unity.SetLineStyle(2)
    unity.Draw()
    keep.append(unity)
    legend2 = make_legend(ROOT, 0.19, 0.19, 0.66, 0.42, text_size=0.031)
    for graph, (_, _, _, _, label) in zip(scatter_graphs, series, strict=True):
        legend2.AddEntry(graph, label, "p")
    legend2.Draw()
    keep.append(legend2)
    keep.extend(draw_flash_label(ROOT, pad, extra="Simulation", right="FLASH vs analog"))
    save_canvas(canvas, output_dir / "full_detector_occupancy")
    return keep


def _scenario_name(pitch_um: float, recovery_ns: float) -> str:
    return FullDetectorScenario(PDE, pitch_um, recovery_ns).name


def _plot_response(ROOT, palette, arrays, pion_beam: float, output_dir: Path):
    target = _scenario_name(TARGET_PITCH_UM, TARGET_RECOVERY_NS)
    coarse = _scenario_name(COARSE_PITCH_UM, TARGET_RECOVERY_NS)
    curves = [
        (
            f"{100 * PDE:g}% PDE only (linear)",
            arrays["pion_pde_only_c_energy"] / pion_beam,
            palette["black"],
            2,
        ),
        ("Analog C-Series proxy", arrays["pion_analog_c_energy"] / pion_beam, palette["purple"], 7),
        (
            f"FLASH {TARGET_PITCH_UM:g} #mum",
            arrays[f"pion_c_{target}"] / pion_beam,
            palette["blue"],
            1,
        ),
        ("FLASH 50 #mum", arrays[f"pion_c_{coarse}"] / pion_beam, palette["red"], 2),
    ]
    dr_curves = [
        (
            f"{100 * PDE:g}% PDE only (linear)",
            arrays["pion_pde_only_dr_energy"] / pion_beam,
            palette["black"],
            2,
        ),
        ("Analog C-Series proxy", arrays["pion_analog_dr_energy"] / pion_beam, palette["purple"], 7),
        (
            f"FLASH {TARGET_PITCH_UM:g} #mum",
            arrays[f"pion_dr_{target}"] / pion_beam,
            palette["blue"],
            1,
        ),
        ("FLASH 50 #mum", arrays[f"pion_dr_{coarse}"] / pion_beam, palette["red"], 2),
    ]
    canvas = make_canvas(ROOT, "c_response", width=1400, height=650)
    canvas.Divide(2, 1, 0.015, 0.0)
    keep = []
    for pad_index, (items, low, high, xlabel, right) in enumerate(
        [
            (curves, 0.20, 1.05, "Electron-calibrated C response / E_{beam}", "Cherenkov channel"),
            (dr_curves, 0.40, 1.50, "Illustrative dual-readout response / E_{beam}", "fixed #chi from photon baseline"),
        ],
        start=1,
    ):
        pad = canvas.cd(pad_index)
        style_pad(pad, left=0.16, right=0.045, top=0.22, bottom=0.16)
        histograms = []
        for curve_index, (label, values, color, style) in enumerate(items):
            hist = _histogram(
                ROOT,
                f"h_response_{pad_index}_{curve_index}",
                values,
                14,
                low,
                high,
                color=color,
                style=style,
            )
            hist.GetXaxis().SetTitle(xlabel)
            histograms.append(hist)
        ymax = 1.35 * max(hist.GetMaximum() for hist in histograms)
        histograms[0].SetMaximum(ymax)
        histograms[0].Draw("HIST")
        for hist in histograms[1:]:
            hist.Draw("HIST SAME")
        reference = ROOT.TLine(1.0, 0.0, 1.0, ymax)
        reference.SetLineColor(palette["gray"])
        reference.SetLineStyle(7)
        reference.SetLineWidth(2)
        reference.Draw()
        legend = make_legend(ROOT, 0.18, 0.48, 0.66, 0.73, text_size=0.034)
        for hist, item in zip(histograms, items, strict=True):
            legend.AddEntry(hist, item[0], "l")
        legend.Draw()
        keep.extend(histograms)
        keep.extend([reference, legend])
        keep.extend(
            draw_flash_label(
                ROOT,
                pad,
                extra="Simulation",
            )
        )
        sample = ROOT.TLatex()
        sample.SetNDC(True)
        sample.SetTextFont(42)
        sample.SetTextSize(0.029)
        sample.SetTextAlign(31)
        sample.DrawLatex(0.95, 0.865, f"90 GeV #pi^{{+}}, N={len(items[0][1])}")
        detail = ROOT.TLatex()
        detail.SetNDC(True)
        detail.SetTextFont(42)
        detail.SetTextColor(palette["gray"])
        detail.SetTextSize(0.030)
        detail.SetTextAlign(31)
        detail.DrawLatex(0.95, 0.815, right)
        keep.extend([sample, detail])
    save_canvas(canvas, output_dir / "full_detector_response")
    return keep


def _metric_lookup(metrics: list[dict[str, object]]) -> dict[tuple[float, float], dict[str, object]]:
    return {
        (float(row["pitch_um"]), float(row["recovery_ns"])): row
        for row in metrics
        if row["pitch_um"] != ""
    }


def _plot_performance(ROOT, palette, metrics, output_dir: Path):
    lookup = _metric_lookup(metrics)
    pde_only = next(
        row
        for row in metrics
        if row["scenario"] == "reference_pde_only_no_recovery"
    )
    analog = next(row for row in metrics if row["scenario"] == ANALOG_SIPM.name)
    canvas = make_canvas(ROOT, "c_performance", width=1400, height=650)
    canvas.Divide(2, 1, 0.015, 0.0)
    keep = []
    specs = [
        (
            "pion_C_mean_shift_vs_pde_only",
            "pion_C_mean_shift_error",
            "Pion C response shift after e calibration [%]",
            (-1.0, 25.0),
        ),
        (
            "pion_C_rms_over_mean",
            "pion_C_rms_over_mean_error",
            "Pion C response RMS / mean [%]",
            (10.0, 26.0),
        ),
    ]
    colors = [palette["green"], palette["blue"], palette["red"]]
    markers = [22, 20, 21]
    for pad_index, (metric_name, error_name, ylabel, yrange) in enumerate(specs, start=1):
        pad = canvas.cd(pad_index)
        style_pad(pad, left=0.18, right=0.045, top=0.135, bottom=0.16)
        frame = ROOT.TH1D(f"frame_{pad_index}", f";SPAD pitch [#mum];{ylabel}", 100, 0.0, 105.0)
        frame.SetMinimum(yrange[0])
        frame.SetMaximum(yrange[1])
        frame.Draw("AXIS")
        style_axes(frame, title_size=0.047, label_size=0.040)
        keep.append(frame)
        study_band = ROOT.TBox(10.0, yrange[0], 20.0, yrange[1])
        study_band.SetFillColorAlpha(palette["green"], 0.10)
        study_band.SetLineColor(0)
        study_band.Draw("SAME")
        keep.append(study_band)
        baseline_value = 100.0 * float(pde_only[metric_name])
        baseline_line = ROOT.TLine(0.0, baseline_value, 105.0, baseline_value)
        baseline_line.SetLineColor(palette["black"])
        baseline_line.SetLineWidth(3)
        baseline_line.Draw()
        keep.append(baseline_line)
        graphs = []
        for recovery, color, marker in zip(RECOVERY_NS, colors, markers, strict=True):
            y = np.asarray(
                [100.0 * float(lookup[(float(pitch), float(recovery))][metric_name]) for pitch in PITCHES_UM]
            )
            ey = np.asarray(
                [100.0 * float(lookup[(float(pitch), float(recovery))][error_name]) for pitch in PITCHES_UM]
            )
            graph = ROOT.TGraphErrors(
                len(PITCHES_UM),
                np.ascontiguousarray(PITCHES_UM, dtype=np.float64),
                np.ascontiguousarray(y, dtype=np.float64),
                np.zeros(len(PITCHES_UM), dtype=np.float64),
                np.ascontiguousarray(ey, dtype=np.float64),
            )
            graph.SetLineColor(color)
            graph.SetMarkerColor(color)
            graph.SetMarkerStyle(marker)
            graph.SetLineWidth(3)
            graph.Draw("LP SAME")
            graphs.append(graph)
            keep.append(graph)
        analog_graph = ROOT.TGraphErrors(
            1,
            np.asarray([ANALOG_SIPM.effective_pitch_um], dtype=np.float64),
            np.asarray([100.0 * float(analog[metric_name])], dtype=np.float64),
            np.zeros(1, dtype=np.float64),
            np.asarray([100.0 * float(analog[error_name])], dtype=np.float64),
        )
        analog_graph.SetMarkerColor(palette["purple"])
        analog_graph.SetLineColor(palette["purple"])
        analog_graph.SetMarkerStyle(29)
        analog_graph.SetMarkerSize(1.5)
        analog_graph.SetLineWidth(2)
        analog_graph.Draw("P SAME")
        keep.append(analog_graph)
        if pad_index == 1:
            legend = make_legend(ROOT, 0.18, 0.47, 0.68, 0.84, text_size=0.029)
            legend.AddEntry(study_band, "10-20 #mum study region", "f")
            legend.AddEntry(
                baseline_line,
                f"{100 * PDE:g}% PDE, no recovery loss",
                "l",
            )
            for graph, recovery in zip(graphs, RECOVERY_NS, strict=True):
                legend.AddEntry(graph, f"FLASH recovery = {recovery:g} ns", "lp")
            legend.AddEntry(analog_graph, "Analog C-Series proxy", "p")
            legend.Draw()
            keep.append(legend)
        keep.extend(draw_flash_label(ROOT, pad, extra="Simulation", right="90 GeV #pi^{+}, N=100"))
        if pad_index == 2:
            caution = ROOT.TLatex()
            caution.SetNDC(True)
            caution.SetTextFont(42)
            caution.SetTextColor(palette["gray"])
            caution.SetTextSize(0.028)
            caution.DrawLatex(0.31, 0.24, "Star: analog, 4761-cell grid, recharge 82 ns")
            caution.DrawLatex(0.45, 0.19, "Narrowing can indicate saturation")
            keep.append(caution)
    save_canvas(canvas, output_dir / "full_detector_performance_scan")
    return keep


def _reread_event(path: Path, event_index: int) -> FullDetectorPhotonEvent:
    for event in iter_full_detector_events(path, max_events=event_index + 1):
        if event.event_index == event_index:
            return event
    raise IndexError(event_index)


def _rod_histogram(
    ROOT,
    name: str,
    event: FullDetectorPhotonEvent,
    mask: np.ndarray,
    *,
    normalize: bool,
    weights: np.ndarray | None = None,
    range_mm: float | None = None,
    normalization_total: float | None = None,
    normalize_by_area: bool = False,
):
    if not np.any(mask):
        raise ValueError("Cannot draw an empty event map")
    ix = event.rod_ix[mask]
    iy = event.rod_iy[mask]
    if range_mm is None:
        low_x, high_x = int(np.min(ix)) - 1, int(np.max(ix)) + 2
        low_y, high_y = int(np.min(iy)) - 1, int(np.max(iy)) + 2
    else:
        if range_mm <= 0:
            raise ValueError("range_mm must be positive")
        half_bins = int(np.ceil(range_mm / 4.0))
        low_x, high_x = -half_bins, half_bins
        low_y, high_y = -half_bins, half_bins
    hist = ROOT.TH2D(
        name,
        ";x [mm];y [mm]",
        high_x - low_x,
        4.0 * low_x,
        4.0 * high_x,
        high_y - low_y,
        4.0 * low_y,
        4.0 * high_y,
    )
    hist.SetDirectory(0)
    if weights is None:
        selected_weights = np.ones(ix.size, dtype=np.float64)
    else:
        selected_weights = np.asarray(weights, dtype=np.float64)
        if selected_weights.shape == event.time_ns.shape:
            selected_weights = selected_weights[mask]
        elif selected_weights.shape != ix.shape:
            raise ValueError(
                "weights must match either the full photon arrays or the selected photons"
            )
    for x_index, y_index, weight in zip(ix, iy, selected_weights, strict=True):
        hist.Fill(4.0 * float(x_index), 4.0 * float(y_index), float(weight))
    if normalize and hist.Integral() > 0:
        denominator = hist.Integral() if normalization_total is None else normalization_total
        if denominator <= 0:
            raise ValueError("normalization_total must be positive")
        area_mm2 = 1.0
        if normalize_by_area:
            area_mm2 = hist.GetXaxis().GetBinWidth(1) * hist.GetYaxis().GetBinWidth(1)
        hist.Scale(100.0 / (denominator * area_mm2))
    return hist


def _rectangular_channel_histogram(
    ROOT,
    name: str,
    layout,
    channel_ids: np.ndarray,
    mask: np.ndarray,
    *,
    range_mm: float,
    normalization_total: float,
):
    """Draw incident-normalized photon fraction per channel footprint area."""

    valid = np.asarray(mask, dtype=bool) & (channel_ids >= 0)
    counts = np.bincount(
        channel_ids[valid],
        minlength=layout.n_channels,
    ).astype(np.float64)
    if counts.sum() <= 0:
        raise ValueError("Cannot draw an empty HG-DREAM channel map")
    if normalization_total <= 0:
        raise ValueError("normalization_total must be positive")
    hist = ROOT.TH2Poly(
        name,
        ";x [mm];y [mm]",
        -range_mm,
        range_mm,
        -range_mm,
        range_mm,
    )
    hist.SetDirectory(0)
    empty_boxes = []
    for channel_id in range(layout.n_channels):
        half_x = 0.5 * layout.width_x_mm[channel_id]
        half_y = 0.5 * layout.width_y_mm[channel_id]
        low_x = layout.center_x_mm[channel_id] - half_x
        low_y = layout.center_y_mm[channel_id] - half_y
        high_x = layout.center_x_mm[channel_id] + half_x
        high_y = layout.center_y_mm[channel_id] + half_y
        bin_id = hist.AddBin(
            low_x,
            low_y,
            high_x,
            high_y,
        )
        channel_area_mm2 = layout.width_x_mm[channel_id] * layout.width_y_mm[channel_id]
        density = 100.0 * counts[channel_id] / (normalization_total * channel_area_mm2)
        # TH2Poly paints exact zeros with the first palette color.  A negative
        # sentinel is outside the logarithmic color domain and leaves an empty
        # channel white while retaining its outlined polygon.
        content = density if counts[channel_id] > 0 else -1.0
        hist.SetBinContent(bin_id, content)
        if counts[channel_id] <= 0:
            box = ROOT.TBox(low_x, low_y, high_x, high_y)
            box.SetFillColor(ROOT.kWhite)
            box.SetFillStyle(1001)
            box.SetLineColor(ROOT.TColor.GetColor("#B8BDC4"))
            box.SetLineWidth(1)
            empty_boxes.append(box)
    return hist, empty_boxes


def _fiber_channel_histogram(
    ROOT,
    name: str,
    event: FullDetectorPhotonEvent,
    mask: np.ndarray,
    *,
    range_mm: float,
    normalization_total: float,
):
    """Draw ideal fiber truth with enlarged markers for slide readability."""

    selected = np.asarray(mask, dtype=bool)
    if not np.any(selected):
        raise ValueError("Cannot draw an empty ideal fiber map")
    if normalization_total <= 0:
        raise ValueError("normalization_total must be positive")
    sensor_ids = event.sensor_id[selected]
    unique_ids, first_indices, inverse = np.unique(
        sensor_ids,
        return_index=True,
        return_inverse=True,
    )
    counts = np.bincount(inverse, minlength=unique_ids.size).astype(np.float64)
    center_x = event.sensor_x_mm[selected][first_indices]
    center_y = event.sensor_y_mm[selected][first_indices]
    display_radius_mm = 0.5 * IDEAL_FIBER_DISPLAY_DIAMETER_MM
    # Four Cherenkov fibers sample each nominal 4 x 4 mm2 rod footprint.
    # Normalize each fiber channel by its assigned 4 mm2 calorimeter area.  The
    # displayed polygon is deliberately larger than the physical 0.8 mm fiber
    # so an individual fiber remains visible at the full-detector slide scale.
    channel_area_mm2 = (4.0 * 4.0) / 4.0
    angles = np.arange(16, dtype=np.float64) * (2.0 * np.pi / 16.0)

    hist = ROOT.TH2Poly(
        name,
        ";x [mm];y [mm]",
        -range_mm,
        range_mm,
        -range_mm,
        range_mm,
    )
    hist.SetDirectory(0)
    for x_mm, y_mm, count in zip(center_x, center_y, counts, strict=True):
        polygon_x = x_mm + display_radius_mm * np.cos(angles)
        polygon_y = y_mm + display_radius_mm * np.sin(angles)
        bin_id = hist.AddBin(16, polygon_x, polygon_y)
        density = 100.0 * count / (normalization_total * channel_area_mm2)
        hist.SetBinContent(bin_id, density)
    return hist


def _pmt_channel_histogram(
    ROOT,
    name: str,
    centers_mm: np.ndarray,
    channel_ids: np.ndarray,
    mask: np.ndarray,
    *,
    range_mm: float,
    normalization_total: float,
):
    """Draw incident-normalized photon fraction per original DREAM tower area."""

    valid = np.asarray(mask, dtype=bool) & (channel_ids >= 0)
    counts = np.bincount(
        channel_ids[valid],
        minlength=centers_mm.shape[0],
    ).astype(np.float64)
    if counts.sum() <= 0:
        raise ValueError("Cannot draw an empty original-DREAM PMT map")
    if normalization_total <= 0:
        raise ValueError("normalization_total must be positive")
    hist = ROOT.TH2Poly(
        name,
        ";x [mm];y [mm]",
        -range_mm,
        range_mm,
        -range_mm,
        range_mm,
    )
    hist.SetDirectory(0)
    tower_radius_mm = ORIGINAL_DREAM_TOWER_RADIUS_MM
    # Each original DREAM PMT sums about 270 nominal 4 x 4 mm2 rod cells.
    # Using the readout footprint rather than the drawn approximate hexagon
    # makes the normalization directly comparable with HG-DREAM and FLASH.
    tower_area_mm2 = 270.0 * 4.0 * 4.0
    angles = np.pi / 6.0 + np.arange(6, dtype=np.float64) * np.pi / 3.0
    for channel_id, (center_x, center_y) in enumerate(centers_mm):
        polygon_x = center_x + tower_radius_mm * np.cos(angles)
        polygon_y = center_y + tower_radius_mm * np.sin(angles)
        bin_id = hist.AddBin(6, polygon_x, polygon_y)
        density = 100.0 * counts[channel_id] / (normalization_total * tower_area_mm2)
        content = density if counts[channel_id] > 0 else -1.0
        hist.SetBinContent(bin_id, content)
    return hist


def _plot_event_display(ROOT, palette, event: FullDetectorPhotonEvent, output_dir: Path):
    scenario = FullDetectorScenario(PDE, TARGET_PITCH_UM, TARGET_RECOVERY_NS)
    digitized = digitize_full_detector_event(
        event,
        scenario,
        seed=WORKFLOW.reference.seed,
    )
    analog = digitize_analog_sipm_event(
        event,
        ANALOG_SIPM,
        seed=WORKFLOW.reference.seed,
    )
    all_photons = np.ones(event.time_ns.size, dtype=bool)
    incident_photons = float(event.selected_photon_count)
    plot_range_mm = 200.0
    accepted = _rod_histogram(
        ROOT,
        "h_event_accepted",
        event,
        digitized.accepted_mask,
        normalize=True,
        range_mm=plot_range_mm,
        normalization_total=incident_photons,
        normalize_by_area=True,
    )
    ideal_fibers = _fiber_channel_histogram(
        ROOT,
        "h_event_ideal_fibers",
        event,
        all_photons,
        range_mm=plot_range_mm,
        normalization_total=incident_photons,
    )
    hgdream_layout = hgdream_cherenkov_channel_layout()
    hgdream_ids = hgdream_cherenkov_channel_ids(event, hgdream_layout)
    hgdream_channels, hgdream_empty_boxes = _rectangular_channel_histogram(
        ROOT,
        "h_event_hgdream_channels",
        hgdream_layout,
        hgdream_ids,
        all_photons,
        range_mm=plot_range_mm,
        normalization_total=incident_photons,
    )
    pmt_centers = original_dream_pmt_centers()
    pmt_ids = original_dream_pmt_channel_ids(event, pmt_centers)
    pmt_channels = _pmt_channel_histogram(
        ROOT,
        "h_event_original_pmt_channels",
        pmt_centers,
        pmt_ids,
        all_photons,
        range_mm=plot_range_mm,
        normalization_total=incident_photons,
    )
    event_hists = (pmt_channels, hgdream_channels, accepted, ideal_fibers)
    for hist in event_hists:
        positive_contents = np.asarray(
            [
            hist.GetBinContent(bin_id)
            for bin_id in range(hist.GetNcells())
            if hist.GetBinContent(bin_id) > 0
            ],
            dtype=np.float64,
        )
        if positive_contents.size >= 20:
            # A few hot or barely occupied channels should not consume most of
            # the logarithmic palette.  Clipping only the 2% tails exposes the
            # spatial structure while preserving the underlying bin contents.
            low, high = np.quantile(positive_contents, [0.02, 0.98])
            panel_minimum = 0.50 * low
            panel_maximum = 1.20 * high
        else:
            panel_minimum = 0.80 * np.min(positive_contents)
            panel_maximum = 1.25 * np.max(positive_contents)
        panel_minimum = max(panel_minimum, panel_maximum * 1.0e-5)
        hist.GetZaxis().SetTitle("Photon fraction / area [% mm^{-2}]")
        hist.SetMinimum(panel_minimum)
        hist.SetMaximum(panel_maximum)

    # Retain the familiar dark viridis low end but replace its yellow high end
    # with a saturated red that remains legible on a white slide background.
    stops = np.asarray([0.00, 0.28, 0.55, 0.78, 1.00], dtype=np.float64)
    red = np.asarray([0.267, 0.231, 0.129, 0.369, 0.843], dtype=np.float64)
    green = np.asarray([0.004, 0.322, 0.569, 0.788, 0.188], dtype=np.float64)
    blue = np.asarray([0.329, 0.545, 0.549, 0.384, 0.153], dtype=np.float64)
    first_color = ROOT.TColor.CreateGradientColorTable(
        len(stops),
        stops,
        red,
        green,
        blue,
        255,
    )
    density_palette = np.arange(first_color, first_color + 255, dtype=np.int32)
    ROOT.gStyle.SetPalette(len(density_palette), density_palette)

    canvas = make_canvas(ROOT, "c_event_display", width=2200, height=650)
    canvas.Divide(4, 1, 0.006, 0.0)
    keep = [pmt_channels, hgdream_channels, accepted, ideal_fibers]
    for pad_index, (hist, title, subtitle) in enumerate(
        [
            (pmt_channels, "Original DREAM", "19 C PMTs, 270 rods / PMT"),
            (
                hgdream_channels,
                "HG-DREAM",
                "3 rods / core channel, 12 / outer",
            ),
            (
                accepted,
                "DREAM with FLASH",
                f"{100 * PDE:g}% PDE, 4 mm rod channels",
            ),
            (
                ideal_fibers,
                "Ideal photon truth",
                "0.8 mm fibers, 2.0 mm display markers",
            ),
        ],
        start=1,
    ):
        pad = canvas.cd(pad_index)
        style_pad(pad, left=0.15, right=0.20, top=0.19, bottom=0.16)
        pad.SetLogz(True)
        style_axes(hist, title_size=0.044, label_size=0.036)
        hist.GetZaxis().SetTitleSize(0.036)
        hist.GetZaxis().SetLabelSize(0.032)
        # Retain tower boundaries where they carry useful geometric meaning.
        # HG-DREAM has hundreds of small polygons whose black outlines obscure
        # their color, so draw that panel without polygon borders.
        draw_option = "COLZ L" if hist is pmt_channels else "COLZ"
        hist.Draw(draw_option)
        if hist is hgdream_channels:
            for box in hgdream_empty_boxes:
                box.Draw("SAME")
            keep.extend(hgdream_empty_boxes)
        keep.extend(draw_panel_title(ROOT, pad, title, subtitle=subtitle, y=0.965))
    save_canvas(canvas, output_dir / "full_detector_event_display")
    ROOT.gStyle.SetPalette(ROOT.kViridis)
    return keep, digitized, analog


def _plot_manifest() -> list[dict[str, str]]:
    return [
        {
            "figure": "full_detector_validation.pdf",
            "question": "Do the explicit optical photons track the stored detector-level Cherenkov signal?",
            "takeaway": "Yes, but not identically; correlations support comparative response studies and require an explicit caveat.",
            "talk_role": "validation / backup",
        },
        {
            "figure": "full_detector_occupancy.pdf",
            "question": "Where do conventional and digital SiPMs encounter dynamic-range stress in a complete shower?",
            "takeaway": "The stress is localized in the brightest fiber; analog and digital recovery can now be compared on the same photons.",
            "talk_role": "main",
        },
        {
            "figure": "full_detector_event_display.pdf",
            "question": (
                "How does the same complete shower appear with original DREAM, "
                "HG-DREAM, FLASH, and ideal photon-level readout?"
            ),
            "takeaway": (
                "Normalizing the photon fraction by each readout footprint removes "
                "the HG-DREAM core-to-outer channel-area bias. Robust independent "
                "logarithmic ranges and enlarged ideal-fiber display markers expose "
                "the spatial structure within each readout view."
            ),
            "talk_role": "main",
        },
        {
            "figure": "full_detector_response.pdf",
            "question": "How does the sensor model change the event-level Cherenkov response and an illustrative dual-readout estimator?",
            "takeaway": "Electron calibration removes the overall scale, but analog microcell recharge and digital holdoff produce different particle-dependent response transfer.",
            "talk_role": "main",
        },
        {
            "figure": "full_detector_performance_scan.pdf",
            "question": "How does the FLASH pitch/recovery surface compare with a representative conventional SiPM?",
            "takeaway": "The analog benchmark and coarse digital pixels can both narrow distributions by saturation while creating particle-dependent response shifts.",
            "talk_role": "main",
        },
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--electron", type=Path, default=DEFAULT_ELECTRON)
    parser.add_argument("--pion", type=Path, default=DEFAULT_PION)
    parser.add_argument("--output", type=Path, default=REPO_DIR / "results" / "full-detector-v1")
    parser.add_argument("--max-electron-events", type=int)
    parser.add_argument("--max-pion-events", type=int)
    parser.add_argument("--bootstrap-samples", type=int, default=500)
    parser.add_argument("--seed", type=int, default=WORKFLOW.reference.seed)
    parser.add_argument(
        "--reuse-event-tables",
        action="store_true",
        help="Regenerate plots from cached electron_events.csv and pion_events.csv",
    )
    args = parser.parse_args()
    output_dir = args.output.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.reuse_event_tables:
        _validate_cached_workflow(output_dir)
        electron_rows = _read_csv(output_dir / "electron_events.csv")
        pion_rows = _read_csv(output_dir / "pion_events.csv")
        electron_accepted = _accepted_cube(electron_rows)
        pion_accepted = _accepted_cube(pion_rows)
    else:
        electron_rows, electron_accepted = _read_and_scan(
            args.electron,
            max_events=args.max_electron_events,
            seed=args.seed,
        )
        pion_rows, pion_accepted = _read_and_scan(
            args.pion,
            max_events=args.max_pion_events,
            seed=args.seed,
        )
        _write_csv(output_dir / "electron_events.csv", electron_rows)
        _write_csv(output_dir / "pion_events.csv", pion_rows)

    metrics, arrays, calibration = _calibrated_responses(
        electron_rows,
        pion_rows,
        electron_accepted,
        pion_accepted,
        bootstrap_samples=args.bootstrap_samples,
    )
    _write_csv(output_dir / "full_detector_metrics.csv", metrics)
    manifest = _plot_manifest()
    _write_csv(output_dir / "plot_manifest.csv", manifest)

    pion_selected = _column(pion_rows, "selected_cherenkov_photons")
    representative_position = int(np.argmin(np.abs(pion_selected - np.median(pion_selected))))
    representative_index = int(pion_rows[representative_position]["event"])
    representative_event = _reread_event(args.pion, representative_index)

    import ROOT

    palette = set_flash_style(ROOT)
    ROOT.gStyle.SetPalette(ROOT.kViridis)
    keep = []
    keep.extend(_plot_validation(ROOT, palette, electron_rows, pion_rows, output_dir))
    keep.extend(_plot_occupancy(ROOT, palette, electron_rows, pion_rows, output_dir))
    keep.extend(
        _plot_response(
            ROOT,
            palette,
            arrays,
            float(calibration["pion_beam_energy_gev_mean"]),
            output_dir,
        )
    )
    keep.extend(_plot_performance(ROOT, palette, metrics, output_dir))
    event_keep, event_digitization, event_analog_digitization = _plot_event_display(
        ROOT,
        palette,
        representative_event,
        output_dir,
    )
    keep.extend(event_keep)
    event_hgdream_layout = hgdream_cherenkov_channel_layout()
    event_hgdream_ids = hgdream_cherenkov_channel_ids(
        representative_event,
        event_hgdream_layout,
    )
    event_dream_ids = original_dream_pmt_channel_ids(representative_event)
    rod_centroid_mm = np.asarray(
        [
            4.0 * np.mean(representative_event.rod_ix),
            4.0 * np.mean(representative_event.rod_iy),
        ],
        dtype=np.float64,
    )
    ideal_fiber_centroid_mm = np.asarray(
        [
            np.mean(representative_event.sensor_x_mm),
            np.mean(representative_event.sensor_y_mm),
        ],
        dtype=np.float64,
    )

    target_row = next(
        row
        for row in metrics
        if row["pitch_um"] == TARGET_PITCH_UM
        and row["recovery_ns"] == TARGET_RECOVERY_NS
    )
    current_row = metrics[0]
    ideal_row = next(row for row in metrics if row["scenario"] == "ideal_explicit_photon_C_plus_S")
    pde_only_row = next(
        row
        for row in metrics
        if row["scenario"] == "reference_pde_only_no_recovery"
    )
    analog_row = next(row for row in metrics if row["scenario"] == ANALOG_SIPM.name)
    summary = {
        "study": "HG-DREAM full-detector Cherenkov readout replacement",
        "sensor_workflow": {
            "path": str(WORKFLOW.source_path),
            "sha256": WORKFLOW.sha256,
            "schema_version": WORKFLOW.schema_version,
            "reference_scenario": WORKFLOW.reference.name,
        },
        "inputs": {
            "electron": str(args.electron.expanduser().resolve()),
            "electron_events": len(electron_rows),
            "pion": str(args.pion.expanduser().resolve()),
            "pion_events": len(pion_rows),
            "gate_ns": [0.0, 40.0],
        },
        "scan": {
            "pde": PDE,
            "pitches_um": PITCHES_UM.tolist(),
            "recovery_times_ns": RECOVERY_NS.tolist(),
            "analog_sipm": {
                "model": ANALOG_SIPM.name,
                "pde_at_420nm": ANALOG_SIPM.pde,
                "overvoltage_v": ANALOG_SIPM.overvoltage_v,
                "nominal_microcell_size_um": ANALOG_SIPM.microcell_size_um,
                "datasheet_microcells": ANALOG_SIPM.n_microcells,
                "modeled_microcells": ANALOG_SIPM.modeled_microcells,
                "effective_grid_pitch_um": ANALOG_SIPM.effective_pitch_um,
                "recovery_tau_ns": ANALOG_SIPM.recovery_tau_ns,
                "fast_output_fwhm_ns": ANALOG_SIPM.fast_output_fwhm_ns,
                "source_url": ANALOG_SIPM.source_url,
            },
        },
        "calibration": calibration,
        "event_display_geometry": {
            "coordinate_frame": "rod-centered detector x-y frame",
            "rod_truth_centroid_mm": rod_centroid_mm.tolist(),
            "ideal_fiber_centroid_mm": ideal_fiber_centroid_mm.tolist(),
            "ideal_minus_rod_centroid_mm": (
                ideal_fiber_centroid_mm - rod_centroid_mm
            ).tolist(),
            "hgdream_cherenkov_channels": event_hgdream_layout.n_channels,
            "hgdream_core_channels": int(np.count_nonzero(event_hgdream_layout.is_core)),
            "hgdream_outer_channels": int(np.count_nonzero(~event_hgdream_layout.is_core)),
            "hgdream_mapped_photons": int(np.count_nonzero(event_hgdream_ids >= 0)),
            "hgdream_occupied_channels": int(
                np.unique(event_hgdream_ids[event_hgdream_ids >= 0]).size
            ),
            "original_dream_cherenkov_pmts": 19,
            "original_dream_rods_per_tower": 270,
            "original_dream_tower_apex_to_apex_mm": 80.0,
            "ideal_cherenkov_fiber_diameter_mm": ORIGINAL_DREAM_FIBER_DIAMETER_MM,
            "ideal_fiber_display_diameter_mm": IDEAL_FIBER_DISPLAY_DIAMETER_MM,
            "original_dream_mapped_photons": int(np.count_nonzero(event_dream_ids >= 0)),
            "original_dream_occupied_pmts": int(
                np.unique(event_dream_ids[event_dream_ids >= 0]).size
            ),
            "color_scale": (
                "independent robust logarithmic ranges, with the 2 percent tails "
                "clipped, of incident-normalized photon fraction per assigned "
                "calorimeter channel area in percent per square millimeter; white "
                "means zero and colors are compared only within a panel"
            ),
        },
        "headline": {
            "current_DR_mean_bias_percent": 100.0 * float(current_row["pion_DR_mean_bias"]),
            "current_DR_rms_over_mean_percent": 100.0 * float(current_row["pion_DR_rms_over_mean"]),
            "ideal_photon_DR_mean_bias_percent": 100.0 * float(ideal_row["pion_DR_mean_bias"]),
            "ideal_photon_DR_rms_over_mean_percent": 100.0 * float(ideal_row["pion_DR_rms_over_mean"]),
            "reference_pde_only_C_rms_over_mean_percent": 100.0 * float(
                pde_only_row["pion_C_rms_over_mean"]
            ),
            "analog_C_mean_shift_vs_pde_only_percent": 100.0 * float(analog_row["pion_C_mean_shift_vs_pde_only"]),
            "analog_C_rms_over_mean_percent": 100.0 * float(analog_row["pion_C_rms_over_mean"]),
            "analog_DR_mean_bias_percent": 100.0 * float(analog_row["pion_DR_mean_bias"]),
            "analog_DR_rms_over_mean_percent": 100.0 * float(analog_row["pion_DR_rms_over_mean"]),
            "analog_electron_charge_survival_percent": 100.0 * float(analog_row["electron_recovery_survival"]),
            "analog_pion_charge_survival_percent": 100.0 * float(analog_row["pion_recovery_survival"]),
            "FLASH_reference_DR_mean_bias_percent": 100.0 * float(
                target_row["pion_DR_mean_bias"]
            ),
            "FLASH_reference_DR_rms_over_mean_percent": 100.0 * float(
                target_row["pion_DR_rms_over_mean"]
            ),
            "FLASH_reference_C_mean_shift_vs_pde_only_percent": 100.0 * float(
                target_row["pion_C_mean_shift_vs_pde_only"]
            ),
            "FLASH_reference_C_rms_over_mean_percent": 100.0 * float(
                target_row["pion_C_rms_over_mean"]
            ),
            "FLASH_reference_electron_recovery_survival_percent": 100.0 * float(
                target_row["electron_recovery_survival"]
            ),
            "FLASH_reference_pion_recovery_survival_percent": 100.0 * float(
                target_row["pion_recovery_survival"]
            ),
            "representative_pion_event": representative_index,
            "representative_event_selected_photons": representative_event.selected_photon_count,
            "representative_event_detected_photons": event_digitization.n_detected,
            "representative_event_accepted_photons": event_digitization.n_accepted,
            "representative_event_recovery_survival": event_digitization.recovery_survival,
            "representative_event_analog_detected_photons": event_analog_digitization.n_detected,
            "representative_event_analog_charge_pe": event_analog_digitization.charge_pe,
            "representative_event_analog_charge_survival": event_analog_digitization.charge_survival,
        },
        "scope_limits": [
            "Only Cherenkov optical photons are explicit in these ROOT files; sum3dSS is kept unchanged.",
            "The electron and pion samples are single-energy samples and were produced with different run configurations.",
            "The dSiPM scan uses a chi derived from the internally consistent explicit-photon baseline and holds it fixed for every sensor scenario.",
            "The conventional-SiPM comparison is an Onsemi MicroFC-30035 proxy using 31% peak PDE, 4,774 nominal microcells, an 82 ns recharge constant, and an optimistic partial-charge model.",
            "The analog proxy does not yet include wavelength-dependent PDE, recharge-dependent trigger probability, crosstalk, afterpulsing, electronics noise, or waveform integration.",
            "The current stored S+C response is shown as context with its own chi; it is not treated as the same photon observable.",
            "The model currently excludes dark counts, optical cross-talk, afterpulsing, channel masking, electronics noise, and link-level throughput.",
            "The explicit photon branch is correlated with, but not identical to, the stored sum3dCC signal.",
        ],
    }
    with (output_dir / "full_detector_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)
    print(json.dumps(summary["headline"], indent=2))


if __name__ == "__main__":
    main()
