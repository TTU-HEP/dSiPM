#!/usr/bin/env python3
"""Run the physically scaled FLASH/HG-DREAM sensor-requirement study.

The 2025 slides are reproduced elsewhere as a provenance/closure product.  In
that historical code, a 3 mm sensor segmentation was drawn across a 1 mm-wide
fiber-local coordinate window.  The resulting cell width in DREAMSim
coordinates was the labelled pitch divided by three.  This study removes that
display-coordinate scaling: photons remain in physical millimeters and are
assigned to a centered 3 x 3 mm sensor with the requested physical SPAD pitch.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import replace
import json
from pathlib import Path
import time

import numpy as np

from dsipm.characterization import PDERequirementPoint, scan_pde_pitch_recovery
from dsipm.config import load_device_catalog
from dsipm.digitizer import (
    _apply_nonparalyzable_dead_time,
    _cell_indices,
    _quantize_time,
)
from dsipm.requirements import accepted_counts_by_recovery, scan_pitch_and_recovery
from dsipm.root_io import inspect_root_file, iter_root_events, transform_coordinates


REPO_DIR = Path(__file__).resolve().parent
SAMPLES = {
    "100 GeV positron": REPO_DIR / "mc_dreamsim_e+_100_101.root",
    "1 TeV pion": REPO_DIR / "mc_dreamsim_pi+_1000_1001.root",
}
SENSOR_MIN_MM = -1.5
SENSOR_MAX_MM = 1.5
PITCHES_UM = (5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0)
PDES = (0.10, 0.20, 0.30, 0.50, 1.00)
RECOVERY_NS = (1.0, 4.3, 40.0)
EVENT_PITCHES_UM = (10.0, 20.0, 25.0)
EVENT_PDE = 0.20
EVENT_RECOVERY_NS = (1.0, 4.3, 40.0)
TIMING_SCENARIOS = (
    ("ideal", 0.0, None),
    ("10 ps response + 10 ps TDC", 10.0, 10.0),
    ("50 ps response + 50 ps TDC", 50.0, 50.0),
    ("100 ps response + 50 ps TDC", 100.0, 50.0),
    ("250 ps response + 50 ps TDC", 250.0, 50.0),
    ("1145 ps response + 37 ps TDC", 1145.0, 37.0),
)


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Cannot write empty table: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _point_row(sample: str, point: PDERequirementPoint) -> dict[str, object]:
    return {
        "sample": sample,
        "events": point.events,
        "pde": point.pde,
        "pitch_um": point.pitch_um,
        "recovery_ns": point.recovery_ns,
        "incident_photons": point.incident_photons,
        "detected_photons": point.detected_photons,
        "accepted_photons": point.accepted_photons,
        "usable_fraction": point.usable_fraction,
        "usable_fraction_low68": point.usable_fraction_low68,
        "usable_fraction_high68": point.usable_fraction_high68,
        "recovery_survival": point.recovery_retention,
        "recovery_survival_low68": point.recovery_retention_low68,
        "recovery_survival_high68": point.recovery_retention_high68,
    }


def _root_style():
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTitleFont(42, "XYZ")
    ROOT.gStyle.SetLabelFont(42, "XYZ")
    return ROOT


def _audit_inputs(*, max_events: int | None) -> dict[str, object]:
    """Record an explicit sequential cutflow and the selected coordinate extent."""

    import ROOT

    audit: dict[str, object] = {}
    for sample, path in SAMPLES.items():
        root_file = ROOT.TFile.Open(str(path), "READ")
        tree = root_file.Get("tree")
        cutflow = {
            "raw_photons": 0,
            "valid_production_fiber": 0,
            "in_cherenkov_fiber_core": 0,
            "exits_positive_z": 0,
            "inside_0_40ns_gate": 0,
            "production_fiber_0": 0,
            "inside_centered_3mm_sensor": 0,
        }
        x_chunks = []
        y_chunks = []
        event_counts = []
        processed = 0
        beam = None
        for event_index, event in enumerate(tree):
            if max_events is not None and event_index >= max_events:
                break
            if beam is None:
                beam = {
                    "beam_type": bytes(event.beamType).decode(errors="replace"),
                    "beam_id": int(event.beamID),
                    "beam_min_energy_GeV": float(event.beamMinE),
                    "beam_max_energy_GeV": float(event.beamMaxE),
                }
            time_ns = np.asarray(event.OP_time_final)
            fiber = np.asarray(event.OP_productionFiber).astype(np.int64)
            core = np.asarray(event.OP_isCoreC).astype(bool)
            final_z = np.asarray(event.OP_pos_final_z)
            valid = (fiber >= 0) & (fiber < 4)
            after_core = valid & core
            after_exit = after_core & (final_z > 0)
            after_gate = after_exit & (time_ns > 0.0) & (time_ns < 40.0)
            fiber_zero = after_gate & (fiber == 0)
            x_mm, y_mm = transform_coordinates(
                np.asarray(event.OP_pos_final_x)[fiber_zero],
                np.asarray(event.OP_pos_final_y)[fiber_zero],
                fiber[fiber_zero],
                "old_2025",
            )
            in_sensor = (
                (x_mm >= SENSOR_MIN_MM)
                & (x_mm < SENSOR_MAX_MM)
                & (y_mm >= SENSOR_MIN_MM)
                & (y_mm < SENSOR_MAX_MM)
            )
            cutflow["raw_photons"] += int(time_ns.size)
            cutflow["valid_production_fiber"] += int(np.count_nonzero(valid))
            cutflow["in_cherenkov_fiber_core"] += int(
                np.count_nonzero(after_core)
            )
            cutflow["exits_positive_z"] += int(np.count_nonzero(after_exit))
            cutflow["inside_0_40ns_gate"] += int(np.count_nonzero(after_gate))
            cutflow["production_fiber_0"] += int(np.count_nonzero(fiber_zero))
            cutflow["inside_centered_3mm_sensor"] += int(
                np.count_nonzero(in_sensor)
            )
            x_chunks.append(x_mm[in_sensor])
            y_chunks.append(y_mm[in_sensor])
            event_counts.append(int(np.count_nonzero(in_sensor)))
            processed += 1
        root_file.Close()
        x_all = np.concatenate(x_chunks) if x_chunks else np.empty(0)
        y_all = np.concatenate(y_chunks) if y_chunks else np.empty(0)
        counts = np.asarray(event_counts, dtype=np.float64)
        audit[sample] = {
            "source": inspect_root_file(path),
            "processed_events": processed,
            "beam": beam,
            "sequential_cutflow": cutflow,
            "selected_x_mm_range": [float(np.min(x_all)), float(np.max(x_all))],
            "selected_y_mm_range": [float(np.min(y_all)), float(np.max(y_all))],
            "selected_photons_per_event_quantiles": {
                label: float(value)
                for label, value in zip(
                    ("min", "q10", "median", "q90", "q99", "max"),
                    np.quantile(counts, (0.0, 0.10, 0.50, 0.90, 0.99, 1.0)),
                    strict=True,
                )
            },
        }
    return audit


def _graph(ROOT, x, y, *, color, marker=20, style=1, width=3):
    x_values = np.ascontiguousarray(x, dtype=np.float64)
    y_values = np.ascontiguousarray(y, dtype=np.float64)
    graph = ROOT.TGraph(len(x_values), x_values, y_values)
    graph.SetLineColor(color)
    graph.SetMarkerColor(color)
    graph.SetMarkerStyle(marker)
    graph.SetMarkerSize(1.0)
    graph.SetLineStyle(style)
    graph.SetLineWidth(width)
    return graph


def _save_canvas(canvas, stem: Path) -> list[Path]:
    outputs = [stem.with_suffix(".png"), stem.with_suffix(".pdf")]
    for output in outputs:
        canvas.SaveAs(str(output))
    return outputs


def _scan_physical(
    *, max_events: int | None, bootstrap_samples: int
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for sample_index, (sample, path) in enumerate(SAMPLES.items()):
        events = iter_root_events(
            path,
            max_events=max_events,
            geometry_mode="old_2025",
            gate_ns=(0.0, 40.0),
            detector_min_mm=SENSOR_MIN_MM,
            detector_max_mm=SENSOR_MAX_MM,
            fiber_index=0,
        )
        points = scan_pde_pitch_recovery(
            events,
            pdes=PDES,
            pitches_um=PITCHES_UM,
            recovery_times_ns=RECOVERY_NS,
            seed=4307,
            coordinate_scale=1.0,
            detector_min_mm=SENSOR_MIN_MM,
            detector_max_mm=SENSOR_MAX_MM,
            bootstrap_samples=bootstrap_samples,
            bootstrap_seed=20261010 + sample_index,
        )
        rows.extend(_point_row(sample, point) for point in points)
    return rows


def _scan_mapping_comparison(
    physical_rows: list[dict[str, object]], *, max_events: int | None
) -> list[dict[str, object]]:
    physical = {
        (
            str(row["sample"]),
            float(row["pitch_um"]),
            float(row["recovery_ns"]),
        ): row
        for row in physical_rows
        if np.isclose(float(row["pde"]), 1.0)
    }
    rows: list[dict[str, object]] = []
    for sample, path in SAMPLES.items():
        events = iter_root_events(
            path,
            max_events=max_events,
            geometry_mode="old_2025",
            gate_ns=(0.0, 40.0),
            detector_min_mm=SENSOR_MIN_MM,
            detector_max_mm=SENSOR_MAX_MM,
            fiber_index=0,
        )
        legacy_points = scan_pitch_and_recovery(
            events,
            pitches_um=PITCHES_UM,
            recovery_times_ns=RECOVERY_NS,
            coordinate_scale=3.0,
            detector_min_mm=SENSOR_MIN_MM,
            detector_max_mm=SENSOR_MAX_MM,
        )
        for legacy in legacy_points:
            key = (sample, legacy.pitch_um, legacy.recovery_ns)
            corrected = physical[key]
            rows.append(
                {
                    "sample": sample,
                    "pitch_um": legacy.pitch_um,
                    "recovery_ns": legacy.recovery_ns,
                    "historical_display_scaled_survival": legacy.retained_fraction,
                    "physical_mm_survival": corrected["recovery_survival"],
                    "absolute_difference": float(
                        legacy.retained_fraction
                        - float(corrected["recovery_survival"])
                    ),
                    "historical_effective_pitch_um": legacy.pitch_um / 3.0,
                }
            )
    return rows


def _event_response_rows(*, max_events: int | None) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for sample, path in SAMPLES.items():
        events = iter_root_events(
            path,
            max_events=max_events,
            geometry_mode="old_2025",
            gate_ns=(0.0, 40.0),
            detector_min_mm=SENSOR_MIN_MM,
            detector_max_mm=SENSOR_MAX_MM,
            fiber_index=0,
        )
        for event in events:
            event = event.sorted_by_time()
            incident = event.selected_photon_count
            rng = np.random.default_rng(
                np.random.SeedSequence([4307, event.event_index, 3409])
            )
            detected_mask = rng.random(incident) < EVENT_PDE
            detected = int(np.count_nonzero(detected_mask))
            for pitch_um in EVENT_PITCHES_UM:
                if detected:
                    _, _, cell_ids, _ = _cell_indices(
                        event.x_mm[detected_mask],
                        event.y_mm[detected_mask],
                        pitch_um / 1000.0,
                        SENSOR_MIN_MM,
                        SENSOR_MAX_MM,
                    )
                    accepted = accepted_counts_by_recovery(
                        cell_ids,
                        event.time_ns[detected_mask],
                        EVENT_RECOVERY_NS,
                    )
                else:
                    accepted = np.zeros(len(EVENT_RECOVERY_NS), dtype=np.int64)
                for recovery_ns, accepted_count in zip(
                    EVENT_RECOVERY_NS, accepted, strict=True
                ):
                    rows.append(
                        {
                            "sample": sample,
                            "event": event.event_index,
                            "pde": EVENT_PDE,
                            "pitch_um": pitch_um,
                            "recovery_ns": recovery_ns,
                            "incident_photons": incident,
                            "detected_photons": detected,
                            "accepted_photons": int(accepted_count),
                            "recovery_survival": (
                                float(accepted_count / detected) if detected else ""
                            ),
                        }
                    )
    return rows


def _physical_pulse_rows(*, max_events: int | None) -> list[dict[str, object]]:
    """Aggregate the pion arrival profile before and after physical-cell recovery."""

    edges = np.linspace(0.0, 40.0, 801)
    detected_hist = np.zeros(edges.size - 1, dtype=np.int64)
    accepted_hist = {
        pitch_um: np.zeros(edges.size - 1, dtype=np.int64)
        for pitch_um in EVENT_PITCHES_UM
    }
    events = iter_root_events(
        SAMPLES["1 TeV pion"],
        max_events=max_events,
        geometry_mode="old_2025",
        gate_ns=(0.0, 40.0),
        detector_min_mm=SENSOR_MIN_MM,
        detector_max_mm=SENSOR_MAX_MM,
        fiber_index=0,
    )
    for event in events:
        event = event.sorted_by_time()
        rng = np.random.default_rng(
            np.random.SeedSequence([4307, event.event_index, 3409])
        )
        detected_mask = rng.random(event.selected_photon_count) < EVENT_PDE
        detected_time = event.time_ns[detected_mask]
        detected_hist += np.histogram(detected_time, bins=edges)[0]
        for pitch_um in EVENT_PITCHES_UM:
            _, _, cell_ids, _ = _cell_indices(
                event.x_mm[detected_mask],
                event.y_mm[detected_mask],
                pitch_um / 1000.0,
                SENSOR_MIN_MM,
                SENSOR_MAX_MM,
            )
            accepted_mask = _apply_nonparalyzable_dead_time(
                cell_ids,
                detected_time,
                4.3,
            )
            accepted_hist[pitch_um] += np.histogram(
                detected_time[accepted_mask], bins=edges
            )[0]

    rows = []
    for bin_index, (low, high, detected) in enumerate(
        zip(edges[:-1], edges[1:], detected_hist, strict=True)
    ):
        for pitch_um in EVENT_PITCHES_UM:
            accepted = int(accepted_hist[pitch_um][bin_index])
            rows.append(
                {
                    "time_low_ns": float(low),
                    "time_high_ns": float(high),
                    "time_center_ns": float(0.5 * (low + high)),
                    "pde": EVENT_PDE,
                    "pitch_um": pitch_um,
                    "recovery_ns": 4.3,
                    "detected_photons": int(detected),
                    "accepted_photons": accepted,
                    "recovery_survival": (
                        float(accepted / detected) if detected else ""
                    ),
                }
            )
    return rows


def _bootstrap_rmse_interval(
    residuals: np.ndarray, *, seed: int, replicas: int = 1000
) -> tuple[float, float]:
    if residuals.size < 2:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, residuals.size, size=(replicas, residuals.size))
    rms = np.sqrt(np.mean(residuals[indices] ** 2, axis=1))
    low, high = np.quantile(rms, [0.16, 0.84])
    return float(low), float(high)


def _bootstrap_rmse_delta_interval(
    residuals: np.ndarray,
    baseline_residuals: np.ndarray,
    *,
    seed: int,
    replicas: int = 1000,
) -> tuple[float, float]:
    """Return a paired interval for the RMSE change relative to a baseline.

    Absolute RMSE is sensitive to a small number of shower-tail events in this
    proof-of-concept observable.  Pairing the same resampled events between a
    timing scenario and the ideal baseline isolates the incremental detector
    timing contribution without pretending that the proxy itself is precise.
    """

    if residuals.shape != baseline_residuals.shape or residuals.size < 2:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, residuals.size, size=(replicas, residuals.size))
    scenario_rmse = np.sqrt(np.mean(residuals[indices] ** 2, axis=1))
    baseline_rmse = np.sqrt(np.mean(baseline_residuals[indices] ** 2, axis=1))
    low, high = np.quantile(scenario_rmse - baseline_rmse, [0.16, 0.84])
    return float(low), float(high)


def _timing_depth_rows(*, max_events: int | None) -> list[dict[str, object]]:
    """Calibrate a simple pion longitudinal-depth proxy from the 10% time quantile."""

    observations = {name: [] for name, _, _ in TIMING_SCENARIOS}
    truth_depth = []
    event_ids = []
    path = SAMPLES["1 TeV pion"]
    events = iter_root_events(
        path,
        max_events=max_events,
        geometry_mode="old_2025",
        gate_ns=(0.0, 40.0),
        detector_min_mm=SENSOR_MIN_MM,
        detector_max_mm=SENSOR_MAX_MM,
        fiber_index=0,
    )
    for event in events:
        event = event.sorted_by_time()
        depth_mask = (
            (event.production_z_mm >= 0.0) & (event.production_z_mm <= 2000.0)
        )
        if np.count_nonzero(depth_mask) < 5:
            continue
        n_incident = event.selected_photon_count
        rng = np.random.default_rng(
            np.random.SeedSequence([4307, event.event_index, 3409])
        )
        detection_uniform = rng.random(n_incident)
        standard_normal = rng.standard_normal(n_incident)
        detected_mask = detection_uniform < 0.20
        if np.count_nonzero(detected_mask) < 5:
            continue
        _, _, cell_ids, _ = _cell_indices(
            event.x_mm[detected_mask],
            event.y_mm[detected_mask],
            20.0 / 1000.0,
            SENSOR_MIN_MM,
            SENSOR_MAX_MM,
        )
        detected_time = event.time_ns[detected_mask]
        accepted_mask = _apply_nonparalyzable_dead_time(
            cell_ids, detected_time, 4.3
        )
        if np.count_nonzero(accepted_mask) < 5:
            continue
        accepted_time = detected_time[accepted_mask]
        accepted_jitter = standard_normal[detected_mask][accepted_mask]
        for name, response_fwhm_ps, tdc_bin_ps in TIMING_SCENARIOS:
            if response_fwhm_ps > 0:
                sigma_ns = response_fwhm_ps / 2354.820045
                observed = accepted_time + sigma_ns * accepted_jitter
            else:
                observed = accepted_time.copy()
            observed = _quantize_time(observed, tdc_bin_ps)
            observations[name].append(float(np.quantile(observed, 0.10)))
        truth_depth.append(float(np.mean(event.production_z_mm[depth_mask])))
        event_ids.append(event.event_index)

    truth = np.asarray(truth_depth, dtype=np.float64)
    ids = np.asarray(event_ids, dtype=np.int64)
    train = ids % 2 == 0
    test = ~train
    fits = []
    for scenario_index, (name, response_fwhm_ps, tdc_bin_ps) in enumerate(
        TIMING_SCENARIOS
    ):
        observable = np.asarray(observations[name], dtype=np.float64)
        slope, intercept = np.polyfit(observable[train], truth[train], 1)
        prediction = intercept + slope * observable[test]
        residual = prediction - truth[test]
        rmse = float(np.sqrt(np.mean(residual**2)))
        low, high = _bootstrap_rmse_interval(
            residual,
            seed=20261010 + scenario_index,
        )
        fits.append(
            (
                residual,
                {
                    "scenario": name,
                    "response_fwhm_ps": response_fwhm_ps,
                    "tdc_bin_ps": "" if tdc_bin_ps is None else tdc_bin_ps,
                    "events_total": int(truth.size),
                    "events_train": int(np.count_nonzero(train)),
                    "events_test": int(np.count_nonzero(test)),
                    "calibration_slope_mm_per_ns": float(slope),
                    "calibration_intercept_mm": float(intercept),
                    "test_bias_mm": float(np.mean(residual)),
                    "test_rmse_mm": rmse,
                    "test_rmse_low68_mm": low,
                    "test_rmse_high68_mm": high,
                    "test_median_absolute_error_mm": float(
                        np.median(np.abs(residual))
                    ),
                    "test_absolute_error_q68_mm": float(
                        np.quantile(np.abs(residual), 0.68)
                    ),
                    "test_absolute_error_q90_mm": float(
                        np.quantile(np.abs(residual), 0.90)
                    ),
                    "test_correlation_time_depth": float(
                        np.corrcoef(observable[test], truth[test])[0, 1]
                    ),
                    "observable": "10th percentile of accepted photon arrival times",
                    "truth_target": "mean production depth of incident selected photons",
                },
            )
        )
    baseline_residual = fits[0][0]
    baseline_rmse = float(fits[0][1]["test_rmse_mm"])
    rows = []
    for scenario_index, (residual, row) in enumerate(fits):
        delta = float(row["test_rmse_mm"]) - baseline_rmse
        delta_low, delta_high = _bootstrap_rmse_delta_interval(
            residual,
            baseline_residual,
            seed=20262010 + scenario_index,
        )
        row["test_rmse_delta_vs_ideal_mm"] = delta
        row["test_rmse_delta_low68_mm"] = delta_low
        row["test_rmse_delta_high68_mm"] = delta_high
        rows.append(row)
    return rows


def _draw_mapping_correction(
    output_dir: Path, rows: list[dict[str, object]]
) -> list[Path]:
    ROOT = _root_style()
    canvas = ROOT.TCanvas("mapping_correction", "Mapping correction", 1600, 720)
    canvas.Divide(2, 1, 0.006, 0.006)
    retained = []
    for pad_index, sample in enumerate(SAMPLES, 1):
        canvas.cd(pad_index)
        ROOT.gPad.SetLeftMargin(0.13)
        ROOT.gPad.SetRightMargin(0.04)
        ROOT.gPad.SetTopMargin(0.13)
        ROOT.gPad.SetBottomMargin(0.13)
        ROOT.gPad.SetGridy()
        frame = ROOT.TH1D(
            f"mapping_frame_{pad_index}",
            ";Labelled SPAD pitch [#mum];Accepted / detected photons",
            100,
            0.0,
            105.0,
        )
        frame.SetDirectory(0)
        frame.SetMinimum(0.0)
        frame.SetMaximum(1.03)
        frame.Draw("AXIS")
        retained.append(frame)
        selected = sorted(
            (
                row
                for row in rows
                if row["sample"] == sample
                and np.isclose(float(row["recovery_ns"]), 4.3)
            ),
            key=lambda row: float(row["pitch_um"]),
        )
        x = [float(row["pitch_um"]) for row in selected]
        physical = [float(row["physical_mm_survival"]) for row in selected]
        legacy = [
            float(row["historical_display_scaled_survival"]) for row in selected
        ]
        physical_graph = _graph(
            ROOT, x, physical, color=ROOT.TColor.GetColor("#C83E4D"), marker=20
        )
        legacy_graph = _graph(
            ROOT,
            x,
            legacy,
            color=ROOT.TColor.GetColor("#2B6CB0"),
            marker=24,
            style=2,
        )
        physical_graph.Draw("LP SAME")
        legacy_graph.Draw("LP SAME")
        retained.extend([physical_graph, legacy_graph])
        legend = ROOT.TLegend(0.16, 0.17, 0.76, 0.34)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.AddEntry(physical_graph, "Physical millimeters (corrected)", "lp")
        legend.AddEntry(legacy_graph, "2025 display-scaled reproduction", "lp")
        legend.Draw()
        retained.append(legend)
        title = ROOT.TLatex()
        title.SetNDC(True)
        title.SetTextFont(62)
        title.SetTextSize(0.045)
        title.DrawLatex(0.13, 0.94, sample)
        subtitle = ROOT.TLatex()
        subtitle.SetNDC(True)
        subtitle.SetTextFont(42)
        subtitle.SetTextSize(0.030)
        subtitle.DrawLatex(0.13, 0.895, "Unity PDE, 4.3 ns recovery, one fiber")
        retained.extend([title, subtitle])
    return _save_canvas(canvas, output_dir / "legacy_vs_physical_mapping")


def _select_scan(
    rows: list[dict[str, object]],
    *,
    sample: str,
    pde: float | None = None,
    pitch_um: float | None = None,
    recovery_ns: float | None = None,
) -> list[dict[str, object]]:
    selected = []
    for row in rows:
        if row["sample"] != sample:
            continue
        if pde is not None and not np.isclose(float(row["pde"]), pde):
            continue
        if pitch_um is not None and not np.isclose(
            float(row["pitch_um"]), pitch_um
        ):
            continue
        if recovery_ns is not None and not np.isclose(
            float(row["recovery_ns"]), recovery_ns
        ):
            continue
        selected.append(row)
    return selected


def _draw_requirement_surface(
    output_dir: Path, rows: list[dict[str, object]]
) -> list[Path]:
    ROOT = _root_style()
    canvas = ROOT.TCanvas("requirements", "Physical sensor requirements", 1600, 720)
    canvas.Divide(2, 1, 0.006, 0.006)
    retained = []
    colors = {
        1.0: ROOT.TColor.GetColor("#159D66"),
        4.3: ROOT.TColor.GetColor("#7A5195"),
        40.0: ROOT.TColor.GetColor("#C83E4D"),
    }
    canvas.cd(1)
    ROOT.gPad.SetLeftMargin(0.13)
    ROOT.gPad.SetRightMargin(0.04)
    ROOT.gPad.SetTopMargin(0.13)
    ROOT.gPad.SetBottomMargin(0.13)
    ROOT.gPad.SetGridy()
    frame = ROOT.TH1D(
        "requirement_pitch_frame",
        ";Physical SPAD pitch [#mum];Accepted / detected photons",
        100,
        0.0,
        105.0,
    )
    frame.SetDirectory(0)
    frame.SetMinimum(0.45)
    frame.SetMaximum(1.02)
    frame.Draw("AXIS")
    retained.append(frame)
    legend = ROOT.TLegend(0.54, 0.66, 0.94, 0.87)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.035)
    for recovery_ns in RECOVERY_NS:
        selected = sorted(
            _select_scan(
                rows,
                sample="1 TeV pion",
                pde=0.20,
                recovery_ns=recovery_ns,
            ),
            key=lambda row: float(row["pitch_um"]),
        )
        graph = _graph(
            ROOT,
            [float(row["pitch_um"]) for row in selected],
            [float(row["recovery_survival"]) for row in selected],
            color=colors[recovery_ns],
            marker=20 if recovery_ns != 40.0 else 24,
        )
        graph.Draw("LP SAME")
        legend.AddEntry(
            graph,
            "40 ns (one hit/event)" if recovery_ns == 40.0 else f"{recovery_ns:g} ns",
            "lp",
        )
        retained.append(graph)
    legend.Draw()
    retained.append(legend)
    for target in (0.90, 0.95):
        line = ROOT.TLine(0.0, target, 105.0, target)
        line.SetLineStyle(2 if target == 0.90 else 3)
        line.SetLineColor(ROOT.kGray + 1)
        line.Draw("SAME")
        retained.append(line)
    title = ROOT.TLatex()
    title.SetNDC(True)
    title.SetTextFont(62)
    title.SetTextSize(0.040)
    title.DrawLatex(0.13, 0.94, "PDE = 20%: pitch #times recovery")
    retained.append(title)

    canvas.cd(2)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.04)
    ROOT.gPad.SetTopMargin(0.13)
    ROOT.gPad.SetBottomMargin(0.13)
    ROOT.gPad.SetGridy()
    frame2 = ROOT.TH1D(
        "requirement_pde_frame",
        ";Photon-detection efficiency;Accepted / detected photons",
        100,
        0.05,
        1.05,
    )
    frame2.SetDirectory(0)
    frame2.SetMinimum(0.45)
    frame2.SetMaximum(1.02)
    frame2.Draw("AXIS")
    retained.append(frame2)
    pitch_colors = {
        10.0: ROOT.TColor.GetColor("#2B6CB0"),
        20.0: ROOT.TColor.GetColor("#159D66"),
        25.0: ROOT.TColor.GetColor("#C83E4D"),
    }
    legend2 = ROOT.TLegend(0.59, 0.66, 0.94, 0.87)
    legend2.SetBorderSize(0)
    legend2.SetFillStyle(0)
    legend2.SetTextSize(0.035)
    for pitch_um in EVENT_PITCHES_UM:
        selected = sorted(
            _select_scan(
                rows,
                sample="1 TeV pion",
                pitch_um=pitch_um,
                recovery_ns=4.3,
            ),
            key=lambda row: float(row["pde"]),
        )
        graph = _graph(
            ROOT,
            [float(row["pde"]) for row in selected],
            [float(row["recovery_survival"]) for row in selected],
            color=pitch_colors[pitch_um],
            marker=20,
        )
        graph.Draw("LP SAME")
        legend2.AddEntry(graph, f"{pitch_um:g} #mum pitch", "lp")
        retained.append(graph)
    legend2.Draw()
    retained.append(legend2)
    for target in (0.90, 0.95):
        line = ROOT.TLine(0.05, target, 1.05, target)
        line.SetLineStyle(2 if target == 0.90 else 3)
        line.SetLineColor(ROOT.kGray + 1)
        line.Draw("SAME")
        retained.append(line)
    title2 = ROOT.TLatex()
    title2.SetNDC(True)
    title2.SetTextFont(62)
    title2.SetTextSize(0.040)
    title2.DrawLatex(0.12, 0.94, "4.3 ns: PDE #times pitch")
    retained.append(title2)
    return _save_canvas(canvas, output_dir / "physical_sensor_requirements")


def _draw_event_tail(
    output_dir: Path, rows: list[dict[str, object]]
) -> list[Path]:
    ROOT = _root_style()
    canvas = ROOT.TCanvas("event_tail", "Event-level tail risk", 950, 720)
    canvas.SetLeftMargin(0.13)
    canvas.SetRightMargin(0.04)
    canvas.SetTopMargin(0.13)
    canvas.SetBottomMargin(0.13)
    canvas.SetGridx()
    frame = ROOT.TH1D(
        "event_tail_frame",
        ";Event recovery survival;Fraction of events below value",
        100,
        0.0,
        1.01,
    )
    frame.SetDirectory(0)
    frame.SetMinimum(0.0)
    frame.SetMaximum(1.02)
    frame.Draw("AXIS")
    retained = [frame]
    colors = {
        10.0: ROOT.TColor.GetColor("#2B6CB0"),
        20.0: ROOT.TColor.GetColor("#159D66"),
        25.0: ROOT.TColor.GetColor("#C83E4D"),
    }
    legend = ROOT.TLegend(0.17, 0.67, 0.55, 0.87)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    for pitch_um in EVENT_PITCHES_UM:
        values = sorted(
            float(row["recovery_survival"])
            for row in rows
            if row["sample"] == "1 TeV pion"
            and np.isclose(float(row["pitch_um"]), pitch_um)
            and np.isclose(float(row["recovery_ns"]), 4.3)
            and row["recovery_survival"] != ""
        )
        cdf = (np.arange(len(values), dtype=np.float64) + 1) / len(values)
        graph = _graph(
            ROOT, values, cdf, color=colors[pitch_um], marker=1, width=3
        )
        graph.Draw("L SAME")
        tenth = float(np.quantile(values, 0.10))
        legend.AddEntry(graph, f"{pitch_um:g} #mum (10th pct {tenth:.1%})", "l")
        retained.append(graph)
    legend.Draw()
    retained.append(legend)
    title = ROOT.TLatex()
    title.SetNDC(True)
    title.SetTextFont(62)
    title.SetTextSize(0.044)
    title.DrawLatex(0.13, 0.94, "Aggregate retention hides bright-event tails")
    subtitle = ROOT.TLatex()
    subtitle.SetNDC(True)
    subtitle.SetTextFont(42)
    subtitle.SetTextSize(0.030)
    subtitle.DrawLatex(0.13, 0.895, "1 TeV pion; 20% PDE; 4.3 ns recovery")
    retained.extend([title, subtitle])
    return _save_canvas(canvas, output_dir / "event_level_tail_risk")


def _draw_timing_depth(
    output_dir: Path, rows: list[dict[str, object]]
) -> list[Path]:
    ROOT = _root_style()
    canvas = ROOT.TCanvas("timing_depth", "Timing depth proxy", 980, 720)
    canvas.SetLeftMargin(0.13)
    canvas.SetRightMargin(0.04)
    canvas.SetTopMargin(0.15)
    canvas.SetBottomMargin(0.13)
    canvas.SetGridy()
    frame = ROOT.TH1D(
        "timing_depth_frame",
        ";Single-photon response FWHM [ps];Added depth-proxy RMSE vs. ideal [mm]",
        100,
        -25.0,
        1200.0,
    )
    frame.SetDirectory(0)
    values = [float(row["test_rmse_delta_vs_ideal_mm"]) for row in rows]
    lows = [float(row["test_rmse_delta_low68_mm"]) for row in rows]
    highs = [float(row["test_rmse_delta_high68_mm"]) for row in rows]
    frame.SetMinimum(min(-0.5, min(lows) - 0.5))
    frame.SetMaximum(max(values + highs) + 1.5)
    frame.Draw("AXIS")
    retained = [frame]
    graph = ROOT.TGraphAsymmErrors(len(rows))
    for index, row in enumerate(rows):
        x = float(row["response_fwhm_ps"])
        y = float(row["test_rmse_delta_vs_ideal_mm"])
        graph.SetPoint(index, x, y)
        graph.SetPointError(
            index,
            0.0,
            0.0,
            y - float(row["test_rmse_delta_low68_mm"]),
            float(row["test_rmse_delta_high68_mm"]) - y,
        )
    graph.SetLineColor(ROOT.TColor.GetColor("#7A5195"))
    graph.SetMarkerColor(ROOT.TColor.GetColor("#7A5195"))
    graph.SetLineWidth(3)
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(1.2)
    graph.Draw("LP SAME")
    retained.append(graph)
    title = ROOT.TLatex()
    title.SetNDC(True)
    title.SetTextFont(62)
    title.SetTextSize(0.043)
    title.DrawLatex(0.13, 0.945, "Sub-100 ps smearing is not the limiting term for this proxy")
    subtitle = ROOT.TLatex()
    subtitle.SetNDC(True)
    subtitle.SetTextFont(42)
    subtitle.SetTextSize(0.028)
    subtitle.DrawLatex(
        0.13,
        0.900,
        "Paired event bootstrap; ideal test RMSE = 236 mm (tail-sensitive proof of concept)",
    )
    retained.extend([title, subtitle])
    return _save_canvas(canvas, output_dir / "timing_depth_proxy")


def _draw_physical_pulse(
    output_dir: Path, rows: list[dict[str, object]]
) -> list[Path]:
    ROOT = _root_style()
    canvas = ROOT.TCanvas("physical_pulse", "Physical pulse shape", 1050, 780)
    canvas.Divide(1, 2, 0.0, 0.0)
    colors = {
        10.0: ROOT.TColor.GetColor("#2B6CB0"),
        20.0: ROOT.TColor.GetColor("#159D66"),
        25.0: ROOT.TColor.GetColor("#C83E4D"),
    }
    retained = []
    selected_20 = [row for row in rows if np.isclose(float(row["pitch_um"]), 20.0)]
    time_ns = [float(row["time_center_ns"]) for row in selected_20]
    detected_k = [float(row["detected_photons"]) / 1000.0 for row in selected_20]
    accepted_k = [float(row["accepted_photons"]) / 1000.0 for row in selected_20]

    canvas.cd(1)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.04)
    ROOT.gPad.SetTopMargin(0.19)
    ROOT.gPad.SetBottomMargin(0.02)
    frame = ROOT.TH1D(
        "physical_pulse_frame_top",
        ";;Photons / 50 ps [10^{3}]",
        100,
        5.0,
        20.0,
    )
    frame.SetDirectory(0)
    frame.SetMinimum(0.0)
    frame.SetMaximum(1.12 * max(detected_k))
    frame.GetXaxis().SetLabelSize(0.0)
    frame.Draw("AXIS")
    retained.append(frame)
    detected_graph = _graph(
        ROOT, time_ns, detected_k, color=ROOT.TColor.GetColor("#4B5563"), width=3
    )
    accepted_graph = _graph(ROOT, time_ns, accepted_k, color=colors[20.0], width=3)
    detected_graph.Draw("L SAME")
    accepted_graph.Draw("L SAME")
    retained.extend([detected_graph, accepted_graph])
    legend = ROOT.TLegend(0.57, 0.68, 0.94, 0.84)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(detected_graph, "After 20% PDE", "l")
    legend.AddEntry(accepted_graph, "20 #mum + 4.3 ns recovery", "l")
    legend.Draw()
    retained.append(legend)
    title = ROOT.TLatex()
    title.SetNDC(True)
    title.SetTextFont(62)
    title.SetTextSize(0.052)
    title.DrawLatex(0.12, 0.945, "Recovery distorts the prompt photon-arrival peak")
    subtitle = ROOT.TLatex()
    subtitle.SetNDC(True)
    subtitle.SetTextFont(42)
    subtitle.SetTextSize(0.032)
    subtitle.DrawLatex(0.12, 0.895, "1 TeV pion; one Cherenkov fiber; physical 3 #times 3 mm^{2} sensor")
    retained.extend([title, subtitle])

    canvas.cd(2)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.04)
    ROOT.gPad.SetTopMargin(0.03)
    ROOT.gPad.SetBottomMargin(0.24)
    ROOT.gPad.SetGridy()
    frame2 = ROOT.TH1D(
        "physical_pulse_frame_bottom",
        ";Photon arrival time [ns];Accepted / detected",
        100,
        5.0,
        20.0,
    )
    frame2.SetDirectory(0)
    frame2.SetMinimum(0.45)
    frame2.SetMaximum(1.03)
    frame2.Draw("AXIS")
    retained.append(frame2)
    legend2 = ROOT.TLegend(0.66, 0.30, 0.94, 0.58)
    legend2.SetBorderSize(0)
    legend2.SetFillStyle(0)
    for pitch_um in EVENT_PITCHES_UM:
        selected = [
            row
            for row in rows
            if np.isclose(float(row["pitch_um"]), pitch_um)
            and row["recovery_survival"] != ""
            and int(row["detected_photons"]) >= 100
        ]
        graph = _graph(
            ROOT,
            [float(row["time_center_ns"]) for row in selected],
            [float(row["recovery_survival"]) for row in selected],
            color=colors[pitch_um],
            width=3,
        )
        graph.Draw("L SAME")
        legend2.AddEntry(graph, f"{pitch_um:g} #mum pitch", "l")
        retained.append(graph)
    legend2.Draw()
    retained.append(legend2)
    return _save_canvas(canvas, output_dir / "physical_pulse_shape")


def _find_scan_point(
    rows: list[dict[str, object]],
    *,
    sample: str,
    pde: float,
    pitch_um: float,
    recovery_ns: float,
) -> dict[str, object]:
    matches = _select_scan(
        rows,
        sample=sample,
        pde=pde,
        pitch_um=pitch_um,
        recovery_ns=recovery_ns,
    )
    if len(matches) != 1:
        raise ValueError(f"Expected one scan point, found {len(matches)}")
    return matches[0]


def run_study(
    output_dir: Path,
    *,
    max_events: int | None,
    bootstrap_samples: int,
) -> list[Path]:
    output_dir = output_dir.expanduser().resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    started = time.time()

    input_audit = _audit_inputs(max_events=max_events)
    input_audit_path = output_dir / "input_selection_audit.json"
    with input_audit_path.open("w") as handle:
        json.dump(input_audit, handle, indent=2)
        handle.write("\n")

    physical_rows = _scan_physical(
        max_events=max_events,
        bootstrap_samples=bootstrap_samples,
    )
    physical_csv = output_dir / "physical_pde_pitch_recovery.csv"
    _write_csv(physical_csv, physical_rows)

    mapping_rows = _scan_mapping_comparison(
        physical_rows,
        max_events=max_events,
    )
    mapping_csv = output_dir / "legacy_coordinate_comparison.csv"
    _write_csv(mapping_csv, mapping_rows)

    event_rows = _event_response_rows(max_events=max_events)
    event_csv = output_dir / "event_response.csv"
    _write_csv(event_csv, event_rows)

    timing_rows = _timing_depth_rows(max_events=max_events)
    timing_csv = output_dir / "timing_depth_proxy.csv"
    _write_csv(timing_csv, timing_rows)

    pulse_rows = _physical_pulse_rows(max_events=max_events)
    pulse_csv = output_dir / "physical_pulse_shape.csv"
    _write_csv(pulse_csv, pulse_rows)

    inputs = {name: inspect_root_file(path) for name, path in SAMPLES.items()}
    point_20 = _find_scan_point(
        physical_rows,
        sample="1 TeV pion",
        pde=0.20,
        pitch_um=20.0,
        recovery_ns=4.3,
    )
    point_25 = _find_scan_point(
        physical_rows,
        sample="1 TeV pion",
        pde=0.20,
        pitch_um=25.0,
        recovery_ns=4.3,
    )
    mapping_25 = next(
        row
        for row in mapping_rows
        if row["sample"] == "1 TeV pion"
        and np.isclose(float(row["pitch_um"]), 25.0)
        and np.isclose(float(row["recovery_ns"]), 4.3)
    )
    summary = {
        "central_correction": (
            "The historical plot is retained as an exact software closure, but its "
            "3 mm segmentation was applied across a 1 mm fiber-local display range. "
            "Physical sensor requirements therefore use coordinate_scale=1."
        ),
        "pion_unity_pde_25um_4p3ns": {
            "historical_display_scaled_survival": mapping_25[
                "historical_display_scaled_survival"
            ],
            "physical_mm_survival": mapping_25["physical_mm_survival"],
        },
        "pion_20percent_pde_4p3ns": {
            "20um_survival": point_20["recovery_survival"],
            "20um_usable_fraction": point_20["usable_fraction"],
            "25um_survival": point_25["recovery_survival"],
            "25um_usable_fraction": point_25["usable_fraction"],
        },
        "timing_depth_proxy": timing_rows,
        "interpretation": [
            "No universal pitch working point exists without a PDE and recovery assumption.",
            "The 100 GeV positron and 1 TeV pion files are unmatched-energy stress samples, not a species comparison.",
            "The timing-depth observable is a proof of information retention, not a calorimeter longitudinal-resolution measurement.",
            "The ideal-to-50 ps change is small for this many-photon quantile, so a 10--50 ps device target needs additional observables or an independent system requirement.",
        ],
    }
    summary_path = output_dir / "study_summary.json"
    with summary_path.open("w") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")

    plot_paths = [
        *_draw_mapping_correction(output_dir, mapping_rows),
        *_draw_requirement_surface(output_dir, physical_rows),
        *_draw_event_tail(output_dir, event_rows),
        *_draw_physical_pulse(output_dir, pulse_rows),
        *_draw_timing_depth(output_dir, timing_rows),
    ]

    manifest = {
        "study": "FLASH/HG-DREAM physical sensor requirements v2",
        "inputs": inputs,
        "selection": {
            "geometry_mode": "old_2025 fiber-center shifts",
            "production_fiber": 0,
            "photon_type": "Cherenkov photons exiting the positive-z end",
            "gate_ns": [0.0, 40.0],
            "sensor_bounds_mm": [SENSOR_MIN_MM, SENSOR_MAX_MM],
            "coordinate_scale": 1.0,
        },
        "scan": {
            "pde": list(PDES),
            "pitch_um": list(PITCHES_UM),
            "recovery_ns": list(RECOVERY_NS),
            "seed": 4307,
            "bootstrap_samples": bootstrap_samples,
            "max_events_per_sample": max_events,
        },
        "timing_depth_proxy": {
            "sample": "1 TeV pion",
            "front_end": "20% PDE, 20 um pitch, 4.3 ns recovery",
            "observable": "10th percentile of accepted arrival times",
            "truth": "mean selected-photon production depth in 0--2000 mm",
            "split": "even event numbers train; odd event numbers test",
            "scenarios": [
                {
                    "name": name,
                    "response_fwhm_ps": response,
                    "tdc_bin_ps": tdc,
                }
                for name, response, tdc in TIMING_SCENARIOS
            ],
        },
        "limitations": [
            "One Cherenkov fiber is analyzed; this is a local occupancy benchmark, not a full calorimeter response.",
            "The two large samples have different beam energies.",
            "PDE is wavelength-independent in the requirement grid.",
            "No dark counts, optical crosstalk, afterpulsing, hot-cell masking, optical coupling, or bandwidth/power model is included.",
            "The timing-depth proxy uses photon-production truth and needs validation with matched beam-test observables.",
        ],
        "elapsed_seconds": time.time() - started,
    }
    manifest_path = output_dir / "run_manifest.json"
    with manifest_path.open("w") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")

    return [
        physical_csv,
        input_audit_path,
        mapping_csv,
        event_csv,
        timing_csv,
        pulse_csv,
        summary_path,
        manifest_path,
        *plot_paths,
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_DIR / "results" / "flash-study-v2",
    )
    parser.add_argument("--max-events", type=int)
    parser.add_argument("--bootstrap-samples", type=int, default=400)
    args = parser.parse_args()
    outputs = run_study(
        args.output,
        max_events=args.max_events,
        bootstrap_samples=args.bootstrap_samples,
    )
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
