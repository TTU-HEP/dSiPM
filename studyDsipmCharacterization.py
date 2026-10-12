#!/usr/bin/env python3
"""Run standard dSiPM characterization studies for FLASH and HG-DREAM."""

from __future__ import annotations

import argparse
from array import array
import csv
import json
from pathlib import Path
import time

import numpy as np

from dsipm.characterization import (
    PDERequirementPoint,
    expected_unique_cells,
    scan_pde_pitch_recovery,
)
from dsipm.root_io import inspect_root_file, iter_root_events


REPO_DIR = Path(__file__).resolve().parent
DEFAULT_SAMPLES = {
    "100 GeV positron": REPO_DIR / "mc_dreamsim_e+_100_101.root",
    "1 TeV pion": REPO_DIR / "mc_dreamsim_pi+_1000_1001.root",
}
DEFAULT_PDES = (0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00)
DEFAULT_PITCHES_UM = (10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0)
DEFAULT_RECOVERY_NS = (1.0, 4.3, 40.0)
LINEARITY_CELL_COUNTS = (1024, 4096, 16384)
LINEARITY_OCCUPANCIES = (0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1.0, 2.0, 4.0, 8.0)


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Cannot write empty table: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _validate_unity_pde_against_historical_scan(
    rows: list[dict[str, object]],
    reference_path: Path,
) -> dict[str, object]:
    """Check the new scan against the already validated unity-PDE baseline."""

    with reference_path.open() as handle:
        reference_rows = list(csv.DictReader(handle))
    reference = {
        (
            row["sample"],
            float(row["pitch_um"]),
            float(row["recovery_ns"]),
        ): row
        for row in reference_rows
    }
    comparisons = []
    for row in rows:
        if not np.isclose(float(row["pde"]), 1.0):
            continue
        key = (
            str(row["sample"]),
            float(row["pitch_um"]),
            float(row["recovery_ns"]),
        )
        if key not in reference:
            continue
        old = reference[key]
        comparisons.append(
            {
                "sample": key[0],
                "pitch_um": key[1],
                "recovery_ns": key[2],
                "incident_difference": int(row["incident_photons"])
                - int(old["incident_photons"]),
                "accepted_difference": int(row["accepted_photons"])
                - int(old["accepted_photons"]),
                "retained_fraction_difference": float(
                    row["recovery_retention_accepted_over_detected"]
                )
                - float(old["retained_fraction"]),
            }
        )
    if not comparisons:
        raise ValueError("No overlapping unity-PDE points found for closure test")
    return {
        "reference": str(reference_path.resolve()),
        "checked_points": len(comparisons),
        "all_integer_counts_exact": all(
            comparison["incident_difference"] == 0
            and comparison["accepted_difference"] == 0
            for comparison in comparisons
        ),
        "maximum_absolute_retained_fraction_difference": max(
            abs(float(comparison["retained_fraction_difference"]))
            for comparison in comparisons
        ),
        "comparisons": comparisons,
    }


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
        "usable_fraction_accepted_over_incident": point.usable_fraction,
        "usable_fraction_low68": point.usable_fraction_low68,
        "usable_fraction_high68": point.usable_fraction_high68,
        "recovery_retention_accepted_over_detected": point.recovery_retention,
        "recovery_retention_low68": point.recovery_retention_low68,
        "recovery_retention_high68": point.recovery_retention_high68,
    }


def _matching_rows(
    rows: list[dict[str, object]],
    *,
    sample: str,
    pitch_um: float | None = None,
    recovery_ns: float | None = None,
) -> list[dict[str, object]]:
    selected = []
    for row in rows:
        if row["sample"] != sample:
            continue
        if pitch_um is not None and not np.isclose(float(row["pitch_um"]), pitch_um):
            continue
        if recovery_ns is not None and not np.isclose(
            float(row["recovery_ns"]), recovery_ns
        ):
            continue
        selected.append(row)
    return selected


def _root_style() -> object:
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTitleFont(42, "XYZ")
    ROOT.gStyle.SetLabelFont(42, "XYZ")
    return ROOT


def _graph_with_intervals(
    ROOT: object,
    selected: list[dict[str, object]],
    value_key: str,
    low_key: str,
    high_key: str,
) -> object:
    selected = sorted(selected, key=lambda row: float(row["pde"]))
    graph = ROOT.TGraphAsymmErrors(len(selected))
    for index, row in enumerate(selected):
        x = float(row["pde"])
        y = float(row[value_key])
        low = float(row[low_key])
        high = float(row[high_key])
        graph.SetPoint(index, x, y)
        graph.SetPointError(index, 0.0, 0.0, max(0.0, y - low), max(0.0, high - y))
    return graph


def _draw_pde_tradeoff(
    output_dir: Path,
    rows: list[dict[str, object]],
    samples: list[str],
    *,
    pitch_um: float = 25.0,
) -> Path:
    ROOT = _root_style()
    colors = {
        1.0: ROOT.TColor.GetColor("#159D66"),
        4.3: ROOT.TColor.GetColor("#7A5195"),
        40.0: ROOT.TColor.GetColor("#C83E4D"),
    }
    retained: list[object] = []
    n_columns = len(samples)
    canvas = ROOT.TCanvas(
        "pde_tradeoff",
        "PDE pitch recovery tradeoff",
        820 * n_columns,
        1120,
    )
    canvas.Divide(n_columns, 2, 0.006, 0.006)

    for column, sample in enumerate(samples):
        sample_rows = _matching_rows(rows, sample=sample, pitch_um=pitch_um)
        for row_index, metric in enumerate(("usable", "survival")):
            pad_index = column + 1 + row_index * n_columns
            canvas.cd(pad_index)
            ROOT.gPad.SetLeftMargin(0.13 if column == 0 else 0.11)
            ROOT.gPad.SetRightMargin(0.035)
            ROOT.gPad.SetTopMargin(0.13)
            ROOT.gPad.SetBottomMargin(0.14)
            ROOT.gPad.SetGridy()
            ROOT.gPad.SetTicks(1, 1)

            if metric == "usable":
                title = ";Photon-detection efficiency;Accepted / incident photons"
                y_min = 0.0
                value_key = "usable_fraction_accepted_over_incident"
                low_key = "usable_fraction_low68"
                high_key = "usable_fraction_high68"
            else:
                title = ";Photon-detection efficiency;Accepted / detected photons"
                values = [
                    float(row["recovery_retention_accepted_over_detected"])
                    for row in sample_rows
                ]
                y_min = max(0.0, np.floor((min(values) - 0.02) * 20.0) / 20.0)
                value_key = "recovery_retention_accepted_over_detected"
                low_key = "recovery_retention_low68"
                high_key = "recovery_retention_high68"

            frame = ROOT.TH1D(
                f"pde_frame_{column}_{row_index}",
                title,
                100,
                0.05,
                1.05,
            )
            frame.SetDirectory(0)
            frame.SetMinimum(float(y_min))
            frame.SetMaximum(1.02)
            frame.GetXaxis().SetTitleSize(0.047)
            frame.GetYaxis().SetTitleSize(0.047)
            frame.GetXaxis().SetLabelSize(0.039)
            frame.GetYaxis().SetLabelSize(0.039)
            frame.GetXaxis().SetTitleOffset(1.12)
            frame.GetYaxis().SetTitleOffset(1.18 if column == 0 else 1.0)
            frame.Draw("AXIS")
            retained.append(frame)

            if metric == "usable":
                upper_bound = ROOT.TLine(0.05, 0.05, 1.0, 1.0)
                upper_bound.SetLineColor(ROOT.TColor.GetColor("#707880"))
                upper_bound.SetLineStyle(2)
                upper_bound.SetLineWidth(2)
                upper_bound.Draw("SAME")
                retained.append(upper_bound)

            legend = ROOT.TLegend(0.16, 0.69, 0.60, 0.88)
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            legend.SetTextSize(0.035)
            for recovery_ns in DEFAULT_RECOVERY_NS:
                selected = [
                    row
                    for row in sample_rows
                    if np.isclose(float(row["recovery_ns"]), recovery_ns)
                ]
                graph = _graph_with_intervals(
                    ROOT,
                    selected,
                    value_key,
                    low_key,
                    high_key,
                )
                graph.SetName(
                    f"pde_{metric}_{column}_{str(recovery_ns).replace('.', 'p')}"
                )
                graph.SetLineColor(colors[recovery_ns])
                graph.SetMarkerColor(colors[recovery_ns])
                graph.SetLineWidth(3)
                graph.SetMarkerStyle(20 if recovery_ns != 40.0 else 24)
                graph.SetMarkerSize(1.0)
                graph.Draw("LP SAME")
                label = "40 ns (one hit / event)" if recovery_ns == 40.0 else f"{recovery_ns:g} ns"
                legend.AddEntry(graph, label, "lp")
                retained.append(graph)
            if row_index == 0:
                legend.Draw()
                retained.append(legend)

            title_text = ROOT.TLatex()
            title_text.SetNDC(True)
            title_text.SetTextFont(62)
            title_text.SetTextSize(0.040)
            title_text.DrawLatex(0.13, 0.945, sample)
            subtitle = ROOT.TLatex()
            subtitle.SetNDC(True)
            subtitle.SetTextFont(42)
            subtitle.SetTextSize(0.029)
            subtitle.SetTextColor(ROOT.TColor.GetColor("#4B5563"))
            description = (
                f"{pitch_um:g} #mum pitch; signal yield"
                if metric == "usable"
                else f"{pitch_um:g} #mum pitch; recovery survival"
            )
            subtitle.DrawLatex(0.13, 0.900, description)
            retained.extend([title_text, subtitle])

    output = output_dir / "pde_tradeoff_25um.png"
    canvas.SaveAs(str(output))
    return output


def _bin_edges(values: tuple[float, ...]) -> array:
    points = np.asarray(values, dtype=np.float64)
    if points.size == 1:
        return array("d", [points[0] - 0.5, points[0] + 0.5])
    edges = np.empty(points.size + 1, dtype=np.float64)
    edges[1:-1] = 0.5 * (points[:-1] + points[1:])
    edges[0] = points[0] - 0.5 * (points[1] - points[0])
    edges[-1] = points[-1] + 0.5 * (points[-1] - points[-2])
    return array("d", edges)


def _draw_design_maps(
    output_dir: Path,
    rows: list[dict[str, object]],
    sample: str,
) -> Path:
    ROOT = _root_style()
    ROOT.gStyle.SetPaintTextFormat(".3f")
    canvas = ROOT.TCanvas("design_maps", "PDE pitch design maps", 2250, 1200)
    canvas.Divide(3, 2, 0.006, 0.006)
    retained: list[object] = []
    pde_edges = _bin_edges(DEFAULT_PDES)
    pitch_edges = _bin_edges(DEFAULT_PITCHES_UM)

    for column, recovery_ns in enumerate(DEFAULT_RECOVERY_NS):
        for row_index, metric in enumerate(("usable", "survival")):
            canvas.cd(column + 1 + row_index * 3)
            ROOT.gPad.SetLeftMargin(0.12)
            ROOT.gPad.SetRightMargin(0.15)
            ROOT.gPad.SetTopMargin(0.14)
            ROOT.gPad.SetBottomMargin(0.13)
            ROOT.gPad.SetTicks(1, 1)
            key = (
                "usable_fraction_accepted_over_incident"
                if metric == "usable"
                else "recovery_retention_accepted_over_detected"
            )
            z_title = (
                "Accepted / incident" if metric == "usable" else "Accepted / detected"
            )
            histogram = ROOT.TH2D(
                f"map_{column}_{row_index}",
                f";Photon-detection efficiency;Physical SPAD pitch [#mum];{z_title}",
                len(DEFAULT_PDES),
                pde_edges,
                len(DEFAULT_PITCHES_UM),
                pitch_edges,
            )
            histogram.SetDirectory(0)
            selected = _matching_rows(
                rows,
                sample=sample,
                recovery_ns=recovery_ns,
            )
            for point in selected:
                histogram.Fill(
                    float(point["pde"]),
                    float(point["pitch_um"]),
                    float(point[key]),
                )
            histogram.SetMinimum(0.0 if metric == "usable" else 0.80)
            histogram.SetMaximum(1.0)
            histogram.GetXaxis().SetTitleSize(0.045)
            histogram.GetYaxis().SetTitleSize(0.045)
            histogram.GetZaxis().SetTitleSize(0.040)
            histogram.GetXaxis().SetLabelSize(0.036)
            histogram.GetYaxis().SetLabelSize(0.036)
            histogram.GetZaxis().SetLabelSize(0.032)
            histogram.GetXaxis().SetTitleOffset(1.15)
            histogram.GetYaxis().SetTitleOffset(1.15)
            histogram.GetZaxis().SetTitleOffset(1.25)
            histogram.SetMarkerSize(1.15)
            histogram.Draw("COLZ TEXT")
            retained.append(histogram)

            title = ROOT.TLatex()
            title.SetNDC(True)
            title.SetTextFont(62)
            title.SetTextSize(0.038)
            title.DrawLatex(0.12, 0.955, f"{recovery_ns:g} ns recovery")
            subtitle = ROOT.TLatex()
            subtitle.SetNDC(True)
            subtitle.SetTextFont(42)
            subtitle.SetTextSize(0.027)
            subtitle.SetTextColor(ROOT.TColor.GetColor("#4B5563"))
            subtitle.DrawLatex(0.12, 0.915, f"{sample}; {z_title.lower()}")
            retained.extend([title, subtitle])

    output = output_dir / "pion_pde_pitch_recovery_maps.png"
    canvas.SaveAs(str(output))
    return output


def _simulate_uniform_flash(
    *,
    trials: int,
    seed: int,
) -> list[dict[str, object]]:
    if trials <= 0:
        raise ValueError("linearity trials must be positive")
    rng = np.random.default_rng(seed)
    rows: list[dict[str, object]] = []
    for n_cells in LINEARITY_CELL_COUNTS:
        for occupancy in LINEARITY_OCCUPANCIES:
            incident = max(1, int(round(occupancy * n_cells)))
            observed = np.empty(trials, dtype=np.float64)
            for trial in range(trials):
                cells = rng.integers(0, n_cells, size=incident)
                observed[trial] = np.unique(cells).size
            exact = float(expected_unique_cells(incident, n_cells))
            mean = float(observed.mean())
            sem = float(observed.std(ddof=1) / np.sqrt(trials)) if trials > 1 else 0.0
            rows.append(
                {
                    "n_cells": n_cells,
                    "trials": trials,
                    "incident_photons": incident,
                    "incident_photons_per_cell": incident / n_cells,
                    "exact_mean_fired_cells": exact,
                    "mc_mean_fired_cells": mean,
                    "mc_sem_fired_cells": sem,
                    "exact_fired_fraction": exact / n_cells,
                    "mc_fired_fraction": mean / n_cells,
                    "mc_minus_exact_per_mille_of_cells": 1000.0
                    * (mean - exact)
                    / n_cells,
                }
            )
    return rows


def _draw_linearity_closure(
    output_dir: Path,
    rows: list[dict[str, object]],
) -> Path:
    ROOT = _root_style()
    colors = {
        1024: ROOT.TColor.GetColor("#2F6690"),
        4096: ROOT.TColor.GetColor("#7A5195"),
        16384: ROOT.TColor.GetColor("#E07A1F"),
    }
    canvas = ROOT.TCanvas("linearity_closure", "Uniform flash closure", 1100, 1050)
    upper = ROOT.TPad("linearity_upper", "", 0.0, 0.31, 1.0, 1.0)
    lower = ROOT.TPad("linearity_lower", "", 0.0, 0.0, 1.0, 0.31)
    upper.SetLeftMargin(0.13)
    upper.SetRightMargin(0.04)
    upper.SetTopMargin(0.12)
    upper.SetBottomMargin(0.02)
    lower.SetLeftMargin(0.13)
    lower.SetRightMargin(0.04)
    lower.SetTopMargin(0.02)
    lower.SetBottomMargin(0.31)
    upper.Draw()
    lower.Draw()
    retained: list[object] = [upper, lower]

    upper.cd()
    upper.SetTicks(1, 1)
    upper.SetGridy()
    frame = ROOT.TH1D(
        "linearity_frame",
        ";;Fired cells / available cells",
        100,
        0.0,
        8.2,
    )
    frame.SetDirectory(0)
    frame.SetMinimum(0.0)
    frame.SetMaximum(1.05)
    frame.GetYaxis().SetTitleSize(0.052)
    frame.GetYaxis().SetLabelSize(0.043)
    frame.GetYaxis().SetTitleOffset(1.16)
    frame.GetXaxis().SetLabelSize(0.0)
    frame.Draw("AXIS")
    retained.append(frame)

    ideal = ROOT.TLine(0.0, 0.0, 1.0, 1.0)
    ideal.SetLineColor(ROOT.TColor.GetColor("#80878E"))
    ideal.SetLineStyle(2)
    ideal.SetLineWidth(2)
    ideal.Draw("SAME")
    retained.append(ideal)

    legend = ROOT.TLegend(0.58, 0.22, 0.92, 0.48)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.037)
    first_exact_graph = None
    for n_cells in LINEARITY_CELL_COUNTS:
        selected = [row for row in rows if int(row["n_cells"]) == n_cells]
        graph = ROOT.TGraphErrors(len(selected))
        for index, row in enumerate(selected):
            graph.SetPoint(
                index,
                float(row["incident_photons_per_cell"]),
                float(row["mc_fired_fraction"]),
            )
            graph.SetPointError(
                index,
                0.0,
                float(row["mc_sem_fired_cells"]) / n_cells,
            )
        graph.SetLineColor(colors[n_cells])
        graph.SetMarkerColor(colors[n_cells])
        graph.SetMarkerStyle(20)
        graph.SetMarkerSize(1.05)
        graph.SetLineWidth(2)
        graph.Draw("P SAME")
        legend.AddEntry(graph, f"{n_cells:,} cells: Monte Carlo", "p")
        retained.append(graph)

        exact_x = np.linspace(0.0, 8.0, 400)
        exact_n = exact_x * n_cells
        exact_y = expected_unique_cells(exact_n, n_cells) / n_cells
        exact_graph = ROOT.TGraph(len(exact_x), exact_x, exact_y)
        exact_graph.SetLineColor(colors[n_cells])
        exact_graph.SetLineWidth(3)
        exact_graph.Draw("L SAME")
        if first_exact_graph is None:
            first_exact_graph = exact_graph
        retained.append(exact_graph)
    if first_exact_graph is not None:
        legend.AddEntry(first_exact_graph, "Exact finite-cell model", "l")
    legend.Draw()
    retained.append(legend)

    title = ROOT.TLatex()
    title.SetNDC(True)
    title.SetTextFont(62)
    title.SetTextSize(0.044)
    title.DrawLatex(0.13, 0.945, "Standard prompt-flash saturation closure")
    subtitle = ROOT.TLatex()
    subtitle.SetNDC(True)
    subtitle.SetTextFont(42)
    subtitle.SetTextSize(0.031)
    subtitle.SetTextColor(ROOT.TColor.GetColor("#4B5563"))
    subtitle.DrawLatex(0.13, 0.900, "Uniform illumination; one count per cell; no recovery during pulse")
    retained.extend([title, subtitle])

    lower.cd()
    lower.SetTicks(1, 1)
    lower.SetGridy()
    residual_values = [
        abs(float(row["mc_minus_exact_per_mille_of_cells"])) for row in rows
    ]
    residual_limit = max(0.30, 1.35 * max(residual_values))
    residual_frame = ROOT.TH1D(
        "linearity_residual_frame",
        ";Incident photons / available cells;(MC - exact) / cells [#permil]",
        100,
        0.0,
        8.2,
    )
    residual_frame.SetDirectory(0)
    residual_frame.SetMinimum(-residual_limit)
    residual_frame.SetMaximum(residual_limit)
    residual_frame.GetXaxis().SetTitleSize(0.105)
    residual_frame.GetYaxis().SetTitleSize(0.085)
    residual_frame.GetXaxis().SetLabelSize(0.082)
    residual_frame.GetYaxis().SetLabelSize(0.072)
    residual_frame.GetXaxis().SetTitleOffset(1.15)
    residual_frame.GetYaxis().SetTitleOffset(0.64)
    residual_frame.GetYaxis().SetNdivisions(505)
    residual_frame.Draw("AXIS")
    retained.append(residual_frame)
    zero = ROOT.TLine(0.0, 0.0, 8.2, 0.0)
    zero.SetLineColor(ROOT.TColor.GetColor("#4B5563"))
    zero.SetLineStyle(2)
    zero.Draw("SAME")
    retained.append(zero)

    for n_cells in LINEARITY_CELL_COUNTS:
        selected = [row for row in rows if int(row["n_cells"]) == n_cells]
        graph = ROOT.TGraphErrors(len(selected))
        for index, row in enumerate(selected):
            graph.SetPoint(
                index,
                float(row["incident_photons_per_cell"]),
                float(row["mc_minus_exact_per_mille_of_cells"]),
            )
            graph.SetPointError(
                index,
                0.0,
                1000.0 * float(row["mc_sem_fired_cells"]) / n_cells,
            )
        graph.SetMarkerColor(colors[n_cells])
        graph.SetLineColor(colors[n_cells])
        graph.SetMarkerStyle(20)
        graph.SetMarkerSize(0.95)
        graph.Draw("P SAME")
        retained.append(graph)

    output = output_dir / "standard_prompt_flash_linearity.png"
    canvas.SaveAs(str(output))
    return output


def _nonlinearity_threshold(target_response_fraction: float) -> float:
    x = np.linspace(1.0e-8, 2.0, 2_000_000)
    response = -np.expm1(-x) / x
    index = int(np.argmin(np.abs(response - target_response_fraction)))
    return float(x[index])


def run_study(
    output_dir: Path,
    *,
    sample_names: list[str],
    max_events: int | None,
    bootstrap_samples: int,
    linearity_trials: int,
) -> list[Path]:
    output_dir = output_dir.expanduser().resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    started = time.time()
    scan_rows: list[dict[str, object]] = []
    inputs: dict[str, object] = {}
    for sample_index, sample in enumerate(sample_names):
        path = DEFAULT_SAMPLES[sample]
        inputs[sample] = inspect_root_file(path)
        events = iter_root_events(
            path,
            max_events=max_events,
            geometry_mode="old_2025",
            gate_ns=(0.0, 40.0),
            detector_min_mm=-100.0,
            detector_max_mm=100.0,
            fiber_index=0,
        )
        points = scan_pde_pitch_recovery(
            events,
            pdes=DEFAULT_PDES,
            pitches_um=DEFAULT_PITCHES_UM,
            recovery_times_ns=DEFAULT_RECOVERY_NS,
            seed=4307,
            coordinate_scale=3.0,
            bootstrap_samples=bootstrap_samples,
            bootstrap_seed=20261010 + sample_index,
        )
        scan_rows.extend(_point_row(sample, point) for point in points)

    scan_csv = output_dir / "pde_pitch_recovery_scan.csv"
    _write_csv(scan_csv, scan_rows)

    linearity_rows = _simulate_uniform_flash(
        trials=linearity_trials,
        seed=20261010,
    )
    linearity_csv = output_dir / "standard_prompt_flash_linearity.csv"
    _write_csv(linearity_csv, linearity_rows)

    closure_reference = (
        REPO_DIR / "results" / "flash-requirements-2026" / "flash_pitch_recovery_scan.csv"
    )
    closure_path = output_dir / "unity_pde_historical_closure.json"
    if max_events is None and closure_reference.exists():
        unity_closure = _validate_unity_pde_against_historical_scan(
            scan_rows,
            closure_reference,
        )
        with closure_path.open("w") as handle:
            json.dump(unity_closure, handle, indent=2)
            handle.write("\n")
    else:
        unity_closure = None

    key_results = []
    for row in scan_rows:
        if (
            np.isclose(float(row["pitch_um"]), 25.0)
            and float(row["pde"]) in (0.2, 0.5, 1.0)
            and float(row["recovery_ns"]) in (4.3, 40.0)
        ):
            key_results.append(row)

    summary = {
        "interpretation": {
            "usable_fraction": "accepted photons / incident photons",
            "recovery_retention": "accepted photons / PDE-selected photons",
            "why_both_are_required": (
                "Lower PDE reduces cell collisions and can raise recovery retention "
                "while simultaneously reducing usable photon yield."
            ),
        },
        "prompt_flash_dynamic_range": {
            "five_percent_nonlinearity_incident_photons_per_cell": _nonlinearity_threshold(0.95),
            "ten_percent_nonlinearity_incident_photons_per_cell": _nonlinearity_threshold(0.90),
            "model": "uniform simultaneous photons, one count per cell, no recovery",
            "maximum_absolute_mc_pull": max(
                abs(
                    (float(row["mc_mean_fired_cells"]) - float(row["exact_mean_fired_cells"]))
                    / float(row["mc_sem_fired_cells"])
                )
                for row in linearity_rows
                if float(row["mc_sem_fired_cells"]) > 0
            ),
        },
        "unity_pde_historical_closure": unity_closure,
        "key_25um_results": key_results,
    }
    summary_path = output_dir / "characterization_summary.json"
    with summary_path.open("w") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")

    plot_paths = [
        _draw_pde_tradeoff(output_dir, scan_rows, sample_names),
        _draw_design_maps(
            output_dir,
            scan_rows,
            "1 TeV pion" if "1 TeV pion" in sample_names else sample_names[0],
        ),
        _draw_linearity_closure(output_dir, linearity_rows),
    ]

    manifest = {
        "purpose": (
            "Add standard dSiPM photon-counting characterization to the validated "
            "HG-DREAM occupancy study."
        ),
        "inputs": inputs,
        "selection": {
            "geometry_mode": "old_2025",
            "fiber_index": 0,
            "gate_ns": [0.0, 40.0],
            "coordinate_scale": 3.0,
        },
        "pde_pitch_recovery_scan": {
            "pdes": list(DEFAULT_PDES),
            "pitches_um": list(DEFAULT_PITCHES_UM),
            "recovery_times_ns": list(DEFAULT_RECOVERY_NS),
            "seed": 4307,
            "max_events_per_sample": max_events,
            "bootstrap_samples": bootstrap_samples,
            "bootstrap_interval": "event-resampled 16th--84th percentile",
            "common_random_numbers": True,
        },
        "linearity_closure": {
            "cell_counts": list(LINEARITY_CELL_COUNTS),
            "incident_photons_per_cell": list(LINEARITY_OCCUPANCIES),
            "trials_per_point": linearity_trials,
            "seed": 20261010,
            "analytic_model": "M * (1 - (1 - 1/M)^N)",
        },
        "limitations": [
            "PDE is wavelength-independent in this scan.",
            "No dark counts, optical crosstalk, afterpulsing, or hot-cell masking.",
            "No timing-group, buffer, arbitration, bandwidth, or power model in these plots.",
            "The 100 GeV positron and 1 TeV pion samples are unmatched-energy stress cases.",
            "The uniform-flash plot validates the finite-cell transfer function; it is not an HG-DREAM event model.",
        ],
        "outputs": [
            path.name
            for path in [
                scan_csv,
                linearity_csv,
                *([closure_path] if unity_closure is not None else []),
                summary_path,
                *plot_paths,
            ]
        ],
        "elapsed_seconds": time.time() - started,
    }
    manifest_path = output_dir / "run_manifest.json"
    with manifest_path.open("w") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")

    return [
        scan_csv,
        linearity_csv,
        *([closure_path] if unity_closure is not None else []),
        summary_path,
        manifest_path,
        *plot_paths,
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_DIR / "results" / "dsipm-characterization-2026",
    )
    parser.add_argument(
        "--sample",
        choices=("both", "electron", "pion"),
        default="both",
    )
    parser.add_argument("--max-events", type=int)
    parser.add_argument("--bootstrap-samples", type=int, default=400)
    parser.add_argument("--linearity-trials", type=int, default=200)
    args = parser.parse_args()

    sample_names = {
        "both": list(DEFAULT_SAMPLES),
        "electron": ["100 GeV positron"],
        "pion": ["1 TeV pion"],
    }[args.sample]
    outputs = run_study(
        args.output,
        sample_names=sample_names,
        max_events=args.max_events,
        bootstrap_samples=args.bootstrap_samples,
        linearity_trials=args.linearity_trials,
    )
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
