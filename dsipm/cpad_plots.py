"""Slide-ready plots that map the historical study onto public dSiPMs.

The occupancy panel deliberately uses only the old unity-PDE, full-event
one-hit curves.  Device markers therefore answer one narrow question: where
does a public SPAD pitch sit on Chris Madrid's original occupancy study?  They
must not be read as full device-performance predictions.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np


CSV_FIELDS = [
    "name",
    "label",
    "role",
    "maturity",
    "pitch_um",
    "pitch_basis",
    "electron_pitch_only_efficiency",
    "pion_pitch_only_efficiency",
    "timing_low_ps",
    "timing_high_ps",
    "timing_value_ps",
    "timing_metric",
    "timing_status",
    "source_url",
]


def _load_historical_curves(
    output_root: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    summary_path = output_root / "electron_pion_efficiency_summary.csv"
    if not summary_path.exists():
        raise FileNotFoundError(
            f"Historical electron/pion summary does not exist: {summary_path}"
        )
    with summary_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"Historical summary has no rows: {summary_path}")

    rows.sort(key=lambda row: float(row["pitch_um"]))
    pitch = np.asarray([float(row["pitch_um"]) for row in rows])
    electron = np.asarray(
        [float(row["electron_one_hit_efficiency"]) for row in rows]
    )
    pion = np.asarray([float(row["pion_one_hit_efficiency"]) for row in rows])
    return pitch, electron, pion


def pitch_only_projection(
    pitch_um: float,
    historical_pitch_um: np.ndarray,
    historical_efficiency: np.ndarray,
) -> float:
    """Linearly interpolate a device pitch on the historical scan."""

    if pitch_um < historical_pitch_um[0] or pitch_um > historical_pitch_um[-1]:
        raise ValueError(
            f"Pitch {pitch_um:g} um is outside the historical scan "
            f"[{historical_pitch_um[0]:g}, {historical_pitch_um[-1]:g}] um"
        )
    return float(np.interp(pitch_um, historical_pitch_um, historical_efficiency))


def _comparison_rows(
    benchmarks: dict[str, Any],
    pitch: np.ndarray,
    electron: np.ndarray,
    pion: np.ndarray,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    flash = benchmarks["study_context"]["flash"]
    flash_pitch = float(flash["pitch_um"])
    common = {
        "name": "flash_working_point",
        "label": flash["label"],
        "role": "FLASH",
        "maturity": flash["pitch_status"],
        "pitch_um": flash_pitch,
        "pitch_basis": flash["pitch_basis"],
        "electron_pitch_only_efficiency": pitch_only_projection(
            flash_pitch, pitch, electron
        ),
        "pion_pitch_only_efficiency": pitch_only_projection(
            flash_pitch, pitch, pion
        ),
        "timing_low_ps": flash["timing_low_ps"],
        "timing_high_ps": flash["timing_high_ps"],
        "source_url": flash["source_url"],
    }
    rows.append(
        {
            **common,
            "timing_value_ps": "",
            "timing_metric": flash["timing_metric"],
            "timing_status": flash["timing_status"],
        }
    )

    for device in benchmarks["competitors"]:
        device_pitch = float(device["pitch_um"])
        common = {
            "name": device["name"],
            "label": device["short_label"],
            "role": "competitor",
            "maturity": device["maturity"],
            "pitch_um": device_pitch,
            "pitch_basis": device["pitch_basis"],
            "electron_pitch_only_efficiency": pitch_only_projection(
                device_pitch, pitch, electron
            ),
            "pion_pitch_only_efficiency": pitch_only_projection(
                device_pitch, pitch, pion
            ),
            "timing_low_ps": "",
            "timing_high_ps": "",
            "source_url": device["source_url"],
        }
        for timing in device["timing"]:
            rows.append(
                {
                    **common,
                    "timing_value_ps": timing["value_ps"],
                    "timing_metric": timing["metric"],
                    "timing_status": timing["status"],
                }
            )
    return rows


def _write_comparison_csv(output_root: Path, rows: list[dict[str, object]]) -> Path:
    destination = output_root / "cpad_device_landscape.csv"
    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    return destination


def _root():
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    return ROOT


def _root_color(ROOT, value: str) -> int:
    return int(ROOT.TColor.GetColor(value))


def _save_canvas(canvas, output_root: Path, stem: str) -> list[Path]:
    destinations = [output_root / f"{stem}.png", output_root / f"{stem}.pdf"]
    for destination in destinations:
        canvas.SaveAs(str(destination))
    return destinations


def _draw_occupancy_landscape(
    output_root: Path,
    benchmarks: dict[str, Any],
    rows: list[dict[str, object]],
    pitch: np.ndarray,
    electron: np.ndarray,
    pion: np.ndarray,
) -> list[Path]:
    ROOT = _root()
    electron_color = _root_color(ROOT, "#276FBF")
    pion_color = _root_color(ROOT, "#C83E4D")
    flash_color = _root_color(ROOT, "#159D66")
    measured_color = _root_color(ROOT, "#455A64")
    design_color = _root_color(ROOT, "#E07A1F")
    gray = _root_color(ROOT, "#6B7280")

    canvas = ROOT.TCanvas(
        "canvas_cpad_occupancy_landscape",
        "CPAD occupancy landscape",
        1600,
        900,
    )
    canvas.SetLeftMargin(0.105)
    canvas.SetRightMargin(0.035)
    canvas.SetTopMargin(0.11)
    canvas.SetBottomMargin(0.17)
    canvas.SetGridy()
    canvas.SetTicks(1, 1)

    frame = ROOT.TH1D(
        "frame_cpad_occupancy",
        ";Physical or effective SPAD pitch [#mum];Accepted / incident photons",
        103,
        0.0,
        103.0,
    )
    frame.SetDirectory(0)
    frame.SetMinimum(0.45)
    frame.SetMaximum(1.025)
    frame.GetXaxis().SetTitleSize(0.048)
    frame.GetYaxis().SetTitleSize(0.048)
    frame.GetXaxis().SetLabelSize(0.038)
    frame.GetYaxis().SetLabelSize(0.038)
    frame.GetYaxis().SetTitleOffset(0.95)
    frame.GetXaxis().SetTitleOffset(1.12)
    frame.Draw("AXIS")

    flash = benchmarks["study_context"]["flash"]
    flash_pitch = float(flash["pitch_um"])
    retained_objects: list[object] = [frame]
    flash_region = ROOT.TBox(0.0, 0.45, flash_pitch, 1.025)
    flash_region.SetFillColorAlpha(flash_color, 0.11)
    flash_region.SetLineColor(0)
    flash_region.Draw("SAME")
    retained_objects.append(flash_region)

    flash_line = ROOT.TLine(flash_pitch, 0.45, flash_pitch, 1.025)
    flash_line.SetLineColor(flash_color)
    flash_line.SetLineWidth(3)
    flash_line.SetLineStyle(7)
    flash_line.Draw("SAME")
    retained_objects.append(flash_line)
    retention_line = ROOT.TLine(0.0, 0.90, 103.0, 0.90)
    retention_line.SetLineColor(gray)
    retention_line.SetLineStyle(3)
    retention_line.Draw("SAME")
    retained_objects.append(retention_line)

    electron_graph = ROOT.TGraph(int(pitch.size), pitch, electron)
    electron_graph.SetName("cpad_historical_electron")
    electron_graph.SetLineColor(electron_color)
    electron_graph.SetMarkerColor(electron_color)
    electron_graph.SetLineWidth(4)
    electron_graph.SetMarkerStyle(20)
    electron_graph.SetMarkerSize(1.25)
    electron_graph.Draw("LP SAME")
    pion_graph = ROOT.TGraph(int(pitch.size), pitch, pion)
    pion_graph.SetName("cpad_historical_pion")
    pion_graph.SetLineColor(pion_color)
    pion_graph.SetMarkerColor(pion_color)
    pion_graph.SetLineWidth(4)
    pion_graph.SetMarkerStyle(20)
    pion_graph.SetMarkerSize(1.25)
    pion_graph.Draw("LP SAME")
    retained_objects.extend([electron_graph, pion_graph])

    flash_row = next(row for row in rows if row["role"] == "FLASH")
    flash_x = np.asarray([flash_pitch, flash_pitch], dtype=np.float64)
    flash_y = np.asarray(
        [
            float(flash_row["electron_pitch_only_efficiency"]),
            float(flash_row["pion_pitch_only_efficiency"]),
        ],
        dtype=np.float64,
    )
    flash_graph = ROOT.TGraph(2, flash_x, flash_y)
    flash_graph.SetName("cpad_flash_study_point")
    flash_graph.SetMarkerStyle(29)
    flash_graph.SetMarkerSize(2.5)
    flash_graph.SetMarkerColor(flash_color)
    flash_graph.Draw("P SAME")
    retained_objects.append(flash_graph)

    competitor_rows: dict[str, dict[str, object]] = {}
    for row in rows:
        if row["role"] == "competitor":
            competitor_rows.setdefault(str(row["name"]), row)
    label_positions = {
        "litavis_2026": (12.5, 0.952),
        "aspides_2025": (27.0, 0.825),
        "desy_measured": (37.0, 0.790),
        "digilog_prompt": (51.5, 0.705),
        "philips_dpc3200": (63.0, 0.625),
        "pdc_2026": (79.0, 0.515),
    }
    for name, row in competitor_rows.items():
        x_value = float(row["pitch_um"])
        y_value = float(row["pion_pitch_only_efficiency"])
        is_design = "design" in str(row["maturity"])
        marker = ROOT.TMarker(x_value, y_value, 22 if is_design else 33)
        marker.SetMarkerColor(design_color if is_design else measured_color)
        marker.SetMarkerSize(1.65)
        marker.Draw("SAME")
        retained_objects.append(marker)
        label_x, label_y = label_positions.get(name, (x_value + 2.0, y_value - 0.04))
        connector = ROOT.TLine(x_value, y_value, label_x, label_y + 0.014)
        connector.SetLineColor(_root_color(ROOT, "#9AA0A6"))
        connector.SetLineWidth(1)
        connector.Draw("SAME")
        label = ROOT.TLatex(
            label_x,
            label_y,
            f"#splitline{{{row['label']}}}{{{x_value:g} #mum}}",
        )
        label.SetTextFont(42)
        label.SetTextSize(0.031)
        label.SetTextColor(_root_color(ROOT, "#263238"))
        label.Draw("SAME")
        retained_objects.extend([connector, label])

    text = ROOT.TLatex()
    text.SetTextFont(42)
    text.SetTextColor(_root_color(ROOT, "#0B6F49"))
    text.SetTextSize(0.035)
    text.DrawLatex(2.5, 0.655, "#splitline{FLASH study region}{#leq 25 #mum physical pitch}")
    retained_objects.append(text)
    retention_text = ROOT.TLatex(98.0, 0.907, "90% retention")
    retention_text.SetTextFont(42)
    retention_text.SetTextSize(0.029)
    retention_text.SetTextColor(gray)
    retention_text.SetTextAlign(31)
    retention_text.Draw("SAME")
    retained_objects.append(retention_text)

    title = ROOT.TLatex()
    title.SetNDC(True)
    title.SetTextFont(62)
    title.SetTextSize(0.040)
    title.DrawLatex(
        0.105,
        0.955,
        "Old HG-DREAM occupancy study translated to public dSiPM pitches",
    )
    retained_objects.append(title)

    legend = ROOT.TLegend(0.105, 0.18, 0.66, 0.235)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.029)
    legend.SetNColumns(3)
    legend.AddEntry(electron_graph, "100 GeV positron", "lp")
    legend.AddEntry(pion_graph, "1 TeV pion", "lp")
    legend.AddEntry(flash_graph, "FLASH 25 #mum study point", "p")
    legend.Draw()
    retained_objects.append(legend)

    footer = ROOT.TLatex()
    footer.SetNDC(True)
    footer.SetTextFont(42)
    footer.SetTextSize(0.023)
    footer.SetTextColor(_root_color(ROOT, "#4B5563"))
    footer.DrawLatex(
        0.105,
        0.045,
        "Diamond: measured/public device; orange triangle: design stage. Pitch-only projection; finite recovery, PDE, grouping, and readout are excluded.",
    )
    retained_objects.append(footer)
    canvas.RedrawAxis()
    destinations = _save_canvas(
        canvas, output_root, "deadtime_efficiency_electron_vs_pion_cpad"
    )
    return destinations


def _timing_marker(metric: str) -> int:
    metric_lower = metric.lower()
    if "tdc" in metric_lower:
        return 21
    if "mip" in metric_lower:
        return 22
    return 20


def _timing_is_design(status: str) -> bool:
    return "design" in status or "target" in status


def _draw_timing_landscape(
    output_root: Path,
    benchmarks: dict[str, Any],
) -> list[Path]:
    ROOT = _root()
    competitors = benchmarks["competitors"]
    labels = ["FLASH"] + [device["short_label"] for device in competitors]
    n_rows = len(labels)
    measured_color = _root_color(ROOT, "#2F6690")
    design_color = _root_color(ROOT, "#E07A1F")
    flash_color = _root_color(ROOT, "#159D66")
    gray = _root_color(ROOT, "#B8C0C8")

    canvas = ROOT.TCanvas(
        "canvas_cpad_timing_landscape",
        "CPAD timing landscape",
        1600,
        900,
    )
    canvas.SetLogx()
    canvas.SetLeftMargin(0.15)
    canvas.SetRightMargin(0.04)
    canvas.SetTopMargin(0.11)
    canvas.SetBottomMargin(0.21)
    canvas.SetGridx()
    canvas.SetTicks(1, 1)
    frame = ROOT.TH2D(
        "frame_cpad_timing",
        ";Published timing scale [ps];",
        200,
        8.0,
        1800.0,
        n_rows,
        -0.5,
        n_rows - 0.5,
    )
    frame.SetDirectory(0)
    frame.GetXaxis().SetTitleSize(0.048)
    frame.GetXaxis().SetLabelSize(0.038)
    frame.GetXaxis().SetTitleOffset(1.12)
    frame.GetYaxis().SetLabelSize(0.041)
    frame.GetYaxis().SetTickLength(0.0)
    for index, label in enumerate(labels):
        y_value = n_rows - 1 - index
        y_bin = frame.GetYaxis().FindBin(y_value)
        frame.GetYaxis().SetBinLabel(y_bin, label)
    frame.Draw("AXIS")
    retained_objects: list[object] = [frame]

    flash = benchmarks["study_context"]["flash"]
    timing_low = float(flash["timing_low_ps"])
    timing_high = float(flash["timing_high_ps"])
    flash_y = float(n_rows - 1)
    flash_box = ROOT.TBox(timing_low, flash_y - 0.23, timing_high, flash_y + 0.23)
    flash_box.SetFillColorAlpha(flash_color, 0.86)
    flash_box.SetLineColor(0)
    flash_box.Draw("SAME")
    flash_text = ROOT.TLatex(
        np.sqrt(timing_low * timing_high), flash_y, "10-50 ps target"
    )
    flash_text.SetTextFont(62)
    flash_text.SetTextSize(0.030)
    flash_text.SetTextColor(ROOT.kWhite)
    flash_text.SetTextAlign(22)
    flash_text.Draw("SAME")
    retained_objects.extend([flash_box, flash_text])

    for row_index, device in enumerate(competitors, start=1):
        y_value = float(n_rows - 1 - row_index)
        values = [float(item["value_ps"]) for item in device["timing"]]
        if len(values) > 1:
            connector = ROOT.TLine(min(values), y_value, max(values), y_value)
            connector.SetLineColor(gray)
            connector.SetLineWidth(2)
            connector.Draw("SAME")
            retained_objects.append(connector)
        for timing in device["timing"]:
            value = float(timing["value_ps"])
            marker_style = _timing_marker(str(timing["metric"]))
            color = (
                design_color
                if _timing_is_design(str(timing["status"]))
                else measured_color
            )
            marker = ROOT.TMarker(value, y_value, marker_style)
            marker.SetMarkerColor(color)
            marker.SetMarkerSize(1.55)
            marker.Draw("SAME")
            value_text = ROOT.TLatex(
                value * 1.10,
                y_value + (-0.17 if marker_style == 21 else 0.14),
                f"{value:g} ps",
            )
            value_text.SetTextFont(42)
            value_text.SetTextSize(0.029)
            value_text.SetTextColor(_root_color(ROOT, "#263238"))
            value_text.Draw("SAME")
            retained_objects.extend([marker, value_text])

    title = ROOT.TLatex()
    title.SetNDC(True)
    title.SetTextFont(62)
    title.SetTextSize(0.040)
    title.DrawLatex(0.15, 0.955, "FLASH timing target and public device benchmarks")
    retained_objects.append(title)

    footer = ROOT.TLatex()
    footer.SetNDC(True)
    footer.SetTextFont(42)
    footer.SetTextSize(0.023)
    footer.SetTextColor(_root_color(ROOT, "#4B5563"))
    footer.DrawLatex(
        0.15,
        0.073,
        "Circle: photon/sensor response FWHM; square: TDC bin/resolution; triangle: MIP system timing; orange: design-stage value.",
    )
    footer.DrawLatex(
        0.15,
        0.038,
        "TDC binning, SPTR/IRF, and full-system time resolution are not interchangeable; Litavis shows the distinction at 450 nm.",
    )
    retained_objects.append(footer)
    canvas.RedrawAxis()
    destinations = _save_canvas(canvas, output_root, "cpad_timing_landscape")
    return destinations


def draw_cpad_landscape(
    output_root: str | Path,
    benchmark_path: str | Path,
) -> list[Path]:
    """Build the CPAD occupancy/timing plots and their source table."""

    output_root = Path(output_root).expanduser().resolve()
    benchmark_path = Path(benchmark_path).expanduser().resolve()
    with benchmark_path.open() as handle:
        benchmarks = json.load(handle)
    pitch, electron, pion = _load_historical_curves(output_root)
    rows = _comparison_rows(benchmarks, pitch, electron, pion)
    outputs = [_write_comparison_csv(output_root, rows)]
    outputs.extend(
        _draw_occupancy_landscape(
            output_root, benchmarks, rows, pitch, electron, pion
        )
    )
    outputs.extend(_draw_timing_landscape(output_root, benchmarks))
    return outputs
