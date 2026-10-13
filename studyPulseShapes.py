#!/usr/bin/env python3
"""Compare one-fiber analog-SiPM pulses with dSiPM hit streams.

All response curves begin from the same DREAMSim Cherenkov photons and the
same deterministic 20% PDE selection.  The public-device curves isolate pitch
and recovery/holdoff; unknown wavelength-dependent PDE and output-policy
effects are intentionally not folded in.  dSiPM curves labelled
"waveform-equivalent" are time-tagged hits convolved with the same display
kernel as the conventional SiPM.  They are not predicted analog voltages.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np

from dsipm.digitizer import _apply_nonparalyzable_dead_time, _cell_indices
from dsipm.pulse_shape import (
    histogram_observables,
    partial_recharge_weights,
    scope_waveform,
)
from dsipm.root_io import iter_root_events


REPO = Path(__file__).resolve().parent
DEFAULT_INPUT = REPO / "mc_dreamsim_pi+_1000_1001.root"
SENSOR_MIN_MM = -1.5
SENSOR_MAX_MM = 1.5
COMMON_PDE = 0.20
SEED = 4307
TIME_LOW_NS = 5.0
TIME_HIGH_NS = 40.0
BIN_WIDTH_NS = 0.05
TIME_EDGES_NS = np.arange(
    TIME_LOW_NS,
    TIME_HIGH_NS + 0.5 * BIN_WIDTH_NS,
    BIN_WIDTH_NS,
    dtype=np.float64,
)
PROMPT_END_NS = 12.0
ANALOG_MICROCELLS = 4774
ANALOG_GRID_SIDE = int(round(np.sqrt(ANALOG_MICROCELLS)))
ANALOG_EFFECTIVE_PITCH_UM = 3000.0 / ANALOG_GRID_SIDE


@dataclass(frozen=True)
class PulseScenario:
    key: str
    label: str
    pitch_um: float
    recovery_ns: float
    source_status: str
    color: str
    line_style: int = 1


DIGITAL_SCENARIOS = (
    PulseScenario(
        "litavis",
        "Litavis-like 10.17 #mum, 4.3 ns",
        10.17,
        4.3,
        "fabricated pitch and published 4.3 ns minimum dead time",
        "#2B6CB0",
    ),
    PulseScenario(
        "flash",
        "FLASH study point 20 #mum, 4.3 ns",
        20.0,
        4.3,
        "FLASH simulation working point",
        "#159D66",
    ),
    PulseScenario(
        "dpc3200",
        "Philips DPC3200 61.7 #mum, 40 ns",
        61.7,
        40.0,
        "fabricated effective cell pitch; 40 ns recovery proxy",
        "#E07A1F",
    ),
    PulseScenario(
        "pdc78",
        "Sherbrooke PDC 78 #mum, 32 ns",
        78.0,
        32.0,
        "largest fabricated dSiPM-array pitch found; published minimum holdoff",
        "#C83E4D",
    ),
    PulseScenario(
        "stress200",
        "200 #mum stress point, 4.3 ns",
        200.0,
        4.3,
        "simulation stress point from the capstone scan; not a device claim",
        "#263238",
        2,
    ),
)

PITCH_SCAN_UM = np.asarray(
    [5.0, 10.0, 20.0, 25.0, 35.0, 50.0, 61.7, 78.0, 100.0, 150.0, 200.0]
)
RECOVERY_SCAN_NS = np.asarray([1.0, 4.3, 32.0, 40.0])


def _histogram(times: np.ndarray, weights: np.ndarray | None = None) -> np.ndarray:
    return np.histogram(times, bins=TIME_EDGES_NS, weights=weights)[0].astype(
        np.float64
    )


def _cell_ids(x_mm: np.ndarray, y_mm: np.ndarray, pitch_um: float) -> np.ndarray:
    return _cell_indices(
        x_mm,
        y_mm,
        pitch_um / 1000.0,
        SENSOR_MIN_MM,
        SENSOR_MAX_MM,
    )[2]


def _common_detection(event) -> np.ndarray:
    random = np.random.default_rng(
        np.random.SeedSequence([SEED, event.event_index, 3409])
    )
    return random.random(event.selected_photon_count) < COMMON_PDE


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _accumulate(path: Path, max_events: int | None) -> dict[str, object]:
    detected_hist = np.zeros(TIME_EDGES_NS.size - 1, dtype=np.float64)
    analog_hist = np.zeros_like(detected_hist)
    digital_hist = {
        scenario.key: np.zeros_like(detected_hist) for scenario in DIGITAL_SCENARIOS
    }
    scan_hist = {
        (float(pitch), float(recovery)): np.zeros_like(detected_hist)
        for pitch in PITCH_SCAN_UM
        for recovery in RECOVERY_SCAN_NS
    }
    event_products: list[dict[str, object]] = []
    incident_total = 0
    detected_total = 0

    source = iter_root_events(
        path,
        max_events=max_events,
        geometry_mode="old_2025",
        gate_ns=(TIME_LOW_NS, TIME_HIGH_NS),
        detector_min_mm=SENSOR_MIN_MM,
        detector_max_mm=SENSOR_MAX_MM,
        fiber_index=0,
    )
    for event in source:
        event = event.sorted_by_time()
        incident_total += event.selected_photon_count
        detected_mask = _common_detection(event)
        x = event.x_mm[detected_mask]
        y = event.y_mm[detected_mask]
        time = event.time_ns[detected_mask]
        detected_total += int(time.size)
        event_detected_hist = _histogram(time)
        detected_hist += event_detected_hist

        analog_ids = _cell_ids(x, y, ANALOG_EFFECTIVE_PITCH_UM)
        analog_weights = partial_recharge_weights(analog_ids, time, 82.0)
        event_analog_hist = _histogram(time, analog_weights)
        analog_hist += event_analog_hist

        event_digital: dict[str, np.ndarray] = {}
        for scenario in DIGITAL_SCENARIOS:
            ids = _cell_ids(x, y, scenario.pitch_um)
            accepted = _apply_nonparalyzable_dead_time(
                ids,
                time,
                scenario.recovery_ns,
            )
            response = _histogram(time[accepted])
            digital_hist[scenario.key] += response
            event_digital[scenario.key] = response

        for pitch in PITCH_SCAN_UM:
            ids = _cell_ids(x, y, float(pitch))
            for recovery in RECOVERY_SCAN_NS:
                accepted = _apply_nonparalyzable_dead_time(
                    ids,
                    time,
                    float(recovery),
                )
                scan_hist[(float(pitch), float(recovery))] += _histogram(
                    time[accepted]
                )

        event_products.append(
            {
                "event_index": event.event_index,
                "incident": event.selected_photon_count,
                "detected": int(time.size),
                "detected_hist": event_detected_hist,
                "analog_hist": event_analog_hist,
                "digital_hist": event_digital,
            }
        )

    if not event_products:
        raise ValueError(f"No selected events found in {path}")
    counts = np.asarray([row["detected"] for row in event_products], dtype=np.int64)
    median_count = float(np.median(counts))
    representative = min(
        event_products,
        key=lambda row: (abs(int(row["detected"]) - median_count), row["event_index"]),
    )
    return {
        "events": len(event_products),
        "incident_total": incident_total,
        "detected_total": detected_total,
        "detected_hist": detected_hist,
        "analog_hist": analog_hist,
        "digital_hist": digital_hist,
        "scan_hist": scan_hist,
        "representative": representative,
        "median_detected": median_count,
    }


def _root():
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTitleFont(42, "XYZ")
    ROOT.gStyle.SetLabelFont(42, "XYZ")
    return ROOT


def _color(ROOT, value: str) -> int:
    return int(ROOT.TColor.GetColor(value))


def _graph(ROOT, x: np.ndarray, y: np.ndarray, *, color: str, style: int = 1, width: int = 3):
    graph = ROOT.TGraph(
        int(x.size),
        np.ascontiguousarray(x, dtype=np.float64),
        np.ascontiguousarray(y, dtype=np.float64),
    )
    graph.SetLineColor(_color(ROOT, color))
    graph.SetMarkerColor(_color(ROOT, color))
    graph.SetLineStyle(style)
    graph.SetLineWidth(width)
    return graph


def _frame(ROOT, name: str, x_title: str, y_title: str, x_min: float, x_max: float, y_min: float, y_max: float):
    frame = ROOT.TH1D(name, f";{x_title};{y_title}", 100, x_min, x_max)
    frame.SetDirectory(0)
    frame.SetMinimum(y_min)
    frame.SetMaximum(y_max)
    frame.GetXaxis().SetTitleSize(0.050)
    frame.GetYaxis().SetTitleSize(0.050)
    frame.GetXaxis().SetLabelSize(0.040)
    frame.GetYaxis().SetLabelSize(0.040)
    frame.GetXaxis().SetTitleOffset(1.05)
    frame.GetYaxis().SetTitleOffset(0.88)
    frame.Draw("AXIS")
    return frame


def _save(canvas, output_dir: Path, stem: str) -> list[Path]:
    paths = [output_dir / f"{stem}.png", output_dir / f"{stem}.pdf"]
    for path in paths:
        canvas.SaveAs(str(path))
    return paths


def _draw_header(ROOT, title_text: str, subtitle_text: str, retained: list[object]) -> None:
    title = ROOT.TLatex()
    title.SetNDC(True)
    title.SetTextFont(62)
    title.SetTextSize(0.038)
    title.DrawLatex(0.08, 0.965, title_text)
    subtitle = ROOT.TLatex()
    subtitle.SetNDC(True)
    subtitle.SetTextFont(42)
    subtitle.SetTextSize(0.024)
    subtitle.SetTextColor(_color(ROOT, "#4B5563"))
    subtitle.DrawLatex(0.08, 0.927, subtitle_text)
    retained.extend([title, subtitle])


def _draw_single_event(ROOT, output_dir: Path, products: dict[str, object]) -> list[Path]:
    event = products["representative"]
    centers = 0.5 * (TIME_EDGES_NS[:-1] + TIME_EDGES_NS[1:])
    detected = np.asarray(event["detected_hist"], dtype=np.float64)
    analog = np.asarray(event["analog_hist"], dtype=np.float64)
    digital = event["digital_hist"]
    ideal_waveform = scope_waveform(detected, bin_width_ns=BIN_WIDTH_NS)
    analog_waveform = scope_waveform(analog, bin_width_ns=BIN_WIDTH_NS)
    digital_waveforms = {
        key: scope_waveform(np.asarray(value), bin_width_ns=BIN_WIDTH_NS)
        for key, value in digital.items()
    }
    scale = max(float(np.max(ideal_waveform)), 1.0)

    canvas = ROOT.TCanvas("pulse_scope_single_event", "Pulse comparison", 1500, 1050)
    retained: list[object] = []
    pads = []
    for index, (bottom, top) in enumerate(((0.66, 0.89), (0.37, 0.65), (0.08, 0.36))):
        pad = ROOT.TPad(f"scope_event_pad_{index}", "", 0.07, bottom, 0.98, top)
        pad.SetLeftMargin(0.10)
        pad.SetRightMargin(0.025)
        pad.SetTopMargin(0.03)
        pad.SetBottomMargin(0.04 if index < 2 else 0.20)
        pad.SetTicks(1, 1)
        pad.SetGridy()
        pad.Draw()
        pads.append(pad)
        retained.append(pad)

    x_min, x_max = 8.5, 25.0
    pads[0].cd()
    first = _frame(
        ROOT,
        "scope_event_detected_frame",
        "",
        "Detected / 50 ps",
        x_min,
        x_max,
        0.0,
        1.16 * max(float(np.max(detected)), 1.0),
    )
    first.GetXaxis().SetLabelSize(0.0)
    detected_graph = _graph(ROOT, centers, detected, color="#18212B", width=3)
    detected_graph.Draw("L SAME")
    label = ROOT.TLatex()
    label.SetNDC(True)
    label.SetTextFont(62)
    label.SetTextSize(0.055)
    label.DrawLatex(0.14, 0.82, "Same detected photons")
    retained.extend([first, detected_graph, label])

    pads[1].cd()
    second = _frame(
        ROOT,
        "scope_event_analog_frame",
        "",
        "Fast-output proxy [a.u.]",
        x_min,
        x_max,
        0.0,
        1.12,
    )
    second.GetXaxis().SetLabelSize(0.0)
    ideal_graph = _graph(ROOT, centers, ideal_waveform / scale, color="#8B949E", style=2, width=2)
    analog_graph = _graph(ROOT, centers, analog_waveform / scale, color="#7A5195", width=4)
    ideal_graph.Draw("L SAME")
    analog_graph.Draw("L SAME")
    legend = ROOT.TLegend(0.57, 0.58, 0.965, 0.90)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.048)
    legend.AddEntry(ideal_graph, "Fully charged reference", "l")
    legend.AddEntry(analog_graph, "Conventional SiPM: partial recharge", "l")
    legend.Draw()
    analog_label = ROOT.TLatex()
    analog_label.SetNDC(True)
    analog_label.SetTextFont(62)
    analog_label.SetTextSize(0.055)
    analog_label.DrawLatex(0.14, 0.82, "What an oscilloscope-like fast output retains")
    retained.extend([second, ideal_graph, analog_graph, legend, analog_label])

    pads[2].cd()
    third = _frame(
        ROOT,
        "scope_event_digital_frame",
        "Photon arrival time [ns]",
        "Waveform-equivalent [a.u.]",
        x_min,
        x_max,
        0.0,
        1.12,
    )
    digital_legend = ROOT.TLegend(0.48, 0.52, 0.965, 0.90)
    digital_legend.SetBorderSize(0)
    digital_legend.SetFillStyle(0)
    digital_legend.SetTextSize(0.044)
    shown = {"flash", "pdc78", "stress200"}
    for scenario in DIGITAL_SCENARIOS:
        if scenario.key not in shown:
            continue
        graph = _graph(
            ROOT,
            centers,
            digital_waveforms[scenario.key] / scale,
            color=scenario.color,
            style=scenario.line_style,
            width=4,
        )
        graph.Draw("L SAME")
        digital_legend.AddEntry(graph, scenario.label, "l")
        retained.append(graph)
    digital_legend.Draw()
    digital_label = ROOT.TLatex()
    digital_label.SetNDC(True)
    digital_label.SetTextFont(62)
    digital_label.SetTextSize(0.055)
    digital_label.DrawLatex(0.14, 0.82, "What sparse dSiPM hits can reconstruct")
    retained.extend([third, digital_legend, digital_label])

    canvas.cd()
    _draw_header(
        ROOT,
        "One fiber, one shower: from detected photons to readout pulse",
        (
            f"1 TeV pion event {event['event_index']}; {event['detected']} photons after common 20% PDE; "
            "50 ps bins"
        ),
        retained,
    )
    footer = ROOT.TLatex()
    footer.SetNDC(True)
    footer.SetTextFont(42)
    footer.SetTextSize(0.019)
    footer.SetTextColor(_color(ROOT, "#4B5563"))
    footer.DrawLatex(
        0.08,
        0.025,
        "dSiPM curves are hit timestamps convolved with one common display kernel, not predicted device voltages. 200 #mum is a stress test, not a fabricated-device claim.",
    )
    retained.append(footer)
    return _save(canvas, output_dir, "pulse_scope_single_event")


def _draw_ensemble(ROOT, output_dir: Path, products: dict[str, object]) -> list[Path]:
    centers = 0.5 * (TIME_EDGES_NS[:-1] + TIME_EDGES_NS[1:])
    events = int(products["events"])
    detected = np.asarray(products["detected_hist"], dtype=np.float64) / events
    analog = np.asarray(products["analog_hist"], dtype=np.float64) / events
    digital = {
        key: np.asarray(value, dtype=np.float64) / events
        for key, value in products["digital_hist"].items()
    }
    ideal_waveform = scope_waveform(detected, bin_width_ns=BIN_WIDTH_NS)
    scale = max(float(np.max(ideal_waveform)), 1.0)
    waveforms = {
        "analog": scope_waveform(analog, bin_width_ns=BIN_WIDTH_NS) / scale,
        **{
            key: scope_waveform(value, bin_width_ns=BIN_WIDTH_NS) / scale
            for key, value in digital.items()
        },
    }

    canvas = ROOT.TCanvas("pulse_architecture_comparison", "Pulse architectures", 1600, 900)
    upper = ROOT.TPad("pulse_architecture_upper", "", 0.07, 0.39, 0.98, 0.89)
    lower = ROOT.TPad("pulse_architecture_lower", "", 0.07, 0.08, 0.98, 0.38)
    retained: list[object] = [upper, lower]
    for pad in (upper, lower):
        pad.SetLeftMargin(0.10)
        pad.SetRightMargin(0.025)
        pad.SetTicks(1, 1)
        pad.SetGridy()
        pad.Draw()
    upper.SetBottomMargin(0.03)
    lower.SetTopMargin(0.02)
    lower.SetBottomMargin(0.22)

    x_min, x_max = 8.5, 30.0
    upper.cd()
    frame = _frame(
        ROOT,
        "pulse_architecture_frame",
        "",
        "Mean pulse / ideal peak",
        x_min,
        x_max,
        0.0,
        1.08,
    )
    frame.GetXaxis().SetLabelSize(0.0)
    ideal_graph = _graph(ROOT, centers, ideal_waveform / scale, color="#8B949E", style=2, width=2)
    ideal_graph.Draw("L SAME")
    analog_graph = _graph(ROOT, centers, waveforms["analog"], color="#7A5195", width=4)
    analog_graph.Draw("L SAME")
    legend = ROOT.TLegend(0.50, 0.45, 0.965, 0.91)
    legend.SetNColumns(2)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.034)
    legend.AddEntry(ideal_graph, "Detected-photon reference", "l")
    legend.AddEntry(analog_graph, "Conventional SiPM 35 #mum, #tau=82 ns", "l")
    retained.extend([frame, ideal_graph, analog_graph, legend])
    for scenario in DIGITAL_SCENARIOS:
        graph = _graph(
            ROOT,
            centers,
            waveforms[scenario.key],
            color=scenario.color,
            style=scenario.line_style,
            width=3,
        )
        graph.Draw("L SAME")
        legend.AddEntry(graph, scenario.label, "l")
        retained.append(graph)
    legend.Draw()

    lower.cd()
    lower_frame = _frame(
        ROOT,
        "pulse_cumulative_frame",
        "Photon arrival time [ns]",
        "Cumulative output / detected",
        x_min,
        x_max,
        0.0,
        1.03,
    )
    denominator = max(float(np.sum(detected)), 1.0)
    analog_cumulative = np.cumsum(analog) / denominator
    cumulative_analog_graph = _graph(ROOT, centers, analog_cumulative, color="#7A5195", width=4)
    cumulative_analog_graph.Draw("L SAME")
    retained.extend([lower_frame, cumulative_analog_graph])
    for scenario in DIGITAL_SCENARIOS:
        cumulative = np.cumsum(digital[scenario.key]) / denominator
        graph = _graph(
            ROOT,
            centers,
            cumulative,
            color=scenario.color,
            style=scenario.line_style,
            width=3,
        )
        graph.Draw("L SAME")
        retained.append(graph)

    canvas.cd()
    _draw_header(
        ROOT,
        "The same photons look different after analog recharge or digital holdoff",
        f"Mean of {events} one-fiber 1 TeV pion events; common 20% PDE isolates pitch and recovery",
        retained,
    )
    footer = ROOT.TLatex()
    footer.SetNDC(True)
    footer.SetTextFont(42)
    footer.SetTextSize(0.020)
    footer.SetTextColor(_color(ROOT, "#4B5563"))
    footer.DrawLatex(
        0.08,
        0.022,
        "Scope rendering uses the MicroFC-30035-SMT fast-output timing (0.6 ns rise, 1.5 ns FWHM). dSiPM PDE, jitter, grouping, and exported timestamp policy are not included here.",
    )
    retained.append(footer)
    return _save(canvas, output_dir, "pulse_architecture_comparison")


def _scan_rows(products: dict[str, object]) -> list[dict[str, object]]:
    detected = np.asarray(products["detected_hist"], dtype=np.float64)
    reference = histogram_observables(
        detected,
        TIME_EDGES_NS,
        prompt_end_ns=PROMPT_END_NS,
    )
    rows: list[dict[str, object]] = []
    for pitch in PITCH_SCAN_UM:
        for recovery in RECOVERY_SCAN_NS:
            response = np.asarray(
                products["scan_hist"][(float(pitch), float(recovery))],
                dtype=np.float64,
            )
            values = histogram_observables(
                response,
                TIME_EDGES_NS,
                prompt_end_ns=PROMPT_END_NS,
            )
            rows.append(
                {
                    "response": "digital hard holdoff",
                    "pitch_um": float(pitch),
                    "recovery_ns": float(recovery),
                    "retained_fraction": values["integral"] / reference["integral"],
                    "prompt_fraction": values["prompt_fraction"],
                    "prompt_fraction_shift_percentage_points": 100.0
                    * (values["prompt_fraction"] - reference["prompt_fraction"]),
                    "mean_time_shift_ps": 1000.0
                    * (values["mean_time_ns"] - reference["mean_time_ns"]),
                    "central80_width_shift_ps": 1000.0
                    * (
                        values["central80_width_ns"]
                        - reference["central80_width_ns"]
                    ),
                }
            )
    analog = histogram_observables(
        np.asarray(products["analog_hist"], dtype=np.float64),
        TIME_EDGES_NS,
        prompt_end_ns=PROMPT_END_NS,
    )
    rows.append(
        {
            "response": "conventional SiPM partial recharge",
            "pitch_um": ANALOG_EFFECTIVE_PITCH_UM,
            "recovery_ns": 82.0,
            "retained_fraction": analog["integral"] / reference["integral"],
            "prompt_fraction": analog["prompt_fraction"],
            "prompt_fraction_shift_percentage_points": 100.0
            * (analog["prompt_fraction"] - reference["prompt_fraction"]),
            "mean_time_shift_ps": 1000.0
            * (analog["mean_time_ns"] - reference["mean_time_ns"]),
            "central80_width_shift_ps": 1000.0
            * (analog["central80_width_ns"] - reference["central80_width_ns"]),
        }
    )
    return rows


def _draw_pitch_scan(ROOT, output_dir: Path, rows: list[dict[str, object]]) -> list[Path]:
    colors = {1.0: "#2B6CB0", 4.3: "#159D66", 32.0: "#C83E4D", 40.0: "#E07A1F"}
    line_styles = {1.0: 1, 4.3: 1, 32.0: 2, 40.0: 3}
    marker_styles = {1.0: 20, 4.3: 21, 32.0: 24, 40.0: 25}
    canvas = ROOT.TCanvas("pulse_pitch_scan", "Pulse pitch scan", 1600, 850)
    retained: list[object] = []
    pads = []
    for index in range(2):
        pad = ROOT.TPad(
            f"pulse_pitch_pad_{index}",
            "",
            0.055 + 0.475 * index,
            0.10,
            0.515 + 0.475 * index,
            0.88,
        )
        pad.SetLeftMargin(0.15)
        pad.SetRightMargin(0.04)
        pad.SetBottomMargin(0.15)
        pad.SetTopMargin(0.05)
        pad.SetLogx()
        pad.SetGridy()
        pad.SetTicks(1, 1)
        pad.Draw()
        pads.append(pad)
        retained.append(pad)

    digital_rows = [row for row in rows if row["response"] == "digital hard holdoff"]
    analog_row = next(row for row in rows if row["response"].startswith("conventional"))
    shift_values = [float(row["mean_time_shift_ps"]) for row in digital_rows]
    shift_values.append(float(analog_row["mean_time_shift_ps"]))
    shift_min = min(shift_values)
    shift_max = max(shift_values)
    shift_padding = max(8.0, 0.12 * (shift_max - shift_min))

    pads[0].cd()
    left = _frame(
        ROOT,
        "pulse_pitch_retention_frame",
        "SPAD pitch [#mum]",
        "Integrated output / detected",
        4.5,
        220.0,
        0.0,
        1.02,
    )
    retained.append(left)
    legend = ROOT.TLegend(0.17, 0.18, 0.62, 0.43)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.036)
    for recovery in RECOVERY_SCAN_NS:
        selected = [row for row in digital_rows if np.isclose(row["recovery_ns"], recovery)]
        graph = _graph(
            ROOT,
            np.asarray([row["pitch_um"] for row in selected], dtype=np.float64),
            np.asarray([row["retained_fraction"] for row in selected], dtype=np.float64),
            color=colors[float(recovery)],
            style=line_styles[float(recovery)],
            width=4,
        )
        graph.SetMarkerStyle(marker_styles[float(recovery)])
        graph.SetMarkerSize(1.0)
        graph.Draw("LP SAME")
        legend.AddEntry(graph, f"Hard holdoff {recovery:g} ns", "lp")
        retained.append(graph)
    analog_marker = ROOT.TMarker(
        float(analog_row["pitch_um"]),
        float(analog_row["retained_fraction"]),
        29,
    )
    analog_marker.SetMarkerColor(_color(ROOT, "#7A5195"))
    analog_marker.SetMarkerSize(2.3)
    analog_marker.Draw("SAME")
    legend.AddEntry(analog_marker, "Conventional SiPM partial recharge", "p")
    legend.Draw()
    retained.extend([analog_marker, legend])

    pads[1].cd()
    right = _frame(
        ROOT,
        "pulse_pitch_shift_frame",
        "SPAD pitch [#mum]",
        "Mean arrival-time shift [ps]",
        4.5,
        220.0,
        shift_min - shift_padding,
        shift_max + shift_padding,
    )
    right.GetYaxis().SetTitleOffset(1.10)
    retained.append(right)
    zero = ROOT.TLine(4.5, 0.0, 220.0, 0.0)
    zero.SetLineColor(_color(ROOT, "#8B949E"))
    zero.SetLineStyle(2)
    zero.Draw("SAME")
    retained.append(zero)
    for recovery in RECOVERY_SCAN_NS:
        selected = [row for row in digital_rows if np.isclose(row["recovery_ns"], recovery)]
        graph = _graph(
            ROOT,
            np.asarray([row["pitch_um"] for row in selected], dtype=np.float64),
            np.asarray([row["mean_time_shift_ps"] for row in selected], dtype=np.float64),
            color=colors[float(recovery)],
            style=line_styles[float(recovery)],
            width=4,
        )
        graph.SetMarkerStyle(marker_styles[float(recovery)])
        graph.SetMarkerSize(1.0)
        graph.Draw("LP SAME")
        retained.append(graph)
    analog_shift_marker = ROOT.TMarker(
        float(analog_row["pitch_um"]),
        float(analog_row["mean_time_shift_ps"]),
        29,
    )
    analog_shift_marker.SetMarkerColor(_color(ROOT, "#7A5195"))
    analog_shift_marker.SetMarkerSize(2.3)
    analog_shift_marker.Draw("SAME")
    retained.append(analog_shift_marker)

    canvas.cd()
    _draw_header(
        ROOT,
        "Pitch and recovery jointly determine whether the pulse survives",
        "1 TeV pion one-fiber sample; common 20% PDE; scan extended to 200 #mum",
        retained,
    )
    footer = ROOT.TLatex()
    footer.SetNDC(True)
    footer.SetTextFont(42)
    footer.SetTextSize(0.020)
    footer.SetTextColor(_color(ROOT, "#4B5563"))
    footer.DrawLatex(
        0.07,
        0.035,
        "78 #mum is the largest pitch of a fabricated dSiPM-class array identified here (Sherbrooke PDC). 100--200 #mum points are deliberately labelled simulation stress tests.",
    )
    retained.append(footer)
    return _save(canvas, output_dir, "pulse_pitch_recovery_scan")


def _time_rows(products: dict[str, object]) -> list[dict[str, object]]:
    centers = 0.5 * (TIME_EDGES_NS[:-1] + TIME_EDGES_NS[1:])
    rows: list[dict[str, object]] = []
    series = {
        "detected_reference": np.asarray(products["detected_hist"]),
        "conventional_sipm": np.asarray(products["analog_hist"]),
        **{
            scenario.key: np.asarray(products["digital_hist"][scenario.key])
            for scenario in DIGITAL_SCENARIOS
        },
    }
    for name, histogram in series.items():
        waveform = scope_waveform(histogram, bin_width_ns=BIN_WIDTH_NS)
        cumulative = np.cumsum(histogram)
        for index, center in enumerate(centers):
            rows.append(
                {
                    "response": name,
                    "time_low_ns": TIME_EDGES_NS[index],
                    "time_high_ns": TIME_EDGES_NS[index + 1],
                    "time_center_ns": center,
                    "avalanche_or_charge_weight": histogram[index],
                    "scope_equivalent": waveform[index],
                    "cumulative_output": cumulative[index],
                }
            )
    return rows


def run(input_path: Path, output_dir: Path, *, max_events: int | None) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    products = _accumulate(input_path, max_events)
    ROOT = _root()
    outputs: list[Path] = []
    outputs.extend(_draw_single_event(ROOT, output_dir, products))
    outputs.extend(_draw_ensemble(ROOT, output_dir, products))
    scan_rows = _scan_rows(products)
    outputs.extend(_draw_pitch_scan(ROOT, output_dir, scan_rows))
    time_rows = _time_rows(products)
    time_csv = output_dir / "pulse_architecture_time_series.csv"
    scan_csv = output_dir / "pulse_pitch_recovery_scan.csv"
    _write_csv(time_csv, time_rows)
    _write_csv(scan_csv, scan_rows)
    outputs.extend([time_csv, scan_csv])

    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "input": str(input_path.resolve()),
        "max_events": max_events,
        "events_processed": products["events"],
        "common_pde": COMMON_PDE,
        "time_gate_ns": [TIME_LOW_NS, TIME_HIGH_NS],
        "time_bin_ns": BIN_WIDTH_NS,
        "prompt_end_ns": PROMPT_END_NS,
        "representative_event": {
            "event_index": products["representative"]["event_index"],
            "selected_photons": products["representative"]["incident"],
            "detected_photons": products["representative"]["detected"],
            "selection": "detected count nearest the sample median, then earliest event",
            "sample_median_detected_photons": products["median_detected"],
        },
        "conventional_sipm": {
            "model": "MicroFC-30035-SMT sensor-only proxy",
            "active_area_mm2": [3.0, 3.0],
            "nominal_microcell_um": 35.0,
            "microcells": ANALOG_MICROCELLS,
            "effective_square_grid_pitch_um": ANALOG_EFFECTIVE_PITCH_UM,
            "recovery_tau_ns": 82.0,
            "fast_output_rise_10_90_ns": 0.6,
            "fast_output_fwhm_ns": 1.5,
            "source": "https://www.onsemi.com/pdf/datasheet/microc-series-d.pdf",
        },
        "digital_scenarios": [scenario.__dict__ for scenario in DIGITAL_SCENARIOS],
        "interpretation": (
            "Digital curves are accepted avalanche times convolved with the common "
            "SiPM display kernel. They are not predicted analog voltages. Public-device "
            "PDE, timing jitter, readout grouping, and output policy are excluded to "
            "isolate pitch and recovery."
        ),
        "largest_pitch_statement": (
            "78 um is the largest pitch of a fabricated dSiPM-class array found in the "
            "reviewed public primary material. The 100, 150, and 200 um points are "
            "simulation stress tests; the 200 um point also reproduces the capstone range."
        ),
        "outputs": [str(path.resolve()) for path in outputs],
    }
    manifest_path = output_dir / "pulse_study_manifest.json"
    with manifest_path.open("w") as handle:
        json.dump(manifest, handle, indent=2)
    outputs.append(manifest_path)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO / "results" / "pulse-shape-v1",
    )
    parser.add_argument("--max-events", type=int, default=None)
    args = parser.parse_args()
    outputs = run(args.input, args.output, max_events=args.max_events)
    print(json.dumps([str(path) for path in outputs], indent=2))


if __name__ == "__main__":
    main()
