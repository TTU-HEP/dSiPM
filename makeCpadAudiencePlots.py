#!/usr/bin/env python3
"""Make audience-facing HG-DREAM photon-image and arrival-profile plots.

The plots deliberately stay on the photon side of the simulation.  The time
curves are photon-arrival distributions after an idealized SPAD holdoff, not
analog voltage waveforms.  This keeps the comparison tied to the existing
DREAMSim truth while making the next electronics-convolution step explicit.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from dsipm.digitizer import _apply_nonparalyzable_dead_time, _cell_indices
from dsipm.root_io import inspect_root_file, iter_root_events


REPO_DIR = Path(__file__).resolve().parent
DEFAULT_SAMPLES = {
    "100 GeV positron": REPO_DIR / "mc_dreamsim_e+_100_101.root",
    "1 TeV pion": REPO_DIR / "mc_dreamsim_pi+_1000_1001.root",
}
RECOVERY_TIMES_NS = (1.0, 4.3, 40.0)
PITCH_UM = 25.0
COORDINATE_SCALE = 3.0
TIME_RANGE_NS = (5.0, 40.0)
TIME_BINS = 700
XY_RANGE_MM = (-0.5, 0.5)
XY_BINS = 120


def _fill_1d(histogram, values: np.ndarray) -> None:
    if values.size == 0:
        return
    values = np.ascontiguousarray(values, dtype=np.float64)
    weights = np.ones(values.size, dtype=np.float64)
    histogram.FillN(int(values.size), values, weights)


def _fill_2d(histogram, x: np.ndarray, y: np.ndarray) -> None:
    if x.size == 0:
        return
    x = np.ascontiguousarray(x, dtype=np.float64)
    y = np.ascontiguousarray(y, dtype=np.float64)
    weights = np.ones(x.size, dtype=np.float64)
    histogram.FillN(int(x.size), x, y, weights)


def _safe_name(label: str) -> str:
    return label.lower().replace(" ", "_").replace("+", "plus")


def _make_histograms(ROOT, label: str) -> dict[str, object]:
    stem = _safe_name(label)
    xy_low, xy_high = XY_RANGE_MM
    time_low, time_high = TIME_RANGE_NS
    histograms: dict[str, object] = {
        "xy": ROOT.TH2D(
            f"incident_xy_{stem}",
            f"{label};Sensor x [mm];Sensor y [mm];Mean photons / event / map bin",
            XY_BINS,
            xy_low,
            xy_high,
            XY_BINS,
            xy_low,
            xy_high,
        ),
        "event_xy": ROOT.TH2D(
            f"representative_event_xy_{stem}",
            f"{label};Sensor x [mm];Sensor y [mm];Photons / map bin",
            XY_BINS,
            xy_low,
            xy_high,
            XY_BINS,
            xy_low,
            xy_high,
        ),
        "incident_time": ROOT.TH1D(
            f"incident_time_{stem}",
            f"{label};Photon arrival time [ns];Fraction of incident photons / 50 ps",
            TIME_BINS,
            time_low,
            time_high,
        ),
    }
    for recovery_ns in RECOVERY_TIMES_NS:
        tag = str(recovery_ns).replace(".", "p")
        histograms[f"accepted_{recovery_ns:g}"] = ROOT.TH1D(
            f"accepted_time_{tag}ns_{stem}",
            f"{label};Photon arrival time [ns];Fraction of incident photons / 50 ps",
            TIME_BINS,
            time_low,
            time_high,
        )
    for histogram in histograms.values():
        histogram.SetDirectory(0)
        histogram.SetStats(0)
    return histograms


def _accumulate_sample(
    ROOT,
    label: str,
    path: Path,
    *,
    max_events: int | None,
) -> tuple[dict[str, object], dict[str, object]]:
    histograms = _make_histograms(ROOT, label)
    accepted_counts = {recovery: 0 for recovery in RECOVERY_TIMES_NS}
    event_indices: list[int] = []
    event_counts: list[int] = []
    events = 0
    incident = 0

    source = iter_root_events(
        path,
        max_events=max_events,
        geometry_mode="old_2025",
        gate_ns=(0.0, 40.0),
        detector_min_mm=-100.0,
        detector_max_mm=100.0,
        fiber_index=0,
    )
    for event in source:
        event = event.sorted_by_time()
        events += 1
        incident += event.selected_photon_count
        event_indices.append(event.event_index)
        event_counts.append(event.selected_photon_count)
        _fill_2d(histograms["xy"], event.x_mm, event.y_mm)
        _fill_1d(histograms["incident_time"], event.time_ns)

        _, _, cell_ids, _ = _cell_indices(
            event.x_mm * COORDINATE_SCALE,
            event.y_mm * COORDINATE_SCALE,
            PITCH_UM / 1000.0,
            -100.0 * COORDINATE_SCALE,
            100.0 * COORDINATE_SCALE,
        )
        for recovery_ns in RECOVERY_TIMES_NS:
            accepted = _apply_nonparalyzable_dead_time(
                cell_ids,
                event.time_ns,
                recovery_ns,
            )
            accepted_counts[recovery_ns] += int(np.count_nonzero(accepted))
            _fill_1d(
                histograms[f"accepted_{recovery_ns:g}"],
                event.time_ns[accepted],
            )

    representative: dict[str, object] | None = None
    if event_counts:
        counts = np.asarray(event_counts, dtype=np.int64)
        median_count = float(np.median(counts))
        representative_position = int(np.argmin(np.abs(counts - median_count)))
        representative_index = int(event_indices[representative_position])
        representative_count = int(counts[representative_position])

        representative_source = iter_root_events(
            path,
            max_events=max_events,
            geometry_mode="old_2025",
            gate_ns=(0.0, 40.0),
            detector_min_mm=-100.0,
            detector_max_mm=100.0,
            fiber_index=0,
        )
        for event in representative_source:
            if event.event_index != representative_index:
                continue
            _fill_2d(histograms["event_xy"], event.x_mm, event.y_mm)
            break

        representative = {
            "selection": (
                "Earliest event with selected-photon count nearest the sample median"
            ),
            "event_index": representative_index,
            "selected_photons": representative_count,
            "sample_median_selected_photons": median_count,
        }

    summary = {
        "sample": label,
        "events": events,
        "incident_photons": incident,
        "mean_incident_photons_per_event": incident / events if events else 0.0,
        "representative_event": representative,
        "accepted": {
            f"{recovery:g}_ns": {
                "photons": accepted_counts[recovery],
                "fraction": accepted_counts[recovery] / incident if incident else 0.0,
            }
            for recovery in RECOVERY_TIMES_NS
        },
    }
    return histograms, summary


def _style_axis(histogram) -> None:
    for axis in (histogram.GetXaxis(), histogram.GetYaxis()):
        axis.SetTitleFont(42)
        axis.SetLabelFont(42)
        axis.SetTitleSize(0.047)
        axis.SetLabelSize(0.038)
    histogram.GetXaxis().SetTitleOffset(1.08)
    histogram.GetYaxis().SetTitleOffset(1.18)


def _draw_xy_pair(
    ROOT,
    output_dir: Path,
    products,
    summaries,
    *,
    histogram_key: str,
    output_stem: str,
    event_level: bool,
) -> list[Path]:
    ROOT.gStyle.SetPalette(ROOT.kViridis)
    ROOT.gStyle.SetNumberContours(255)
    canvas = ROOT.TCanvas(
        f"cpad_{output_stem}",
        "Photon images",
        1600,
        900,
    )
    retained: list[object] = []

    display = []
    for label in DEFAULT_SAMPLES:
        plot = products[label][histogram_key].Clone(
            f"display_{histogram_key}_{_safe_name(label)}"
        )
        plot.SetDirectory(0)
        events = int(summaries[label]["events"])
        if not event_level and events:
            plot.Scale(1.0 / events)
        display.append(plot)
    shared_max = max(plot.GetMaximum() for plot in display) * 1.08
    shared_min = (
        0.8 if event_level else max(shared_max * 2.0e-4, 1.0e-5)
    )

    for column, (label, plot) in enumerate(zip(DEFAULT_SAMPLES, display)):
        canvas.cd()
        plot_pad = ROOT.TPad(
            f"{output_stem}_pad_{column}",
            label,
            0.5 * column,
            0.0,
            0.5 * (column + 1),
            0.83,
        )
        plot_pad.SetLeftMargin(0.13)
        plot_pad.SetRightMargin(0.18)
        plot_pad.SetTopMargin(0.04)
        plot_pad.SetBottomMargin(0.20)
        plot_pad.SetTicks(1, 1)
        plot_pad.SetLogz(True)
        plot_pad.Draw()
        plot_pad.cd()
        z_title = (
            "Photons / map bin"
            if event_level
            else "Mean photons / event / map bin"
        )
        plot.SetTitle(f";Sensor x [mm];Sensor y [mm];{z_title}")
        plot.SetMinimum(shared_min)
        plot.SetMaximum(shared_max)
        _style_axis(plot)
        plot.GetZaxis().SetTitleFont(42)
        plot.GetZaxis().SetLabelFont(42)
        plot.GetZaxis().SetTitleSize(0.040)
        plot.GetZaxis().SetLabelSize(0.033)
        plot.GetZaxis().SetTitleOffset(1.18)
        plot.Draw("COLZ")

        fiber = ROOT.TEllipse(0.0, 0.0, 0.40, 0.40)
        fiber.SetFillStyle(0)
        fiber.SetLineColor(
            ROOT.TColor.GetColor("#4B5563") if event_level else ROOT.kWhite
        )
        fiber.SetLineStyle(2)
        fiber.SetLineWidth(3)
        fiber.Draw("SAME")

        retained.extend([plot_pad, plot, fiber])

    canvas.cd()
    for column, label in enumerate(DEFAULT_SAMPLES):
        title = ROOT.TLatex()
        title.SetNDC(True)
        title.SetTextFont(62)
        title.SetTextSize(0.040)
        title_text = title.DrawLatex(0.065 + 0.5 * column, 0.945, label)
        subtitle = ROOT.TLatex()
        subtitle.SetNDC(True)
        subtitle.SetTextFont(42)
        subtitle.SetTextSize(0.026)
        subtitle.SetTextColor(ROOT.TColor.GetColor("#4B5563"))
        if event_level:
            representative = summaries[label]["representative_event"]
            subtitle_text = (
                f"Event {representative['event_index']}; "
                f"{representative['selected_photons']} selected photons "
                "(nearest median)"
            )
        else:
            mean = float(summaries[label]["mean_incident_photons_per_event"])
            events = int(summaries[label]["events"])
            subtitle_text = (
                f"{mean:.1f} selected photons/event; mean of {events:,} events"
            )
        subtitle_text_object = subtitle.DrawLatex(
            0.065 + 0.5 * column,
            0.895,
            subtitle_text,
        )
        retained.extend(
            [
                title,
                title_text,
                subtitle,
                subtitle_text_object,
            ]
        )

    outputs = [
        output_dir / f"{output_stem}.png",
        output_dir / f"{output_stem}.pdf",
    ]
    for path in outputs:
        canvas.SaveAs(str(path))
    return outputs


def _draw_xy(ROOT, output_dir: Path, products, summaries) -> list[Path]:
    outputs = _draw_xy_pair(
        ROOT,
        output_dir,
        products,
        summaries,
        histogram_key="xy",
        output_stem="photon_xy_comparison",
        event_level=False,
    )
    outputs.extend(
        _draw_xy_pair(
            ROOT,
            output_dir,
            products,
            summaries,
            histogram_key="event_xy",
            output_stem="photon_xy_single_event",
            event_level=True,
        )
    )
    return outputs


def _normalized_time_histogram(source, name: str, incident_total: float):
    plot = source.Clone(name)
    plot.SetDirectory(0)
    if incident_total:
        plot.Scale(1.0 / incident_total)
    plot.SetStats(0)
    return plot


def _draw_time_profiles(ROOT, output_dir: Path, products, summaries) -> list[Path]:
    colors = {
        "incident": ROOT.TColor.GetColor("#18212B"),
        1.0: ROOT.TColor.GetColor("#159D66"),
        4.3: ROOT.TColor.GetColor("#E07A1F"),
        40.0: ROOT.TColor.GetColor("#CC0000"),
    }
    canvas = ROOT.TCanvas(
        "cpad_photon_time_profiles",
        "Photon arrival profiles",
        1600,
        820,
    )
    retained: list[object] = []
    display: dict[str, dict[str | float, object]] = {}
    shared_max = 0.0

    for label in DEFAULT_SAMPLES:
        incident_total = float(summaries[label]["incident_photons"])
        sample_display: dict[str | float, object] = {
            "incident": _normalized_time_histogram(
                products[label]["incident_time"],
                f"display_incident_{_safe_name(label)}",
                incident_total,
            )
        }
        for recovery_ns in RECOVERY_TIMES_NS:
            sample_display[recovery_ns] = _normalized_time_histogram(
                products[label][f"accepted_{recovery_ns:g}"],
                f"display_{recovery_ns:g}ns_{_safe_name(label)}",
                incident_total,
            )
        display[label] = sample_display
        shared_max = max(shared_max, sample_display["incident"].GetMaximum())

    for column, label in enumerate(DEFAULT_SAMPLES):
        canvas.cd()
        x_low = 0.5 * column
        x_high = x_low + 0.5
        upper = ROOT.TPad(
            f"time_upper_{column}",
            label,
            x_low,
            0.31,
            x_high,
            1.0,
        )
        lower = ROOT.TPad(
            f"time_lower_{column}",
            label,
            x_low,
            0.0,
            x_high,
            0.31,
        )
        for pad in (upper, lower):
            pad.SetLeftMargin(0.13)
            pad.SetRightMargin(0.035)
            pad.SetTicks(1, 1)
            pad.Draw()
        upper.SetTopMargin(0.20)
        upper.SetBottomMargin(0.02)
        lower.SetTopMargin(0.03)
        lower.SetBottomMargin(0.32)
        upper.cd()

        curves = display[label]
        incident = curves["incident"]
        incident.SetTitle(
            ";Photon arrival time [ns];Fraction of incident photons / 50 ps"
        )
        incident.SetMinimum(0.0)
        incident.SetMaximum(shared_max * 1.18)
        incident.SetLineColor(colors["incident"])
        incident.SetLineWidth(4)
        _style_axis(incident)
        incident.GetXaxis().SetTitleSize(0.0)
        incident.GetXaxis().SetLabelSize(0.0)
        incident.Draw("HIST")

        legend = ROOT.TLegend(0.55, 0.46, 0.95, 0.77)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextFont(42)
        legend.SetTextSize(0.030)
        legend.AddEntry(incident, "Incident photons", "l")
        for recovery_ns in RECOVERY_TIMES_NS:
            curve = curves[recovery_ns]
            curve.SetLineColor(colors[recovery_ns])
            curve.SetLineWidth(3)
            if recovery_ns == 40.0:
                curve.SetLineStyle(2)
            curve.Draw("HIST SAME")
            if recovery_ns == 40.0:
                text = "One hit / event (old study)"
            elif recovery_ns == 4.3:
                text = "4.3 ns recovery (Litavis)"
            else:
                text = "1 ns recovery target"
            legend.AddEntry(curve, text, "l")
        legend.Draw()

        lower.cd()
        ratio_curves = []
        for recovery_ns in RECOVERY_TIMES_NS:
            ratio = curves[recovery_ns].Clone(
                f"ratio_{recovery_ns:g}ns_{_safe_name(label)}"
            )
            ratio.SetDirectory(0)
            ratio.Divide(curves["incident"])
            ratio.SetTitle(";Photon arrival time [ns];Accepted / incident")
            ratio.GetXaxis().SetRangeUser(8.5, 18.0)
            ratio.SetMinimum(0.74)
            ratio.SetMaximum(1.015)
            ratio.SetLineColor(colors[recovery_ns])
            ratio.SetLineWidth(3)
            if recovery_ns == 40.0:
                ratio.SetLineStyle(2)
            ratio.GetXaxis().SetTitleFont(42)
            ratio.GetYaxis().SetTitleFont(42)
            ratio.GetXaxis().SetLabelFont(42)
            ratio.GetYaxis().SetLabelFont(42)
            ratio.GetXaxis().SetTitleSize(0.11)
            ratio.GetYaxis().SetTitleSize(0.085)
            ratio.GetXaxis().SetLabelSize(0.085)
            ratio.GetYaxis().SetLabelSize(0.075)
            ratio.GetXaxis().SetTitleOffset(1.05)
            ratio.GetYaxis().SetTitleOffset(0.65)
            ratio.GetYaxis().SetNdivisions(505)
            ratio.Draw("HIST" if not ratio_curves else "HIST SAME")
            ratio_curves.append(ratio)

        unity = ROOT.TLine(8.5, 1.0, 18.0, 1.0)
        unity.SetLineColor(ROOT.TColor.GetColor("#7A838C"))
        unity.SetLineStyle(3)
        unity.Draw("SAME")
        retained.extend(
            [upper, lower, *curves.values(), *ratio_curves, unity, legend]
        )

    canvas.cd()
    for column, label in enumerate(DEFAULT_SAMPLES):
        title = ROOT.TLatex()
        title.SetNDC(True)
        title.SetTextFont(62)
        title.SetTextSize(0.038)
        title_text = title.DrawLatex(0.065 + 0.5 * column, 0.965, label)
        subtitle = ROOT.TLatex()
        subtitle.SetNDC(True)
        subtitle.SetTextFont(42)
        subtitle.SetTextSize(0.024)
        subtitle.SetTextColor(ROOT.TColor.GetColor("#4B5563"))
        subtitle_text = subtitle.DrawLatex(
            0.065 + 0.5 * column,
            0.925,
            "25 #mum pitch, unity PDE; arrival profile, not analog voltage",
        )
        retained.extend([title, title_text, subtitle, subtitle_text])

    outputs = [
        output_dir / "photon_time_profiles.png",
        output_dir / "photon_time_profiles.pdf",
    ]
    for path in outputs:
        canvas.SaveAs(str(path))
    return outputs


def make_plots(output_dir: Path, *, max_events: int | None) -> list[Path]:
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    output_dir = output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    products: dict[str, dict[str, object]] = {}
    summaries: dict[str, dict[str, object]] = {}
    inputs: dict[str, dict[str, object]] = {}
    for label, path in DEFAULT_SAMPLES.items():
        inputs[label] = inspect_root_file(path)
        products[label], summaries[label] = _accumulate_sample(
            ROOT,
            label,
            path,
            max_events=max_events,
        )

    root_path = output_dir / "audience_plot_histograms.root"
    root_file = ROOT.TFile(str(root_path), "RECREATE")
    for sample_products in products.values():
        for histogram in sample_products.values():
            histogram.Write()
    root_file.Close()

    summary_path = output_dir / "audience_plot_summary.json"
    with summary_path.open("w") as handle:
        json.dump(
            {
                "purpose": "Audience-facing spatial and temporal photon comparisons for CPAD 2026.",
                "inputs": inputs,
                "selection": {
                    "geometry_mode": "old_2025",
                    "fiber_index": 0,
                    "gate_ns": [0.0, 40.0],
                    "time_histogram_ns": list(TIME_RANGE_NS),
                    "time_bin_ps": 50.0,
                    "xy_display_mm": list(XY_RANGE_MM),
                    "xy_bins": XY_BINS,
                    "single_event_selection": (
                        "Earliest event with selected-photon count nearest each "
                        "sample's median"
                    ),
                    "pitch_um": PITCH_UM,
                    "coordinate_scale": COORDINATE_SCALE,
                    "pde": "unity",
                    "recovery_ns": list(RECOVERY_TIMES_NS),
                },
                "samples": summaries,
                "interpretation": (
                    "Time curves are photon-arrival distributions after a "
                    "non-paralyzable SPAD holdoff; they are not analog voltage waveforms."
                ),
            },
            handle,
            indent=2,
        )
        handle.write("\n")

    outputs = [root_path, summary_path]
    outputs.extend(_draw_xy(ROOT, output_dir, products, summaries))
    outputs.extend(_draw_time_profiles(ROOT, output_dir, products, summaries))
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_DIR / "results" / "cpad2026-audience",
    )
    parser.add_argument(
        "--max-events",
        type=int,
        default=None,
        help="Optional event limit per sample; default uses all events.",
    )
    args = parser.parse_args()
    outputs = make_plots(args.output, max_events=args.max_events)
    print("Wrote:")
    for path in outputs:
        print(f"  {path}")


if __name__ == "__main__":
    main()
