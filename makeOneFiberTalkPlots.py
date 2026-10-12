#!/usr/bin/env python3
"""Create the one-fiber segmentation plot used in the CPAD narrative.

The saved 2025 multiplicity histograms are exact software closure for the
original occupancy study.  Their displayed pitch was three times the physical
pitch because the 3 mm sensor segmentation was applied to a 1 mm fiber-local
coordinate window.  This plot keeps the validated histogram contents and
labels the curves with the corresponding physical pitch.
"""

from __future__ import annotations

import argparse
import csv
from array import array
from pathlib import Path

from dsipm.cms_style import (
    draw_flash_label,
    make_legend,
    save_canvas,
    set_flash_style,
    style_axes,
    style_pad,
)


REPO_DIR = Path(__file__).resolve().parent
DEFAULT_ELECTRON = REPO_DIR / "results/2025-both/100GeV_eplus/comparison.root"
DEFAULT_PION = REPO_DIR / "results/2025-both/1TeV_piplus/comparison.root"
DEFAULT_OUTPUT = REPO_DIR / "results/cpad2026-one-fiber"
LEGACY_LABELS_UM = (10, 25, 50, 75, 100)
DISPLAY_TO_PHYSICAL = 1.0 / 3.0


def _format_pitch(value: float) -> str:
    if abs(value - round(value)) < 1e-9:
        return f"{value:.0f}"
    return f"{value:.1f}"


def _occupied_distribution(ROOT, source, name: str, *, maximum: int = 16):
    """Return P(N photons | N >= 1), with the final bin inclusive."""

    plot = ROOT.TH1D(
        name,
        ";Photons in an occupied SPAD / event;Fraction of occupied SPADs / bin",
        maximum,
        0.5,
        maximum + 0.5,
    )
    plot.SetDirectory(0)
    # Source bin 1 is N=0, bin 2 is N=1, and so on.
    occupied = source.Integral(2, source.GetNbinsX() + 1)
    if occupied <= 0:
        return plot, 0.0
    for multiplicity in range(1, maximum):
        source_bin = multiplicity + 1
        content = source.GetBinContent(source_bin)
        error = source.GetBinError(source_bin)
        plot.SetBinContent(multiplicity, content / occupied)
        plot.SetBinError(multiplicity, error / occupied)
    tail = source.Integral(maximum + 1, source.GetNbinsX() + 1)
    plot.SetBinContent(maximum, tail / occupied)
    multi_hit = 1.0 - source.GetBinContent(2) / occupied
    return plot, multi_hit


def _draw_panel(ROOT, pad, root_file, sample_label: str, colors: list[int]):
    pad.cd()
    pad.SetLogy(True)
    style_pad(pad, left=0.15, right=0.035, top=0.25, bottom=0.16)

    frame = ROOT.TH1D(
        f"frame_{sample_label.replace(' ', '_')}",
        ";Photons in an occupied SPAD / event;Fraction of occupied SPADs / bin",
        16,
        0.5,
        16.5,
    )
    frame.SetDirectory(0)
    frame.SetMinimum(1.0e-6)
    frame.SetMaximum(1.7)
    frame.GetXaxis().SetNdivisions(516)
    style_axes(frame)
    frame.Draw("AXIS")

    legend = make_legend(ROOT, 0.61, 0.40, 0.94, 0.73, text_size=0.037)
    legend.SetHeader("Physical SPAD pitch", "C")
    legend.SetMargin(0.24)

    keep = [frame, legend]
    rows: list[dict[str, float | str]] = []
    for index, legacy_pitch in enumerate(LEGACY_LABELS_UM):
        source = root_file.Get(f"nPhotonsPerChannel_{legacy_pitch}x{legacy_pitch}")
        if not source:
            raise RuntimeError(f"Missing nPhotonsPerChannel_{legacy_pitch}x{legacy_pitch}")
        plot, multi_hit = _occupied_distribution(
            ROOT,
            source,
            f"plot_{sample_label}_{legacy_pitch}",
        )
        plot.SetLineColor(colors[index])
        plot.SetLineWidth(3)
        plot.Draw("HIST SAME")
        physical_pitch = legacy_pitch * DISPLAY_TO_PHYSICAL
        legend.AddEntry(plot, f"{_format_pitch(physical_pitch)} #mum", "l")
        keep.append(plot)
        all_spads = source.Integral(1, source.GetNbinsX() + 1)
        occupied_spads = source.Integral(2, source.GetNbinsX() + 1)
        for bin_index in range(1, plot.GetNbinsX() + 1):
            multiplicity = int(round(plot.GetBinCenter(bin_index)))
            if multiplicity < plot.GetNbinsX():
                all_fraction = source.GetBinContent(multiplicity + 1) / all_spads
            else:
                all_fraction = source.Integral(
                    multiplicity + 1,
                    source.GetNbinsX() + 1,
                ) / all_spads
            rows.append(
                {
                    "sample": sample_label,
                    "legacy_display_pitch_um": legacy_pitch,
                    "physical_pitch_um": physical_pitch,
                    "photons_per_spad": multiplicity,
                    "fraction_of_all_spads": all_fraction,
                    "fraction_of_occupied_spads": plot.GetBinContent(bin_index),
                    "multi_hit_fraction_of_occupied_spads": multi_hit,
                }
            )

    legend.Draw()
    keep.extend(draw_flash_label(ROOT, pad, right=sample_label, y=0.955))

    note = ROOT.TLatex()
    note.SetNDC(True)
    note.SetTextFont(42)
    note.SetTextSize(0.030)
    note.SetTextColor(ROOT.TColor.GetColor("#4B5563"))
    note.DrawLatex(0.15, 0.84, "Unity PDE; before recovery")
    note.DrawLatex(0.15, 0.79, "Occupied SPADs only; final bin is #geq 16 photons")
    keep.append(note)
    return keep, rows


def _draw_multi_hit_panel(ROOT, pad, files, palette):
    """Update the CPAD 2025 fill-fraction scan with corrected physical pitch."""

    pad.cd()
    pad.SetLogy(True)
    style_pad(pad, left=0.15, right=0.035, top=0.15, bottom=0.16)
    frame = ROOT.TH1D(
        "one_fiber_multi_hit_frame",
        ";Physical SPAD pitch [#mum];Occupied SPADs with #geq 2 photons [%]",
        100,
        0.0,
        36.0,
    )
    frame.SetDirectory(0)
    frame.SetMinimum(0.08)
    frame.SetMaximum(85.0)
    style_axes(frame)
    frame.Draw("AXIS")
    keep = [frame]

    band = ROOT.TBox(10.0, 0.08, 20.0, 85.0)
    band.SetFillColorAlpha(palette["green"], 0.10)
    band.SetLineColor(0)
    band.Draw("SAME")
    keep.append(band)

    for threshold in (1.0, 10.0):
        line = ROOT.TLine(0.0, threshold, 36.0, threshold)
        line.SetLineColor(palette["gray"])
        line.SetLineStyle(2)
        line.Draw()
        keep.append(line)

    rows: list[dict[str, float | str]] = []
    graphs = []
    for root_file, sample, color, marker in (
        (files[0], "100 GeV positron", palette["blue"], 20),
        (files[1], "1 TeV pion", palette["red"], 21),
    ):
        x_values = []
        y_values = []
        for legacy_pitch in LEGACY_LABELS_UM:
            source = root_file.Get(f"nPhotonsPerChannel_{legacy_pitch}x{legacy_pitch}")
            if not source:
                raise RuntimeError(f"Missing nPhotonsPerChannel_{legacy_pitch}x{legacy_pitch}")
            occupied = source.Integral(2, source.GetNbinsX() + 1)
            multiple = source.Integral(3, source.GetNbinsX() + 1)
            physical_pitch = legacy_pitch * DISPLAY_TO_PHYSICAL
            fraction = 100.0 * multiple / occupied if occupied else 0.0
            x_values.append(physical_pitch)
            y_values.append(fraction)
            rows.append(
                {
                    "sample": sample,
                    "legacy_display_pitch_um": legacy_pitch,
                    "physical_pitch_um": physical_pitch,
                    "photons_per_spad": ">=2",
                    "fraction_of_all_spads": 0.0,
                    "fraction_of_occupied_spads": fraction / 100.0,
                    "multi_hit_fraction_of_occupied_spads": fraction / 100.0,
                }
            )
        graph = ROOT.TGraph(
            len(x_values),
            array("d", x_values),
            array("d", y_values),
        )
        graph.SetLineColor(color)
        graph.SetMarkerColor(color)
        graph.SetLineWidth(3)
        graph.SetMarkerStyle(marker)
        graph.SetMarkerSize(1.05)
        graph.Draw("LP SAME")
        graphs.append(graph)
        keep.append(graph)

    legend = make_legend(ROOT, 0.19, 0.66, 0.52, 0.81, text_size=0.037)
    legend.AddEntry(graphs[0], "100 GeV e^{+}", "lp")
    legend.AddEntry(graphs[1], "1 TeV #pi^{+}", "lp")
    legend.Draw()
    keep.append(legend)
    keep.extend(draw_flash_label(ROOT, pad, right="one Cherenkov fiber", y=0.955))

    note = ROOT.TLatex()
    note.SetNDC(True)
    note.SetTextFont(42)
    note.SetTextColor(palette["gray"])
    note.SetTextSize(0.030)
    note.DrawLatex(0.54, 0.27, "Shaded band: 10-20 #mum")
    note.DrawLatex(0.54, 0.22, "Unity PDE; before recovery")
    keep.append(note)
    return keep, rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--electron", type=Path, default=DEFAULT_ELECTRON)
    parser.add_argument("--pion", type=Path, default=DEFAULT_PION)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    import ROOT  # type: ignore

    palette = set_flash_style(ROOT)
    args.output.mkdir(parents=True, exist_ok=True)
    files = [ROOT.TFile.Open(str(args.electron)), ROOT.TFile.Open(str(args.pion))]
    if any(not item or item.IsZombie() for item in files):
        raise RuntimeError("Could not open one or more validated comparison ROOT files")

    colors = [
        palette["red"],
        palette["blue"],
        palette["green"],
        palette["purple"],
        palette["orange"],
    ]
    canvas = ROOT.TCanvas("one_fiber_segmentation", "One-fiber segmentation", 1600, 800)
    canvas.Divide(2, 1, 0.01, 0.0)
    keep: list[object] = [canvas]
    rows: list[dict[str, float | str]] = []
    panel_keep, panel_rows = _draw_panel(
        ROOT,
        canvas.cd(1),
        files[0],
        "100 GeV positron",
        colors,
    )
    keep.extend(panel_keep)
    rows.extend(panel_rows)
    scan_keep, scan_rows = _draw_multi_hit_panel(ROOT, canvas.cd(2), files, palette)
    keep.extend(scan_keep)
    rows.extend(scan_rows)

    canvas.cd()
    save_canvas(canvas, args.output / "one_fiber_segmentation")

    with (args.output / "one_fiber_segmentation.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    for root_file in files:
        root_file.Close()
    _ = keep


if __name__ == "__main__":
    main()
