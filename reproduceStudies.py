#!/usr/bin/env python3
"""Convenient entry points for the local historical dSiPM samples."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

from dsipm.cpad_plots import draw_cpad_landscape
from dsipm.outputs import PLOT_STYLES
from runDeviceComparison import DEFAULT_CATALOG, run_comparison


REPO_DIR = Path(__file__).resolve().parent
DEFAULT_CPAD_BENCHMARKS = REPO_DIR / "cpad_benchmarks.json"

STUDIES = {
    "quick": {
        "samples": [
            (
                "10GeV_eplus_2evt",
                "mc_testjob_run001_003_Test_2evt_e+_10.0_10.0.root",
            )
        ],
        "geometry": "old_2025",
        "fiber_index": 0,
        "profile": "past_study",
        "plot_style": "slides",
        "default_max_events": 2,
    },
    "2025-electron": {
        "samples": [("100GeV_eplus", "mc_dreamsim_e+_100_101.root")],
        "geometry": "old_2025",
        "fiber_index": 0,
        "profile": "past_study",
        "plot_style": "slides",
        "default_max_events": None,
    },
    "2025-pion": {
        "samples": [("1TeV_piplus", "mc_dreamsim_pi+_1000_1001.root")],
        "geometry": "old_2025",
        "fiber_index": 0,
        "profile": "past_study",
        "plot_style": "slides",
        "default_max_events": None,
    },
    "2025-both": {
        "samples": [
            ("100GeV_eplus", "mc_dreamsim_e+_100_101.root"),
            ("1TeV_piplus", "mc_dreamsim_pi+_1000_1001.root"),
        ],
        "geometry": "old_2025",
        "fiber_index": 0,
        "profile": "past_study",
        "plot_style": "slides",
        "default_max_events": None,
    },
    "rod90": {
        "samples": [("90GeV_piplus_rod", "sim_200evt_pi+_90+rod.root")],
        "geometry": "student_2026",
        "fiber_index": None,
        "profile": "default",
        "plot_style": "comparison",
        "default_max_events": None,
    },
}


def validate_saved_electron_reference(output_dir: Path) -> bool:
    """Compare the regenerated slide projections with the saved 2025 ROOT file."""

    import ROOT

    reference_path = REPO_DIR / "outfile_100_e+.root"
    generated_path = output_dir / "comparison.root"
    if not reference_path.exists() or not generated_path.exists():
        return False

    reference = ROOT.TFile.Open(str(reference_path), "READ")
    generated = ROOT.TFile.Open(str(generated_path), "READ")
    comparisons = []
    try:
        for pitch in (10, 25, 50, 75, 100):
            for stem in (
                "nPhotons_xyt_all",
                "nPhotons_xyt_oneHit",
                "nPhotonsPerChannel",
            ):
                name = f"{stem}_{pitch}x{pitch}"
                old_histogram = reference.Get(name)
                new_histogram = generated.Get(name)
                if not old_histogram or not new_histogram:
                    comparisons.append({"histogram": name, "present": False})
                    continue
                old_integral = float(old_histogram.Integral())
                new_integral = float(new_histogram.Integral())
                comparisons.append(
                    {
                        "histogram": name,
                        "present": True,
                        "xy_binning_match": (
                            old_histogram.GetNbinsX() == new_histogram.GetNbinsX()
                            and old_histogram.GetNbinsY()
                            == new_histogram.GetNbinsY()
                        ),
                        "reference_integral": old_integral,
                        "generated_integral": new_integral,
                        "integral_delta": new_integral - old_integral,
                    }
                )
    finally:
        reference.Close()
        generated.Close()

    exact = all(
        item.get("present")
        and item.get("xy_binning_match")
        and item.get("integral_delta") == 0.0
        for item in comparisons
    )
    report = {
        "reference": str(reference_path),
        "generated": str(generated_path),
        "exact_match": exact,
        "scope": (
            "XY bin counts and integrals for incident, one-hit, and "
            "photons-per-channel histograms at 10, 25, 50, 75, and 100 um"
        ),
        "comparisons": comparisons,
    }
    with (output_dir / "historical_reference_validation.json").open("w") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")
    print(f"Historical electron reference match: {'PASS' if exact else 'FAIL'}")
    return exact


def draw_electron_pion_multiplicity(output_root: Path) -> None:
    """Recreate the two-panel photons-per-channel comparison from slide 10."""

    import ROOT

    ROOT.gROOT.SetBatch(True)
    sources = (
        ("100GeV_eplus", "100 GeV positron"),
        ("1TeV_piplus", "1 TeV pion"),
    )
    canvas = ROOT.TCanvas(
        "canvas_electron_pion_multiplicity",
        "Number of photons per channel",
        1600,
        700,
    )
    canvas.Divide(2, 1)
    colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2, ROOT.kMagenta, ROOT.kOrange + 7]
    retained = []
    files = []
    for pad_index, (sample_name, title) in enumerate(sources, 1):
        source = ROOT.TFile.Open(
            str(output_root / sample_name / "comparison.root"), "READ"
        )
        files.append(source)
        canvas.cd(pad_index)
        ROOT.gPad.SetLogy()
        ROOT.gPad.SetTicks(1, 1)
        ROOT.gPad.SetLeftMargin(0.12)
        ROOT.gPad.SetRightMargin(0.25)
        legend = ROOT.TLegend(0.76, 0.50, 0.99, 0.88)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        for index, pitch in enumerate((10, 25, 50, 75, 100)):
            original = source.Get(f"nPhotonsPerChannel_{pitch}x{pitch}")
            plot = original.Clone(f"combined_{sample_name}_{pitch}um")
            plot.SetDirectory(0)
            integral = plot.Integral(0, plot.GetNbinsX() + 1)
            if integral:
                plot.Scale(1.0 / integral)
            plot.SetLineColor(colors[index])
            plot.SetLineWidth(2)
            plot.SetStats(0)
            plot.SetTitle(f"{title};NPhotons per channel;Normalized")
            plot.GetXaxis().SetRangeUser(0.0, 15.0)
            plot.SetMinimum(1.0e-8)
            plot.SetMaximum(1.5)
            plot.Draw("HIST" if index == 0 else "HIST SAME")
            legend.AddEntry(plot, f"{pitch}#times{pitch} #mum^{{2}}", "l")
            retained.append(plot)
        legend.Draw()
        retained.append(legend)
    destination = output_root / "photons_per_channel_electron_vs_pion.png"
    canvas.SaveAs(str(destination))
    for source in files:
        source.Close()
    print(f"Wrote combined slide comparison to {destination}")


def draw_electron_pion_efficiency(output_root: Path) -> None:
    """Summarize the historical dead-time loss for both slide samples."""

    import ROOT
    import numpy as np

    sample_files = {
        "electron": output_root / "100GeV_eplus" / "device_summary.csv",
        "pion": output_root / "1TeV_piplus" / "device_summary.csv",
    }
    by_sample: dict[str, dict[float, dict[str, str]]] = {}
    for sample, path in sample_files.items():
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        by_sample[sample] = {
            float(row["spad_pitch_um"]): row
            for row in rows
            if row["category"] == "historical idealization"
        }

    common_pitches = sorted(
        set(by_sample["electron"]) & set(by_sample["pion"])
    )
    summary_path = output_root / "electron_pion_efficiency_summary.csv"
    with summary_path.open("w", newline="") as handle:
        fieldnames = [
            "pitch_um",
            "electron_mean_incident",
            "electron_one_hit_efficiency",
            "pion_mean_incident",
            "pion_one_hit_efficiency",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for pitch in common_pitches:
            electron = by_sample["electron"][pitch]
            pion = by_sample["pion"][pitch]
            writer.writerow(
                {
                    "pitch_um": pitch,
                    "electron_mean_incident": electron["mean_selected"],
                    "electron_one_hit_efficiency": electron[
                        "aggregate_retained_fraction"
                    ],
                    "pion_mean_incident": pion["mean_selected"],
                    "pion_one_hit_efficiency": pion[
                        "aggregate_retained_fraction"
                    ],
                }
            )

    ROOT.gROOT.SetBatch(True)
    x = np.asarray(common_pitches, dtype=np.float64)
    electron_y = np.asarray(
        [
            float(by_sample["electron"][pitch]["aggregate_retained_fraction"])
            for pitch in common_pitches
        ],
        dtype=np.float64,
    )
    pion_y = np.asarray(
        [
            float(by_sample["pion"][pitch]["aggregate_retained_fraction"])
            for pitch in common_pitches
        ],
        dtype=np.float64,
    )
    electron_graph = ROOT.TGraph(int(x.size), x, electron_y)
    pion_graph = ROOT.TGraph(int(x.size), x, pion_y)
    electron_graph.SetLineColor(ROOT.kBlue + 1)
    electron_graph.SetMarkerColor(ROOT.kBlue + 1)
    pion_graph.SetLineColor(ROOT.kRed + 1)
    pion_graph.SetMarkerColor(ROOT.kRed + 1)
    for graph in (electron_graph, pion_graph):
        graph.SetLineWidth(3)
        graph.SetMarkerStyle(20)
        graph.SetMarkerSize(1.3)

    multigraph = ROOT.TMultiGraph()
    multigraph.SetTitle(
        "Full-event one-hit efficiency;Channel pitch [#mum];Accepted / incident"
    )
    multigraph.Add(electron_graph, "LP")
    multigraph.Add(pion_graph, "LP")
    multigraph.SetMinimum(0.0)
    multigraph.SetMaximum(1.05)
    canvas = ROOT.TCanvas(
        "canvas_electron_pion_efficiency",
        "Historical efficiency comparison",
        950,
        700,
    )
    canvas.SetGridy()
    canvas.SetTicks(1, 1)
    canvas.SetLeftMargin(0.12)
    canvas.SetBottomMargin(0.12)
    multigraph.Draw("A")
    legend = ROOT.TLegend(0.16, 0.18, 0.43, 0.34)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(electron_graph, "100 GeV positron", "lp")
    legend.AddEntry(pion_graph, "1 TeV pion", "lp")
    legend.Draw()
    destination = output_root / "deadtime_efficiency_electron_vs_pion.png"
    canvas.SaveAs(str(destination))
    print(f"Wrote combined efficiency summary to {destination}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Reproduce a named local dSiPM study."
    )
    parser.add_argument(
        "study",
        choices=sorted(STUDIES),
        nargs="?",
        default="quick",
        help="Study to run (default: %(default)s)",
    )
    parser.add_argument(
        "--output-root",
        default=str(REPO_DIR / "results"),
        help="Parent directory for study outputs (default: %(default)s)",
    )
    parser.add_argument(
        "--max-events",
        type=int,
        help="Override the study event count; useful for commissioning",
    )
    parser.add_argument("--seed", type=int, default=4307)
    parser.add_argument("--catalog", default=str(DEFAULT_CATALOG))
    parser.add_argument(
        "--plot-style",
        choices=sorted(PLOT_STYLES),
        help="Override the study's default plotting style",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    study = STUDIES[args.study]
    max_events = (
        args.max_events
        if args.max_events is not None
        else study["default_max_events"]
    )
    output_root = Path(args.output_root).expanduser().resolve() / args.study
    plot_style = args.plot_style or str(study["plot_style"])
    completed_outputs: dict[str, Path] = {}

    try:
        for sample_name, filename in study["samples"]:
            sample_output = output_root / sample_name
            run_comparison(
                input_path=REPO_DIR / filename,
                output_dir=sample_output,
                profile=str(study["profile"]),
                device_names=None,
                catalog_path=args.catalog,
                max_events=max_events,
                seed=args.seed,
                geometry_mode=str(study["geometry"]),
                detector_min_mm=-100.0,
                detector_max_mm=100.0,
                gate_ns=(0.0, 40.0),
                fiber_index=study["fiber_index"],
                plot_style=plot_style,
            )
            completed_outputs[sample_name] = sample_output
            if (
                sample_name == "100GeV_eplus"
                and max_events is None
                and plot_style in {"slides", "both"}
            ):
                validate_saved_electron_reference(sample_output)
        if (
            args.study == "2025-both"
            and plot_style in {"slides", "both"}
            and {"100GeV_eplus", "1TeV_piplus"} <= completed_outputs.keys()
        ):
            draw_electron_pion_multiplicity(output_root)
            draw_electron_pion_efficiency(output_root)
            for path in draw_cpad_landscape(
                output_root,
                DEFAULT_CPAD_BENCHMARKS,
            ):
                print(f"Wrote CPAD comparison product to {path}")
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
