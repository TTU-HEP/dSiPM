#!/usr/bin/env python3
"""Scan FLASH physical pitch and SPAD recovery on HG-DREAM photon samples."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import time

import numpy as np

from dsipm.requirements import interpolated_pitch_limit, scan_pitch_and_recovery
from dsipm.root_io import inspect_root_file, iter_root_events


REPO_DIR = Path(__file__).resolve().parent
DEFAULT_SAMPLES = {
    "100 GeV positron": REPO_DIR / "mc_dreamsim_e+_100_101.root",
    "1 TeV pion": REPO_DIR / "mc_dreamsim_pi+_1000_1001.root",
}
DEFAULT_PITCHES_UM = (10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0)
DEFAULT_RECOVERY_NS = (1.0, 2.0, 4.3, 10.0, 20.0, 40.0)
TARGETS = (0.90, 0.95, 0.99)


def _write_scan_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _draw_scan(output_dir: Path, rows: list[dict[str, object]]) -> list[Path]:
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    colors = {
        1.0: ROOT.TColor.GetColor("#159D66"),
        2.0: ROOT.TColor.GetColor("#2F6690"),
        4.3: ROOT.TColor.GetColor("#7A5195"),
        10.0: ROOT.TColor.GetColor("#E07A1F"),
        20.0: ROOT.TColor.GetColor("#C83E4D"),
        40.0: ROOT.TColor.GetColor("#4B5563"),
    }
    samples = ("100 GeV positron", "1 TeV pion")
    retained: list[object] = []
    canvas = ROOT.TCanvas("flash_requirement_scan", "FLASH requirements", 1600, 760)
    canvas.Divide(2, 1, 0.008, 0.0)

    for pad_index, sample in enumerate(samples, 1):
        canvas.cd(pad_index)
        ROOT.gPad.SetLeftMargin(0.12 if pad_index == 1 else 0.105)
        ROOT.gPad.SetRightMargin(0.035)
        ROOT.gPad.SetTopMargin(0.12)
        ROOT.gPad.SetBottomMargin(0.15)
        ROOT.gPad.SetGridy()
        ROOT.gPad.SetTicks(1, 1)
        frame = ROOT.TH1D(
            f"frame_{pad_index}",
            ";Physical SPAD pitch [#mum];Accepted / incident photons",
            100,
            5.0,
            105.0,
        )
        frame.SetDirectory(0)
        frame.SetMinimum(0.45)
        frame.SetMaximum(1.015)
        frame.GetXaxis().SetTitleSize(0.047)
        frame.GetYaxis().SetTitleSize(0.047)
        frame.GetXaxis().SetLabelSize(0.038)
        frame.GetYaxis().SetLabelSize(0.038)
        frame.GetXaxis().SetTitleOffset(1.15)
        frame.GetYaxis().SetTitleOffset(1.20 if pad_index == 1 else 1.03)
        frame.Draw("AXIS")
        retained.append(frame)

        for target, style in ((0.90, 3), (0.95, 2)):
            line = ROOT.TLine(5.0, target, 105.0, target)
            line.SetLineColor(ROOT.TColor.GetColor("#B8C0C8"))
            line.SetLineStyle(style)
            line.Draw("SAME")
            retained.append(line)
            target_label = ROOT.TLatex(102.0, target + 0.006, f"{target:.0%}")
            target_label.SetTextAlign(31)
            target_label.SetTextFont(42)
            target_label.SetTextSize(0.028)
            target_label.SetTextColor(ROOT.TColor.GetColor("#7A838C"))
            target_label.Draw("SAME")
            retained.append(target_label)

        legend = ROOT.TLegend(0.14, 0.20, 0.63, 0.46)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.033)
        legend.SetNColumns(2)
        for recovery_ns in DEFAULT_RECOVERY_NS:
            selected = sorted(
                (
                    row
                    for row in rows
                    if row["sample"] == sample
                    and np.isclose(float(row["recovery_ns"]), recovery_ns)
                ),
                key=lambda row: float(row["pitch_um"]),
            )
            x = np.asarray([float(row["pitch_um"]) for row in selected])
            y = np.asarray([float(row["retained_fraction"]) for row in selected])
            graph = ROOT.TGraph(len(selected), x, y)
            graph.SetName(f"retention_{pad_index}_{recovery_ns:g}ns")
            graph.SetLineColor(colors[recovery_ns])
            graph.SetMarkerColor(colors[recovery_ns])
            graph.SetLineWidth(3 if recovery_ns != 40.0 else 4)
            graph.SetMarkerStyle(20 if recovery_ns != 40.0 else 24)
            graph.SetMarkerSize(0.95)
            graph.Draw("LP SAME")
            label = (
                "40 ns (one hit / event)"
                if recovery_ns == 40.0
                else f"{recovery_ns:g} ns"
            )
            legend.AddEntry(graph, label, "lp")
            retained.append(graph)
        legend.Draw()
        retained.append(legend)

        sample_title = ROOT.TLatex()
        sample_title.SetNDC(True)
        sample_title.SetTextFont(62)
        sample_title.SetTextSize(0.039)
        sample_title.DrawLatex(0.13, 0.945, sample)
        subtitle = ROOT.TLatex()
        subtitle.SetNDC(True)
        subtitle.SetTextFont(42)
        subtitle.SetTextSize(0.028)
        subtitle.SetTextColor(ROOT.TColor.GetColor("#4B5563"))
        subtitle.DrawLatex(0.13, 0.905, "Unity PDE, 0--40 ns gate, historical geometry")
        retained.extend([sample_title, subtitle])
    outputs = [
        output_dir / "flash_pitch_recovery_scan.png",
        output_dir / "flash_pitch_recovery_scan.pdf",
    ]
    for path in outputs:
        canvas.SaveAs(str(path))
    return outputs


def run_scan(
    output_dir: Path,
    *,
    max_events: int | None,
    pitches_um: tuple[float, ...] = DEFAULT_PITCHES_UM,
    recovery_times_ns: tuple[float, ...] = DEFAULT_RECOVERY_NS,
) -> list[Path]:
    output_dir = output_dir.expanduser().resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    started = time.time()
    csv_rows: list[dict[str, object]] = []
    all_points = {}
    inputs = {}
    for sample, path in DEFAULT_SAMPLES.items():
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
        points = scan_pitch_and_recovery(
            events,
            pitches_um=pitches_um,
            recovery_times_ns=recovery_times_ns,
            coordinate_scale=3.0,
        )
        all_points[sample] = points
        for point in points:
            csv_rows.append(
                {
                    "sample": sample,
                    "events": point.events,
                    "pitch_um": point.pitch_um,
                    "recovery_ns": point.recovery_ns,
                    "incident_photons": point.incident_photons,
                    "accepted_photons": point.accepted_photons,
                    "retained_fraction": point.retained_fraction,
                }
            )

    scan_csv = output_dir / "flash_pitch_recovery_scan.csv"
    _write_scan_csv(scan_csv, csv_rows)

    threshold_rows = []
    for sample, points in all_points.items():
        for recovery_ns in recovery_times_ns:
            for target in TARGETS:
                limit = interpolated_pitch_limit(
                    points,
                    recovery_ns=recovery_ns,
                    target_retention=target,
                )
                threshold_rows.append(
                    {
                        "sample": sample,
                        "recovery_ns": recovery_ns,
                        "target_retention": target,
                        "interpolated_max_pitch_um": "" if limit is None else limit,
                    }
                )
    threshold_csv = output_dir / "flash_pitch_requirement_thresholds.csv"
    _write_scan_csv(threshold_csv, threshold_rows)

    manifest = {
        "purpose": "Separate the historical physical-pitch and full-event dead-time assumptions.",
        "inputs": inputs,
        "selection": {
            "geometry_mode": "old_2025",
            "fiber_index": 0,
            "gate_ns": [0.0, 40.0],
            "coordinate_scale": 3.0,
            "pde": "unity",
        },
        "scan": {
            "pitches_um": list(pitches_um),
            "recovery_times_ns": list(recovery_times_ns),
            "max_events_per_sample": max_events,
        },
        "limitations": [
            "No PDE loss, crosstalk, afterpulsing, dark counts, or hot-cell masking.",
            "No readout grouping, buffering, arbitration, bandwidth, or power model.",
            "The 100 GeV positron and 1 TeV pion samples are stress tests, not a collider spectrum.",
            "Pitch limits are linearly interpolated between the simulated grid points.",
        ],
        "elapsed_seconds": time.time() - started,
    }
    manifest_path = output_dir / "run_manifest.json"
    with manifest_path.open("w") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")

    plot_paths = _draw_scan(output_dir, csv_rows)
    return [scan_csv, threshold_csv, manifest_path, *plot_paths]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_DIR / "results" / "flash-requirements-2026",
    )
    parser.add_argument(
        "--max-events",
        type=int,
        default=None,
        help="limit each sample for commissioning; default uses every event",
    )
    args = parser.parse_args()
    outputs = run_scan(args.output, max_events=args.max_events)
    for output in outputs:
        print(f"Wrote {output}")


if __name__ == "__main__":
    main()
