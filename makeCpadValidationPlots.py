#!/usr/bin/env python3
"""Build validation plots for the CPAD FLASH/HG-DREAM simulation story.

The suite checks four distinct links in the chain:

1. a deterministic event-level view of incident, accepted, and lost photons;
2. optical propagation time and production-depth-dependent arrival profiles;
3. event-to-event photon-count and occupancy-loss distributions;
4. numerical closure between the reproduced historical scan and the new
   pitch/recovery implementation.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from dsipm.digitizer import _apply_nonparalyzable_dead_time, _cell_indices
from dsipm.root_io import iter_root_events


REPO_DIR = Path(__file__).resolve().parent
DEFAULT_SAMPLES = {
    "100 GeV positron": REPO_DIR / "mc_dreamsim_e+_100_101.root",
    "1 TeV pion": REPO_DIR / "mc_dreamsim_pi+_1000_1001.root",
}
PITCH_UM = 25.0
COORDINATE_SCALE = 3.0
RECOVERY_TIMES_NS = (1.0, 4.3, 40.0)
GATE_NS = (0.0, 40.0)
C_CM_PER_NS = 29.9792458
FUSED_SILICA_INDEX = 1.49


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


def _cell_response(event, recovery_ns: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    event = event.sorted_by_time()
    ix, iy, cell_ids, _ = _cell_indices(
        event.x_mm * COORDINATE_SCALE,
        event.y_mm * COORDINATE_SCALE,
        PITCH_UM / 1000.0,
        -100.0 * COORDINATE_SCALE,
        100.0 * COORDINATE_SCALE,
    )
    accepted = _apply_nonparalyzable_dead_time(
        cell_ids,
        event.time_ns,
        recovery_ns,
    )
    return accepted, ix, iy


def _collect_event_statistics(ROOT, *, max_events: int | None):
    counts: dict[str, np.ndarray] = {}
    loss: dict[str, dict[float, np.ndarray]] = {}
    depth_time = ROOT.TH2D(
        "validation_pion_depth_time",
        ";Photon production depth [mm];Arrival time [ns];Photons",
        40,
        0.0,
        2000.0,
        240,
        8.0,
        20.0,
    )
    depth_time.SetDirectory(0)
    depth_slices = []
    for index in range(4):
        histogram = ROOT.TH1D(
            f"validation_depth_slice_{index}",
            ";Arrival time [ns];Normalized photons / 50 ps",
            240,
            8.0,
            20.0,
        )
        histogram.SetDirectory(0)
        depth_slices.append(histogram)

    for label, path in DEFAULT_SAMPLES.items():
        event_counts: list[int] = []
        event_loss = {recovery: [] for recovery in RECOVERY_TIMES_NS}
        events = iter_root_events(
            path,
            max_events=max_events,
            geometry_mode="old_2025",
            gate_ns=GATE_NS,
            detector_min_mm=-100.0,
            detector_max_mm=100.0,
            fiber_index=0,
        )
        for event in events:
            event = event.sorted_by_time()
            n_photons = event.selected_photon_count
            event_counts.append(n_photons)
            for recovery_ns in RECOVERY_TIMES_NS:
                accepted, _, _ = _cell_response(event, recovery_ns)
                retained = np.count_nonzero(accepted) / n_photons if n_photons else 1.0
                event_loss[recovery_ns].append(1.0 - retained)

            if label == "1 TeV pion":
                _fill_2d(depth_time, event.production_z_mm, event.time_ns)
                for index in range(4):
                    low = 500.0 * index
                    high = low + 500.0
                    selected = (
                        (event.production_z_mm >= low)
                        & (event.production_z_mm < high)
                    )
                    _fill_1d(depth_slices[index], event.time_ns[selected])

        counts[label] = np.asarray(event_counts, dtype=np.int64)
        loss[label] = {
            recovery: np.asarray(values, dtype=np.float64)
            for recovery, values in event_loss.items()
        }
    return counts, loss, depth_time, depth_slices


def _representative_pion_event(counts: np.ndarray, *, max_events: int | None):
    target = float(np.quantile(counts, 0.90))
    event_index = int(np.argmin(np.abs(counts.astype(float) - target)))
    events = iter_root_events(
        DEFAULT_SAMPLES["1 TeV pion"],
        max_events=max_events,
        geometry_mode="old_2025",
        gate_ns=GATE_NS,
        detector_min_mm=-100.0,
        detector_max_mm=100.0,
        fiber_index=0,
    )
    for event in events:
        if event.event_index == event_index:
            return event.sorted_by_time(), target
    raise RuntimeError(f"Could not recover representative event {event_index}")


def _style_axes(histogram, *, size: float = 0.045) -> None:
    for axis in (histogram.GetXaxis(), histogram.GetYaxis()):
        axis.SetTitleFont(42)
        axis.SetLabelFont(42)
        axis.SetTitleSize(size)
        axis.SetLabelSize(size * 0.82)
    histogram.GetXaxis().SetTitleOffset(1.08)
    histogram.GetYaxis().SetTitleOffset(1.22)


def _draw_event_occupancy(ROOT, output_dir: Path, event, percentile_target: float):
    accepted, ix, iy = _cell_response(event, 40.0)
    lost = ~accepted
    # Use a collision-safe key independent of the large detector indexing.
    packed = (ix.astype(np.int64) << 32) | iy.astype(np.int64)
    unique, multiplicity = np.unique(packed, return_counts=True)
    hottest = unique[int(np.argmax(multiplicity))]
    hot_index = int(np.flatnonzero(packed == hottest)[0])
    display_pitch_mm = (PITCH_UM / 1000.0) / COORDINATE_SCALE
    scaled_min = -100.0 * COORDINATE_SCALE
    hot_center_x = (
        scaled_min + (float(ix[hot_index]) + 0.5) * (PITCH_UM / 1000.0)
    ) / COORDINATE_SCALE
    hot_center_y = (
        scaled_min + (float(iy[hot_index]) + 0.5) * (PITCH_UM / 1000.0)
    ) / COORDINATE_SCALE

    ROOT.gStyle.SetPalette(ROOT.kViridis)
    ROOT.gStyle.SetPaintTextFormat("1.0f")
    canvas = ROOT.TCanvas("validation_event_occupancy", "Event occupancy", 1800, 680)
    canvas.Divide(3, 1, 0.008, 0.0)
    retained: list[object] = []

    canvas.cd(1)
    ROOT.gPad.SetLeftMargin(0.13)
    ROOT.gPad.SetRightMargin(0.04)
    ROOT.gPad.SetTopMargin(0.17)
    ROOT.gPad.SetBottomMargin(0.14)
    ROOT.gPad.SetTicks(1, 1)
    frame = ROOT.TH2D(
        "event_xy_frame",
        ";Sensor x [mm];Sensor y [mm]",
        10,
        -0.5,
        0.5,
        10,
        -0.5,
        0.5,
    )
    frame.SetDirectory(0)
    frame.SetStats(0)
    _style_axes(frame)
    frame.Draw("AXIS")
    accepted_graph = ROOT.TGraph(
        int(np.count_nonzero(accepted)),
        np.ascontiguousarray(event.x_mm[accepted], dtype=np.float64),
        np.ascontiguousarray(event.y_mm[accepted], dtype=np.float64),
    )
    accepted_graph.SetMarkerStyle(20)
    accepted_graph.SetMarkerSize(0.50)
    accepted_graph.SetMarkerColor(ROOT.TColor.GetColor("#2F6690"))
    accepted_graph.Draw("P SAME")
    lost_graph = ROOT.TGraph(
        int(np.count_nonzero(lost)),
        np.ascontiguousarray(event.x_mm[lost], dtype=np.float64),
        np.ascontiguousarray(event.y_mm[lost], dtype=np.float64),
    )
    lost_graph.SetMarkerStyle(34)
    lost_graph.SetMarkerSize(0.95)
    lost_graph.SetMarkerColor(ROOT.TColor.GetColor("#CC0000"))
    lost_graph.Draw("P SAME")
    legend = ROOT.TLegend(0.14, 0.73, 0.56, 0.84)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.034)
    legend.AddEntry(accepted_graph, "Accepted photon", "p")
    legend.AddEntry(lost_graph, "Repeated hit lost", "p")
    legend.Draw()

    canvas.cd(2)
    ROOT.gPad.SetLeftMargin(0.13)
    ROOT.gPad.SetRightMargin(0.18)
    ROOT.gPad.SetTopMargin(0.17)
    ROOT.gPad.SetBottomMargin(0.14)
    ROOT.gPad.SetTicks(1, 1)
    half_width = 10.0 * display_pitch_mm
    occupancy = ROOT.TH2D(
        "event_occupancy_zoom",
        ";Sensor x [mm];Sensor y [mm];Photons in one SPAD",
        20,
        hot_center_x - half_width,
        hot_center_x + half_width,
        20,
        hot_center_y - half_width,
        hot_center_y + half_width,
    )
    occupancy.SetDirectory(0)
    occupancy.SetStats(0)
    cell_center_x = (
        scaled_min + (ix.astype(np.float64) + 0.5) * (PITCH_UM / 1000.0)
    ) / COORDINATE_SCALE
    cell_center_y = (
        scaled_min + (iy.astype(np.float64) + 0.5) * (PITCH_UM / 1000.0)
    ) / COORDINATE_SCALE
    _fill_2d(occupancy, cell_center_x, cell_center_y)
    occupancy.SetMinimum(0.0)
    occupancy.SetMaximum(max(3.0, occupancy.GetMaximum()))
    _style_axes(occupancy, size=0.042)
    occupancy.GetZaxis().SetTitleSize(0.038)
    occupancy.GetZaxis().SetLabelSize(0.032)
    occupancy.Draw("COLZ TEXT")

    canvas.cd(3)
    ROOT.gPad.SetLeftMargin(0.13)
    ROOT.gPad.SetRightMargin(0.04)
    ROOT.gPad.SetTopMargin(0.17)
    ROOT.gPad.SetBottomMargin(0.14)
    ROOT.gPad.SetTicks(1, 1)
    time_all = ROOT.TH1D(
        "event_time_all",
        ";Photon arrival time [ns];Photons / 50 ps",
        700,
        5.0,
        40.0,
    )
    time_accepted = time_all.Clone("event_time_accepted")
    time_lost = time_all.Clone("event_time_lost")
    for histogram in (time_all, time_accepted, time_lost):
        histogram.SetDirectory(0)
        histogram.SetStats(0)
    _fill_1d(time_all, event.time_ns)
    _fill_1d(time_accepted, event.time_ns[accepted])
    _fill_1d(time_lost, event.time_ns[lost])
    time_all.SetLineColor(ROOT.TColor.GetColor("#18212B"))
    time_all.SetLineWidth(3)
    time_accepted.SetLineColor(ROOT.TColor.GetColor("#159D66"))
    time_accepted.SetLineWidth(2)
    time_lost.SetLineColor(ROOT.TColor.GetColor("#CC0000"))
    time_lost.SetLineWidth(2)
    _style_axes(time_all)
    time_all.SetMaximum(time_all.GetMaximum() * 1.22)
    time_all.Draw("HIST")
    time_accepted.Draw("HIST SAME")
    time_lost.Draw("HIST SAME")
    time_legend = ROOT.TLegend(0.56, 0.66, 0.94, 0.84)
    time_legend.SetBorderSize(0)
    time_legend.SetFillStyle(0)
    time_legend.SetTextSize(0.032)
    time_legend.AddEntry(time_all, "Incident", "l")
    time_legend.AddEntry(time_accepted, "Accepted", "l")
    time_legend.AddEntry(time_lost, "Lost repeat", "l")
    time_legend.Draw()

    labels = (
        "Full fiber face",
        "20#times20 SPAD zoom around densest cell",
        "Same event in time",
    )
    for pad_index, label in enumerate(labels, 1):
        canvas.cd(pad_index)
        title = ROOT.TLatex()
        title.SetNDC(True)
        title.SetTextFont(62)
        title.SetTextSize(0.041)
        title.DrawLatex(0.13, 0.93, label)
        retained.append(title)
    canvas.cd()
    summary = ROOT.TLatex()
    summary.SetNDC(True)
    summary.SetTextFont(42)
    summary.SetTextSize(0.025)
    summary.SetTextColor(ROOT.TColor.GetColor("#4B5563"))
    summary.DrawLatex(
        0.035,
        0.965,
        (
            f"1 TeV pion event {event.event_index}, selected objectively near the "
            f"90th-percentile photon count ({percentile_target:.0f}); "
            f"{event.selected_photon_count} incident, {int(np.count_nonzero(lost))} lost"
        ),
    )
    retained.extend(
        [
            frame,
            accepted_graph,
            lost_graph,
            legend,
            occupancy,
            time_all,
            time_accepted,
            time_lost,
            time_legend,
            summary,
        ]
    )
    outputs = [
        output_dir / "validation_event_occupancy.png",
        output_dir / "validation_event_occupancy.pdf",
    ]
    for path in outputs:
        canvas.SaveAs(str(path))
    metadata = {
        "event_index": event.event_index,
        "selection": "closest event to the 90th percentile of selected photons/event",
        "target_photons": percentile_target,
        "incident_photons": event.selected_photon_count,
        "accepted_photons": int(np.count_nonzero(accepted)),
        "lost_photons": int(np.count_nonzero(lost)),
        "active_spads": int(unique.size),
        "maximum_cell_multiplicity": int(multiplicity.max()),
    }
    return outputs, metadata


def _fill_transport_histogram(ROOT, path: Path, *, max_events: int | None):
    root_file = ROOT.TFile.Open(str(path), "READ")
    tree = root_file.Get("tree")
    histogram = ROOT.TH2D(
        "validation_transport_time",
        ";Axial distance to sensor [cm];Photon propagation time [ns];Photons",
        100,
        100.0,
        205.0,
        160,
        4.0,
        12.0,
    )
    histogram.SetDirectory(0)
    selected_photons = 0
    for event_index, event in enumerate(tree):
        if max_events is not None and event_index >= max_events:
            break
        final_time = np.asarray(event.OP_time_final, dtype=np.float64)
        produced_time = np.asarray(event.OP_time_produced, dtype=np.float64)
        produced_z = np.asarray(event.OP_pos_produced_z, dtype=np.float64)
        final_z = np.asarray(event.OP_pos_final_z, dtype=np.float64)
        production_fiber = np.asarray(event.OP_productionFiber, dtype=np.int64)
        core = np.asarray(event.OP_isCoreC, dtype=bool)
        mask = (
            (production_fiber == 0)
            & core
            & (final_z > 0.0)
            & (final_time > GATE_NS[0])
            & (final_time < GATE_NS[1])
        )
        distance_cm = final_z[mask] - produced_z[mask]
        propagation_ns = final_time[mask] - produced_time[mask]
        physical = (distance_cm > 0.0) & (propagation_ns > 0.0)
        _fill_2d(histogram, distance_cm[physical], propagation_ns[physical])
        selected_photons += int(np.count_nonzero(physical))
    root_file.Close()
    return histogram, selected_photons


def _draw_transport_depth(ROOT, output_dir: Path, transport, depth_slices):
    ROOT.gStyle.SetPalette(ROOT.kViridis)
    ROOT.gStyle.SetNumberContours(255)
    canvas = ROOT.TCanvas(
        "validation_transport_depth",
        "Optical transport and depth timing",
        1600,
        720,
    )
    canvas.Divide(2, 1, 0.008, 0.0)
    retained: list[object] = []

    canvas.cd(1)
    ROOT.gPad.SetLeftMargin(0.13)
    ROOT.gPad.SetRightMargin(0.17)
    ROOT.gPad.SetTopMargin(0.16)
    ROOT.gPad.SetBottomMargin(0.14)
    ROOT.gPad.SetLogz(True)
    ROOT.gPad.SetTicks(1, 1)
    transport.SetStats(0)
    transport.SetMinimum(1.0)
    _style_axes(transport)
    transport.GetZaxis().SetTitleSize(0.040)
    transport.GetZaxis().SetLabelSize(0.032)
    transport.Draw("COLZ")
    x = np.asarray([100.0, 205.0], dtype=np.float64)
    y = FUSED_SILICA_INDEX * x / C_CM_PER_NS
    direct = ROOT.TGraph(2, x, y)
    direct.SetLineColor(ROOT.kWhite)
    direct.SetLineWidth(4)
    direct.SetLineStyle(2)
    direct.Draw("L SAME")
    direct_label = ROOT.TLatex(108.0, FUSED_SILICA_INDEX * 108.0 / C_CM_PER_NS + 0.55, "nL/c, n=1.49")
    direct_label.SetTextColor(ROOT.kWhite)
    direct_label.SetTextFont(62)
    direct_label.SetTextSize(0.034)
    direct_label.Draw()

    canvas.cd(2)
    ROOT.gPad.SetLeftMargin(0.13)
    ROOT.gPad.SetRightMargin(0.04)
    ROOT.gPad.SetTopMargin(0.16)
    ROOT.gPad.SetBottomMargin(0.14)
    ROOT.gPad.SetTicks(1, 1)
    colors = ["#2F6690", "#159D66", "#E07A1F", "#CC0000"]
    legend = ROOT.TLegend(0.57, 0.55, 0.94, 0.83)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.032)
    for index, (histogram, color) in enumerate(zip(depth_slices, colors)):
        plot = histogram.Clone(f"display_depth_slice_{index}")
        plot.SetDirectory(0)
        integral = plot.Integral(0, plot.GetNbinsX() + 1)
        if integral:
            plot.Scale(1.0 / integral)
        plot.SetLineColor(ROOT.TColor.GetColor(color))
        plot.SetLineWidth(3)
        plot.SetStats(0)
        plot.SetTitle(";Photon arrival time [ns];Normalized photons / 50 ps")
        plot.GetXaxis().SetRangeUser(8.0, 20.0)
        plot.SetMaximum(0.09)
        _style_axes(plot)
        plot.Draw("HIST" if index == 0 else "HIST SAME")
        legend.AddEntry(plot, f"{500 * index}-{500 * (index + 1)} mm", "l")
        retained.append(plot)
    legend.Draw()

    titles = (
        "Transport time follows the optical path",
        "Production depth changes the arrival profile",
    )
    for pad_index, label in enumerate(titles, 1):
        canvas.cd(pad_index)
        title = ROOT.TLatex()
        title.SetNDC(True)
        title.SetTextFont(62)
        title.SetTextSize(0.040)
        title.DrawLatex(0.13, 0.93, label)
        retained.append(title)
    retained.extend([transport, direct, direct_label, legend])
    outputs = [
        output_dir / "validation_transport_depth.png",
        output_dir / "validation_transport_depth.pdf",
    ]
    for path in outputs:
        canvas.SaveAs(str(path))
    return outputs


def _draw_event_distributions(ROOT, output_dir: Path, counts, loss):
    canvas = ROOT.TCanvas(
        "validation_event_distributions",
        "Event distributions",
        1600,
        720,
    )
    canvas.Divide(2, 1, 0.008, 0.0)
    retained: list[object] = []
    colors = {
        "100 GeV positron": ROOT.TColor.GetColor("#2F6690"),
        "1 TeV pion": ROOT.TColor.GetColor("#CC0000"),
        1.0: ROOT.TColor.GetColor("#159D66"),
        4.3: ROOT.TColor.GetColor("#E07A1F"),
        40.0: ROOT.TColor.GetColor("#CC0000"),
    }
    upper = max(float(np.quantile(values, 0.995)) for values in counts.values())
    x_max = max(500.0, 50.0 * np.ceil(1.08 * upper / 50.0))

    canvas.cd(1)
    ROOT.gPad.SetLeftMargin(0.13)
    ROOT.gPad.SetRightMargin(0.04)
    ROOT.gPad.SetTopMargin(0.16)
    ROOT.gPad.SetBottomMargin(0.14)
    ROOT.gPad.SetLogy(True)
    ROOT.gPad.SetTicks(1, 1)
    legend = ROOT.TLegend(0.56, 0.68, 0.94, 0.83)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.033)
    for index, label in enumerate(DEFAULT_SAMPLES):
        histogram = ROOT.TH1D(
            f"event_count_{index}",
            ";Selected photons / event;Fraction of events",
            80,
            0.0,
            x_max,
        )
        histogram.SetDirectory(0)
        histogram.SetStats(0)
        _fill_1d(histogram, counts[label].astype(np.float64))
        integral = histogram.Integral(0, histogram.GetNbinsX() + 1)
        if integral:
            histogram.Scale(1.0 / integral)
        histogram.SetLineColor(colors[label])
        histogram.SetLineWidth(3)
        histogram.SetMinimum(2.0e-5)
        histogram.SetMaximum(0.35)
        _style_axes(histogram)
        histogram.Draw("HIST" if index == 0 else "HIST SAME")
        legend.AddEntry(histogram, label, "l")
        retained.append(histogram)
    legend.Draw()

    canvas.cd(2)
    ROOT.gPad.SetLeftMargin(0.13)
    ROOT.gPad.SetRightMargin(0.04)
    ROOT.gPad.SetTopMargin(0.16)
    ROOT.gPad.SetBottomMargin(0.14)
    ROOT.gPad.SetTicks(1, 1)
    frame = ROOT.TH1D(
        "loss_profile_frame",
        ";Selected photons / event;Mean fraction lost to repeated hits",
        80,
        0.0,
        x_max,
    )
    frame.SetDirectory(0)
    frame.SetStats(0)
    frame.SetMinimum(0.0)
    frame.SetMaximum(0.26)
    _style_axes(frame)
    frame.Draw("AXIS")
    loss_legend = ROOT.TLegend(0.58, 0.65, 0.94, 0.83)
    loss_legend.SetBorderSize(0)
    loss_legend.SetFillStyle(0)
    loss_legend.SetTextSize(0.033)
    pion_counts = counts["1 TeV pion"].astype(np.float64)
    for recovery_ns in RECOVERY_TIMES_NS:
        profile = ROOT.TProfile(
            f"pion_loss_{recovery_ns:g}",
            "",
            28,
            0.0,
            x_max,
        )
        profile.SetDirectory(0)
        profile.SetStats(0)
        for n_photons, lost_fraction in zip(
            pion_counts,
            loss["1 TeV pion"][recovery_ns],
            strict=True,
        ):
            profile.Fill(float(n_photons), float(lost_fraction))
        profile.SetLineColor(colors[recovery_ns])
        profile.SetMarkerColor(colors[recovery_ns])
        profile.SetLineWidth(3)
        profile.SetMarkerStyle(20 if recovery_ns != 40.0 else 24)
        profile.SetMarkerSize(0.8)
        profile.Draw("LP SAME")
        label = (
            "One hit / event"
            if recovery_ns == 40.0
            else f"{recovery_ns:g} ns recovery"
        )
        loss_legend.AddEntry(profile, label, "lp")
        retained.append(profile)
    loss_legend.Draw()

    titles = (
        "Photon counts fluctuate event to event",
        "Brighter events lose a larger fraction",
    )
    for pad_index, label in enumerate(titles, 1):
        canvas.cd(pad_index)
        title = ROOT.TLatex()
        title.SetNDC(True)
        title.SetTextFont(62)
        title.SetTextSize(0.041)
        title.DrawLatex(0.13, 0.93, label)
        retained.append(title)
    retained.extend([legend, frame, loss_legend])
    outputs = [
        output_dir / "validation_event_distributions.png",
        output_dir / "validation_event_distributions.pdf",
    ]
    for path in outputs:
        canvas.SaveAs(str(path))
    return outputs


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _draw_closure(ROOT, output_dir: Path):
    historical = _read_csv(
        REPO_DIR / "results" / "2025-both" / "electron_pion_efficiency_summary.csv"
    )
    current = _read_csv(
        REPO_DIR / "results" / "flash-requirements-2026" / "flash_pitch_recovery_scan.csv"
    )
    historical_lookup = {
        ("100 GeV positron", float(row["pitch_um"])): float(
            row["electron_one_hit_efficiency"]
        )
        for row in historical
    }
    historical_lookup.update(
        {
            ("1 TeV pion", float(row["pitch_um"])): float(
                row["pion_one_hit_efficiency"]
            )
            for row in historical
        }
    )
    current_lookup = {
        (row["sample"], float(row["pitch_um"])): float(row["retained_fraction"])
        for row in current
        if np.isclose(float(row["recovery_ns"]), 40.0)
    }
    pitches = sorted(
        {
            pitch
            for sample, pitch in historical_lookup
            if pitch >= 10.0 and (sample, pitch) in current_lookup
        }
    )
    canvas = ROOT.TCanvas("validation_closure", "Historical closure", 1600, 720)
    canvas.Divide(2, 1, 0.008, 0.0)
    retained: list[object] = []
    colors = {
        "100 GeV positron": ROOT.TColor.GetColor("#2F6690"),
        "1 TeV pion": ROOT.TColor.GetColor("#CC0000"),
    }

    canvas.cd(1)
    ROOT.gPad.SetLeftMargin(0.13)
    ROOT.gPad.SetRightMargin(0.04)
    ROOT.gPad.SetTopMargin(0.16)
    ROOT.gPad.SetBottomMargin(0.14)
    ROOT.gPad.SetTicks(1, 1)
    frame = ROOT.TH1D(
        "closure_frame",
        ";Physical SPAD pitch [#mum];Accepted / incident photons",
        100,
        5.0,
        105.0,
    )
    frame.SetDirectory(0)
    frame.SetStats(0)
    frame.SetMinimum(0.45)
    frame.SetMaximum(1.02)
    _style_axes(frame)
    frame.Draw("AXIS")
    legend = ROOT.TLegend(0.14, 0.19, 0.64, 0.42)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.030)
    max_difference = 0.0
    residual_graphs = []
    for sample in DEFAULT_SAMPLES:
        x = np.asarray(pitches, dtype=np.float64)
        old = np.asarray(
            [historical_lookup[(sample, pitch)] for pitch in pitches],
            dtype=np.float64,
        )
        new = np.asarray(
            [current_lookup[(sample, pitch)] for pitch in pitches],
            dtype=np.float64,
        )
        old_graph = ROOT.TGraph(len(x), x, old)
        old_graph.SetLineColor(colors[sample])
        old_graph.SetLineWidth(4)
        old_graph.Draw("L SAME")
        new_graph = ROOT.TGraph(len(x), x, new)
        new_graph.SetMarkerColor(colors[sample])
        new_graph.SetMarkerStyle(24)
        new_graph.SetMarkerSize(1.2)
        new_graph.SetLineColor(colors[sample])
        new_graph.Draw("P SAME")
        legend.AddEntry(old_graph, f"{sample}: historical output", "l")
        legend.AddEntry(new_graph, f"{sample}: new 40 ns scan", "p")
        residual = 100.0 * (new - old)
        max_difference = max(max_difference, float(np.max(np.abs(residual))))
        residual_graphs.append((sample, ROOT.TGraph(len(x), x, residual)))
        retained.extend([old_graph, new_graph])
    legend.Draw()

    canvas.cd(2)
    ROOT.gPad.SetLeftMargin(0.13)
    ROOT.gPad.SetRightMargin(0.04)
    ROOT.gPad.SetTopMargin(0.16)
    ROOT.gPad.SetBottomMargin(0.14)
    ROOT.gPad.SetTicks(1, 1)
    residual_frame = ROOT.TH1D(
        "closure_residual_frame",
        ";Physical SPAD pitch [#mum];New - historical [percentage points]",
        100,
        5.0,
        105.0,
    )
    residual_frame.SetDirectory(0)
    residual_frame.SetStats(0)
    residual_frame.SetMinimum(-0.002)
    residual_frame.SetMaximum(0.002)
    _style_axes(residual_frame)
    residual_frame.GetYaxis().SetTitleOffset(1.55)
    residual_frame.Draw("AXIS")
    zero = ROOT.TLine(5.0, 0.0, 105.0, 0.0)
    zero.SetLineColor(ROOT.TColor.GetColor("#7A838C"))
    zero.SetLineStyle(2)
    zero.Draw("SAME")
    for sample, graph in residual_graphs:
        graph.SetMarkerColor(colors[sample])
        graph.SetMarkerStyle(20 if sample.startswith("100") else 21)
        graph.SetMarkerSize(1.1)
        graph.Draw("P SAME")
        retained.append(graph)
    result = ROOT.TLatex()
    result.SetNDC(True)
    result.SetTextFont(62)
    result.SetTextSize(0.046)
    result.SetTextColor(ROOT.TColor.GetColor("#159D66"))
    result.DrawLatex(0.20, 0.75, f"Maximum difference = {max_difference:.3g} points")

    titles = (
        "Two independent paths give the same efficiency",
        "Numerical closure",
    )
    for pad_index, label in enumerate(titles, 1):
        canvas.cd(pad_index)
        title = ROOT.TLatex()
        title.SetNDC(True)
        title.SetTextFont(62)
        title.SetTextSize(0.039)
        title.DrawLatex(0.13, 0.93, label)
        retained.append(title)
    retained.extend([frame, legend, residual_frame, zero, result])
    outputs = [
        output_dir / "validation_historical_closure.png",
        output_dir / "validation_historical_closure.pdf",
    ]
    for path in outputs:
        canvas.SaveAs(str(path))
    return outputs, max_difference


def make_plots(output_dir: Path, *, max_events: int | None) -> list[Path]:
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    output_dir = output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    counts, loss, depth_time, depth_slices = _collect_event_statistics(
        ROOT,
        max_events=max_events,
    )
    event, percentile_target = _representative_pion_event(
        counts["1 TeV pion"],
        max_events=max_events,
    )
    transport, transport_photons = _fill_transport_histogram(
        ROOT,
        DEFAULT_SAMPLES["1 TeV pion"],
        max_events=max_events,
    )

    outputs: list[Path] = []
    event_outputs, event_metadata = _draw_event_occupancy(
        ROOT,
        output_dir,
        event,
        percentile_target,
    )
    outputs.extend(event_outputs)
    outputs.extend(_draw_transport_depth(ROOT, output_dir, transport, depth_slices))
    outputs.extend(_draw_event_distributions(ROOT, output_dir, counts, loss))
    closure_outputs, max_difference = _draw_closure(ROOT, output_dir)
    outputs.extend(closure_outputs)

    root_path = output_dir / "validation_histograms.root"
    root_file = ROOT.TFile(str(root_path), "RECREATE")
    transport.Write()
    depth_time.Write()
    for histogram in depth_slices:
        histogram.Write()
    root_file.Close()
    outputs.append(root_path)

    summary = {
        "purpose": "Validation suite for the CPAD 2026 FLASH/HG-DREAM simulation narrative.",
        "selection": {
            "geometry_mode": "old_2025",
            "fiber_index": 0,
            "gate_ns": list(GATE_NS),
            "pitch_um": PITCH_UM,
            "coordinate_scale": COORDINATE_SCALE,
            "pde": "unity",
            "recovery_ns": list(RECOVERY_TIMES_NS),
        },
        "event_example": event_metadata,
        "event_count_quantiles": {
            label: {
                str(quantile): float(np.quantile(values, quantile))
                for quantile in (0.1, 0.5, 0.9, 0.99)
            }
            for label, values in counts.items()
        },
        "transport_validation": {
            "photons": transport_photons,
            "direct_path_model": "propagation time = n * axial distance / c",
            "refractive_index": FUSED_SILICA_INDEX,
            "c_cm_per_ns": C_CM_PER_NS,
        },
        "historical_closure": {
            "maximum_absolute_difference_percentage_points": max_difference,
            "comparison": "reproduced historical one-hit/event output versus new 40 ns recovery scan",
        },
        "caution": (
            "The 100 GeV positron and 1 TeV pion samples have different energies; "
            "they define two occupancy stress cases and are not a controlled particle-species comparison."
        ),
    }
    summary_path = output_dir / "validation_summary.json"
    with summary_path.open("w") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")
    outputs.append(summary_path)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_DIR / "results" / "cpad2026-validation",
    )
    parser.add_argument(
        "--max-events",
        type=int,
        default=None,
        help="Optional event limit per sample; default uses all 5,000 events.",
    )
    args = parser.parse_args()
    outputs = make_plots(args.output, max_events=args.max_events)
    print("Wrote:")
    for path in outputs:
        print(f"  {path}")


if __name__ == "__main__":
    main()
