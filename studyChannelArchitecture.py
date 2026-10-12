#!/usr/bin/env python3
"""Scan SPAD pitch and shared readout-channel size on one-fiber samples.

This study deliberately separates three response stages:

1. constant-PDE photon detection;
2. recovery of each physical SPAD;
3. finite resources shared by a square readout/timing channel.

The first scan fixes PDE and recovery so the new channel-size dependence is
visible.  A channel exports one finite count word and up to a fixed number of
timestamps after its own resolving time.  The model is an architecture stress
test, not a finished FLASH electronics specification.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass, field
import hashlib
import json
from pathlib import Path

import numpy as np

from dsipm.architecture import ChannelReadoutConfig, channel_readout
from dsipm.cms_style import (
    draw_panel_title,
    make_legend,
    save_canvas,
    set_flash_style,
    style_axes,
    style_pad,
)
from dsipm.digitizer import _apply_nonparalyzable_dead_time, _cell_indices
from dsipm.root_io import iter_root_events
from dsipm.workflow import load_sensor_workflow


REPO_DIR = Path(__file__).resolve().parent
WORKFLOW = load_sensor_workflow()
SAMPLES = {
    "100 GeV positron": REPO_DIR / "mc_dreamsim_e+_100_101.root",
    "1 TeV pion": REPO_DIR / "mc_dreamsim_pi+_1000_1001.root",
}
PITCHES_UM = WORKFLOW.scans.one_fiber_spad_pitch_um
# Include powers of two plus device-relevant intermediate sharing factors.
# In particular, 6 x 6 represents the published DIGILOG micro-SiPM topology.
SPADS_PER_CHANNEL_SIDE = WORKFLOW.scans.spads_per_readout_channel_side
SENSOR_MIN_MM = -1.5
SENSOR_MAX_MM = 1.5
PDE = WORKFLOW.reference.effective_pde
SPAD_RECOVERY_NS = WORKFLOW.reference.spad_recovery_ns
SEED = WORKFLOW.reference.seed
REFERENCE_GROUP_SIDE = max(
    1,
    int(
        round(
            WORKFLOW.reference.target_readout_channel_side_um
            / WORKFLOW.reference.spad_pitch_um
        )
    ),
)
MIN_TIMESTAMP_SURVIVAL = 0.95
MAX_OVERFLOW_EVENT_FRACTION = 0.05
MIN_COUNT_SURVIVAL = 0.99
TOPOLOGY_PROJECTIONS = (
    {
        "label": "FLASH example (not frozen)",
        "scan_pitch_um": WORKFLOW.reference.spad_pitch_um,
        "published_pitch_um": "",
        "spads_per_channel_side": REFERENCE_GROUP_SIDE,
        "note": (
            f"Illustrative {WORKFLOW.reference.spad_pitch_um:g} um pitch with "
            f"{REFERENCE_GROUP_SIDE} x {REFERENCE_GROUP_SIDE} shared readout; "
            "not an LDRD specification."
        ),
    },
    {
        "label": "Litavis topology",
        "scan_pitch_um": 10.0,
        "published_pitch_um": 10.17,
        "spads_per_channel_side": 4,
        "note": "Nearest scanned pitch; the commercial 16-bin macropixel histogram is not modeled here.",
    },
    {
        "label": "DIGILOG topology",
        "scan_pitch_um": 50.0,
        "published_pitch_um": 50.0,
        "spads_per_channel_side": 6,
        "note": "Physical grouping match; its first-time hybrid readout is not modeled here.",
    },
    {
        "label": "ASPIDES topology",
        "scan_pitch_um": 30.0,
        "published_pitch_um": 30.0,
        "spads_per_channel_side": 32,
        "note": "Physical grouping match; its count plus first/last-time policy is not modeled here.",
    },
)


@dataclass
class _Accumulator:
    events: int = 0
    incident_photons: int = 0
    detected_photons: int = 0
    accepted_photons: int = 0
    counted_hits: int = 0
    timestamp_candidates: int = 0
    timestamp_records: int = 0
    resolving_lost_hits: int = 0
    buffer_lost_records: int = 0
    overflow_events: int = 0
    count_saturation_events: int = 0
    active_channels: list[int] = field(default_factory=list)
    max_channel_hits: list[int] = field(default_factory=list)
    payload_bits: list[int] = field(default_factory=list)


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Cannot write empty table: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _safe_ratio(numerator: int, denominator: int, *, empty: float = 1.0) -> float:
    return numerator / denominator if denominator else empty


def _scan_sample(
    sample_index: int,
    sample: str,
    path: Path,
    *,
    max_events: int | None,
    readout_config: ChannelReadoutConfig,
) -> list[dict[str, object]]:
    accumulators = {
        (pitch, group_side): _Accumulator()
        for pitch in PITCHES_UM
        for group_side in SPADS_PER_CHANNEL_SIDE
    }
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
        random = np.random.default_rng(
            np.random.SeedSequence([SEED, sample_index, event.event_index, 6421])
        )
        detected_mask = random.random(event.selected_photon_count) < PDE
        detected_x = event.x_mm[detected_mask]
        detected_y = event.y_mm[detected_mask]
        detected_time = event.time_ns[detected_mask]

        for pitch_um in PITCHES_UM:
            ix, iy, cell_ids, n_spads_axis = _cell_indices(
                detected_x,
                detected_y,
                pitch_um / 1000.0,
                SENSOR_MIN_MM,
                SENSOR_MAX_MM,
            )
            accepted_mask = _apply_nonparalyzable_dead_time(
                cell_ids,
                detected_time,
                SPAD_RECOVERY_NS,
            )
            accepted_ix = ix[accepted_mask]
            accepted_iy = iy[accepted_mask]
            accepted_time = detected_time[accepted_mask]

            for group_side in SPADS_PER_CHANNEL_SIDE:
                key = (pitch_um, group_side)
                accumulator = accumulators[key]
                config = ChannelReadoutConfig(
                    **{
                        **asdict(readout_config),
                        "spads_per_channel_side": group_side,
                    }
                )
                result = channel_readout(
                    accepted_ix,
                    accepted_iy,
                    accepted_time,
                    n_spads_axis=n_spads_axis,
                    config=config,
                )
                accumulator.events += 1
                accumulator.incident_photons += event.selected_photon_count
                accumulator.detected_photons += int(np.count_nonzero(detected_mask))
                accumulator.accepted_photons += result.accepted_hits
                accumulator.counted_hits += result.counted_hits
                accumulator.timestamp_candidates += result.timestamp_candidates
                accumulator.timestamp_records += result.timestamp_records
                accumulator.resolving_lost_hits += result.resolving_lost_hits
                accumulator.buffer_lost_records += result.buffer_lost_records
                accumulator.overflow_events += int(result.any_buffer_overflow)
                accumulator.count_saturation_events += int(
                    result.count_saturated_channels > 0
                )
                accumulator.active_channels.append(result.active_channels)
                accumulator.max_channel_hits.append(result.max_hits_in_channel)
                accumulator.payload_bits.append(result.payload_bits)

    rows: list[dict[str, object]] = []
    sensor_width_um = 1000.0 * (SENSOR_MAX_MM - SENSOR_MIN_MM)
    for pitch_um in PITCHES_UM:
        n_spads_axis = int(np.ceil(sensor_width_um / pitch_um))
        for group_side in SPADS_PER_CHANNEL_SIDE:
            accumulator = accumulators[(pitch_um, group_side)]
            active = np.asarray(accumulator.active_channels, dtype=np.float64)
            maximum = np.asarray(accumulator.max_channel_hits, dtype=np.float64)
            payload = np.asarray(accumulator.payload_bits, dtype=np.float64)
            channels_axis = int(np.ceil(n_spads_axis / group_side))
            rows.append(
                {
                    "sample": sample,
                    "events": accumulator.events,
                    "pde": PDE,
                    "spad_pitch_um": pitch_um,
                    "spad_recovery_ns": SPAD_RECOVERY_NS,
                    "spads_per_channel_side": group_side,
                    "nominal_channel_side_um": pitch_um * group_side,
                    "readout_channels_per_sensor": channels_axis**2,
                    "timestamp_buffer_depth": readout_config.timestamp_buffer_depth,
                    "channel_resolving_time_ns": readout_config.channel_resolving_time_ns,
                    "count_bits": readout_config.count_bits,
                    "incident_photons": accumulator.incident_photons,
                    "detected_photons": accumulator.detected_photons,
                    "accepted_spad_hits": accumulator.accepted_photons,
                    "recovery_survival": _safe_ratio(
                        accumulator.accepted_photons,
                        accumulator.detected_photons,
                    ),
                    "counted_hits": accumulator.counted_hits,
                    "count_survival": _safe_ratio(
                        accumulator.counted_hits,
                        accumulator.accepted_photons,
                    ),
                    "timestamp_candidates": accumulator.timestamp_candidates,
                    "timestamp_records": accumulator.timestamp_records,
                    "timestamp_survival": _safe_ratio(
                        accumulator.timestamp_records,
                        accumulator.accepted_photons,
                    ),
                    "resolving_survival": _safe_ratio(
                        accumulator.timestamp_candidates,
                        accumulator.accepted_photons,
                    ),
                    "buffer_survival_after_resolving": _safe_ratio(
                        accumulator.timestamp_records,
                        accumulator.timestamp_candidates,
                    ),
                    "events_with_buffer_overflow_fraction": _safe_ratio(
                        accumulator.overflow_events,
                        accumulator.events,
                        empty=0.0,
                    ),
                    "events_with_count_saturation_fraction": _safe_ratio(
                        accumulator.count_saturation_events,
                        accumulator.events,
                        empty=0.0,
                    ),
                    "mean_active_channels": float(np.mean(active)) if active.size else 0.0,
                    "median_max_hits_per_channel": float(np.median(maximum)) if maximum.size else 0.0,
                    "p90_max_hits_per_channel": float(np.quantile(maximum, 0.90)) if maximum.size else 0.0,
                    "mean_payload_kbit": float(np.mean(payload) / 1000.0) if payload.size else 0.0,
                    "p90_payload_kbit": float(np.quantile(payload, 0.90) / 1000.0) if payload.size else 0.0,
                }
            )
    return rows


def _lookup(rows: list[dict[str, object]], sample: str) -> dict[tuple[float, int], dict[str, object]]:
    return {
        (float(row["spad_pitch_um"]), int(row["spads_per_channel_side"])): row
        for row in rows
        if row["sample"] == sample
    }


def _select_working_points(
    rows: list[dict[str, object]],
    sample: str,
) -> list[dict[str, object]]:
    """Choose the lowest-payload architecture satisfying explicit requirements.

    The criteria are deliberately simple and visible.  They are not claimed as
    final FLASH requirements; they make the consequence of the assumed
    eight-word timing buffer testable and give the next electronics iteration a
    concrete target.
    """

    selected_rows: list[dict[str, object]] = []
    for pitch_um in PITCHES_UM:
        candidates = [
            row
            for row in rows
            if row["sample"] == sample
            and float(row["spad_pitch_um"]) == pitch_um
            and float(row["timestamp_survival"]) >= MIN_TIMESTAMP_SURVIVAL
            and float(row["events_with_buffer_overflow_fraction"])
            <= MAX_OVERFLOW_EVENT_FRACTION
            and float(row["count_survival"]) >= MIN_COUNT_SURVIVAL
        ]
        if not candidates:
            selected_rows.append(
                {
                    "sample": sample,
                    "spad_pitch_um": pitch_um,
                    "status": "no passing point",
                    "spads_per_channel_side": "",
                    "nominal_channel_side_um": "",
                    "readout_channels_per_sensor": "",
                    "timestamp_survival": "",
                    "events_with_buffer_overflow_fraction": "",
                    "count_survival": "",
                    "mean_payload_kbit": "",
                    "requirement_timestamp_survival_min": MIN_TIMESTAMP_SURVIVAL,
                    "requirement_overflow_event_fraction_max": MAX_OVERFLOW_EVENT_FRACTION,
                    "requirement_count_survival_min": MIN_COUNT_SURVIVAL,
                }
            )
            continue
        best = min(
            candidates,
            key=lambda row: (
                float(row["mean_payload_kbit"]),
                -float(row["timestamp_survival"]),
            ),
        )
        selected_rows.append(
            {
                "sample": sample,
                "spad_pitch_um": pitch_um,
                "status": "passes",
                "spads_per_channel_side": best["spads_per_channel_side"],
                "nominal_channel_side_um": best["nominal_channel_side_um"],
                "readout_channels_per_sensor": best["readout_channels_per_sensor"],
                "timestamp_survival": best["timestamp_survival"],
                "events_with_buffer_overflow_fraction": best[
                    "events_with_buffer_overflow_fraction"
                ],
                "count_survival": best["count_survival"],
                "mean_payload_kbit": best["mean_payload_kbit"],
                "requirement_timestamp_survival_min": MIN_TIMESTAMP_SURVIVAL,
                "requirement_overflow_event_fraction_max": MAX_OVERFLOW_EVENT_FRACTION,
                "requirement_count_survival_min": MIN_COUNT_SURVIVAL,
            }
        )
    return selected_rows


def _project_device_topologies(
    rows: list[dict[str, object]],
    sample: str,
) -> list[dict[str, object]]:
    """Report real-device groupings under one common electronics assumption."""

    lookup = _lookup(rows, sample)
    projections: list[dict[str, object]] = []
    for topology in TOPOLOGY_PROJECTIONS:
        pitch = float(topology["scan_pitch_um"])
        group_side = int(topology["spads_per_channel_side"])
        row = lookup[(pitch, group_side)]
        projections.append(
            {
                "sample": sample,
                "label": topology["label"],
                "scan_pitch_um": pitch,
                "published_pitch_um": topology["published_pitch_um"],
                "spads_per_channel_side": group_side,
                "nominal_channel_side_um": row["nominal_channel_side_um"],
                "timestamp_survival": row["timestamp_survival"],
                "events_with_buffer_overflow_fraction": row[
                    "events_with_buffer_overflow_fraction"
                ],
                "count_survival": row["count_survival"],
                "mean_payload_kbit": row["mean_payload_kbit"],
                "interpretation": "common 8-word-buffer projection; not a complete device prediction",
                "note": topology["note"],
            }
        )
    return projections


def _category_heatmap(ROOT, name: str, title: str):
    histogram = ROOT.TH2D(
        name,
        title,
        len(PITCHES_UM),
        0.5,
        len(PITCHES_UM) + 0.5,
        len(SPADS_PER_CHANNEL_SIDE),
        0.5,
        len(SPADS_PER_CHANNEL_SIDE) + 0.5,
    )
    histogram.SetDirectory(0)
    for index, pitch in enumerate(PITCHES_UM, 1):
        histogram.GetXaxis().SetBinLabel(index, f"{pitch:g}")
    for index, group_side in enumerate(SPADS_PER_CHANNEL_SIDE, 1):
        histogram.GetYaxis().SetBinLabel(index, f"{group_side}#times{group_side}")
    return histogram


def _draw_heatmaps(
    output_dir: Path,
    rows: list[dict[str, object]],
    sample: str,
) -> list[Path]:
    import ROOT

    palette = set_flash_style(ROOT)
    ROOT.gStyle.SetPalette(ROOT.kViridis)
    ROOT.gStyle.SetPaintTextFormat(".1f")
    lookup = _lookup(rows, sample)
    metrics = (
        ("timestamp_survival", "Timestamp retention [%]", 100.0, 0.0, 100.0),
        (
            "events_with_buffer_overflow_fraction",
            "Events with any timestamp overflow [%]",
            100.0,
            0.0,
            100.0,
        ),
        ("mean_payload_kbit", "Mean payload [kbit / event]", 1.0, 0.0, None),
        ("count_survival", "8-bit count retention [%]", 100.0, 90.0, 100.0),
    )
    canvas = ROOT.TCanvas("channel_architecture_matrix", "channel architecture matrix", 1800, 1120)
    canvas.Divide(2, 2, 0.004, 0.004)
    keep: list[object] = []
    for pad_index, (key, heading, scale, minimum, maximum) in enumerate(metrics, 1):
        pad = canvas.cd(pad_index)
        style_pad(pad, left=0.16, right=0.17, top=0.17, bottom=0.16)
        histogram = _category_heatmap(
            ROOT,
            f"matrix_{pad_index}",
            ";SPAD pitch [#mum];SPADs per readout-channel side",
        )
        for x_index, pitch in enumerate(PITCHES_UM, 1):
            for y_index, group_side in enumerate(SPADS_PER_CHANNEL_SIDE, 1):
                value = scale * float(lookup[(pitch, group_side)][key])
                histogram.SetBinContent(x_index, y_index, value)
        histogram.SetMinimum(minimum)
        if maximum is not None:
            histogram.SetMaximum(maximum)
        style_axes(histogram, title_size=0.047, label_size=0.038)
        histogram.GetZaxis().SetTitle(heading)
        histogram.GetZaxis().SetTitleOffset(1.30)
        histogram.SetMarkerColor(ROOT.kBlack)
        histogram.SetMarkerSize(0.62)
        histogram.Draw("COLZ TEXT")
        keep.append(histogram)
        keep.extend(
            draw_panel_title(
                ROOT,
                pad,
                heading,
                subtitle=(
                    f"{100 * PDE:g}% PDE, {SPAD_RECOVERY_NS:g} ns SPAD recovery, "
                    f"{WORKFLOW.reference.channel_resolving_time_ps:g} ps channel resolving time"
                ),
            )
        )
    canvas._keep = keep
    return save_canvas(canvas, output_dir / "channel_architecture_matrix")


def _graph(ROOT, x, y, *, color, marker):
    x = np.ascontiguousarray(x, dtype=np.float64)
    y = np.ascontiguousarray(y, dtype=np.float64)
    graph = ROOT.TGraph(len(x), x, y)
    graph.SetLineColor(color)
    graph.SetMarkerColor(color)
    graph.SetLineWidth(3)
    graph.SetMarkerStyle(marker)
    graph.SetMarkerSize(1.05)
    return graph


def _draw_channel_tradeoff(
    output_dir: Path,
    rows: list[dict[str, object]],
    sample: str,
) -> list[Path]:
    import ROOT

    palette = set_flash_style(ROOT)
    lookup = _lookup(rows, sample)
    selected_pitches = (10.0, 20.0, 30.0, 50.0)
    colors = (palette["blue"], palette["red"], palette["green"], palette["purple"])
    metrics = (
        ("timestamp_survival", "Timestamp retention [%]", 100.0, 0.0, 103.0),
        ("events_with_buffer_overflow_fraction", "Events with any buffer overflow [%]", 100.0, 0.0, 103.0),
        ("mean_payload_kbit", "Mean output payload [kbit / event]", 1.0, 0.0, None),
    )
    canvas = ROOT.TCanvas("channel_size_tradeoff", "channel size tradeoff", 1800, 620)
    canvas.Divide(3, 1, 0.004, 0.004)
    keep: list[object] = []
    legend = None
    for pad_index, (key, y_title, scale, y_min, y_max) in enumerate(metrics, 1):
        pad = canvas.cd(pad_index)
        style_pad(pad, left=0.15, right=0.04, top=0.17, bottom=0.17, grid_y=True)
        pad.SetLogx(True)
        values = [scale * float(row[key]) for row in lookup.values()]
        upper = y_max if y_max is not None else 1.20 * max(values)
        frame = ROOT.TH1D(
            f"channel_tradeoff_frame_{pad_index}",
            f";Nominal readout-channel side [#mum];{y_title}",
            100,
            8.0,
            4000.0,
        )
        frame.SetDirectory(0)
        frame.SetMinimum(y_min)
        frame.SetMaximum(upper)
        style_axes(frame, title_size=0.047, label_size=0.038)
        frame.Draw("AXIS")
        keep.append(frame)
        if pad_index == 1:
            legend = make_legend(ROOT, 0.18, 0.21, 0.64, 0.43, text_size=0.032, columns=2)
        for curve_index, (pitch, color) in enumerate(zip(selected_pitches, colors, strict=True)):
            selected = [lookup[(pitch, group)] for group in SPADS_PER_CHANNEL_SIDE]
            x = [float(row["nominal_channel_side_um"]) for row in selected]
            y = [scale * float(row[key]) for row in selected]
            graph = _graph(ROOT, x, y, color=color, marker=20 + curve_index)
            graph.Draw("LP SAME")
            keep.append(graph)
            if legend is not None and pad_index == 1:
                legend.AddEntry(graph, f"{pitch:g} #mum", "lp")
        if legend is not None and pad_index == 1:
            legend.Draw()
            keep.append(legend)
        keep.extend(
            draw_panel_title(
                ROOT,
                pad,
                ("Timing retention", "Overflow risk", "Payload cost")[pad_index - 1],
                subtitle="8 timestamp words and one 8-bit count per active channel",
            )
        )
    canvas._keep = keep
    return save_canvas(canvas, output_dir / "channel_size_tradeoff")


def run_study(
    output_dir: Path,
    *,
    max_events: int | None,
    timestamp_buffer_depth: int,
    channel_resolving_time_ns: float,
    count_bits: int,
) -> list[Path]:
    output_dir = output_dir.expanduser().resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    base_config = ChannelReadoutConfig(
        spads_per_channel_side=1,
        timestamp_buffer_depth=timestamp_buffer_depth,
        channel_resolving_time_ns=channel_resolving_time_ns,
        count_bits=count_bits,
    )
    base_config.validate()

    rows: list[dict[str, object]] = []
    for sample_index, (sample, path) in enumerate(SAMPLES.items()):
        rows.extend(
            _scan_sample(
                sample_index,
                sample,
                path,
                max_events=max_events,
                readout_config=base_config,
            )
        )
    table_path = output_dir / "channel_architecture_scan.csv"
    _write_csv(table_path, rows)
    focus_sample = "1 TeV pion"
    working_points_path = output_dir / "architecture_working_points.csv"
    _write_csv(
        working_points_path,
        _select_working_points(rows, focus_sample),
    )
    topology_path = output_dir / "common_resource_topology_projections.csv"
    _write_csv(
        topology_path,
        _project_device_topologies(rows, focus_sample),
    )

    plot_manifest = [
        {
            "figure": "channel_architecture_matrix.pdf",
            "question": "How do SPAD pitch and shared channel size jointly control counting, timing, overflow, and payload?",
            "takeaway": "SPAD recovery is set locally, while channel grouping primarily controls timestamp survival, overflow, and data volume.",
        },
        {
            "figure": "channel_size_tradeoff.pdf",
            "question": "What is gained and lost as the physical readout channel grows for several SPAD pitches?",
            "takeaway": "Larger shared channels reduce payload but create timing-buffer pressure in bright events.",
        },
    ]
    manifest_path = output_dir / "plot_manifest.csv"
    _write_csv(manifest_path, plot_manifest)

    output_paths = [
        table_path,
        working_points_path,
        topology_path,
        manifest_path,
    ]
    output_paths.extend(_draw_heatmaps(output_dir, rows, focus_sample))
    output_paths.extend(_draw_channel_tradeoff(output_dir, rows, focus_sample))

    summary = {
        "status": "architecture stress test",
        "scope": "one Cherenkov fiber; existing unmatched-energy occupancy samples",
        "sensor_workflow": {
            "path": str(WORKFLOW.source_path),
            "sha256": WORKFLOW.sha256,
            "schema_version": WORKFLOW.schema_version,
            "reference_scenario": WORKFLOW.reference.name,
        },
        "inputs": {
            sample: {
                "path": str(path),
                "bytes": path.stat().st_size,
                "sha256": _sha256(path),
            }
            for sample, path in SAMPLES.items()
        },
        "scan": {
            "spad_pitch_um": list(PITCHES_UM),
            "spads_per_channel_side": list(SPADS_PER_CHANNEL_SIDE),
            "nominal_channel_side_definition": "SPAD pitch times SPADs per channel side",
            "pde": PDE,
            "spad_recovery_ns": SPAD_RECOVERY_NS,
            "channel_readout": asdict(base_config),
            "seed": SEED,
            "max_events": max_events,
        },
        "interpretation": {
            "count_path": "Accepted SPAD avalanches are accumulated in one finite count word per active readout channel.",
            "timing_path": "Accepted avalanches first encounter the shared channel resolving time, then a finite timestamp buffer.",
            "payload": "One address and count word per active channel plus one word per retained timestamp.",
            "working_point_selection": {
                "sample": focus_sample,
                "objective": "lowest mean payload among scanned points that pass all requirements",
                "minimum_timestamp_survival": MIN_TIMESTAMP_SURVIVAL,
                "maximum_events_with_any_timestamp_overflow_fraction": MAX_OVERFLOW_EVENT_FRACTION,
                "minimum_count_survival": MIN_COUNT_SURVIVAL,
                "status": "study criteria, not frozen FLASH requirements",
            },
        },
        "limitations": [
            "The SPAD pitch is an effective square-cell pitch; junction diameter and fill factor are not yet independent parameters.",
            "PDE and recovery are fixed to expose the pitch-by-channel-size dependence.",
            (
                f"The {WORKFLOW.reference.timestamp_buffer_depth}-word timestamp "
                f"buffer, {WORKFLOW.reference.channel_resolving_time_ps:g} ps "
                f"resolving time, and {WORKFLOW.reference.count_bits}-bit counter "
                "are explicit stress-test assumptions, not a frozen FLASH design."
            ),
            "No arbitration, link bandwidth, power, DCR, crosstalk, or afterpulsing is modeled.",
            "The electron and pion files have different energies and are not a particle-species comparison.",
        ],
    }
    summary_path = output_dir / "study_summary.json"
    with summary_path.open("w") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")
    output_paths.append(summary_path)
    return output_paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_DIR / "results" / "channel-architecture-v2",
    )
    parser.add_argument("--max-events", type=int)
    parser.add_argument(
        "--timestamp-buffer-depth",
        type=int,
        default=WORKFLOW.reference.timestamp_buffer_depth,
    )
    parser.add_argument(
        "--channel-resolving-time-ns",
        type=float,
        default=WORKFLOW.reference.channel_resolving_time_ps / 1000.0,
    )
    parser.add_argument(
        "--count-bits",
        type=int,
        default=WORKFLOW.reference.count_bits,
    )
    args = parser.parse_args()
    outputs = run_study(
        args.output,
        max_events=args.max_events,
        timestamp_buffer_depth=args.timestamp_buffer_depth,
        channel_resolving_time_ns=args.channel_resolving_time_ns,
        count_bits=args.count_bits,
    )
    print("Wrote:")
    for output in outputs:
        print(f"  {output}")


if __name__ == "__main__":
    main()
