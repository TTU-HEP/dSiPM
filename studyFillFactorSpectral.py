#!/usr/bin/env python3
"""Scan SPAD pitch, active-area fill factor, and spectral response.

This study factorizes the front end into four explicit stages:

1. a photon lands inside the photosensitive region of its pitch cell;
2. the junction detects it with a wavelength-dependent probability;
3. the physical SPAD is available after its recovery time;
4. a shared readout channel retains the timestamp.

The spectral curve and square active region are requirement proxies.  They are
not measurements of the FLASH prototype.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, field
import hashlib
import json
from pathlib import Path

import numpy as np

from dsipm.active_area import (
    active_area_mask,
    active_width_um,
    junction_pde_probability,
)
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
FILL_FACTORS = WORKFLOW.scans.fill_factor
SENSOR_MIN_MM = -1.5
SENSOR_MAX_MM = 1.5
ACTIVE_SHAPE = WORKFLOW.factorized_proxy.active_area_shape
JUNCTION_PEAK_PDE = WORKFLOW.factorized_proxy.junction_peak_probability
SPECTRAL_MODEL = WORKFLOW.factorized_proxy.spectral_model
SPAD_RECOVERY_NS = WORKFLOW.reference.spad_recovery_ns
TARGET_CHANNEL_SIDE_UM = WORKFLOW.reference.target_readout_channel_side_um
SEED = WORKFLOW.reference.seed
WAVELENGTH_EDGES_NM = np.linspace(250.0, 850.0, 61)
FOCUS_SAMPLE = "1 TeV pion"
FOCUS_PITCH_UM = WORKFLOW.reference.spad_pitch_um
FOCUS_FILL_FACTOR = WORKFLOW.factorized_proxy.fill_factor


@dataclass
class _Accumulator:
    events: int = 0
    incident: int = 0
    active_area: int = 0
    detected: int = 0
    accepted: int = 0
    timestamps: int = 0
    junction_probability_sum: float = 0.0
    overflow_events: int = 0
    payload_kbit: list[float] = field(default_factory=list)


@dataclass
class _SpectralAccumulator:
    incident: np.ndarray = field(
        default_factory=lambda: np.zeros(len(WAVELENGTH_EDGES_NM) - 1, dtype=np.int64)
    )
    active_area: np.ndarray = field(
        default_factory=lambda: np.zeros(len(WAVELENGTH_EDGES_NM) - 1, dtype=np.int64)
    )
    detected: np.ndarray = field(
        default_factory=lambda: np.zeros(len(WAVELENGTH_EDGES_NM) - 1, dtype=np.int64)
    )
    accepted: np.ndarray = field(
        default_factory=lambda: np.zeros(len(WAVELENGTH_EDGES_NM) - 1, dtype=np.int64)
    )


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


def _ratio(numerator: int | float, denominator: int | float) -> float:
    return float(numerator / denominator) if denominator else 0.0


def _histogram(values: np.ndarray) -> np.ndarray:
    finite = values[np.isfinite(values)]
    return np.histogram(finite, bins=WAVELENGTH_EDGES_NM)[0].astype(np.int64)


def _scan_sample(
    sample_index: int,
    sample: str,
    path: Path,
    *,
    max_events: int | None,
) -> tuple[list[dict[str, object]], _SpectralAccumulator | None]:
    accumulators = {
        (pitch, fill_factor): _Accumulator()
        for pitch in PITCHES_UM
        for fill_factor in FILL_FACTORS
    }
    spectral = _SpectralAccumulator() if sample == FOCUS_SAMPLE else None
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
            np.random.SeedSequence([SEED, sample_index, event.event_index, 9187])
        )
        detection_uniform = random.random(event.selected_photon_count)
        junction_probability = junction_pde_probability(
            event.wavelength_nm,
            peak_probability=JUNCTION_PEAK_PDE,
            spectral_model=SPECTRAL_MODEL,
        )
        for pitch_um in PITCHES_UM:
            pitch_mm = pitch_um / 1000.0
            for fill_factor in FILL_FACTORS:
                active_mask = active_area_mask(
                    event.x_mm,
                    event.y_mm,
                    pitch_um=pitch_um,
                    fill_factor=fill_factor,
                    detector_min_mm=SENSOR_MIN_MM,
                    detector_max_mm=SENSOR_MAX_MM,
                    shape=ACTIVE_SHAPE,
                )
                detected_mask = active_mask & (
                    detection_uniform < junction_probability
                )
                detected_x = event.x_mm[detected_mask]
                detected_y = event.y_mm[detected_mask]
                detected_time = event.time_ns[detected_mask]
                detected_wavelength = event.wavelength_nm[detected_mask]
                ix, iy, cell_ids, n_spads_axis = _cell_indices(
                    detected_x,
                    detected_y,
                    pitch_mm,
                    SENSOR_MIN_MM,
                    SENSOR_MAX_MM,
                )
                accepted_mask = _apply_nonparalyzable_dead_time(
                    cell_ids,
                    detected_time,
                    SPAD_RECOVERY_NS,
                )
                group_side = max(1, int(round(TARGET_CHANNEL_SIDE_UM / pitch_um)))
                readout = channel_readout(
                    ix[accepted_mask],
                    iy[accepted_mask],
                    detected_time[accepted_mask],
                    n_spads_axis=n_spads_axis,
                    config=ChannelReadoutConfig(
                        spads_per_channel_side=group_side,
                        timestamp_buffer_depth=WORKFLOW.reference.timestamp_buffer_depth,
                        channel_resolving_time_ns=(
                            WORKFLOW.reference.channel_resolving_time_ps / 1000.0
                        ),
                        count_bits=WORKFLOW.reference.count_bits,
                    ),
                )

                accumulator = accumulators[(pitch_um, fill_factor)]
                accumulator.events += 1
                accumulator.incident += event.selected_photon_count
                accumulator.active_area += int(np.count_nonzero(active_mask))
                accumulator.detected += int(np.count_nonzero(detected_mask))
                accumulator.accepted += readout.accepted_hits
                accumulator.timestamps += readout.timestamp_records
                accumulator.junction_probability_sum += float(
                    np.sum(junction_probability)
                )
                accumulator.overflow_events += int(readout.any_buffer_overflow)
                accumulator.payload_kbit.append(readout.payload_bits / 1000.0)

                if (
                    spectral is not None
                    and pitch_um == FOCUS_PITCH_UM
                    and fill_factor == FOCUS_FILL_FACTOR
                ):
                    spectral.incident += _histogram(event.wavelength_nm)
                    spectral.active_area += _histogram(
                        event.wavelength_nm[active_mask]
                    )
                    spectral.detected += _histogram(detected_wavelength)
                    spectral.accepted += _histogram(
                        detected_wavelength[accepted_mask]
                    )

    rows: list[dict[str, object]] = []
    sensor_width_um = 1000.0 * (SENSOR_MAX_MM - SENSOR_MIN_MM)
    for pitch_um in PITCHES_UM:
        n_spads_axis = int(np.ceil(sensor_width_um / pitch_um))
        group_side = max(1, int(round(TARGET_CHANNEL_SIDE_UM / pitch_um)))
        channels_axis = int(np.ceil(n_spads_axis / group_side))
        for fill_factor in FILL_FACTORS:
            accumulator = accumulators[(pitch_um, fill_factor)]
            payload = np.asarray(accumulator.payload_kbit, dtype=np.float64)
            rows.append(
                {
                    "sample": sample,
                    "events": accumulator.events,
                    "pitch_um": pitch_um,
                    "active_shape": ACTIVE_SHAPE,
                    "nominal_fill_factor": fill_factor,
                    "active_width_um": active_width_um(
                        pitch_um,
                        fill_factor,
                        shape=ACTIVE_SHAPE,
                    ),
                    "spectral_model": SPECTRAL_MODEL,
                    "junction_peak_pde": JUNCTION_PEAK_PDE,
                    "spad_recovery_ns": SPAD_RECOVERY_NS,
                    "target_channel_side_um": TARGET_CHANNEL_SIDE_UM,
                    "spads_per_channel_side": group_side,
                    "realized_channel_side_um": group_side * pitch_um,
                    "readout_channels_per_sensor": channels_axis**2,
                    "incident_photons": accumulator.incident,
                    "active_area_photons": accumulator.active_area,
                    "detected_photons": accumulator.detected,
                    "accepted_spad_hits": accumulator.accepted,
                    "timestamp_records": accumulator.timestamps,
                    "geometric_acceptance": _ratio(
                        accumulator.active_area,
                        accumulator.incident,
                    ),
                    "mean_junction_pde_over_incident_spectrum": _ratio(
                        accumulator.junction_probability_sum,
                        accumulator.incident,
                    ),
                    "junction_detection_given_active": _ratio(
                        accumulator.detected,
                        accumulator.active_area,
                    ),
                    "effective_pde_detected_over_incident": _ratio(
                        accumulator.detected,
                        accumulator.incident,
                    ),
                    "recovery_survival": _ratio(
                        accumulator.accepted,
                        accumulator.detected,
                    ),
                    "timestamp_survival_after_spad": _ratio(
                        accumulator.timestamps,
                        accumulator.accepted,
                    ),
                    "end_to_end_timestamp_yield": _ratio(
                        accumulator.timestamps,
                        accumulator.incident,
                    ),
                    "events_with_any_timestamp_overflow_fraction": _ratio(
                        accumulator.overflow_events,
                        accumulator.events,
                    ),
                    "mean_payload_kbit": float(np.mean(payload)) if payload.size else 0.0,
                    "p90_payload_kbit": float(np.quantile(payload, 0.90)) if payload.size else 0.0,
                }
            )
    return rows, spectral


def _spectral_rows(spectral: _SpectralAccumulator) -> list[dict[str, object]]:
    rows = []
    centers = 0.5 * (WAVELENGTH_EDGES_NM[:-1] + WAVELENGTH_EDGES_NM[1:])
    probability = junction_pde_probability(
        centers,
        peak_probability=JUNCTION_PEAK_PDE,
        spectral_model=SPECTRAL_MODEL,
    )
    for index, center in enumerate(centers):
        incident = int(spectral.incident[index])
        rows.append(
            {
                "wavelength_low_nm": WAVELENGTH_EDGES_NM[index],
                "wavelength_high_nm": WAVELENGTH_EDGES_NM[index + 1],
                "wavelength_center_nm": center,
                "incident_photons": incident,
                "nominal_junction_pde": probability[index],
                "geometric_acceptance": _ratio(
                    spectral.active_area[index], incident
                ),
                "detected_over_incident": _ratio(
                    spectral.detected[index], incident
                ),
                "accepted_over_incident": _ratio(
                    spectral.accepted[index], incident
                ),
            }
        )
    return rows


def _lookup(rows: list[dict[str, object]]) -> dict[tuple[float, float], dict[str, object]]:
    return {
        (float(row["pitch_um"]), float(row["nominal_fill_factor"])): row
        for row in rows
        if row["sample"] == FOCUS_SAMPLE
    }


def _heatmap(ROOT, name: str):
    histogram = ROOT.TH2D(
        name,
        ";SPAD pitch [#mum];Nominal active-area fill factor",
        len(PITCHES_UM),
        0.5,
        len(PITCHES_UM) + 0.5,
        len(FILL_FACTORS),
        0.5,
        len(FILL_FACTORS) + 0.5,
    )
    histogram.SetDirectory(0)
    for index, pitch in enumerate(PITCHES_UM, 1):
        histogram.GetXaxis().SetBinLabel(index, f"{pitch:g}")
    for index, fill_factor in enumerate(FILL_FACTORS, 1):
        histogram.GetYaxis().SetBinLabel(index, f"{100 * fill_factor:.0f}%")
    return histogram


def _draw_matrix(output_dir: Path, rows: list[dict[str, object]]) -> list[Path]:
    import ROOT

    set_flash_style(ROOT)
    ROOT.gStyle.SetPalette(ROOT.kViridis)
    ROOT.gStyle.SetPaintTextFormat(".1f")
    lookup = _lookup(rows)
    metrics = (
        ("geometric_acceptance", "Photons on active area [%]", 0.0, 100.0),
        ("effective_pde_detected_over_incident", "Detected / incident [%]", 0.0, 30.0),
        ("recovery_survival", "Recovery survival [%]", 50.0, 100.0),
        ("end_to_end_timestamp_yield", "Exported timestamps / incident [%]", 0.0, 30.0),
    )
    canvas = ROOT.TCanvas("fill_factor_matrix", "fill factor matrix", 1800, 1040)
    canvas.Divide(2, 2, 0.004, 0.004)
    keep: list[object] = []
    for pad_index, (key, heading, minimum, maximum) in enumerate(metrics, 1):
        pad = canvas.cd(pad_index)
        style_pad(pad, left=0.16, right=0.17, top=0.17, bottom=0.16)
        histogram = _heatmap(ROOT, f"fill_matrix_{pad_index}")
        for x_index, pitch in enumerate(PITCHES_UM, 1):
            for y_index, fill_factor in enumerate(FILL_FACTORS, 1):
                histogram.SetBinContent(
                    x_index,
                    y_index,
                    100.0 * float(lookup[(pitch, fill_factor)][key]),
                )
        histogram.SetMinimum(minimum)
        histogram.SetMaximum(maximum)
        style_axes(histogram, title_size=0.047, label_size=0.038)
        histogram.GetZaxis().SetTitle(heading)
        histogram.SetMarkerColor(ROOT.kBlack)
        histogram.SetMarkerSize(0.68)
        histogram.Draw("COLZ TEXT")
        keep.append(histogram)
        keep.extend(
            draw_panel_title(
                ROOT,
                pad,
                heading,
                subtitle=(
                    f"{100 * JUNCTION_PEAK_PDE:g}% peak junction PDE; spectral proxy; "
                    f"{SPAD_RECOVERY_NS:g} ns recovery"
                ),
            )
        )
    canvas._keep = keep
    return save_canvas(canvas, output_dir / "fill_factor_pitch_matrix")


def _graph(ROOT, x, y, *, color, marker, line_style=1):
    x_values = np.ascontiguousarray(x, dtype=np.float64)
    y_values = np.ascontiguousarray(y, dtype=np.float64)
    graph = ROOT.TGraph(len(x_values), x_values, y_values)
    graph.SetLineColor(color)
    graph.SetMarkerColor(color)
    graph.SetLineWidth(3)
    graph.SetLineStyle(line_style)
    graph.SetMarkerStyle(marker)
    graph.SetMarkerSize(0.95)
    return graph


def _draw_spectrum(
    output_dir: Path,
    spectral_rows: list[dict[str, object]],
) -> list[Path]:
    import ROOT

    palette = set_flash_style(ROOT)
    selected = [
        row
        for row in spectral_rows
        if 300.0 <= float(row["wavelength_center_nm"]) <= 620.0
        and int(row["incident_photons"]) > 0
    ]
    x = [float(row["wavelength_center_nm"]) for row in selected]
    incident = np.asarray(
        [float(row["incident_photons"]) for row in selected],
        dtype=np.float64,
    )
    incident /= np.sum(incident)
    canvas = ROOT.TCanvas("spectral_response_proxy", "spectral response proxy", 1500, 650)
    canvas.Divide(2, 1, 0.004, 0.004)
    keep: list[object] = []

    pad = canvas.cd(1)
    style_pad(pad, left=0.14, right=0.04, top=0.17, bottom=0.16, grid_y=True)
    frame = ROOT.TH1D(
        "spectrum_frame",
        ";Photon wavelength [nm];Fraction of incident photons / 10 nm",
        100,
        295.0,
        625.0,
    )
    frame.SetDirectory(0)
    frame.SetMinimum(0.0)
    frame.SetMaximum(1.18 * float(np.max(incident)))
    style_axes(frame, title_size=0.047, label_size=0.038)
    frame.Draw("AXIS")
    spectrum = _graph(ROOT, x, incident, color=palette["blue"], marker=20)
    spectrum.Draw("LP SAME")
    keep.extend([frame, spectrum])
    keep.extend(
        draw_panel_title(
            ROOT,
            pad,
            "HG-DREAM photon spectrum",
            subtitle="1 TeV pion one-fiber stress sample",
        )
    )

    pad = canvas.cd(2)
    style_pad(pad, left=0.14, right=0.04, top=0.17, bottom=0.16, grid_y=True)
    frame = ROOT.TH1D(
        "spectral_stage_frame",
        ";Photon wavelength [nm];Fraction of incident photons [%]",
        100,
        295.0,
        625.0,
    )
    frame.SetDirectory(0)
    frame.SetMinimum(0.0)
    frame.SetMaximum(78.0)
    style_axes(frame, title_size=0.047, label_size=0.038)
    frame.Draw("AXIS")
    keep.append(frame)
    curves = (
        ("geometric_acceptance", "on active area", palette["gray"], 24),
        ("nominal_junction_pde", "junction PDE proxy", palette["green"], 25),
        ("detected_over_incident", "detected", palette["blue"], 20),
        ("accepted_over_incident", "after SPAD recovery", palette["red"], 21),
    )
    legend = make_legend(ROOT, 0.19, 0.62, 0.67, 0.84, text_size=0.033)
    for key, label, color, marker in curves:
        y = [100.0 * float(row[key]) for row in selected]
        graph = _graph(ROOT, x, y, color=color, marker=marker)
        graph.Draw("LP SAME")
        legend.AddEntry(graph, label, "lp")
        keep.append(graph)
    legend.Draw()
    keep.append(legend)
    keep.extend(
        draw_panel_title(
            ROOT,
            pad,
            "Factorized response stages",
            subtitle=(
                f"{FOCUS_PITCH_UM:g} #mum pitch, "
                f"{100 * FOCUS_FILL_FACTOR:g}% fill, "
                f"{100 * JUNCTION_PEAK_PDE:g}% peak junction PDE"
            ),
        )
    )
    canvas._keep = keep
    return save_canvas(canvas, output_dir / "spectral_response_proxy")


def run_study(output_dir: Path, *, max_events: int | None) -> list[Path]:
    output_dir = output_dir.expanduser().resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    focus_spectral: _SpectralAccumulator | None = None
    for sample_index, (sample, path) in enumerate(SAMPLES.items()):
        sample_rows, spectral = _scan_sample(
            sample_index,
            sample,
            path,
            max_events=max_events,
        )
        rows.extend(sample_rows)
        if spectral is not None:
            focus_spectral = spectral
    if focus_spectral is None:
        raise RuntimeError("Focus spectral accumulator was not produced")

    scan_path = output_dir / "fill_factor_spectral_scan.csv"
    _write_csv(scan_path, rows)
    spectral_rows = _spectral_rows(focus_spectral)
    spectrum_path = output_dir / "spectral_response_by_wavelength.csv"
    _write_csv(spectrum_path, spectral_rows)

    manifest_rows = [
        {
            "figure": "fill_factor_pitch_matrix.pdf",
            "question": "How do physical pitch and active-area fill factor jointly change yield, recovery, and exported timestamps?",
            "takeaway": "Fill-factor loss and collision loss must be reported separately; lower illumination can increase recovery survival while reducing usable yield.",
        },
        {
            "figure": "spectral_response_proxy.pdf",
            "question": "How does the stored HG-DREAM wavelength distribution pass through a factorized active-area and junction-PDE model?",
            "takeaway": "A single PDE number hides both spectral response and geometric active area; the curve remains a proxy until FLASH measurements exist.",
        },
    ]
    manifest_path = output_dir / "plot_manifest.csv"
    _write_csv(manifest_path, manifest_rows)

    outputs = [scan_path, spectrum_path, manifest_path]
    outputs.extend(_draw_matrix(output_dir, rows))
    outputs.extend(_draw_spectrum(output_dir, spectral_rows))

    summary = {
        "status": "factorized active-area and spectral-response requirement study",
        "sensor_workflow": {
            "path": str(WORKFLOW.source_path),
            "sha256": WORKFLOW.sha256,
            "schema_version": WORKFLOW.schema_version,
            "factorized_proxy": WORKFLOW.factorized_proxy.name,
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
            "pitch_um": list(PITCHES_UM),
            "fill_factor": list(FILL_FACTORS),
            "active_shape": ACTIVE_SHAPE,
            "junction_peak_pde": JUNCTION_PEAK_PDE,
            "spectral_model": SPECTRAL_MODEL,
            "spad_recovery_ns": SPAD_RECOVERY_NS,
            "target_channel_side_um": TARGET_CHANNEL_SIDE_UM,
            "timestamp_buffer_depth": WORKFLOW.reference.timestamp_buffer_depth,
            "channel_resolving_time_ps": WORKFLOW.reference.channel_resolving_time_ps,
            "max_events": max_events,
            "seed": SEED,
        },
        "focus_spectrum": {
            "sample": FOCUS_SAMPLE,
            "pitch_um": FOCUS_PITCH_UM,
            "fill_factor": FOCUS_FILL_FACTOR,
        },
        "limitations": [
            "The centered square active region is a fill-factor proxy, not a fabricated guard-ring or microlens layout.",
            "The normalized capstone polynomial supplies only a traceable spectral shape; it is not a FLASH PDE measurement.",
            (
                f"The {100 * JUNCTION_PEAK_PDE:g}% peak probability is conditional "
                "on landing in active area and is a study assumption."
            ),
            "Wavelength-dependent attenuation, group velocity, instrument response, incidence angle, and optical coupling are not yet modeled.",
            "Noise, crosstalk, afterpulsing, masking, power, and link bandwidth are not included.",
        ],
    }
    summary_path = output_dir / "study_summary.json"
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
        default=REPO_DIR / "results" / "fill-factor-spectral-v1",
    )
    parser.add_argument("--max-events", type=int)
    args = parser.parse_args()
    outputs = run_study(args.output, max_events=args.max_events)
    print("Wrote:")
    for output in outputs:
        print(f"  {output}")


if __name__ == "__main__":
    main()
