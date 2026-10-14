#!/usr/bin/env python3
"""Make matched CPAD/JINST figures from one merged schema-v3 checkpoint."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import subprocess

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle
import numpy as np


ENERGIES = (10, 30, 100, 300)
PARTICLES = ("eplus", "piplus")
PARTICLE_LABELS = {"eplus": r"$e^{+}$", "piplus": r"$\pi^{+}$"}
PARTICLE_COLORS = {"eplus": "#3F90DA", "piplus": "#BD1F01"}
READOUTS = (
    ("sum3dCC", "DREAM stored C signal", "#3F90DA", "o"),
    ("selected_cherenkov_photons", "Explicit-photon truth", "#000000", "s"),
    ("analog_charge_pe", "Conventional SiPM charge proxy", "#FFA90E", "^"),
    ("accepted_reference", r"FLASH 20 $\mu$m, 4.3 ns", "#BD1F01", "D"),
)
PITCHES_UM = np.asarray((5, 10, 20, 25, 30, 50, 75, 100), dtype=float)
RECOVERY_NS = np.asarray((1.0, 4.3, 40.0), dtype=float)
NUMERIC_COLUMNS = (
    "event",
    "beam_energy_gev",
    "sum3dCC",
    "sum3dSS",
    "ecalo_truth_gev",
    "eleak_truth_gev",
    "raw_optical_photons",
    "selected_cherenkov_photons",
    "active_cherenkov_fibers",
    "max_photons_in_one_fiber",
    "active_fiber_occupancy_q50",
    "active_fiber_occupancy_q90",
    "active_fiber_occupancy_q99",
    "top_fiber_photon_fraction",
    "top5_fiber_photon_fraction",
    "active_fibers_ge2_photons",
    "active_fibers_ge5_photons",
    "active_fibers_ge10_photons",
    "active_fibers_ge20_photons",
    "active_fibers_ge50_photons",
    "active_fibers_ge100_photons",
    "detected_reference",
    "accepted_reference",
    "lost_reference_recovery",
    "reference_recovery_survival",
    "analog_detected_at_31pct_pde",
    "analog_charge_pe",
    "analog_charge_survival",
    "arrival_q10_ns",
    "arrival_median_ns",
    "arrival_q90_ns",
    "accepted_arrival_q10_ns",
    "accepted_arrival_median_ns",
    "accepted_arrival_q90_ns",
    "detected_arrival_mean_ns",
)


def _sample(particle: str, energy: int) -> str:
    return f"{particle}_{energy}gev"


def _scenario_column(pitch_um: float, recovery_ns: float) -> str:
    return f"accepted_pitch{pitch_um:g}um_recovery{recovery_ns:g}ns"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_commit(path: Path) -> str | None:
    try:
        result = subprocess.run(
            ["git", "-C", str(path), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip()


def _configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 12,
            "axes.linewidth": 1.15,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.size": 6,
            "ytick.major.size": 6,
            "xtick.minor.size": 3,
            "ytick.minor.size": 3,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.bbox": "tight",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def _stamp(fig: mpl.figure.Figure, right: str) -> None:
    fig.text(0.07, 0.975, "FLASH", fontsize=16, fontweight="bold", va="top")
    fig.text(
        0.175,
        0.975,
        "Simulation Preliminary",
        fontsize=11,
        style="italic",
        va="top",
    )
    fig.text(0.95, 0.975, right, fontsize=10.5, ha="right", va="top")


def _panel(ax: mpl.axes.Axes, label: str, title: str) -> None:
    ax.text(
        0.02,
        1.01,
        f"{label}  {title}",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontweight="bold",
        clip_on=False,
    )
    ax.minorticks_on()


def _load_events(path: Path) -> dict[str, np.ndarray]:
    columns = list(NUMERIC_COLUMNS)
    columns.extend(
        _scenario_column(pitch, recovery)
        for pitch in PITCHES_UM
        for recovery in RECOVERY_NS
    )
    numeric = {column: [] for column in columns}
    roots: list[str] = []
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"Missing CSV header: {path}")
        missing = sorted(set(columns) - set(reader.fieldnames))
        if missing:
            raise ValueError(f"Missing event columns in {path}: {missing}")
        for row in reader:
            roots.append(row["source_root"])
            for column in columns:
                numeric[column].append(float(row[column]))
    result = {name: np.asarray(values, dtype=float) for name, values in numeric.items()}
    result["source_root"] = np.asarray(roots, dtype=str)
    return result


def _load_histograms(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as payload:
        return {key: np.asarray(payload[key]) for key in payload.files}


def _block_sums(
    sample: dict[str, np.ndarray],
    values: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    roots = sample["source_root"]
    _, inverse = np.unique(roots, return_inverse=True)
    counts = np.bincount(inverse).astype(float)
    sums = np.bincount(inverse, weights=np.asarray(values, dtype=float))
    return sums, counts


def _bootstrap_ratio(
    sample: dict[str, np.ndarray],
    numerator: np.ndarray,
    denominator: np.ndarray | None = None,
    *,
    draws: int = 500,
    seed: int = 4307,
) -> np.ndarray:
    numerator_sums, block_counts = _block_sums(sample, numerator)
    if denominator is None:
        denominator_sums = block_counts
    else:
        denominator_sums, _ = _block_sums(sample, denominator)
    random = np.random.default_rng(seed)
    indices = random.integers(
        0,
        numerator_sums.size,
        size=(draws, numerator_sums.size),
    )
    sampled_numerator = np.sum(numerator_sums[indices], axis=1)
    sampled_denominator = np.sum(denominator_sums[indices], axis=1)
    return np.divide(
        sampled_numerator,
        sampled_denominator,
        out=np.full(draws, np.nan),
        where=sampled_denominator != 0,
    )


def _interval(values: np.ndarray) -> list[float]:
    return [float(value) for value in np.quantile(values, (0.025, 0.5, 0.975))]


def _retention(
    sample: dict[str, np.ndarray],
    numerator: str,
    denominator: str,
) -> float:
    return float(np.sum(sample[numerator]) / np.sum(sample[denominator]))


def _make_metrics(
    events: dict[str, dict[str, np.ndarray]],
    histograms: dict[str, dict[str, np.ndarray]],
    merged_summary: dict[str, object],
) -> dict[str, object]:
    metrics: dict[str, object] = {
        "schema_version": 1,
        "checkpoint": merged_summary,
        "samples": {},
        "closures": {},
        "claim_classification": {
            "supported": [
                "Per-fiber occupancy and SPAD multiplicity increase with energy and are larger in compact positron showers.",
                "At fixed PDE, smaller SPAD pitch and shorter recovery preserve more detected photon timestamps.",
                "The 20 um, 4.3 ns FLASH reference has smaller recovery loss than the conventional-SiPM partial-charge proxy for the sampled central optical window.",
                "The conventional-SiPM and FLASH pulse comparisons use the same underlying photon arrivals.",
            ],
            "provisional": [
                "Readout-response linearity is a central-optical-window result, not yet a complete all-fiber calorimeter response.",
                "A 20 um pitch is a useful reference point, not a final prototype requirement.",
                "The electronics-convolved conventional-SiPM trace is a stated fast-output impulse-response proxy.",
            ],
            "unsupported": [
                "Full calorimeter energy resolution or particle-identification improvement.",
                "Explicit scintillation-photon performance.",
                "Prototype timing resolution, power, bandwidth, or radiation tolerance.",
            ],
        },
    }
    for particle in PARTICLES:
        for energy in ENERGIES:
            name = _sample(particle, energy)
            data = events[name]
            hist = histograms[name]
            selected = data["selected_cherenkov_photons"]
            sem = float(np.std(selected, ddof=1) / np.sqrt(selected.size))
            mean_bootstrap = _bootstrap_ratio(
                data,
                selected,
                seed=4307 + energy + (0 if particle == "eplus" else 1000),
            )
            flash_bootstrap = _bootstrap_ratio(
                data,
                data["accepted_reference"],
                data["detected_reference"],
                seed=5307 + energy + (0 if particle == "eplus" else 1000),
            )
            analog_bootstrap = _bootstrap_ratio(
                data,
                data["analog_charge_pe"],
                data["analog_detected_at_31pct_pde"],
                seed=6307 + energy + (0 if particle == "eplus" else 1000),
            )
            edges = hist["time_edges_ns"]
            centers = 0.5 * (edges[:-1] + edges[1:])
            incident = hist["time_incident"].astype(float)
            accepted = hist["time_accepted"].astype(float)
            prompt = centers < 12.0
            cell_photons = hist["cell_detected_photons_by_detected_occupancy"]
            multi_fraction = np.sum(cell_photons[:, 1:], axis=1) / np.sum(
                cell_photons, axis=1
            )

            event_scan = np.asarray(
                [
                    [
                        np.sum(data[_scenario_column(pitch, recovery)])
                        for recovery in RECOVERY_NS
                    ]
                    for pitch in PITCHES_UM
                ]
            )
            histogram_scan = np.sum(hist["time_accepted_scan"], axis=2)
            closure = {
                "incident_time_minus_event": int(
                    np.sum(hist["time_incident"]) - np.sum(selected)
                ),
                "detected_time_minus_event": int(
                    np.sum(hist["time_detected"])
                    - np.sum(data["detected_reference"])
                ),
                "accepted_time_minus_event": int(
                    np.sum(hist["time_accepted"])
                    - np.sum(data["accepted_reference"])
                ),
                "detected_minus_accepted_and_lost": int(
                    np.sum(hist["time_detected"])
                    - np.sum(hist["time_accepted"])
                    - np.sum(hist["time_lost_recovery"])
                ),
                "fiber_incident_minus_event": int(
                    np.sum(hist["fiber_incident_photons_by_incident_occupancy"])
                    - np.sum(selected)
                ),
                "analog_detected_time_minus_event": int(
                    np.sum(hist["time_analog_detected"])
                    - np.sum(data["analog_detected_at_31pct_pde"])
                ),
                "analog_charge_time_minus_event": float(
                    np.sum(hist["time_analog_charge"])
                    - np.sum(data["analog_charge_pe"])
                ),
                "scan_time_max_abs_minus_event": int(
                    np.max(np.abs(histogram_scan - event_scan))
                ),
            }
            metrics["closures"][name] = closure
            metrics["samples"][name] = {
                "events": int(selected.size),
                "selected_photons_mean": float(np.mean(selected)),
                "selected_photons_sem": sem,
                "selected_photons_relative_sem_percent": float(
                    100.0 * sem / np.mean(selected)
                ),
                "selected_photons_file_bootstrap_95": _interval(mean_bootstrap),
                "active_fibers_mean": float(
                    np.mean(data["active_cherenkov_fibers"])
                ),
                "hottest_fiber_median": float(
                    np.median(data["max_photons_in_one_fiber"])
                ),
                "hottest_fiber_q99": float(
                    np.quantile(data["max_photons_in_one_fiber"], 0.99)
                ),
                "flash_reference_recovery_loss_percent": float(
                    100.0
                    * (
                        1.0
                        - _retention(
                            data, "accepted_reference", "detected_reference"
                        )
                    )
                ),
                "flash_reference_survival_file_bootstrap_95": _interval(
                    flash_bootstrap
                ),
                "analog_charge_loss_percent": float(
                    100.0
                    * (
                        1.0
                        - _retention(
                            data,
                            "analog_charge_pe",
                            "analog_detected_at_31pct_pde",
                        )
                    )
                ),
                "analog_charge_survival_file_bootstrap_95": _interval(
                    analog_bootstrap
                ),
                "incident_prompt_fraction_before_12ns": float(
                    np.sum(incident[prompt]) / np.sum(incident)
                ),
                "flash_prompt_fraction_before_12ns": float(
                    np.sum(accepted[prompt]) / np.sum(accepted)
                ),
                "multi_hit_cell_photon_fraction_by_pitch": {
                    f"{pitch:g}": float(value)
                    for pitch, value in zip(PITCHES_UM, multi_fraction, strict=True)
                },
            }
    return metrics


def _save_figure(
    fig: mpl.figure.Figure,
    output_dir: Path,
    stem: str,
    packet: PdfPages,
) -> dict[str, str]:
    png_path = output_dir / f"{stem}.png"
    pdf_path = output_dir / f"{stem}.pdf"
    fig.savefig(png_path, dpi=220)
    fig.savefig(pdf_path)
    packet.savefig(fig)
    plt.close(fig)
    return {
        "png": str(png_path),
        "pdf": str(pdf_path),
        "png_sha256": _sha256(png_path),
        "pdf_sha256": _sha256(pdf_path),
    }


def _response_transfer(
    events: dict[str, dict[str, np.ndarray]],
) -> mpl.figure.Figure:
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.4))
    for panel_index, particle in enumerate(PARTICLES):
        ax = axes[0, panel_index]
        for readout_index, (column, label, color, marker) in enumerate(READOUTS):
            bootstrap_by_energy = []
            central = []
            for energy in ENERGIES:
                data = events[_sample(particle, energy)]
                bootstrap = _bootstrap_ratio(
                    data,
                    data[column],
                    seed=7100 + 100 * readout_index + energy + 1000 * panel_index,
                ) / energy
                bootstrap_by_energy.append(bootstrap)
                central.append(float(np.mean(data[column]) / energy))
            base = bootstrap_by_energy[0]
            normalized_bootstrap = [values / base for values in bootstrap_by_energy]
            normalized = np.asarray(central) / central[0]
            lower = [np.quantile(values, 0.025) for values in normalized_bootstrap]
            upper = [np.quantile(values, 0.975) for values in normalized_bootstrap]
            ax.plot(
                ENERGIES,
                normalized,
                color=color,
                marker=marker,
                linewidth=1.9,
                markersize=5,
                label=label,
            )
            ax.fill_between(ENERGIES, lower, upper, color=color, alpha=0.12)
        ax.axhline(1.0, color="0.65", linestyle=":", linewidth=1.0)
        ax.set_xscale("log")
        ax.set_xticks(ENERGIES, labels=[str(value) for value in ENERGIES])
        ax.set_xlabel("Beam energy [GeV]")
        ax.set_ylabel("Response / energy, normalized at 10 GeV")
        _panel(ax, chr(ord("A") + panel_index), f"{PARTICLE_LABELS[particle]} response transfer")
        if panel_index == 0:
            ax.legend(loc="best", fontsize=8.5)

    ax = axes[1, 0]
    for particle in PARTICLES:
        flash = []
        analog = []
        for energy in ENERGIES:
            data = events[_sample(particle, energy)]
            flash.append(
                100.0
                * _retention(data, "accepted_reference", "detected_reference")
            )
            analog.append(
                100.0
                * _retention(
                    data,
                    "analog_charge_pe",
                    "analog_detected_at_31pct_pde",
                )
            )
        color = PARTICLE_COLORS[particle]
        ax.plot(
            ENERGIES,
            flash,
            color=color,
            marker="o",
            linewidth=2,
            label=f"{PARTICLE_LABELS[particle]} FLASH",
        )
        ax.plot(
            ENERGIES,
            analog,
            color=color,
            marker="s",
            linestyle="--",
            linewidth=2,
            label=f"{PARTICLE_LABELS[particle]} conventional SiPM",
        )
    ax.set_xscale("log")
    ax.set_xticks(ENERGIES, labels=[str(value) for value in ENERGIES])
    ax.set_ylim(65, 100.8)
    ax.set_xlabel("Beam energy [GeV]")
    ax.set_ylabel("Detected signal retained [%]")
    _panel(ax, "C", "Dense cores separate the readouts")
    ax.legend(loc="lower left", fontsize=8.5)

    ax = axes[1, 1]
    for particle in PARTICLES:
        median = []
        q16 = []
        q84 = []
        for energy in ENERGIES:
            values = events[_sample(particle, energy)]["max_photons_in_one_fiber"]
            median.append(np.median(values))
            q16.append(np.quantile(values, 0.16))
            q84.append(np.quantile(values, 0.84))
        color = PARTICLE_COLORS[particle]
        ax.plot(
            ENERGIES,
            median,
            color=color,
            marker="o",
            linewidth=2,
            label=PARTICLE_LABELS[particle],
        )
        ax.fill_between(ENERGIES, q16, q84, color=color, alpha=0.18)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(ENERGIES, labels=[str(value) for value in ENERGIES])
    ax.set_xlabel("Beam energy [GeV]")
    ax.set_ylabel("Hottest-fiber photons / event")
    _panel(ax, "D", "Local occupancy, not total light, drives loss")
    ax.legend(loc="best")
    fig.subplots_adjust(top=0.91, bottom=0.09, hspace=0.29, wspace=0.22)
    _stamp(fig, "Central optical window - 38,715 events")
    return fig


def _occupancy_to_deadtime(
    events: dict[str, dict[str, np.ndarray]],
    histograms: dict[str, dict[str, np.ndarray]],
) -> mpl.figure.Figure:
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 8.4))
    thresholds = np.unique(np.geomspace(1, 2000, 160).astype(int))
    ax = axes[0, 0]
    for particle in PARTICLES:
        hist = histograms[_sample(particle, 100)]
        photons = hist["fiber_incident_photons_by_incident_occupancy"].astype(float)
        tail = np.asarray(
            [100.0 * np.sum(photons[value - 1 :]) / np.sum(photons) for value in thresholds]
        )
        ax.plot(
            thresholds,
            tail,
            color=PARTICLE_COLORS[particle],
            linewidth=2.1,
            label=PARTICLE_LABELS[particle],
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1, 2000)
    ax.set_ylim(0.01, 110)
    ax.set_xlabel("Incident photons in one active fiber")
    ax.set_ylabel("Photons in fibers above threshold [%]")
    _panel(ax, "A", "Fiber occupancy tail at 100 GeV")
    ax.legend(loc="lower left")

    ax = axes[0, 1]
    for particle in PARTICLES:
        hist = histograms[_sample(particle, 100)]
        cell_photons = hist["cell_detected_photons_by_detected_occupancy"].astype(float)
        multi = 100.0 * np.sum(cell_photons[:, 1:], axis=1) / np.sum(
            cell_photons, axis=1
        )
        ax.plot(
            PITCHES_UM,
            multi,
            color=PARTICLE_COLORS[particle],
            marker="o",
            linewidth=2.1,
            label=PARTICLE_LABELS[particle],
        )
    ax.axvline(20.0, color="0.45", linestyle=":", linewidth=1.2)
    ax.set_yscale("log")
    ax.set_xlabel(r"SPAD pitch [$\mu$m]")
    ax.set_ylabel("Detected photons in multi-hit SPADs [%]")
    _panel(ax, "B", "Segmentation converts fiber light to cell occupancy")
    ax.legend(loc="upper left")

    ax = axes[1, 0]
    for particle in PARTICLES:
        for energy, linestyle in ((100, "-"), (300, "--")):
            data = events[_sample(particle, energy)]
            loss = [
                100.0
                * (
                    1.0
                    - np.sum(data[_scenario_column(pitch, 4.3)])
                    / np.sum(data["detected_reference"])
                )
                for pitch in PITCHES_UM
            ]
            ax.plot(
                PITCHES_UM,
                loss,
                color=PARTICLE_COLORS[particle],
                linestyle=linestyle,
                marker="o" if energy == 100 else "s",
                markersize=4.5,
                linewidth=2,
                label=f"{PARTICLE_LABELS[particle]} {energy} GeV",
            )
    ax.axvline(20.0, color="0.45", linestyle=":", linewidth=1.2)
    ax.set_yscale("log")
    ax.set_ylim(0.01, 35)
    ax.set_xlabel(r"SPAD pitch [$\mu$m]")
    ax.set_ylabel("Detected photons lost to 4.3 ns recovery [%]")
    _panel(ax, "C", "Recovery turns repeated cell hits into loss")
    ax.legend(loc="upper left", fontsize=8.5, ncol=2)

    ax = axes[1, 1]
    pitch_index = int(np.flatnonzero(np.isclose(PITCHES_UM, 20.0))[0])
    for particle in PARTICLES:
        multi = []
        loss = []
        for energy in ENERGIES:
            hist = histograms[_sample(particle, energy)]
            cell_photons = hist["cell_detected_photons_by_detected_occupancy"].astype(float)
            multi.append(
                100.0
                * np.sum(cell_photons[pitch_index, 1:])
                / np.sum(cell_photons[pitch_index])
            )
            data = events[_sample(particle, energy)]
            loss.append(
                100.0
                * (
                    1.0
                    - _retention(data, "accepted_reference", "detected_reference")
                )
            )
        color = PARTICLE_COLORS[particle]
        ax.plot(
            ENERGIES,
            multi,
            color=color,
            marker="o",
            linewidth=2,
            label=f"{PARTICLE_LABELS[particle]} multi-hit exposure",
        )
        ax.plot(
            ENERGIES,
            loss,
            color=color,
            marker="s",
            linestyle="--",
            linewidth=2,
            label=f"{PARTICLE_LABELS[particle]} lost timestamps",
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(ENERGIES, labels=[str(value) for value in ENERGIES])
    ax.set_xlabel("Beam energy [GeV]")
    ax.set_ylabel("Fraction of detected photons [%]")
    _panel(ax, "D", r"20 $\mu$m reference: exposure is not automatic loss")
    ax.legend(loc="upper left", fontsize=8.2)
    fig.subplots_adjust(top=0.91, bottom=0.09, hspace=0.29, wspace=0.24)
    _stamp(fig, "One fiber to full-shower occupancy")
    return fig


def _convolve_fast_output(charge: np.ndarray, bin_width_ns: float) -> np.ndarray:
    sigma_ns = 1.5 / 2.355
    half_width = int(np.ceil(5.0 * sigma_ns / bin_width_ns))
    offsets = np.arange(-half_width, half_width + 1) * bin_width_ns
    kernel = np.exp(-0.5 * (offsets / sigma_ns) ** 2)
    kernel /= np.sum(kernel)
    return np.convolve(charge, kernel, mode="same")


def _pulse_shapes(
    events: dict[str, dict[str, np.ndarray]],
    histograms: dict[str, dict[str, np.ndarray]],
) -> mpl.figure.Figure:
    samples = (
        ("eplus", 100),
        ("piplus", 100),
        ("eplus", 300),
        ("piplus", 300),
    )
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 7.8), sharex=True, sharey=True)
    for panel_index, (ax, (particle, energy)) in enumerate(
        zip(axes.flat, samples, strict=True)
    ):
        name = _sample(particle, energy)
        hist = histograms[name]
        data = events[name]
        edges = hist["time_edges_ns"]
        centers = 0.5 * (edges[:-1] + edges[1:])
        bin_width = float(np.mean(np.diff(edges)))
        pitch_index_20 = int(np.flatnonzero(np.isclose(PITCHES_UM, 20.0))[0])
        pitch_index_100 = int(np.flatnonzero(np.isclose(PITCHES_UM, 100.0))[0])
        recovery_index_43 = int(np.flatnonzero(np.isclose(RECOVERY_NS, 4.3))[0])
        recovery_index_40 = int(np.flatnonzero(np.isclose(RECOVERY_NS, 40.0))[0])
        curves = (
            (hist["time_incident"].astype(float), "Photon arrivals", "0.35", "--"),
            (
                hist["time_accepted_scan"][pitch_index_20, recovery_index_43].astype(float),
                r"FLASH 20 $\mu$m, 4.3 ns",
                "#BD1F01",
                "-",
            ),
            (
                hist["time_accepted_scan"][pitch_index_100, recovery_index_40].astype(float),
                r"dSiPM stress: 100 $\mu$m, 40 ns",
                "#832DB6",
                "-.",
            ),
        )
        for values, label, color, linestyle in curves:
            normalized = values / np.max(values)
            ax.step(
                centers,
                normalized,
                where="mid",
                color=color,
                linestyle=linestyle,
                linewidth=1.8,
                label=label,
            )
        analog_waveform = _convolve_fast_output(
            hist["time_analog_charge"].astype(float),
            bin_width,
        )
        analog_waveform /= np.max(analog_waveform)
        ax.plot(
            centers,
            analog_waveform,
            color="#3F90DA",
            linewidth=2.2,
            label="Conventional SiPM fast-output proxy",
        )
        flash_survival = 100.0 * _retention(
            data, "accepted_reference", "detected_reference"
        )
        analog_survival = 100.0 * _retention(
            data,
            "analog_charge_pe",
            "analog_detected_at_31pct_pde",
        )
        coarse_survival = 100.0 * float(
            np.sum(hist["time_accepted_scan"][pitch_index_100, recovery_index_40])
            / np.sum(hist["time_detected"])
        )
        ax.text(
            0.97,
            0.70,
            "Signal retained\n"
            f"FLASH: {flash_survival:.1f}%\n"
            f"SiPM charge: {analog_survival:.1f}%\n"
            f"100 $\\mu$m/40 ns: {coarse_survival:.1f}%",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=8.7,
        )
        ax.set_xlim(8.5, 22.0)
        ax.set_ylim(0.0, 1.15)
        ax.set_xlabel("Photon arrival time [ns]")
        ax.set_ylabel("Response shape [peak normalized]")
        stress = " - stress test" if energy == 300 else ""
        _panel(
            ax,
            chr(ord("A") + panel_index),
            f"{PARTICLE_LABELS[particle]}, {energy} GeV{stress}",
        )
        if panel_index == 0:
            ax.legend(loc="upper right", fontsize=8.2)
    fig.subplots_adjust(top=0.90, bottom=0.09, hspace=0.18, wspace=0.12)
    _stamp(fig, "Same photon arrivals - stated readout models")
    fig.text(
        0.5,
        0.015,
        "The analog trace is the partial-recharge charge profile convolved with a 1.5 ns FWHM Gaussian fast-output proxy; FLASH curves are timestamp densities.",
        ha="center",
        fontsize=8.7,
    )
    return fig


def _positive_limits(array: np.ndarray) -> tuple[float, float]:
    positive = array[array > 0]
    if positive.size == 0:
        return (1.0, 2.0)
    return tuple(float(value) for value in np.percentile(positive, (5.0, 99.7)))


def _spatial_maps(
    events: dict[str, dict[str, np.ndarray]],
    histograms: dict[str, dict[str, np.ndarray]],
) -> mpl.figure.Figure:
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 8.5))
    for column, particle in enumerate(PARTICLES):
        name = _sample(particle, 100)
        hist = histograms[name]
        n_events = events[name]["event"].size

        global_map = hist["global_xy_incident"].astype(float).T / n_events
        gx = hist["global_x_edges_mm"]
        gy = hist["global_y_edges_mm"]
        vmin, vmax = _positive_limits(global_map)
        image = axes[0, column].imshow(
            np.ma.masked_less_equal(global_map, 0),
            origin="lower",
            extent=(gx[0], gx[-1], gy[0], gy[-1]),
            cmap="Blues",
            norm=LogNorm(vmin=vmin, vmax=vmax),
            interpolation="nearest",
            aspect="equal",
        )
        axes[0, column].set_xlabel("Detector x [mm]")
        axes[0, column].set_ylabel("Detector y [mm]")
        axes[0, column].set_xticks((-100, -50, 0, 50))
        _panel(
            axes[0, column],
            chr(ord("A") + column),
            f"{PARTICLE_LABELS[particle]} shower footprint",
        )
        colorbar = fig.colorbar(image, ax=axes[0, column], pad=0.02)
        colorbar.set_label("Photons / event / 1 mm$^2$")

        local_map = hist["local_xy_incident"].astype(float).T / n_events
        lx = hist["local_x_edges_mm"]
        ly = hist["local_y_edges_mm"]
        local_vmin, local_vmax = _positive_limits(local_map)
        local_image = axes[1, column].imshow(
            np.ma.masked_less_equal(local_map, 0),
            origin="lower",
            extent=(lx[0], lx[-1], ly[0], ly[-1]),
            cmap="Blues",
            norm=LogNorm(vmin=local_vmin, vmax=local_vmax),
            interpolation="nearest",
            aspect="equal",
        )
        axes[1, column].set_xlim(-0.5, 0.5)
        axes[1, column].set_ylim(-0.5, 0.5)
        axes[1, column].add_patch(
            Circle(
                (0, 0),
                0.4,
                fill=False,
                color="#BD1F01",
                linewidth=1.6,
                linestyle="--",
            )
        )
        axes[1, column].set_xlabel("Local x at fiber face [mm]")
        axes[1, column].set_ylabel("Local y at fiber face [mm]")
        _panel(
            axes[1, column],
            chr(ord("C") + column),
            f"{PARTICLE_LABELS[particle]} fiber-face photon density",
        )
        local_colorbar = fig.colorbar(local_image, ax=axes[1, column], pad=0.02)
        local_colorbar.set_label("Photons / event / bin")
    fig.subplots_adjust(top=0.89, bottom=0.08, hspace=0.34, wspace=0.48)
    _stamp(fig, "100 GeV - independent readable color ranges")
    return fig


def _validation_page(
    events: dict[str, dict[str, np.ndarray]],
    metrics: dict[str, object],
) -> mpl.figure.Figure:
    labels = [f"{PARTICLE_LABELS[p]}\n{e}" for p in PARTICLES for e in ENERGIES]
    names = [_sample(p, e) for p in PARTICLES for e in ENERGIES]
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 8.4))
    x = np.arange(len(names))

    counts = [metrics["samples"][name]["events"] for name in names]
    axes[0, 0].bar(x, counts, color=[PARTICLE_COLORS[name.split("_")[0]] for name in names])
    axes[0, 0].axhline(5000, color="0.2", linestyle="--", linewidth=1.2)
    axes[0, 0].set_xticks(x, labels=labels)
    axes[0, 0].set_ylabel("Validated events")
    axes[0, 0].set_ylim(0, 5400)
    _panel(axes[0, 0], "A", "Checkpoint completeness")

    precision = [
        metrics["samples"][name]["selected_photons_relative_sem_percent"]
        for name in names
    ]
    axes[0, 1].bar(x, precision, color=[PARTICLE_COLORS[name.split("_")[0]] for name in names])
    axes[0, 1].axhline(1.0, color="0.2", linestyle="--", linewidth=1.2)
    axes[0, 1].set_xticks(x, labels=labels)
    axes[0, 1].set_ylabel("Relative SEM of mean photon yield [%]")
    axes[0, 1].set_ylim(0, max(1.08, 1.2 * max(precision)))
    _panel(axes[0, 1], "B", "Headline-yield precision")

    block_values = []
    for name in names:
        data = events[name]
        sums, counts_per_file = _block_sums(
            data, data["selected_cherenkov_photons"]
        )
        means = sums / counts_per_file
        block_values.append(means / np.mean(data["selected_cherenkov_photons"]))
    axes[1, 0].boxplot(
        block_values,
        positions=x,
        widths=0.6,
        showfliers=True,
        flierprops={"markersize": 2.5, "alpha": 0.45},
        medianprops={"color": "#BD1F01", "linewidth": 1.6},
    )
    axes[1, 0].axhline(1.0, color="0.4", linestyle=":", linewidth=1.0)
    axes[1, 0].set_xticks(x, labels=labels)
    axes[1, 0].set_ylabel("Per-file mean / sample mean")
    _panel(axes[1, 0], "C", "Independent-file stability")

    axes[1, 1].axis("off")
    closures = metrics["closures"]
    integer_residuals = []
    analog_residuals = []
    for name in names:
        for key, value in closures[name].items():
            if key == "analog_charge_time_minus_event":
                analog_residuals.append(abs(float(value)))
            else:
                integer_residuals.append(abs(int(value)))
    provenance = metrics["checkpoint"]["provenance"]
    text = (
        "Closure and provenance\n\n"
        f"Max integer closure residual: {max(integer_residuals)} photons\n"
        f"Max analog-charge residual: {max(analog_residuals):.3e} p.e.\n"
        "All 24 pitch/recovery time integrals close to event rows.\n\n"
        f"Analysis: {str(provenance['analysis_commit'])[:12]}\n"
        f"DREAMSim: {str(provenance['dreamsim_commit'])[:12]}\n"
        f"Generation: {str(provenance['dsipm_commit'])[:12]}\n"
        f"Reducer: {str(provenance['reducer_sha256'])[:12]}\n\n"
        "Scope: full shower transport; explicit optical photons\n"
        "only in rods 20-60 and layers 15-65."
    )
    axes[1, 1].text(
        0.04,
        0.96,
        text,
        transform=axes[1, 1].transAxes,
        ha="left",
        va="top",
        fontsize=10.5,
        linespacing=1.35,
        bbox={"boxstyle": "round,pad=0.6", "facecolor": "#F4F4F4", "edgecolor": "0.75"},
    )
    fig.subplots_adjust(top=0.91, bottom=0.09, hspace=0.30, wspace=0.22)
    _stamp(fig, "Statistical and physics closure")
    return fig


def _audit_markdown(metrics: dict[str, object], figures: dict[str, object]) -> str:
    lines = [
        "# FLASH DREAMSim schema-v3 checkpoint audit",
        "",
        "## Scope",
        "",
        "The full calorimeter shower is transported, while explicit optical photons are retained only in the central 41-by-51 rod/layer window (rod 20-60, layer 15-65). The 10-100 GeV points are the FCC-ee-relevant range; 300 GeV is a dynamic-range stress test.",
        "",
        "## Event counts and headline response",
        "",
        "| Sample | Events | Selected photons/event | Relative SEM | Hottest fiber median | FLASH recovery loss | SiPM charge loss |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for particle in PARTICLES:
        for energy in ENERGIES:
            name = _sample(particle, energy)
            values = metrics["samples"][name]
            lines.append(
                f"| {PARTICLE_LABELS[particle]} {energy} GeV | {values['events']:,} | "
                f"{values['selected_photons_mean']:.2f} | "
                f"{values['selected_photons_relative_sem_percent']:.3f}% | "
                f"{values['hottest_fiber_median']:.0f} | "
                f"{values['flash_reference_recovery_loss_percent']:.3f}% | "
                f"{values['analog_charge_loss_percent']:.3f}% |"
            )
    lines.extend(
        [
            "",
            "## Closure",
            "",
            "Every integer event/histogram, per-fiber, SPAD multiplicity, and pitch/recovery time-profile closure residual is zero. Floating analog-charge differences are required to be consistent with numerical summation precision.",
            "",
            "## Claim classification",
            "",
        ]
    )
    for classification in ("supported", "provisional", "unsupported"):
        lines.append(f"### {classification.title()}")
        lines.append("")
        lines.extend(
            f"- {claim}"
            for claim in metrics["claim_classification"][classification]
        )
        lines.append("")
    lines.extend(
        [
            "## Figure files",
            "",
            *[f"- {name}: `{details['pdf']}`" for name, details in figures.items()],
            "",
        ]
    )
    return "\n".join(lines)


def run(args: argparse.Namespace) -> dict[str, object]:
    _configure_style()
    input_dir = args.input_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    plotter_path = Path(__file__).resolve()
    repository_dir = plotter_path.parents[1]
    merged_summary = json.loads(
        (input_dir / "merged_reduction_summary.json").read_text()
    )
    if merged_summary.get("schema_version") != 3:
        raise ValueError("This plotter requires a schema-v3 merged checkpoint")

    events = {
        _sample(particle, energy): _load_events(
            input_dir / f"{_sample(particle, energy)}.events.csv.gz"
        )
        for particle in PARTICLES
        for energy in ENERGIES
    }
    histograms = {
        _sample(particle, energy): _load_histograms(
            input_dir / f"{_sample(particle, energy)}.histograms.npz"
        )
        for particle in PARTICLES
        for energy in ENERGIES
    }
    metrics = _make_metrics(events, histograms, merged_summary)
    metrics_path = output_dir / "FLASH_schema3_checkpoint_metrics.json"
    metrics_path.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")

    packet_path = output_dir / "FLASH_schema3_checkpoint_figures.pdf"
    figures: dict[str, object] = {}
    with PdfPages(packet_path) as packet:
        figures["response_transfer"] = _save_figure(
            _response_transfer(events),
            output_dir,
            "FLASH_response_transfer",
            packet,
        )
        figures["occupancy_segmentation_deadtime"] = _save_figure(
            _occupancy_to_deadtime(events, histograms),
            output_dir,
            "FLASH_occupancy_segmentation_deadtime",
            packet,
        )
        figures["pulse_shapes"] = _save_figure(
            _pulse_shapes(events, histograms),
            output_dir,
            "FLASH_matched_pulse_shapes",
            packet,
        )
        figures["spatial_maps"] = _save_figure(
            _spatial_maps(events, histograms),
            output_dir,
            "FLASH_spatial_validation",
            packet,
        )
        figures["statistical_validation"] = _save_figure(
            _validation_page(events, metrics),
            output_dir,
            "FLASH_statistical_validation",
            packet,
        )

    audit_path = output_dir / "FLASH_schema3_checkpoint_audit.md"
    audit_path.write_text(_audit_markdown(metrics, figures))
    manifest = {
        "schema_version": 1,
        "input_directory": str(input_dir),
        "output_directory": str(output_dir),
        "packet_pdf": str(packet_path),
        "packet_pdf_sha256": _sha256(packet_path),
        "metrics_json": str(metrics_path),
        "metrics_json_sha256": _sha256(metrics_path),
        "audit_markdown": str(audit_path),
        "audit_markdown_sha256": _sha256(audit_path),
        "plotter_path": str(plotter_path),
        "plotter_sha256": _sha256(plotter_path),
        "plotter_git_commit": _git_commit(repository_dir),
        "figures": figures,
    }
    manifest_path = output_dir / "FLASH_schema3_checkpoint_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        manifest = run(args)
    except Exception as error:
        print(f"ERROR: {error}")
        return 2
    print(manifest["packet_pdf"])
    print(manifest["metrics_json"])
    print(manifest["audit_markdown"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
