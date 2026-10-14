#!/usr/bin/env python3
"""Reduce one validated DREAMSim ROOT file to mergeable plot inputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

REPO_DIR = Path(__file__).resolve().parents[1]
if str(REPO_DIR) not in sys.path:
    sys.path.insert(0, str(REPO_DIR))

from dsipm.full_detector import iter_full_detector_events_uproot
from dsipm.reduction import ANALOG_SIPM, HistogramAccumulator, reduce_event
from dsipm.workflow import DEFAULT_WORKFLOW_PATH, load_sensor_workflow


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_sidecar(
    root_file: Path,
    override: Path | None = None,
) -> tuple[Path, dict[str, object]]:
    sidecar = (
        override.expanduser().resolve()
        if override is not None
        else Path(str(root_file) + ".validation.json")
    )
    if not sidecar.is_file():
        raise FileNotFoundError(f"Missing validation sidecar: {sidecar}")
    with sidecar.open() as handle:
        payload = json.load(handle)
    if payload.get("status") != "valid":
        raise ValueError(f"Validation status is not valid in {sidecar}")
    if int(payload.get("entries", -1)) <= 0:
        raise ValueError(f"Invalid entry count in {sidecar}")
    if int(payload.get("size_bytes", -1)) != root_file.stat().st_size:
        raise ValueError(f"ROOT size no longer matches {sidecar}")
    return sidecar, payload


def _require_value(
    payload: dict[str, object],
    key: str,
    expected: str | None,
) -> None:
    if expected is None:
        return
    observed = payload.get(key)
    if observed != expected:
        raise ValueError(f"Expected {key}={expected!r}, observed {observed!r}")


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("Cannot write an empty event table")
    temporary = path.with_suffix(path.suffix + ".tmp")
    with gzip.open(temporary, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def _write_npz(path: Path, payload: dict[str, np.ndarray]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as handle:
        np.savez_compressed(handle, **payload)
    temporary.replace(path)


def reduce_file(args: argparse.Namespace) -> dict[str, object]:
    started = time.time()
    root_file = args.root_file.expanduser().resolve()
    if not root_file.is_file():
        raise FileNotFoundError(f"Missing ROOT file: {root_file}")
    sidecar_path, sidecar = _load_sidecar(root_file, args.validation_sidecar)
    _require_value(sidecar, "dreamsim_commit", args.expected_dreamsim_commit)
    _require_value(sidecar, "dsipm_commit", args.expected_dsipm_commit)
    _require_value(sidecar, "campaign_sha256", args.expected_campaign_sha256)

    workflow = load_sensor_workflow(args.workflow)
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = root_file.stem
    event_path = output_dir / f"{stem}.events.csv.gz"
    histogram_path = output_dir / f"{stem}.histograms.npz"
    summary_path = output_dir / f"{stem}.reduction.json"
    existing = [path for path in (event_path, histogram_path, summary_path) if path.exists()]
    if existing:
        raise FileExistsError(
            "Reduction output already exists: " + ", ".join(str(path) for path in existing)
        )

    sample = args.sample or root_file.parent.name
    pitches = workflow.scans.full_detector_spad_pitch_um
    recovery = workflow.scans.full_detector_recovery_ns
    rows: list[dict[str, object]] = []
    histograms = HistogramAccumulator(
        pitches_um=tuple(float(value) for value in pitches)
    )
    selected_total = 0
    detected_total = 0
    accepted_total = 0
    for source_entry, event in enumerate(
        iter_full_detector_events_uproot(
            root_file,
            gate_ns=(args.gate_start_ns, args.gate_end_ns),
            event_index_offset=args.event_index_offset,
            step_size=args.step_size,
        )
    ):
        row, detected_mask, accepted_mask, lost_mask = reduce_event(
            event,
            sample=sample,
            source_root=str(root_file),
            source_entry=source_entry,
            pde=workflow.reference.effective_pde,
            pitches_um=pitches,
            recovery_times_ns=recovery,
            reference_pitch_um=workflow.reference.spad_pitch_um,
            reference_recovery_ns=workflow.reference.spad_recovery_ns,
            seed=workflow.reference.seed,
        )
        rows.append(row)
        histograms.add(event, detected_mask, accepted_mask, lost_mask)
        selected_total += event.selected_photon_count
        detected_total += int(np.count_nonzero(detected_mask))
        accepted_total += int(np.count_nonzero(accepted_mask))

    expected_entries = int(sidecar["entries"])
    if len(rows) != expected_entries:
        raise ValueError(
            f"Reduced {len(rows)} events but validation sidecar records {expected_entries}"
        )
    _write_csv(event_path, rows)
    _write_npz(histogram_path, histograms.payload())
    closure = histograms.closure()
    if closure["detected_photons_in_time_histogram"] != detected_total:
        raise AssertionError("Event and histogram detected totals do not close")
    if closure["accepted_photons_in_time_histogram"] != accepted_total:
        raise AssertionError("Event and histogram accepted totals do not close")

    summary: dict[str, object] = {
        "schema_version": 2,
        "status": "valid",
        "created_unix": time.time(),
        "runtime_seconds": time.time() - started,
        "sample": sample,
        "root_file": str(root_file),
        "root_size_bytes": root_file.stat().st_size,
        "validation_sidecar": str(sidecar_path),
        "validation_sidecar_sha256": _sha256(sidecar_path),
        "events": len(rows),
        "event_index_offset": args.event_index_offset,
        "event_index_last": rows[-1]["event"],
        "beam_particle": sidecar.get("observed_particle"),
        "beam_energy_gev": sidecar.get("observed_energy_gev"),
        "analysis_commit": args.analysis_commit,
        "reducer_sha256": _sha256(Path(__file__).resolve()),
        "dreamsim_commit": sidecar.get("dreamsim_commit"),
        "dsipm_commit": sidecar.get("dsipm_commit"),
        "campaign_sha256": sidecar.get("campaign_sha256"),
        "launcher_sha256": sidecar.get("launcher_sha256"),
        "validator_sha256": sidecar.get("validator_sha256"),
        "sim_executable_sha256": sidecar.get("sim_executable_sha256"),
        "macro_sha256": sidecar.get("macro_sha256"),
        "seeds": [sidecar.get("seed_1"), sidecar.get("seed_2")],
        "sensor_workflow": {
            "path": str(workflow.source_path),
            "sha256": workflow.sha256,
            "effective_pde": workflow.reference.effective_pde,
            "reference_pitch_um": workflow.reference.spad_pitch_um,
            "reference_recovery_ns": workflow.reference.spad_recovery_ns,
            "seed": workflow.reference.seed,
            "scan_pitches_um": list(pitches),
            "scan_recovery_ns": list(recovery),
        },
        "analog_sipm": {
            "name": ANALOG_SIPM.name,
            "pde": ANALOG_SIPM.pde,
            "microcell_size_um": ANALOG_SIPM.microcell_size_um,
            "n_microcells": ANALOG_SIPM.n_microcells,
            "recovery_tau_ns": ANALOG_SIPM.recovery_tau_ns,
            "source_url": ANALOG_SIPM.source_url,
        },
        "gate_ns": [args.gate_start_ns, args.gate_end_ns],
        "totals": {
            "selected_cherenkov_photons": selected_total,
            "detected_reference": detected_total,
            "accepted_reference": accepted_total,
            **closure,
        },
        "outputs": {
            "events_csv_gz": str(event_path),
            "events_csv_gz_sha256": _sha256(event_path),
            "histograms_npz": str(histogram_path),
            "histograms_npz_sha256": _sha256(histogram_path),
        },
        "limitations": [
            "Only explicit Cherenkov optical photons are digitized.",
            "sum3dSS is retained as a stored DREAMSim signal; scintillation photons are absent.",
            "The reference digital response uses constant effective PDE.",
            "Dark counts, crosstalk, afterpulsing, and electronics arbitration are not included.",
            "Per-fiber and per-SPAD occupancy products are additive summaries; raw photon vectors remain on HPCC.",
        ],
    }
    temporary_summary = summary_path.with_suffix(summary_path.suffix + ".tmp")
    temporary_summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    temporary_summary.replace(summary_path)
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root_file", type=Path)
    parser.add_argument("--validation-sidecar", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--sample")
    parser.add_argument("--event-index-offset", type=int, default=0)
    parser.add_argument("--workflow", type=Path, default=DEFAULT_WORKFLOW_PATH)
    parser.add_argument("--gate-start-ns", type=float, default=0.0)
    parser.add_argument("--gate-end-ns", type=float, default=40.0)
    parser.add_argument("--step-size", default="64 MB")
    parser.add_argument("--expected-dreamsim-commit")
    parser.add_argument("--expected-dsipm-commit")
    parser.add_argument("--expected-campaign-sha256")
    parser.add_argument("--analysis-commit")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        summary = reduce_file(args)
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(
        f"Reduced {summary['events']} events and "
        f"{summary['totals']['selected_cherenkov_photons']} selected photons"
    )
    print(summary["outputs"]["events_csv_gz"])
    print(summary["outputs"]["histograms_npz"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
