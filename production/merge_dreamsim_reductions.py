#!/usr/bin/env python3
"""Merge validated per-file DREAMSim reductions into plot-ready samples."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np


EDGE_KEYS = {
    "time_edges_ns",
    "local_x_edges_mm",
    "local_y_edges_mm",
    "global_x_edges_mm",
    "global_y_edges_mm",
    "depth_edges_mm",
    "depth_time_edges_ns",
    "fiber_occupancy_edges",
    "cell_occupancy_edges",
    "cell_occupancy_pitches_um",
    "scan_recovery_ns",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _expected_manifest(path: Path) -> dict[tuple[str, int, int], dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    expected: dict[tuple[str, int, int], dict[str, str]] = {}
    for row in rows:
        key = (row["sample"], int(row["seed_1"]), int(row["seed_2"]))
        if key in expected:
            raise ValueError(f"Duplicate sample/seed key in manifest: {key}")
        expected[key] = row
    if not expected:
        raise ValueError(f"Empty production manifest: {path}")
    return expected


def _paths(summary_path: Path) -> tuple[Path, Path]:
    suffix = ".reduction.json"
    if not summary_path.name.endswith(suffix):
        raise ValueError(f"Unexpected reduction-summary name: {summary_path}")
    stem = summary_path.name[: -len(suffix)]
    return (
        summary_path.with_name(f"{stem}.events.csv.gz"),
        summary_path.with_name(f"{stem}.histograms.npz"),
    )


def _load_reduction(summary_path: Path) -> dict[str, object]:
    with summary_path.open() as handle:
        summary = json.load(handle)
    if summary.get("schema_version") not in {1, 2, 3} or summary.get("status") != "valid":
        raise ValueError(f"Invalid reduction summary: {summary_path}")
    event_path, histogram_path = _paths(summary_path)
    if not event_path.is_file() or not histogram_path.is_file():
        raise FileNotFoundError(f"Reduction products are incomplete for {summary_path}")
    outputs = summary.get("outputs", {})
    if outputs.get("events_csv_gz_sha256") != _sha256(event_path):
        raise ValueError(f"Event-table checksum mismatch: {event_path}")
    if outputs.get("histograms_npz_sha256") != _sha256(histogram_path):
        raise ValueError(f"Histogram checksum mismatch: {histogram_path}")
    summary["_summary_path"] = summary_path
    summary["_event_path"] = event_path
    summary["_histogram_path"] = histogram_path
    return summary


def _merge_histograms(paths: list[Path], destination: Path) -> None:
    merged: dict[str, np.ndarray] | None = None
    for path in paths:
        with np.load(path, allow_pickle=False) as payload:
            current = {key: np.asarray(payload[key]) for key in payload.files}
        if merged is None:
            merged = {key: values.copy() for key, values in current.items()}
            continue
        if set(current) != set(merged):
            raise ValueError(f"Histogram keys differ in {path}")
        for key, values in current.items():
            if key in EDGE_KEYS:
                if not np.array_equal(merged[key], values):
                    raise ValueError(f"Histogram edges differ for {key} in {path}")
            else:
                if merged[key].shape != values.shape:
                    raise ValueError(f"Histogram shape differs for {key} in {path}")
                merged[key] += values
    if merged is None:
        raise ValueError("No histogram files to merge")
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    with temporary.open("wb") as handle:
        np.savez_compressed(handle, **merged)
    temporary.replace(destination)


def _merge_event_tables(
    paths: list[Path],
    destination: Path,
    *,
    sample: str,
) -> tuple[int, set[int]]:
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    fieldnames: list[str] | None = None
    event_ids: set[int] = set()
    count = 0
    with gzip.open(temporary, "wt", newline="") as output_handle:
        writer: csv.DictWriter | None = None
        for path in paths:
            with gzip.open(path, "rt", newline="") as input_handle:
                reader = csv.DictReader(input_handle)
                if reader.fieldnames is None:
                    raise ValueError(f"Missing CSV header: {path}")
                if fieldnames is None:
                    fieldnames = list(reader.fieldnames)
                    writer = csv.DictWriter(output_handle, fieldnames=fieldnames)
                    writer.writeheader()
                elif list(reader.fieldnames) != fieldnames:
                    raise ValueError(f"Event-table columns differ in {path}")
                assert writer is not None
                for row in reader:
                    if row["sample"] != sample:
                        raise ValueError(f"Expected sample {sample}, found {row['sample']}")
                    event_id = int(row["event"])
                    if event_id in event_ids:
                        raise ValueError(f"Duplicate event {sample}:{event_id}")
                    event_ids.add(event_id)
                    writer.writerow(row)
                    count += 1
    temporary.replace(destination)
    return count, event_ids


def merge(args: argparse.Namespace) -> dict[str, object]:
    started = time.time()
    input_dir = args.input_dir.expanduser().resolve()
    manifest_path = args.manifest.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    expected = _expected_manifest(manifest_path)
    summary_paths = sorted(input_dir.rglob("*.reduction.json"))
    if not summary_paths:
        raise FileNotFoundError(f"No reduction summaries under {input_dir}")
    reductions = [_load_reduction(path) for path in summary_paths]

    observed: dict[tuple[str, int, int], dict[str, object]] = {}
    provenance_values: dict[str, set[object]] = defaultdict(set)
    by_sample: dict[str, list[dict[str, object]]] = defaultdict(list)
    for reduction in reductions:
        seeds = reduction.get("seeds", [])
        if len(seeds) != 2 or seeds[0] is None or seeds[1] is None:
            raise ValueError(f"Missing production seeds in {reduction['_summary_path']}")
        key = (str(reduction["sample"]), int(seeds[0]), int(seeds[1]))
        if key not in expected:
            raise ValueError(f"Reduction does not appear in production manifest: {key}")
        if key in observed:
            raise ValueError(f"Duplicate reduction for manifest key: {key}")
        observed[key] = reduction
        by_sample[key[0]].append(reduction)
        for name in (
            "analysis_commit",
            "reducer_sha256",
            "dreamsim_commit",
            "dsipm_commit",
            "campaign_sha256",
            "launcher_sha256",
            "validator_sha256",
            "sim_executable_sha256",
        ):
            provenance_values[name].add(reduction.get(name))
    inconsistent = {
        name: sorted(str(value) for value in values)
        for name, values in provenance_values.items()
        if len(values) != 1
    }
    if inconsistent:
        raise ValueError(f"Mixed production provenance: {inconsistent}")

    missing = sorted(set(expected) - set(observed))
    if args.require_complete and missing:
        raise ValueError(f"Missing {len(missing)} reductions from complete merge")

    sample_outputs: dict[str, object] = {}
    index_rows: list[dict[str, object]] = []
    for sample, sample_reductions in sorted(by_sample.items()):
        sample_reductions.sort(key=lambda item: int(item["event_index_offset"]))
        event_destination = output_dir / f"{sample}.events.csv.gz"
        histogram_destination = output_dir / f"{sample}.histograms.npz"
        event_count, event_ids = _merge_event_tables(
            [Path(item["_event_path"]) for item in sample_reductions],
            event_destination,
            sample=sample,
        )
        _merge_histograms(
            [Path(item["_histogram_path"]) for item in sample_reductions],
            histogram_destination,
        )
        expected_events = sum(
            int(expected[key]["events"])
            for key in expected
            if key[0] == sample and key in observed
        )
        if event_count != expected_events:
            raise ValueError(
                f"Sample {sample} has {event_count} rows; expected {expected_events}"
            )
        sample_outputs[sample] = {
            "files_reduced": len(sample_reductions),
            "events": event_count,
            "event_index_min": min(event_ids),
            "event_index_max": max(event_ids),
            "events_csv_gz": str(event_destination),
            "events_csv_gz_sha256": _sha256(event_destination),
            "histograms_npz": str(histogram_destination),
            "histograms_npz_sha256": _sha256(histogram_destination),
        }
        for item in sample_reductions:
            index_rows.append(
                {
                    "sample": sample,
                    "events": item["events"],
                    "event_index_offset": item["event_index_offset"],
                    "root_file": item["root_file"],
                    "reduction_summary": item["_summary_path"],
                    "summary_sha256": _sha256(Path(item["_summary_path"])),
                }
            )

    index_path = output_dir / "reduction_index.csv"
    with index_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(index_rows[0]))
        writer.writeheader()
        writer.writerows(index_rows)
    result: dict[str, object] = {
        "schema_version": 3,
        "status": "complete" if not missing else "partial",
        "created_unix": time.time(),
        "runtime_seconds": time.time() - started,
        "input_directory": str(input_dir),
        "production_manifest": str(manifest_path),
        "production_manifest_sha256": _sha256(manifest_path),
        "expected_files": len(expected),
        "reduced_files": len(observed),
        "missing_files": len(missing),
        "missing_manifest_keys": [
            {"sample": sample, "seed_1": seed_1, "seed_2": seed_2}
            for sample, seed_1, seed_2 in missing
        ],
        "provenance": {
            name: next(iter(values)) for name, values in provenance_values.items()
        },
        "samples": sample_outputs,
        "reduction_index": str(index_path),
        "reduction_index_sha256": _sha256(index_path),
    }
    summary_path = output_dir / "merged_reduction_summary.json"
    summary_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--require-complete", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = merge(args)
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(
        f"Merged {result['reduced_files']} / {result['expected_files']} files "
        f"({result['status']})"
    )
    for sample, details in result["samples"].items():
        print(f"{sample}: {details['events']} events")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
