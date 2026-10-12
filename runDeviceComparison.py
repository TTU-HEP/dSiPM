#!/usr/bin/env python3
"""Run historical pitch scans and real-device comparisons on DREAMSim truth."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

from dsipm.config import DeviceConfig, load_device_catalog
from dsipm.digitizer import digitize_event
from dsipm.outputs import PLOT_STYLES, ComparisonOutputs, print_summary
from dsipm.root_io import GEOMETRY_MODES, inspect_root_file, iter_root_events


REPO_DIR = Path(__file__).resolve().parent
DEFAULT_CATALOG = REPO_DIR / "device_configs.json"


def _parse_device_names(value: str | None) -> list[str] | None:
    if not value:
        return None
    names = [part.strip() for part in value.split(",") if part.strip()]
    return names or None


def _prepare_output(path: str | Path) -> Path:
    output = Path(path).expanduser().resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(
            f"Output directory is not empty: {output}\n"
            "Choose a new directory so an earlier comparison is not overwritten."
        )
    output.mkdir(parents=True, exist_ok=True)
    return output


def _catalog_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_comparison(
    *,
    input_path: str | Path,
    output_dir: str | Path,
    profile: str,
    device_names: list[str] | None,
    catalog_path: str | Path,
    max_events: int | None,
    seed: int,
    geometry_mode: str,
    detector_min_mm: float,
    detector_max_mm: float,
    gate_ns: tuple[float, float],
    fiber_index: int | None = None,
    plot_style: str = "comparison",
) -> list[dict[str, object]]:
    source = Path(input_path).expanduser().resolve()
    if not source.exists():
        raise FileNotFoundError(f"Input ROOT file does not exist: {source}")
    if max_events is not None and max_events <= 0:
        raise ValueError("max_events must be positive")
    if detector_max_mm <= detector_min_mm:
        raise ValueError("detector_max_mm must be greater than detector_min_mm")
    if gate_ns[1] <= gate_ns[0]:
        raise ValueError("The time gate upper edge must be greater than its lower edge")

    output = _prepare_output(output_dir)
    catalog = load_device_catalog(catalog_path)
    devices = catalog.select(profile=profile, device_names=device_names)
    source_info = inspect_root_file(source)
    if source_info["missing_required_branches"]:
        raise ValueError(
            "Input is missing required branches: "
            + ", ".join(source_info["missing_required_branches"])
        )

    manifest = {
        "created_unix": time.time(),
        "input": source_info,
        "output_directory": str(output),
        "profile": profile if device_names is None else None,
        "devices": [device.to_dict() for device in devices],
        "device_catalog": str(catalog.source_path),
        "device_catalog_sha256": _catalog_sha256(catalog.source_path),
        "max_events": max_events,
        "seed": seed,
        "geometry_mode": geometry_mode,
        "detector_range_mm": [detector_min_mm, detector_max_mm],
        "gate_ns": list(gate_ns),
        "fiber_index": fiber_index,
        "plot_style": plot_style,
        "implemented_effects": [
            "wavelength-dependent or constant PDE",
            "common seeded random numbers across devices",
            "non-paralyzable finite SPAD recovery",
            "Gaussian timing response when configured",
            "TDC quantization",
            "timing-group and timestamp-policy reduction",
            "transparent event-payload proxy",
        ],
        "not_yet_implemented": [
            "dark counts",
            "optical crosstalk",
            "afterpulsing",
            "hot-cell masking",
            "wavelength-dependent non-Gaussian timing tails",
            "angular PDE and detailed fiber-to-sensor optics",
            "link arbitration and power",
        ],
    }

    print(f"Input:    {source}")
    print(f"Entries:  {source_info['entries']}")
    print(f"Geometry: {geometry_mode}")
    print(f"Fiber:    {'all' if fiber_index is None else fiber_index}")
    print(f"Seed:     {seed}")
    print("Devices:  " + ", ".join(device.name for device in devices))
    print(f"Output:   {output}")

    output_writer = ComparisonOutputs(output, devices, plot_style=plot_style)
    processed_events = 0
    try:
        for event in iter_root_events(
            source,
            max_events=max_events,
            geometry_mode=geometry_mode,
            gate_ns=gate_ns,
            detector_min_mm=detector_min_mm,
            detector_max_mm=detector_max_mm,
            fiber_index=fiber_index,
        ):
            for device in devices:
                result = digitize_event(
                    event,
                    device,
                    seed=seed,
                    detector_min_mm=detector_min_mm,
                    detector_max_mm=detector_max_mm,
                    gate_ns=gate_ns,
                )
                output_writer.add(result)
            processed_events += 1
            if processed_events <= 3 or processed_events % 100 == 0:
                print(
                    f"Processed event {event.event_index}: "
                    f"{event.selected_photon_count} selected photons"
                )
        if processed_events == 0:
            raise RuntimeError("No events were processed")
        rows = output_writer.close()
    except Exception:
        output_writer._event_handle.close()
        output_writer.root_file.Close()
        raise

    manifest["processed_events"] = processed_events
    manifest["completed_unix"] = time.time()
    manifest["runtime_seconds"] = manifest["completed_unix"] - manifest["created_unix"]
    with (output / "run_manifest.json").open("w") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print_summary(rows)
    print(f"\nWrote comparison products to {output}")
    return rows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Apply historical and real-device response models to one DREAMSim ROOT file."
        )
    )
    parser.add_argument("input", help="DREAMSim ROOT file containing a 'tree'")
    parser.add_argument(
        "--output",
        help="New or empty output directory (default: results/<input>_<profile>)",
    )
    parser.add_argument(
        "--profile",
        default="default",
        help="Named profile from device_configs.json (default: %(default)s)",
    )
    parser.add_argument(
        "--devices",
        help="Comma-separated device configuration names; overrides --profile",
    )
    parser.add_argument(
        "--catalog",
        default=str(DEFAULT_CATALOG),
        help="Device catalog JSON (default: %(default)s)",
    )
    parser.add_argument("--max-events", type=int, help="Process only the first N events")
    parser.add_argument("--seed", type=int, default=4307, help="Recorded RNG seed")
    parser.add_argument(
        "--geometry",
        choices=sorted(GEOMETRY_MODES),
        default="student_2026",
        help="Coordinate mapping (default: %(default)s)",
    )
    parser.add_argument("--detector-min-mm", type=float, default=-100.0)
    parser.add_argument("--detector-max-mm", type=float, default=100.0)
    parser.add_argument("--gate-start-ns", type=float, default=0.0)
    parser.add_argument("--gate-end-ns", type=float, default=40.0)
    parser.add_argument(
        "--fiber",
        type=int,
        choices=range(4),
        help="Restrict the study to one production-fiber index (0--3)",
    )
    parser.add_argument(
        "--plot-style",
        choices=sorted(PLOT_STYLES),
        default="comparison",
        help=(
            "comparison for device-summary plots, slides for the 2025 study "
            "binning, or both (default: %(default)s)"
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    input_path = Path(args.input).expanduser().resolve()
    output = args.output
    if output is None:
        output = REPO_DIR / "results" / f"{input_path.stem}_{args.profile}"

    try:
        run_comparison(
            input_path=input_path,
            output_dir=output,
            profile=args.profile,
            device_names=_parse_device_names(args.devices),
            catalog_path=args.catalog,
            max_events=args.max_events,
            seed=args.seed,
            geometry_mode=args.geometry,
            detector_min_mm=args.detector_min_mm,
            detector_max_mm=args.detector_max_mm,
            gate_ns=(args.gate_start_ns, args.gate_end_ns),
            fiber_index=args.fiber,
            plot_style=args.plot_style,
        )
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
