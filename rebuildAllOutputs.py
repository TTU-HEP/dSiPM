#!/usr/bin/env python3
"""Rebuild the canonical CPAD 2026 and JINST output set.

Run this script with the ROOT-compatible Python after sourcing thisroot.sh:

    source /opt/homebrew/bin/thisroot.sh
    /opt/homebrew/opt/python@3.13/bin/python3.13 rebuildAllOutputs.py

The default path starts from the local ROOT samples rather than cached event
tables.  It reproduces the 2025 study, rebuilds the corrected one-fiber and
full-detector analyses, stages the reviewed figures, compiles the CPAD deck
and JINST draft, and records input/output checksums.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


REPO = Path(__file__).resolve().parent
RESULTS = REPO / "results"
REBUILD = RESULTS / "cpad2026-rebuild"

INPUTS = [
    REPO / "sensor_workflow.json",
    REPO / "mc_dreamsim_e+_100_101.root",
    REPO / "mc_dreamsim_pi+_1000_1001.root",
    REPO / "mc_dreamsim_e+_0_run100_0_FullDetector_10evt_e+_100_101.root",
    REPO / "sim_200evt_pi+_90+rod.root",
]

DECK_ASSETS = {
    RESULTS / "cpad2026-audience" / "photon_xy_comparison.pdf":
        REPO / "cpad2026" / "assets" / "photon_xy_comparison.pdf",
    RESULTS / "cpad2026-audience" / "photon_xy_single_event.pdf":
        REPO / "cpad2026" / "assets" / "photon_xy_single_event.pdf",
    RESULTS / "cpad2026-audience" / "photon_time_profiles.png":
        REPO / "cpad2026" / "assets" / "photon_time_profiles.png",
    RESULTS / "cpad2026-one-fiber" / "one_fiber_segmentation.png":
        REPO / "cpad2026" / "assets" / "one_fiber_segmentation.png",
    RESULTS / "cpad2026-validation" / "validation_event_occupancy.png":
        REPO / "cpad2026" / "assets" / "validation_event_occupancy.png",
    RESULTS / "cpad2026-validation" / "validation_event_distributions.png":
        REPO / "cpad2026" / "assets" / "validation_event_distributions.png",
    RESULTS / "cpad2026-validation" / "validation_transport_depth.png":
        REPO / "cpad2026" / "assets" / "validation_transport_depth.png",
    RESULTS / "cpad2026-validation" / "validation_historical_closure.png":
        REPO / "cpad2026" / "assets" / "validation_historical_closure.png",
    RESULTS / "flash-study-v2-final" / "legacy_vs_physical_mapping.png":
        REPO / "cpad2026" / "assets" / "v2_mapping_correction.png",
    RESULTS / "flash-study-v2-final" / "physical_sensor_requirements.png":
        REPO / "cpad2026" / "assets" / "v2_physical_requirements.png",
    RESULTS / "flash-study-v2-final" / "event_level_tail_risk.png":
        REPO / "cpad2026" / "assets" / "v2_event_tail.png",
    RESULTS / "flash-study-v2-final" / "physical_pulse_shape.png":
        REPO / "cpad2026" / "assets" / "v2_pulse_shape.png",
    RESULTS / "flash-study-v2-final" / "timing_depth_proxy.png":
        REPO / "cpad2026" / "assets" / "v2_timing_depth.png",
    RESULTS / "fill-factor-spectral-v1" / "fill_factor_pitch_matrix.png":
        REPO / "cpad2026" / "assets" / "fill_factor_pitch_matrix.png",
    RESULTS / "channel-architecture-v2" / "channel_architecture_matrix.png":
        REPO / "cpad2026" / "assets" / "channel_architecture_matrix.png",
    RESULTS / "full-detector-v1" / "full_detector_validation.png":
        REPO / "cpad2026" / "assets" / "full_detector_validation.png",
    RESULTS / "full-detector-v1" / "full_detector_event_display.png":
        REPO / "cpad2026" / "assets" / "full_detector_event_display.png",
    RESULTS / "full-detector-v1" / "full_detector_event_display.pdf":
        REPO / "cpad2026" / "assets" / "full_detector_event_display.pdf",
    RESULTS / "full-detector-v1" / "full_detector_occupancy.png":
        REPO / "cpad2026" / "assets" / "full_detector_occupancy.png",
    RESULTS / "full-detector-v1" / "full_detector_performance_scan.png":
        REPO / "cpad2026" / "assets" / "full_detector_performance_scan.png",
    RESULTS / "full-detector-v1" / "full_detector_response.png":
        REPO / "cpad2026" / "assets" / "full_detector_response.png",
}

PAPER_FIGURES = {
    source.with_suffix(".pdf"): REPO / "paper" / "figures" /
    target.with_suffix(".pdf").name
    for source, target in DECK_ASSETS.items()
}


def run(command: list[str], *, cwd: Path = REPO) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, cwd=cwd, check=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def copy_outputs(mapping: dict[Path, Path]) -> None:
    for source, target in mapping.items():
        if not source.exists():
            raise FileNotFoundError(source)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)


def compile_deck() -> Path:
    deck_dir = REPO / "cpad2026"
    source = "CPAD2026_FLASH_rebuild.tex"
    for _ in range(2):
        run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", source],
            cwd=deck_dir)
    return deck_dir / source.replace(".tex", ".pdf")


def compile_paper() -> Path:
    paper_dir = REPO / "paper"
    source = "flash_hgdream_jinst.tex"
    run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", source],
        cwd=paper_dir)
    run(["bibtex", "flash_hgdream_jinst"], cwd=paper_dir)
    for _ in range(2):
        run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", source],
            cwd=paper_dir)
    return paper_dir / "flash_hgdream_jinst.pdf"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    analysis_mode = parser.add_mutually_exclusive_group()
    analysis_mode.add_argument(
        "--skip-analysis",
        action="store_true",
        help="Restage and compile from existing analysis products.",
    )
    analysis_mode.add_argument(
        "--sensor-only",
        action="store_true",
        help=(
            "Rebuild the fill-factor, channel-architecture, and full-detector "
            "studies from the shared sensor workflow, then restage and compile."
        ),
    )
    parser.add_argument(
        "--bootstrap-samples",
        type=int,
        default=500,
        help="Event-bootstrap replicas for the corrected studies.",
    )
    args = parser.parse_args()
    REBUILD.mkdir(parents=True, exist_ok=True)
    manifest_path = REBUILD / "rebuild_manifest.json"
    run_token = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_root = REBUILD / "runs" / run_token
    historical_root = run_root / "historical"
    flash_root = run_root / "flash-study-v2-final"
    fill_factor_root = run_root / "fill-factor-spectral-v1"
    channel_root = run_root / "channel-architecture-v2"
    full_detector_root = run_root / "full-detector-v1"
    commands: list[list[str]] = []
    analysis_run_utc = datetime.now(timezone.utc).isoformat()

    if not args.skip_analysis:
        sensor_commands = [
            [sys.executable, "studyFillFactorSpectral.py",
             "--output", str(fill_factor_root)],
            [sys.executable, "studyChannelArchitecture.py",
             "--output", str(channel_root)],
            [sys.executable, "studyFullDetector.py",
             "--output", str(full_detector_root),
             "--bootstrap-samples", str(args.bootstrap_samples)],
        ]
        commands = [
            [sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"],
        ]
        if args.sensor_only:
            commands.extend(sensor_commands)
        else:
            commands.extend(
                [
                    [sys.executable, "reproduceStudies.py", "2025-both",
                     "--plot-style", "both", "--output-root", str(historical_root)],
                    [sys.executable, "makeCpadAudiencePlots.py"],
                    [sys.executable, "makeOneFiberTalkPlots.py",
                     "--electron", str(historical_root / "2025-both" /
                                        "100GeV_eplus" / "comparison.root"),
                     "--pion", str(historical_root / "2025-both" /
                                    "1TeV_piplus" / "comparison.root")],
                    [sys.executable, "makeCpadValidationPlots.py"],
                    [sys.executable, "studyFlashV2.py",
                     "--output", str(flash_root),
                     "--bootstrap-samples", str(args.bootstrap_samples)],
                ]
            )
            commands.extend(sensor_commands)
        for command in commands:
            run(command)
        if not args.sensor_only:
            shutil.copytree(
                historical_root / "2025-both",
                RESULTS / "2025-both",
                dirs_exist_ok=True,
            )
            shutil.copytree(
                flash_root,
                RESULTS / "flash-study-v2-final",
                dirs_exist_ok=True,
            )
        for source, target in (
            (fill_factor_root, RESULTS / "fill-factor-spectral-v1"),
            (channel_root, RESULTS / "channel-architecture-v2"),
            (full_detector_root, RESULTS / "full-detector-v1"),
        ):
            shutil.copytree(source, target, dirs_exist_ok=True)
    elif manifest_path.exists():
        with manifest_path.open() as handle:
            previous_manifest = json.load(handle)
        commands = previous_manifest.get("commands", [])
        analysis_run_utc = previous_manifest.get(
            "analysis_run_utc",
            previous_manifest.get("generated_utc", analysis_run_utc),
        )

    copy_outputs(DECK_ASSETS)
    copy_outputs(PAPER_FIGURES)
    deck_pdf = compile_deck()
    paper_pdf = compile_paper()

    final_outputs = [deck_pdf, paper_pdf]
    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_run_utc": analysis_run_utc,
        "analysis_rerun_this_invocation": not args.skip_analysis,
        "analysis_scope": (
            "cached products"
            if args.skip_analysis
            else "shared sensor workflow"
            if args.sensor_only
            else "complete CPAD/JINST workflow"
        ),
        "python": sys.executable,
        "commands": commands,
        "inputs": {
            str(path.relative_to(REPO)): {
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in INPUTS
        },
        "final_outputs": {
            str(path.relative_to(REPO)): {
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in final_outputs
        },
        "analysis_products": {
            str(path.relative_to(REPO)): sha256(path)
            for path in sorted(DECK_ASSETS)
        },
    }
    with manifest_path.open("w") as handle:
        json.dump(manifest, handle, indent=2)
    print(json.dumps(manifest["final_outputs"], indent=2))


if __name__ == "__main__":
    main()
