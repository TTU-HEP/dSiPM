#!/usr/bin/env python3
"""Validate a DREAMSim ROOT event file and write a provenance sidecar."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import socket
import sys

try:
    import ROOT
except ModuleNotFoundError:
    ROOT = None


REQUIRED_BRANCHES = (
    "beamType",
    "beamMinE",
    "beamMaxE",
    "eCalotruth",
    "eLeaktruth",
    "sum3dCC",
    "sum3dSS",
    "nOPs",
    "OP_trackid",
    "OP_pos_produced_x",
    "OP_pos_produced_y",
    "OP_pos_produced_z",
    "OP_pos_final_x",
    "OP_pos_final_y",
    "OP_pos_final_z",
    "OP_time_produced",
    "OP_time_final",
    "OP_isCerenkov",
    "OP_isScintillation",
    "OP_productionFiber",
    "OP_finalFiber",
    "OP_isCoreC",
    "OP_isCoreS",
)

VECTOR_BRANCHES = tuple(name for name in REQUIRED_BRANCHES if name.startswith("OP_"))


def _beam_type(value: object) -> str:
    try:
        return bytes(value).decode(errors="replace").rstrip("\x00")
    except TypeError:
        return str(value).rstrip("\x00")


def _count_true(values: object) -> int:
    return sum(bool(value) for value in values)


def _empty_counters() -> dict[str, int]:
    return {
        "events_inspected": 0,
        "events_with_optical_photons": 0,
        "optical_photons": 0,
        "creator_cherenkov": 0,
        "creator_scintillation": 0,
        "produced_in_cherenkov_core": 0,
        "produced_in_scintillation_core": 0,
        "positive_z_cherenkov_core": 0,
        "positive_z_scintillation_core": 0,
    }


def _validate_pyroot(args: argparse.Namespace) -> dict[str, object]:
    ROOT.gROOT.SetBatch(True)
    path = args.root_file.resolve()
    root_file = ROOT.TFile.Open(str(path), "READ")
    if not root_file or root_file.IsZombie():
        raise ValueError(f"ROOT file is unreadable or a zombie: {path}")

    tree = root_file.Get("tree")
    if not tree:
        root_file.Close()
        raise ValueError(f"ROOT file has no 'tree': {path}")

    entries = int(tree.GetEntries())
    if entries != args.expected_events:
        root_file.Close()
        raise ValueError(
            f"Expected {args.expected_events} events but found {entries}: {path}"
        )

    branch_names = {branch.GetName() for branch in tree.GetListOfBranches()}
    missing = sorted(set(REQUIRED_BRANCHES) - branch_names)
    if missing:
        root_file.Close()
        raise ValueError(f"Missing required branches in {path}: {', '.join(missing)}")

    counters = _empty_counters()
    observed_particle = ""
    observed_energy = None
    inspect_events = min(entries, args.max_events_for_counters)

    for index, event in enumerate(tree):
        if index >= inspect_events:
            break
        counters["events_inspected"] += 1
        if index == 0:
            observed_particle = _beam_type(event.beamType)
            observed_energy = 0.5 * (float(event.beamMinE) + float(event.beamMaxE))
        nops = int(event.nOPs)
        counters["optical_photons"] += nops
        counters["events_with_optical_photons"] += int(nops > 0)

        lengths = {name: len(getattr(event, name)) for name in VECTOR_BRANCHES}
        inconsistent = {name: size for name, size in lengths.items() if size != nops}
        if inconsistent:
            root_file.Close()
            raise ValueError(
                f"Event {index} has nOPs={nops} but inconsistent vector sizes: "
                f"{inconsistent}"
            )

        counters["creator_cherenkov"] += _count_true(event.OP_isCerenkov)
        counters["creator_scintillation"] += _count_true(event.OP_isScintillation)
        counters["produced_in_cherenkov_core"] += _count_true(event.OP_isCoreC)
        counters["produced_in_scintillation_core"] += _count_true(event.OP_isCoreS)
        counters["positive_z_cherenkov_core"] += sum(
            bool(core) and float(z) > 0.0
            for core, z in zip(event.OP_isCoreC, event.OP_pos_final_z, strict=True)
        )
        counters["positive_z_scintillation_core"] += sum(
            bool(core) and float(z) > 0.0
            for core, z in zip(event.OP_isCoreS, event.OP_pos_final_z, strict=True)
        )

    root_file.Close()

    return {
        "path": path,
        "entries": entries,
        "branch_names": branch_names,
        "observed_particle": observed_particle,
        "observed_energy": observed_energy,
        "counters": counters,
        "backend": "PyROOT",
    }


def _validate_uproot(args: argparse.Namespace) -> dict[str, object]:
    import awkward as ak
    import uproot

    path = args.root_file.resolve()
    with uproot.open(path) as root_file:
        if "tree" not in root_file:
            raise ValueError(f"ROOT file has no 'tree': {path}")
        tree = root_file["tree"]
        entries = int(tree.num_entries)
        if entries != args.expected_events:
            raise ValueError(
                f"Expected {args.expected_events} events but found {entries}: {path}"
            )

        branch_names = set(tree.keys())
        missing = sorted(set(REQUIRED_BRANCHES) - branch_names)
        if missing:
            raise ValueError(
                f"Missing required branches in {path}: {', '.join(missing)}"
            )

        inspect_events = min(entries, args.max_events_for_counters)
        arrays = tree.arrays(
            list(REQUIRED_BRANCHES),
            entry_stop=inspect_events,
            library="ak",
        )
        counters = _empty_counters()
        observed_particle = ""
        observed_energy = None

        for index in range(inspect_events):
            counters["events_inspected"] += 1
            if index == 0:
                observed_particle = _beam_type(
                    ak.to_list(arrays["beamType"][index])
                )
                observed_energy = 0.5 * (
                    float(arrays["beamMinE"][index])
                    + float(arrays["beamMaxE"][index])
                )

            nops = int(arrays["nOPs"][index])
            counters["optical_photons"] += nops
            counters["events_with_optical_photons"] += int(nops > 0)

            lengths = {
                name: len(arrays[name][index]) for name in VECTOR_BRANCHES
            }
            inconsistent = {
                name: size for name, size in lengths.items() if size != nops
            }
            if inconsistent:
                raise ValueError(
                    f"Event {index} has nOPs={nops} but inconsistent vector sizes: "
                    f"{inconsistent}"
                )

            is_cherenkov = arrays["OP_isCerenkov"][index]
            is_scintillation = arrays["OP_isScintillation"][index]
            is_core_c = arrays["OP_isCoreC"][index]
            is_core_s = arrays["OP_isCoreS"][index]
            final_z = arrays["OP_pos_final_z"][index]
            counters["creator_cherenkov"] += int(ak.sum(is_cherenkov))
            counters["creator_scintillation"] += int(ak.sum(is_scintillation))
            counters["produced_in_cherenkov_core"] += int(ak.sum(is_core_c))
            counters["produced_in_scintillation_core"] += int(ak.sum(is_core_s))
            counters["positive_z_cherenkov_core"] += int(
                ak.sum(is_core_c & (final_z > 0.0))
            )
            counters["positive_z_scintillation_core"] += int(
                ak.sum(is_core_s & (final_z > 0.0))
            )

    return {
        "path": path,
        "entries": entries,
        "branch_names": branch_names,
        "observed_particle": observed_particle,
        "observed_energy": observed_energy,
        "counters": counters,
        "backend": "uproot",
    }


def validate(args: argparse.Namespace) -> dict[str, object]:
    result = _validate_pyroot(args) if ROOT is not None else _validate_uproot(args)
    path = result["path"]
    observed_particle = result["observed_particle"]
    observed_energy = result["observed_energy"]

    if observed_particle != args.expected_particle:
        raise ValueError(
            f"Expected particle {args.expected_particle!r}, observed {observed_particle!r}"
        )
    if observed_energy is None or abs(observed_energy - args.expected_energy_gev) > 1e-6:
        raise ValueError(
            f"Expected {args.expected_energy_gev:g} GeV, observed {observed_energy}"
        )

    return {
        "status": "valid",
        "root_file": str(path),
        "size_bytes": path.stat().st_size,
        "entries": result["entries"],
        "branch_count": len(result["branch_names"]),
        "validation_backend": result["backend"],
        "required_branches": list(REQUIRED_BRANCHES),
        "expected_particle": args.expected_particle,
        "observed_particle": observed_particle,
        "expected_energy_gev": args.expected_energy_gev,
        "observed_energy_gev": observed_energy,
        "dreamsim_commit": args.dreamsim_commit,
        "dsipm_commit": args.dsipm_commit,
        "macro_sha256": args.macro_sha256,
        "campaign_sha256": args.campaign_sha256,
        "launcher_sha256": args.launcher_sha256,
        "validator_sha256": args.validator_sha256,
        "sim_executable_sha256": args.sim_executable_sha256,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "execution_host": socket.gethostname(),
        "seed_1": args.seed_1,
        "seed_2": args.seed_2,
        "optical_branch_counters": result["counters"],
        "scope": (
            "TTU dSiPM DREAMSim full-shower run; positive-z Cherenkov-core "
            "photons retained for rod indices 20--60 and layer indices 15--65."
        ),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("root_file", type=Path)
    parser.add_argument("--expected-events", type=int, required=True)
    parser.add_argument("--expected-particle", required=True)
    parser.add_argument("--expected-energy-gev", type=float, required=True)
    parser.add_argument("--dreamsim-commit", required=True)
    parser.add_argument("--dsipm-commit", required=True)
    parser.add_argument("--macro-sha256", required=True)
    parser.add_argument("--campaign-sha256", required=True)
    parser.add_argument("--launcher-sha256", required=True)
    parser.add_argument("--validator-sha256", required=True)
    parser.add_argument("--sim-executable-sha256", required=True)
    parser.add_argument("--seed-1", type=int, required=True)
    parser.add_argument("--seed-2", type=int, required=True)
    parser.add_argument("--max-events-for-counters", type=int, default=25)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        summary = validate(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    output = args.output or args.root_file.with_suffix(
        args.root_file.suffix + ".validation.json"
    )
    output.write_text(json.dumps(summary, indent=2) + "\n")
    print(
        f"Validated {summary['root_file']}: {summary['entries']} events, "
        f"{summary['optical_branch_counters']['optical_photons']} photons inspected"
    )
    print(f"Wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
