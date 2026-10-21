#!/usr/bin/env python3
"""Calibrate and evaluate immutable compact caches; no event generation.

See docs/DEPTH_RECONSTRUCTION.md. Raw TTU ROOT extraction has a separate
Slurm-only adapter. ``scan`` never tunes on the held-out evaluation input.
"""
import argparse
import json
from pathlib import Path
import subprocess

import numpy as np

from dsipm.depth_io import load_cache, file_sha256
from dsipm.depth_reconstruction import Calibration, fit_calibration, reconstruct
from dsipm.depth_observations import smear_observations
from dsipm.depth_evaluation import evaluate


def calibrate(args):
    obs, contract, truth, meta = load_cache(args.input)
    if not meta["calibration_eligible"] or "kappa_truth" not in truth:
        raise ValueError("Input is not validated calibration truth; do not train on transport proxies")
    model = fit_calibration(obs, truth["kappa_truth"], contract, source_sha256=meta["input_sha256"],
                            min_count=args.min_count)
    model.save(args.output)
    print(f"Calibration saved: {len(model.payload['training_events'])} events; SHA-256 {model.sha256}")


def scan(args):
    obs, contract, truth, meta = load_cache(args.input)
    if meta["provenance"].get("timestamp_sigma_ns") != 0:
        raise ValueError("Matched scan requires an explicitly noiseless timestamp cache")
    calibration = Calibration.load(args.calibration)
    if calibration.payload["training_sigma_ns"] != 0:
        raise ValueError("This CLI scans a fixed ZERO-jitter calibration; use API for retrained studies")
    if not args.sigmas_ps or len(set(args.sigmas_ps)) != len(args.sigmas_ps) or min(args.sigmas_ps) < 0 or args.draws < 1:
        raise ValueError("Invalid scan points/draw count")
    if not np.isfinite(args.sigmas_ps).all() or args.clock_ps < 0 or not np.isfinite(args.clock_ps):
        raise ValueError("Invalid timestamp/clock uncertainty")
    # Refuse overlap/compatibility errors BEFORE creating any output.
    reconstruct(obs, calibration, contract, method="radius")
    target = Path(args.output)
    target.mkdir(parents=True, exist_ok=False)
    try:
        revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=Path(__file__).parent, text=True).strip()
        dirty = bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=Path(__file__).parent, text=True))
    except (OSError, subprocess.CalledProcessError):
        revision, dirty = None, None
    paths = list((Path(__file__).parent/"dsipm").glob("depth_*.py"))+[Path(__file__)]
    code_hashes = {p.name: file_sha256(p) for p in sorted(paths)}
    manifest = dict(status="development_diagnostic_not_certified_physics", source_cache=meta,
        calibration_sha256=calibration.sha256, calibration_mode="fixed_zero_jitter",
        sigma_ps=args.sigmas_ps, entrance_sigma_ps=args.clock_ps, draws=args.draws,
        seed=args.seed, code_revision=revision, dirty_checkout=dirty, code_sha256=code_hashes,
        target="individual accepted-photon creation depth; no photon rejected using truth",
        independent_events=sorted(set(obs.event_key.tolist())), methods=["radius", "radius_delay", "neighbour"],
        claims={"software_regression": "supported_by_tests", "detector_performance": "provisional",
                "known_transport_bug_fixed": "unsupported"}, results=[])
    if args.fixed_n_eff is not None:
        manifest["methods"].insert(0, "fixed")
    for sigma in args.sigmas_ps:
        for draw in range(args.draws):
            noisy = smear_observations(obs, sigma_t_ns=sigma/1000, entrance_sigma_ns=args.clock_ps/1000,
                                       seed=args.seed, draw=draw)
            for method in manifest["methods"]:
                result = reconstruct(noisy, calibration, contract, method=method, sigma_t_ns=sigma/1000,
                    entrance_sigma_ns=args.clock_ps/1000, window=args.window,
                    clean_radius=args.clean_radius, exclude_same_rod=args.leave_rod_out,
                    fixed_n_eff=args.fixed_n_eff)
                tag = f"{method}_{sigma:g}ps_draw{draw}"
                output = target/(tag+".npz")
                np.savez_compressed(output, event_key=obs.event_key, photon_id=obs.photon_id,
                    z_mm=result.z_mm, status=result.status, candidate_z_mm=result.candidate_z_mm,
                    candidate_density=result.candidate_density, prior_seed=result.prior_seed,
                    selected_kappa=result.selected_kappa,
                    photon_timing_sigma_z_mm=result.photon_timing_sigma_z_mm,
                    entrance_timing_sigma_z_mm=result.entrance_timing_sigma_z_mm)
                record = dict(name=tag, options=result.options, sha256=file_sha256(output))
                if "truth_z_mm" in truth:
                    record["evaluation"] = evaluate(noisy, truth["truth_z_mm"], result, contract,
                                                     bootstrap_draws=args.bootstrap_draws, seed=args.seed)
                manifest["results"].append(record)
    with (target/"summary.json").open("x") as f:
        json.dump(manifest, f, indent=2, allow_nan=False); f.write("\n")
    print(f"Saved {len(manifest['results'])} matched reconstructions in {target}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="mode", required=True)
    cal = sub.add_parser("calibrate")
    cal.add_argument("--input", required=True); cal.add_argument("--output", required=True)
    cal.add_argument("--min-count", type=int, default=20); cal.set_defaults(run=calibrate)
    run = sub.add_parser("scan")
    run.add_argument("--input", required=True); run.add_argument("--calibration", required=True)
    run.add_argument("--output", required=True)
    run.add_argument("--sigmas-ps", nargs="+", type=float, default=[0, 50, 200, 500])
    run.add_argument("--clock-ps", type=float, default=0.)
    run.add_argument("--draws", type=int, default=1)
    run.add_argument("--bootstrap-draws", type=int, default=200)
    run.add_argument("--seed", type=int, default=20261021)
    run.add_argument("--window", type=int, default=5)
    run.add_argument("--clean-radius", type=float, default=.7)
    run.add_argument("--leave-rod-out", action="store_true")
    run.add_argument("--fixed-n-eff", type=float, help="Optional frozen historical fixed-path baseline, not refit to test data")
    run.set_defaults(run=scan)
    args = ap.parse_args(); args.run(args)


if __name__ == "__main__":
    main()
