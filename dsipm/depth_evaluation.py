"""Truth-based evaluation, kept outside the reconstruction API.

Inclusive means ALL accepted input photons, not all generated photons or
deposited energy. Gate/PDE/holdoff losses must be reported separately. Truth
is photon creation depth; a channel median is not a resolved energy deposit.
"""
import numpy as np


def residual_summary(truth_mm, reconstructed_mm, *, tolerances_mm=(50., 100., 200.)):
    truth, z = np.asarray(truth_mm, float), np.asarray(reconstructed_mm, float)
    if truth.ndim != 1 or z.shape != truth.shape or not np.isfinite(truth).all():
        raise ValueError("Need finite aligned one-dimensional truth")
    if not np.isfinite(tolerances_mm).all() or np.any(np.asarray(tolerances_mm) <= 0):
        raise ValueError("Invalid tolerances")
    valid = np.isfinite(z)
    delta = z[valid]-truth[valid]
    q = np.percentile(delta, [16, 50, 84]) if len(delta) else [None]*3
    return dict(input=len(z), reconstructed=int(valid.sum()), unreconstructed=int((~valid).sum()),
        failure_fraction=float((~valid).mean()) if len(z) else None,
        conditional_bias_mm=None if not len(delta) else float(q[1]),
        conditional_68_halfwidth_mm=None if not len(delta) else float((q[2]-q[0])/2),
        inclusive_within={str(t): float(np.sum(abs(delta) <= t)/len(z)) if len(z) else None for t in tolerances_mm})


def evaluate(obs, truth_mm, result, contract, *, bootstrap_draws=200, seed=20261021,
             min_channel_photons=5, truth_edges_mm=None):
    truth = np.asarray(truth_mm, float)
    if truth.shape != obs.time_ns.shape or not np.isfinite(truth).all() or np.any((truth < 0) | (truth > contract.length_mm)):
        raise ValueError("Truth depth outside physical detector")
    if len(result.z_mm) != len(truth) or bootstrap_draws < 0 or min_channel_photons < 1:
        raise ValueError("Invalid evaluation options/lengths")
    if np.any(np.isfinite(result.z_mm) & ((result.z_mm < 0) | (result.z_mm > contract.length_mm))):
        raise ValueError("Unphysical depth incorrectly marked reconstructed")
    report = residual_summary(truth, result.z_mm)
    report["target"] = "creation depth of each accepted photon; failures retained"
    report["status_counts"] = {str(k): int(np.sum(result.status == k)) for k in [0, 1, 2]}
    if sum(report["status_counts"].values()) != len(truth):
        raise ValueError("Unknown reconstruction status")
    if not np.array_equal(result.status == 2, ~np.isfinite(result.z_mm)):
        raise ValueError("Failure-status/depth closure mismatch")
    report["excluded_invalid_prior_seeds"] = result.invalid_clean_seed_count
    edges = np.linspace(0, contract.length_mm, 11) if truth_edges_mm is None else np.asarray(truth_edges_mm, float)
    if not np.isfinite(edges).all() or len(edges) < 2 or np.any(np.diff(edges) <= 0):
        raise ValueError("Invalid truth-depth bins")
    bins = np.searchsorted(edges, truth, side="right")-1
    bins[truth == edges[-1]] = len(edges)-2
    report["depth_bins"] = []
    for b in range(len(edges)-1):
        mask = bins == b
        report["depth_bins"].append(dict(low_mm=float(edges[b]), high_mm=float(edges[b+1]),
                                        **residual_summary(truth[mask], result.z_mm[mask])))
    report["truth_underflow"] = int(np.sum(truth < edges[0]))
    report["truth_overflow"] = int(np.sum(truth > edges[-1]))
    if sum(b["input"] for b in report["depth_bins"])+report["truth_underflow"]+report["truth_overflow"] != len(truth):
        raise ValueError("Truth-depth coverage failure")
    events = np.unique(obs.event_key)
    groups = {}
    for j, key in enumerate(zip(obs.event_key.tolist(), obs.rod_x.tolist(), obs.rod_y.tolist())):
        groups.setdefault(key, []).append(j)
    cm_truth, cm_reco = [], []
    for ids in groups.values():
        if len(ids) < min_channel_photons:
            continue
        ids = np.asarray(ids)
        rz = result.z_mm[ids]
        good = np.isfinite(rz)
        cm_truth.append(float(np.median(truth[ids])))
        cm_reco.append(float(np.median(rz[good])) if good.sum() >= min_channel_photons else np.nan)
    report["channel_medians"] = dict(target="median of accepted-photon creation depths; not pre-holdoff or energy truth",
        min_input_and_reconstructed_photons=min_channel_photons, **residual_summary(cm_truth, cm_reco))
    report["events"] = len(events)
    # Event-block bootstrap of PHOTON-WEIGHTED metrics; never resample individual photons.
    draws = []
    if len(events) >= 2 and bootstrap_draws:
        indices = [np.flatnonzero(obs.event_key == e) for e in events]
        rng = np.random.default_rng(seed)
        for _ in range(bootstrap_draws):
            ids = np.concatenate([indices[i] for i in rng.integers(0, len(events), len(events))])
            draws.append(residual_summary(truth[ids], result.z_mm[ids]))
    intervals = {}
    for metric in ("conditional_bias_mm", "conditional_68_halfwidth_mm", "failure_fraction"):
        values = [d[metric] for d in draws if d[metric] is not None]
        intervals[metric] = np.percentile(values, [2.5, 97.5]).tolist() if values else None
    for tolerance in ("50.0", "100.0", "200.0"):
        values = [d["inclusive_within"][tolerance] for d in draws]
        intervals["inclusive_within_"+tolerance+"mm"] = np.percentile(values, [2.5, 97.5]).tolist() if values else None
    report["event_block_bootstrap"] = dict(draws=len(draws), seed=seed, percentile_95=intervals,
        unavailable_reason="need at least two independent events" if len(events) < 2 else None)
    report["closure"] = report["input"] == report["reconstructed"]+report["unreconstructed"]
    return report
