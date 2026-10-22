#!/usr/bin/env python3
"""Calibrate and evaluate immutable compact caches; no event generation.

See README.md. Raw TTU ROOT extraction has a separate
Slurm-only adapter. ``scan`` never tunes on the held-out evaluation input.
"""
import argparse
import json
from pathlib import Path
import subprocess

import numpy as np

from dsipm.io import load_cache, file_sha256
from dsipm.depth import Calibration, fit_calibration, fit_timing_calibration, reconstruct
from dsipm.sensors import smear_observations
from dsipm.depth import evaluate, response_matrix


def calibrate(args):
    obs, contract, truth, meta = load_cache(args.input)
    if meta['provenance'].get('timestamp_sigma_ns')!=0:
        raise ValueError('Calibration requires an explicitly noiseless input cache before controlled jitter')
    noisy=smear_observations(obs,sigma_t_ns=args.sigma_ps/1000,entrance_sigma_ns=args.clock_ps/1000,seed=args.seed)
    options=dict(source_sha256=meta['input_sha256'],min_count=args.min_count,
                 training_sigma_ns=args.sigma_ps/1000,training_entrance_sigma_ns=args.clock_ps/1000)
    if args.factor_kind=='timing_effective':
        if 'truth_z_mm' not in truth:raise ValueError('Timing calibration requires production-Z truth')
        model=fit_timing_calibration(noisy,truth['truth_z_mm'],contract,
            truth_z_validation=meta['provenance'].get('truth_z_validation',''),**options)
    else:
        if not meta["calibration_eligible"] or "kappa_truth" not in truth:
            raise ValueError("Input is not validated calibration truth; do not train on transport proxies")
        model=fit_calibration(noisy,truth['kappa_truth'],contract,**options)
    model=Calibration(dict(model.payload,training_noise_seed=args.seed))
    model.save(args.output)
    print(f"Calibration saved: {len(model.payload['training_events'])} events; SHA-256 {model.sha256}")


def scan(args):
    obs, contract, truth, meta = load_cache(args.input)
    if meta["provenance"].get("timestamp_sigma_ns") != 0:
        raise ValueError("Matched scan requires an explicitly noiseless timestamp cache")
    loaded=[Calibration.load(p) for p in args.calibration]
    calibrations={c.payload['training_sigma_ns']:c for c in loaded}
    if len(calibrations)!=len(loaded):raise ValueError('Duplicate calibration timing sigma')
    if args.allow_jitter_mismatch and len(loaded)!=1:raise ValueError('Fixed-calibration override takes exactly one model')
    if not args.sigmas_ps or len(set(args.sigmas_ps)) != len(args.sigmas_ps) or min(args.sigmas_ps) < 0 or args.draws < 1:
        raise ValueError("Invalid scan points/draw count")
    if not np.isfinite(args.sigmas_ps).all() or args.clock_ps < 0 or not np.isfinite(args.clock_ps):
        raise ValueError("Invalid timestamp/clock uncertainty")
    # Refuse overlap/compatibility errors BEFORE creating any output.
    for sigma in args.sigmas_ps:
        calibration=loaded[0] if args.allow_jitter_mismatch else calibrations.get(sigma/1000)
        if calibration is None:raise ValueError('Missing jitter-matched calibration for '+str(sigma)+' ps')
        reconstruct(obs,calibration,contract,method='radius',sigma_t_ns=sigma/1000,entrance_sigma_ns=args.clock_ps/1000,
                    allow_jitter_mismatch=args.allow_jitter_mismatch)
    target = Path(args.output)
    target.mkdir(parents=True, exist_ok=False)
    try:
        revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=Path(__file__).parent, text=True).strip()
        dirty = bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=Path(__file__).parent, text=True))
    except (OSError, subprocess.CalledProcessError):
        revision, dirty = None, None
    paths = list((Path(__file__).parent/"dsipm").glob("*.py"))+[Path(__file__)]
    code_hashes = {p.name: file_sha256(p) for p in sorted(paths)}
    manifest = dict(status="development_diagnostic_not_certified_physics", source_cache=meta,
        calibration_sha256={str(k):c.sha256 for k,c in calibrations.items()},
        calibration_mode='explicit_fixed_calibration_override' if args.allow_jitter_mismatch else 'matched_per_sigma',
        sigma_ps=args.sigmas_ps, entrance_sigma_ps=args.clock_ps, draws=args.draws,
        seed=args.seed, code_revision=revision, dirty_checkout=dirty, code_sha256=code_hashes,
        target="individual accepted-photon creation depth; no photon rejected using truth",
        independent_events=sorted(set(obs.event_key.tolist())), methods=["radius", "radius_delay", "neighbour"],
        claims={"software_regression": "supported_by_tests", "detector_performance": "provisional",
                "known_transport_bug_fixed": "unsupported"}, results=[])
    if args.fixed_n_eff is not None:
        manifest["methods"].insert(0, "fixed")
    for sigma in args.sigmas_ps:
        calibration=loaded[0] if args.allow_jitter_mismatch else calibrations[sigma/1000]
        for draw in range(args.draws):
            noisy = smear_observations(obs, sigma_t_ns=sigma/1000, entrance_sigma_ns=args.clock_ps/1000,
                                       seed=args.seed, draw=draw)
            for method in manifest["methods"]:
                result = reconstruct(noisy, calibration, contract, method=method, sigma_t_ns=sigma/1000,
                    entrance_sigma_ns=args.clock_ps/1000, window=args.window,
                    clean_radius=args.clean_radius, exclude_same_rod=args.leave_rod_out,
                    fixed_n_eff=args.fixed_n_eff,allow_jitter_mismatch=args.allow_jitter_mismatch)
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


def _plot_axes(label, xlabel, xlim, *, colorbar=False):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "axes.labelweight": "bold", "axes.labelsize": 14, "font.size": 12,
        "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
        "pdf.fonttype": 42, "axes.linewidth": 1.3})
    fig = plt.figure(figsize=(6.9, 6.6))
    ax = fig.add_axes([.16, .15, .68 if colorbar else .77, .72])
    ax.set_box_aspect(1); ax.minorticks_on()
    ax.text(0, 1.04, "FLASH", transform=ax.transAxes, weight="bold", fontsize=15)
    ax.text(.22, 1.04, "Simulation Preliminary", transform=ax.transAxes, style="italic", fontsize=12)
    fig.text(.945, .965, label, ha="right", va="top", fontsize=9)
    ax.set_xlim(*xlim); ax.set_xlabel(xlabel)
    return fig, ax


def render_response(product, target, *, label):
    """Optional scientific plotting; one labeled PDF/PNG per observable."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, LinearSegmentedColormap

    counts = np.asarray(product["counts"], float)
    inputs = np.asarray(product["input_by_truth"], float)
    failed = np.asarray(product["unreconstructed_by_truth"], float)
    edges = np.asarray(product["truth_edges_mm"])/10
    reco_edges = np.asarray(product["reco_edges_mm"])/10
    if not np.array_equal(counts.sum(axis=1)+failed, inputs):
        raise ValueError("Response plot closure failure")
    target = Path(target)
    target.mkdir(parents=True, exist_ok=False)
    cmap = LinearSegmentedColormap.from_list("bird", ["#352a87", "#1062dc", "#108cd2", "#17acb9", "#60bc8b", "#babc62", "#f7c638", "#f9f90f"])
    cmap.set_bad("white")
    def axes(colorbar=False):
        return _plot_axes(label,"True photon-production Z [cm]",(edges[0],edges[-1]),colorbar=colorbar)
    outputs = []
    def save(fig, name):
        for ext in ("pdf", "png"):
            path = target/(name+"."+ext)
            fig.savefig(path, dpi=220, facecolor="white")
            outputs.append(dict(file=path.name, sha256=file_sha256(path)))
        plt.close(fig)
    fig, ax = axes(True)
    probability = np.divide(100*counts, inputs[:, None], out=np.zeros_like(counts), where=inputs[:, None]>0)
    image = ax.pcolormesh(edges, reco_edges, np.ma.masked_equal(probability.T, 0), cmap=cmap,
                         norm=LogNorm(vmin=.01, vmax=100), shading="flat", rasterized=True)
    ax.plot([edges[0], edges[-1]], [edges[0], edges[-1]], "--", color="#ce007a", lw=1.5, label="Perfect reconstruction")
    ax.set_ylim(reco_edges[0], reco_edges[-1]); ax.set_ylabel("Reconstructed Z [cm]")
    ax.legend(loc="lower right", fontsize=9)
    fig.colorbar(image, ax=ax, fraction=.045, pad=.02).set_label("Input fraction per true-Z bin [% / reco bin]", weight="bold", fontsize=11)
    save(fig, "truth_vs_reconstructed")
    centers = (edges[:-1]+edges[1:])/2
    q = np.asarray(product["conditional_residual_q16_q50_q84_mm"], float)/10
    fig, ax = axes()
    ax.axhline(0, color="#ce007a", ls="--", lw=1.3)
    ax.fill_between(centers, q[:, 0], q[:, 2], color="#1676ad", alpha=.2, label="Central 68%, physical returns")
    ax.plot(centers, q[:, 1], color="#1676ad", lw=1.6, label=product.get("median_label", "Median residual"))
    ax.set_ylabel(product.get("residual_ylabel", "Reconstructed Z - true Z [cm]")); ax.legend(loc="lower left", fontsize=9)
    save(fig, "residual_vs_truth")
    fig, ax = axes()
    ax.plot(centers, np.divide(100*failed, inputs, out=np.full_like(inputs, np.nan), where=inputs>0),
            color="#ce007a", label="Unreconstructed")
    if product.get("within100mm_by_truth") is not None:
        within = np.asarray(product["within100mm_by_truth"], float)
        ax.plot(centers, np.divide(100*within, inputs, out=np.full_like(inputs, np.nan), where=inputs>0),
                color="#1676ad", label="Within 10 cm, all inputs")
    ax.set_ylim(0, 100); ax.set_ylabel("Fraction of all input [%]"); ax.legend(loc="lower right", fontsize=9)
    save(fig, "success_and_failure_vs_truth")
    with (target/"plot_manifest.json").open("x") as f:
        json.dump(dict(label=label, product=product, files=outputs,
            scope="development response, not certified detector performance"), f, indent=2, allow_nan=False)


def render_resolution(row, plan, target, *, label):
    """Inclusive all-photon accuracy versus assumed timestamp resolution."""
    import matplotlib.pyplot as plt
    points=np.asarray(plan['sigmas_ps'])
    fig,ax=_plot_axes(label,"Assumed single-photon timing sigma [ps]",(-10,float(points.max())+10))
    for model in plan['baselines']:
        cases=[row['cases'][model['name']+'_'+str(int(s))+'ps'] for s in points]
        y=np.array([r['within100mm_pct'] for r in cases])
        limits=np.array([r['file_bootstrap_95']['within100mm_pct'] for r in cases])
        color='#1676ad' if model['name']=='historical_fixed' else '#777777'
        name=('Frozen-index baseline' if model['name']=='historical_fixed' else 'Axial-path baseline')
        ax.fill_between(points,limits[:,0],limits[:,1],color=color,alpha=.2)
        ax.plot(points,y,'o-',color=color,label=name,markersize=5)
    ax.set_ylim(0,100)
    ax.set_ylabel("Accepted photons within 10 cm of true Z [%]")
    ax.legend(loc='upper left',bbox_to_anchor=(0,.86),fontsize=10)
    paths=[Path(target)/('accuracy_vs_resolution.'+ext) for ext in ('pdf','png')]
    if any(p.exists() for p in paths):
        plt.close(fig); raise FileExistsError('Resolution figure already exists')
    for path in paths: fig.savefig(path,dpi=220,facecolor='white')
    plt.close(fig)
    with (Path(target)/'resolution_manifest.json').open('x') as f:
        json.dump(dict(label=label,cases=row['cases'],assumptions=plan,
            interval='95% independent-source-file bootstrap; statistical only',
            denominator='all scored accepted photons, including failed inversions',
            lines='guides between tested timing assumptions, not an interpolated performance model',
            files=[dict(file=p.name,sha256=file_sha256(p)) for p in paths]),f,indent=2,allow_nan=False)


def plot(args):
    obs, contract, truth, meta = load_cache(args.input)
    if "truth_z_mm" not in truth:
        raise ValueError("Truth-versus-reconstruction plot requires an explicitly aligned true-Z cache")
    run = Path(args.scan_dir)
    summary = json.loads((run/"summary.json").read_text())
    if summary["source_cache"] != meta:
        raise ValueError("Plot cache does not match the recorded scan input")
    records = [r for r in summary["results"] if
        (args.method is None or r["options"]["method"] == args.method) and
        (args.sigma_ps is None or r["options"]["sigma_t_ns"] == args.sigma_ps/1000)]
    if not records:
        raise ValueError("No matching reconstruction results")
    products = []
    for row in records:
        if Path(row["name"]).name != row["name"]:
            raise ValueError("Unsafe reconstruction result name")
        path = run/(row["name"]+".npz")
        if path.parent != run or file_sha256(path) != row["sha256"]:
            raise ValueError("Reconstruction result path/hash mismatch")
        with np.load(path, allow_pickle=False) as f:
            if not np.array_equal(f["event_key"], obs.event_key) or not np.array_equal(f["photon_id"], obs.photon_id):
                raise ValueError("Reconstruction photon identity/alignment mismatch")
            if not np.array_equal(f["status"] == 2, ~np.isfinite(f["z_mm"])):
                raise ValueError("Reconstruction failure status mismatch")
            products.append((row["name"], response_matrix(truth["truth_z_mm"], f["z_mm"], length_mm=contract.length_mm)))
    target = Path(args.output); target.mkdir(parents=True, exist_ok=False)
    for name, product in products:
        render_response(product, target/name, label=args.label+"\n"+name)
    print(f"Saved {len(products)} truth/reconstruction response sets in {target}")


def plot_campaign(args):
    """Only a complete, hash-checked campaign checkpoint may feed paper plots."""
    from reduceDepthCampaign import residual_quantiles
    directory=Path(args.input)
    summary=json.loads((directory/'summary.json').read_text())
    plan=json.loads(Path(args.plan).read_text())
    if file_sha256(args.plan)!=summary['plan_sha256']:
        raise ValueError('Campaign plot reduction-plan checksum mismatch')
    mi=[m['name'] for m in plan['baselines']].index(args.baseline)
    expected={p+'_'+str(e)+'gev' for p in ('eplus','piplus') for e in (10,30,100,300)}
    if (summary['status']!='complete_accounting_physics_provisional' or not summary['closure']
            or summary['files']!=1648 or summary['events']!=40000 or set(summary['samples'])!=expected
            or any(s['events']!=5000 for s in summary['samples'].values())):
        raise ValueError('Plots require the complete separately validated 40,000-event checkpoint')
    products=[]
    for sample,row in sorted(summary['samples'].items()):
        name=row['histograms_file']
        if Path(name).name!=name or file_sha256(directory/name)!=row['histograms_sha256']:
            raise ValueError('Campaign plot histogram path/hash mismatch')
        with np.load(directory/name,allow_pickle=False) as h:
            points=np.flatnonzero(h['sigmas_ps']==args.sigma_ps)
            if len(points)!=1:
                raise ValueError('Unknown timing assumption in checkpoint')
            si=int(points[0])
            counts=h['response'][mi,si]; inputs=h['denominator'][mi,si]; failed=h['failures'][mi,si]
            residual=h['depth_residual'][mi,si]; within=h['within'][mi,si,:,1]
            case=row['cases'][args.baseline+'_'+str(int(args.sigma_ps))+'ps']
            if (not np.array_equal(counts.sum(axis=1)+failed,inputs)
                    or not np.array_equal(residual.sum(axis=1),counts.sum(axis=1))
                    or int(inputs.sum())!=case['input'] or int(failed.sum())!=case['failed']
                    or int(within.sum())!=case['within100mm']):
                raise ValueError('Campaign plot response/event/coverage closure failure')
            product=dict(truth_edges_mm=h['z_edges_mm'].tolist(),reco_edges_mm=h['z_edges_mm'].tolist(),
                counts=counts.tolist(),input_by_truth=inputs.tolist(),unreconstructed_by_truth=failed.tolist(),
                within100mm_by_truth=within.tolist(),conditional_residual_q16_q50_q84_mm=[
                    residual_quantiles(v,h['residual_edges_mm']) for v in residual],
                normalization='all scored accepted individual photons per true-Z bin, including failures',
                source_summary_sha256=file_sha256(directory/'summary.json'),source_histograms_sha256=row['histograms_sha256'],
                reduction_plan_sha256=summary['plan_sha256'],sample=sample,case=case,
                baseline=plan['baselines'][mi],sensor_contract=plan['contract'],
                entrance_sigma_ps=plan['entrance_sigma_ps'],scope=plan['scope'],selection=plan['selection'],
                optical_transport_validated=summary['optical_transport_validated'],
                quantile_bin_precision_mm=summary['residual_quantile_bin_width_mm'])
        particle,energy=sample.split('_')
        label=('e+' if particle=='eplus' else 'pi+')+' '+energy[:-3]+' GeV'
        label+='\n'+('Unfitted axial baseline' if args.baseline=='axial' else 'Frozen historical baseline')
        label+='; '+str(int(args.sigma_ps))+' ps'
        products.append((sample,product,label))
    target=Path(args.output);target.mkdir(parents=True,exist_ok=False)
    for sample,product,label in products:
        render_response(product,target/sample,label=label)
        render_resolution(summary['samples'][sample],plan,target/sample,
            label=label.split('\n')[0]+'\nIndependent-fiber proxy\nExact entrance time')
    import csv
    with (target/'all_resolution_results.csv').open('x',newline='') as f:
        writer=csv.writer(f)
        writer.writerow(['sample','events','baseline','sigma_ps','accepted_inputs','physical_returns','failures',
            'within10cm_pct','within10cm_stat95_low','within10cm_stat95_high',
            'failure_pct','failure_stat95_low','failure_stat95_high',
            'conditional_median_residual_cm','conditional_sigma68_cm'])
        for sample,row in sorted(summary['samples'].items()):
            for name,case in row['cases'].items():
                model,sigma=name.rsplit('_',1)
                q=case['conditional_residual_q16_q50_q84_mm']
                writer.writerow([sample,row['events'],model,int(sigma[:-2]),case['input'],case['reconstructed'],case['failed'],
                    case['within100mm_pct'],*case['file_bootstrap_95']['within100mm_pct'],
                    case['failure_pct'],*case['file_bootstrap_95']['failure_pct'],
                    None if q[1] is None else q[1]/10,
                    None if case['conditional_sigma68_mm'] is None else case['conditional_sigma68_mm']/10])
    print('Saved',len(products),'complete-campaign response sets in',target)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="mode", required=True)
    cal = sub.add_parser("calibrate")
    cal.add_argument("--input", required=True); cal.add_argument("--output", required=True)
    cal.add_argument("--min-count", type=int, default=20)
    cal.add_argument('--sigma-ps',type=float,default=0);cal.add_argument('--clock-ps',type=float,default=0)
    cal.add_argument('--seed',type=int,default=20261021)
    cal.add_argument('--factor-kind',choices=['geometric_path','timing_effective'],default='geometric_path')
    cal.set_defaults(run=calibrate)
    run = sub.add_parser("scan")
    run.add_argument("--input", required=True); run.add_argument("--calibration", required=True,nargs='+')
    run.add_argument('--allow-jitter-mismatch',action='store_true',help='Explicit fixed-calibration robustness study, not matched training')
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
    plots = sub.add_parser("plot", help="True-Z/reconstructed-Z, residual and inclusive failure plots; optional matplotlib")
    plots.add_argument("--input", required=True); plots.add_argument("--run", dest="scan_dir", required=True)
    plots.add_argument("--output", required=True); plots.add_argument("--label", required=True)
    plots.add_argument("--method", choices=["fixed", "radius", "radius_delay", "neighbour"])
    plots.add_argument("--sigma-ps", type=float); plots.set_defaults(run=plot)
    campaign=sub.add_parser("campaign-plot",help="All eight samples from a complete hash-checked campaign checkpoint")
    campaign.add_argument("--input",required=True);campaign.add_argument("--output",required=True)
    campaign.add_argument("--plan",required=True)
    campaign.add_argument("--sigma-ps",type=int,default=50)
    campaign.add_argument("--baseline",choices=["axial","historical_fixed"],default="axial")
    campaign.set_defaults(run=plot_campaign)
    args = ap.parse_args(); args.run(args)


if __name__ == "__main__":
    main()
