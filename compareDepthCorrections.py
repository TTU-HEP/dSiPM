#!/usr/bin/env python3
"""Matched 100 GeV pion depth methods, with whole-source-file cross-fitting.

The primary --method-* study tests dsipm.depth corrections 1-3 using timing-
effective factors and jitter-matched training. Direct-Z lookup models are
separate empirical controls. Legacy --compare checkpoints remain historical
controls, not tests of corrections 1-3. Neither approach validates transport.
Extraction and fitting/scoring require Slurm; plotting uses aggregates only.
"""
import argparse
import csv
import json
import os
from pathlib import Path

import numpy as np

from dsipm.depth import (C_MM_NS, DepthContract, PhotonObservations, channel_delay,
                        fit_timing_calibration_batches, reconstruct, Calibration, digest)
from dsipm.io import file_sha256
from dsipm.sensors import smear_observations
from reduceDepthCampaign import reduce_one

METHODS = ('axial', 'time_calibrated', 'time_radius', 'time_radius_delay',
           'time_radius_delay_neighbours', 'depth_prior')
LABELS = ('Basic timing inversion', 'Calibrated time', '+ Exit radius',
          '+ Same-rod arrival delay', '+ Neighbour timing', 'No-timing depth prior')
COLORS = ('#555555', '#2171b5', '#17a398', '#dd851b', '#a51c7b', '#9b895e')


def feature_spec():
    """Fixed before evaluation; explicit time flows, missing neighbours, zero delay."""
    return dict(time_edges_ns=np.round(np.arange(-2., 42.00001, .1), 6).tolist(),
                radius_edges=[.2, .4, .6, .8, .9, .95, .975, 1.], delay_edges_ns=[0., .1, .3, 1., 3., 10.],
                neighbour_edges_ns=[-1., -.25, .25, 1., 3., 10.],
                target_edges_mm=np.linspace(0., 2000., 101).tolist(),
                min_photons=400, min_sources=30, window_rods=5)


def dimensions(spec):
    return (len(spec['time_edges_ns'])+1, len(spec['radius_edges'])+1,
            len(spec['delay_edges_ns'])+1, len(spec['neighbour_edges_ns'])+2)


def observable_features(obs, contract, spec):
    """No truth inputs. Neighbours exclude the complete target rod, not just self.

    All times are measured after jitter. The context is the median first time
    of other occupied rods in a 5x5 window, with each rod weighted once.
    The known transverse beam-front offset is subtracted before comparison.
    No directions, creation time, or first-step time enter these features.
    """
    if spec['window_rods'] < 1 or spec['window_rods'] % 2 != 1:
        raise ValueError('Need an odd positive neighbourhood')
    # Accepted endpoints may lie in the cladding, not just the core. Radius is
    # the measured pixel centre, with half-diagonal quantization uncertainty.
    # Do not drop this valid observation population or relabel it path truth.
    if np.any(obs.radius > (contract.fiber_outer_radius_mm+1e-5)/contract.core_radius_mm+contract.radius_padding+1e-10):
        raise ValueError('Pixel radius incompatible with sensor contract')
    tau = obs.time_ns-(contract.front_x*obs.sensor_x_mm+contract.front_y*obs.sensor_y_mm
                       +contract.front_offset_mm)/C_MM_NS
    neighbour = np.full(len(tau), np.nan)
    w = spec['window_rods']//2
    for e in np.unique(obs.event_key):
        ids = np.flatnonzero(obs.event_key == e)
        rods = np.column_stack((obs.rod_x[ids], obs.rod_y[ids]))
        keys, inverse = np.unique(rods, axis=0, return_inverse=True)
        first = np.full(len(keys), np.inf)
        np.minimum.at(first, inverse, tau[ids])
        for j, (x, y) in enumerate(keys):
            use = (abs(keys[:,0]-x) <= w) & (abs(keys[:,1]-y) <= w)
            use[j] = False
            if np.any(use):
                at = ids[inverse == j]
                neighbour[at] = tau[at]-np.median(first[use])
    dt = channel_delay(obs)
    ti = np.searchsorted(spec['time_edges_ns'], tau, side='right')
    ri = np.searchsorted(spec['radius_edges'], obs.radius, side='right')
    di = np.searchsorted(spec['delay_edges_ns'], dt, side='right')
    di[dt == 0] = 0
    ni = np.searchsorted(spec['neighbour_edges_ns'], neighbour, side='right')
    ni[~np.isfinite(neighbour)] = len(spec['neighbour_edges_ns'])+1
    return np.column_stack((ti, ri, di, ni)).astype(np.int32)


def fit_tables(hist, source_support, spec):
    """Training-only conditional medians; sparse cells fall back to less detail.

    hist covers *all* training targets in 2 cm bins, including endpoints.
    No truth-dependent cuts or optical-path surrogate are used. Interpolated
    histogram medians have up to 2 cm binning uncertainty, not timing precision.
    """
    expected = dimensions(spec)+(len(spec['target_edges_mm'])-1,)
    if hist.shape != expected or hist.dtype.kind not in 'iu' or np.any(hist < 0) or not hist.sum():
        raise ValueError('Invalid empirical training histogram')
    if len(source_support) != 4:
        raise ValueError('Need independent-source support at every feature level')
    edges = np.asarray(spec['target_edges_mm'])
    tables, supported = [], []
    for level in range(5):
        h = hist.sum(axis=tuple(range(level, 4))) if level < 4 else hist
        count = h.sum(axis=-1)
        cumulative = np.cumsum(h, axis=-1)
        mid = count/2
        b = np.argmax(cumulative >= mid[..., None], axis=-1)
        before = np.take_along_axis(cumulative, b[..., None], axis=-1)[..., 0]
        inside = np.take_along_axis(h, b[..., None], axis=-1)[..., 0]
        before = before-inside
        fraction = np.divide(mid-before, inside, out=np.zeros_like(mid, dtype=float), where=inside > 0)
        value = edges[b]+fraction*(edges[b+1]-edges[b])
        tables.append(value)
        support = count > 0 if level == 0 else ((count >= spec['min_photons']) &
                                               (source_support[level-1] >= spec['min_sources']))
        if level and source_support[level-1].shape != count.shape:
            raise ValueError('Malformed independent-source support')
        supported.append(support)
    return tables, supported


def predict_tables(features, tables, supported):
    """Inference only; the returned level exposes every fallback to a depth prior.

    These conditional-median estimates are physical by construction. Zero
    range failures is NOT evidence of a physically consistent inverse solution.
    The no-timing prior and inclusive errors must accompany their performance.
    """
    n = len(features)
    z, levels = [], []
    current = np.full(n, float(tables[0])); used = np.zeros(n, np.int8)
    for level in range(1, 5):
        index = tuple(features[:, j] for j in range(level))
        valid = supported[level][index]
        current = np.where(valid, tables[level][index], current)
        used = np.where(valid, level, used)
        z.append(current.copy()); levels.append(used.copy())
    return z, levels


def load_observations(directory, row, plan):
    source = directory/'results'/Path(row['root_file']).stem
    report = json.loads((source/'audit.json').read_text())
    if (report['input'] != {k:v for k,v in row.items() if k!='fold'} or report['plan_sha256'] != plan['_sha256'] or
            report['contract_sha256'] != plan['contract_sha256'] or
            file_sha256(source/'accepted.npz') != report['observations_sha256'] or
            file_sha256(source/'histograms.npz') != report['histograms_sha256']):
        raise ValueError('Observation provenance/hash mismatch')
    with np.load(source/'accepted.npz', allow_pickle=False) as f:
        entry = f['entry']
        keys = np.array([row['root_sha256']+':'+str(int(e)) for e in entry])
        obs = PhotonObservations(**{k:keys if k == 'event_key' else f[k] for k in PhotonObservations.__dataclass_fields__})
        truth = f['truth_z_mm']
    if (len(report['records']) != row['events'] or {r['entry'] for r in report['records']} != set(range(row['events'])) or
            len(obs.time_ns) != sum(r['accepted'] for r in report['records'])):
        raise ValueError('Observation event/photon closure failed')
    if (truth.shape != obs.time_ns.shape or not np.isfinite(truth).all() or
            np.any((truth < 0) | (truth > 2000)) or np.any((entry < 0) | (entry >= row['events']))):
        raise ValueError('Invalid cached target or entry')
    for r in report['records']:
        if np.sum(entry == r['entry']) != r['accepted']:
            raise ValueError('Cached per-event count mismatch')
    return obs, truth, report


def score(hist, mi, si, truth, reco, ze, re):
    valid = np.isfinite(reco) & (reco >= 0) & (reco <= 2000)
    error = reco[valid]-truth[valid]
    hist['response'][mi,si] += np.histogram2d(truth[valid], reco[valid], (ze,ze))[0].astype(np.int64)
    hist['failures'][mi,si] += np.histogram(truth[~valid], ze)[0]
    hist['denominator'][mi,si] += np.histogram(truth, ze)[0]
    hist['residual'][mi,si] += np.histogram(error, re)[0]
    for j, tol in enumerate((50,100,200)):
        hist['within'][mi,si,:,j] += np.histogram(truth[valid][abs(error) <= tol], ze)[0]
    return np.array([len(truth), int((~valid).sum()), *[int(np.sum(abs(error) <= t)) for t in (50,100,200)]])


def compare(rows, plan, directory, target):
    if target.exists():
        raise FileExistsError('Immutable comparison output already exists')
    spec, folds = plan['empirical_features'], plan['crossfit_folds']
    contract = DepthContract(**plan['contract'])
    dims = dimensions(spec); nz = len(spec['target_edges_mm'])-1
    train = np.zeros((folds,)+dims+(nz,), np.int32)
    support = [np.zeros((folds,)+dims[:level], np.int32) for level in range(1,5)]
    source_records = []
    for row in rows:
        obs, truth, audit = load_observations(directory, row, plan)
        f = row['fold']; feature = observable_features(obs, contract, spec)
        zbin = np.minimum(np.searchsorted(spec['target_edges_mm'], truth, side='right')-1, nz-1)
        np.add.at(train[f], tuple(feature.T)+(zbin,), 1)
        for level in range(1,5):
            unique = np.unique(feature[:,:level], axis=0)
            np.add.at(support[level-1][f], tuple(unique.T), 1)
        source_records.append(dict(input=row, observations_sha256=audit['observations_sha256'],
            audit_sha256=file_sha256(directory/'results'/Path(row['root_file']).stem/'audit.json')))
    models, calibrations = [], []
    for fold in range(folds):
        # Whole sources, including all events/photons, are excluded at once.
        h = train[np.arange(folds) != fold].sum(axis=0)
        s = [v[np.arange(folds) != fold].sum(axis=0) for v in support]
        models.append(fit_tables(h, s, spec))
        calibrations.append(dict(fold=fold, training_sources=[r['root_sha256'] for r in rows if r['fold'] != fold],
            evaluation_sources=[r['root_sha256'] for r in rows if r['fold'] == fold],
            training_photons=int(h.sum()), supported_cells=[int(v.sum()) for v in models[-1][1]]))
    del train, support
    ze, re = np.array(plan['z_edges_mm']), np.array(plan['residual_edges_mm'])
    nm, ns, nt = len(METHODS), len(plan['sigmas_ps']), len(ze)-1
    hist = dict(response=np.zeros((nm,ns,nt,nt),np.int64), failures=np.zeros((nm,ns,nt),np.int64),
        denominator=np.zeros((nm,ns,nt),np.int64), within=np.zeros((nm,ns,nt,3),np.int64),
        residual=np.zeros((nm,ns,len(re)-1),np.int64), fallback=np.zeros((nm,ns,5),np.int64))
    blocks = np.zeros((len(rows),nm,ns,5), np.int64)
    for fi, row in enumerate(rows):
        obs, truth, audit = load_observations(directory,row,plan)
        tables, supported = models[row['fold']]
        for si, sigma in enumerate(plan['sigmas_ps']):
            noisy = smear_observations(obs, sigma_t_ns=sigma/1000., seed=plan['noise_seed'])
            feature = observable_features(noisy, contract, spec)
            reco, levels = predict_tables(feature,tables,supported)
            results = [contract.inverse(noisy.time_ns,1.,noisy.sensor_x_mm,noisy.sensor_y_mm),
                       *reco,np.full(len(truth),float(tables[0]))]
            for mi, z in enumerate(results):
                blocks[fi,mi,si] = score(hist,mi,si,truth,z,ze,re)
                if 1 <= mi <= 4:
                    hist['fallback'][mi,si] += np.bincount(levels[mi-1],minlength=5)
                elif mi == 5:
                    hist['fallback'][mi,si,0] += len(truth)
        # Exact unfitted-baseline replay on every source, at every timing point.
        with np.load(directory/'results'/Path(row['root_file']).stem/'histograms.npz',allow_pickle=False) as old:
            if not np.array_equal(blocks[fi,0,:,0],old['denominator'][0].sum(axis=-1)):
                raise ValueError('Old/new accepted population mismatch')
            if (not np.array_equal(blocks[fi,0,:,1],old['failures'][0].sum(axis=-1)) or
                    not np.array_equal(blocks[fi,0,:,3],old['within'][0,:,:,1].sum(axis=-1))):
                raise ValueError('Old/new axial inversion replay mismatch')
        print('SCORED',fi+1,len(rows),len(truth),flush=True)
    if (not np.array_equal(hist['response'].sum(axis=-1)+hist['failures'],hist['denominator']) or
            not np.array_equal(hist['residual'].sum(axis=-1),hist['response'].sum(axis=(-1,-2))) or
            np.any(hist['denominator'] != hist['denominator'][0,0])):
        raise ValueError('Matched response/failure/residual closure failed')
    rng = np.random.default_rng(plan['noise_seed'])
    weights = rng.multinomial(len(rows),np.ones(len(rows))/len(rows),size=500)
    boot = np.einsum('bf,fmsp->bmsp',weights.astype(float),blocks.astype(float))
    rates = 100*boot[:,:,:,1:]/boot[:,:,:,0,None]
    central = blocks.sum(axis=0)
    observed = 100*central[:,:,1:]/central[:,:,0,None]
    cases = {}
    for mi, name in enumerate(METHODS):
        for si, sigma in enumerate(plan['sigmas_ps']):
            cases[name+'_'+str(sigma)+'ps'] = dict(input=int(central[mi,si,0]), failure_pct=float(observed[mi,si,0]),
                within5cm_pct=float(observed[mi,si,1]),within10cm_pct=float(observed[mi,si,2]),within20cm_pct=float(observed[mi,si,3]),
                stat95=np.percentile(rates[:,mi,si],[2.5,97.5],axis=0).T.tolist(),
                paired_gain_vs_axial_pp=float(observed[mi,si,2]-observed[0,si,2]),
                paired_gain_stat95=np.percentile(rates[:,mi,si,2]-rates[:,0,si,2],[2.5,97.5]).tolist(),
                paired_gain_vs_time_pp=float(observed[mi,si,2]-observed[1,si,2]),
                paired_gain_vs_time_stat95=np.percentile(rates[:,mi,si,2]-rates[:,1,si,2],[2.5,97.5]).tolist(),
                used_feature_level_counts=hist['fallback'][mi,si].tolist())
    stage = target.with_name(target.name+'.staging_'+os.environ['SLURM_JOB_ID']); stage.mkdir()
    np.savez_compressed(stage/'comparison.npz',**hist,blocks=blocks,z_edges_mm=ze,residual_edges_mm=re)
    for fold,(tables,supported) in enumerate(models):
        np.savez_compressed(stage/('calibration_fold'+str(fold)+'.npz'),
            **{f'table{i}':v for i,v in enumerate(tables)},**{f'support{i}':v for i,v in enumerate(supported)})
    summary = dict(schema='pion-depth-corrections-v1',status='complete_crossfit_physics_provisional',
        plan_sha256=plan['_sha256'],events=sum(r['events'] for r in rows),files=len(rows),
        accepted_photons=int(central[0,0,0]),methods=list(METHODS),sigmas_ps=plan['sigmas_ps'],
        cases=cases,calibrations=calibrations,source_records=source_records,closure=True,
        code_revision=plan.get('analysis_code_revision',plan['code_revision']),
        extraction_code_revision=plan['code_revision'],
        aggregate_source_manifest_sha256=plan.get('aggregate_source_manifest_sha256'),optical_transport_validated=False,
        training='fixed zero-jitter conditional production-Z medians; not kappa/path truth; no fit on evaluation files',
        bootstrap='500 paired source-file draws; conditional on the fitted crossfit maps; no refit or transport systematic',
        physical_range='learned estimates are physical by construction, not by validated timing inversion',
        files_sha256={p.name:file_sha256(p) for p in stage.glob('*.npz')})
    (stage/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    stage.rename(target)
    print('COMPARISON COMPLETE',summary['events'],summary['accepted_photons'],target,flush=True)


def plot(directory, target, plan):
    from studyDepthReconstruction import _plot_axes
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, LinearSegmentedColormap
    s = json.loads((directory/'summary.json').read_text())
    if (s['plan_sha256'] != plan['_sha256'] or not s['closure'] or s['events'] != 5000 or
            s['files'] != len(plan['sources']) or s['methods'] != list(METHODS) or
            s['status'] != 'complete_crossfit_physics_provisional'):
        raise ValueError('Need a complete matched pion comparison')
    for name, sha in s['files_sha256'].items():
        if Path(name).name != name or file_sha256(directory/name) != sha:
            raise ValueError('Comparison output checksum mismatch')
    with np.load(directory/'comparison.npz',allow_pickle=False) as f:
        h = {k:f[k] for k in f.files}
    if not np.array_equal(h['response'].sum(axis=-1)+h['failures'],h['denominator']):
        raise ValueError('Plot response closure failed')
    target.mkdir(parents=True,exist_ok=False)
    # Initialize the same typography/header convention as the all-photon plots.
    temp,_ = _plot_axes('', '', (0,200));plt.close(temp)
    cmap = LinearSegmentedColormap.from_list('bird',['#352a87','#1062dc','#108cd2','#17acb9','#60bc8b','#babc62','#f7c638','#f9f90f'])
    cmap.set_bad('white'); ze = h['z_edges_mm']/10
    si = plan['sigmas_ps'].index(50)
    fig, axes = plt.subplots(2,3,figsize=(14.5,10),layout=None)
    fig.subplots_adjust(left=.065,right=.90,bottom=.08,top=.855,wspace=.27,hspace=.32)
    fig.text(.065,.962,'FLASH',weight='bold',fontsize=20)
    fig.text(.145,.962,'Simulation Preliminary',style='italic',fontsize=16)
    fig.text(.94,.962,'100 GeV pion; 50 ps/photon; exact entrance clock',ha='right',fontsize=12)
    fig.text(.065,.917,'Same 5,000 events and accepted photons; five-fold source-file holdout',fontsize=12)
    fig.text(.065,.888,'Learned panels predict an in-range depth median; this is not a validated optical-path inversion.',fontsize=11)
    for mi, ax in enumerate(axes.flat):
        p = np.divide(100*h['response'][mi,si],h['denominator'][mi,si,:,None],
                      out=np.zeros_like(h['response'][mi,si],float),where=h['denominator'][mi,si,:,None]>0)
        im = ax.pcolormesh(ze,ze,np.ma.masked_equal(p.T,0),cmap=cmap,norm=LogNorm(.01,100),rasterized=True)
        ax.plot([0,200],[0,200],'--',color='#ce007a',lw=1.1)
        ax.set(xlim=(0,200),ylim=(0,200),xlabel='True photon-production Z [cm]',ylabel='Reconstructed Z [cm]')
        ax.set_box_aspect(1); ax.minorticks_on(); ax.tick_params(labelsize=10)
        case = s['cases'][METHODS[mi]+'_50ps']
        note = f"range failures: {case['failure_pct']:.1f}%" if mi==0 else 'in-range median'
        ax.set_title(LABELS[mi]+'\n'+f"Within 10 cm: {case['within10cm_pct']:.1f}%; "+note,fontsize=11,pad=8)
    cax = fig.add_axes([.925,.18,.014,.61]);fig.colorbar(im,cax=cax).set_label('All inputs per true-Z bin [% / 1 cm reco bin]',weight='bold')
    def save(fig,name):
        for ext in ('pdf','png'): fig.savefig(target/(name+'.'+ext),dpi=180,facecolor='white')
        plt.close(fig)
    save(fig,'Z_methods_side_by_side_50ps')
    fig,ax = _plot_axes('100 GeV pion\n5,000 events; five-fold file holdout\nExact entrance time',
                       'Assumed single-photon timing sigma [ps]',(-10,510))
    for mi,name in enumerate(METHODS):
        cases = [s['cases'][name+'_'+str(t)+'ps'] for t in plan['sigmas_ps']]
        limits = np.array([c['stat95'][2] for c in cases])
        ax.fill_between(plan['sigmas_ps'],limits[:,0],limits[:,1],color=COLORS[mi],alpha=.13)
        ax.plot(plan['sigmas_ps'],[c['within10cm_pct'] for c in cases],marker='o',color=COLORS[mi],
                label=LABELS[mi],ls='--' if mi==5 else '-',ms=4)
    ax.set(ylim=(0,100),ylabel='Accepted photons within 10 cm of true Z [%]')
    ax.legend(fontsize=9,loc='upper right');save(fig,'accuracy_vs_timing_all_methods')
    fig,ax = _plot_axes('100 GeV pion\nRelative to calibrated time only\nPaired 95% file-bootstrap intervals',
                       'Assumed single-photon timing sigma [ps]',(-10,510))
    ax.axhline(0,color='#666666',ls='--',lw=1.2)
    for mi in (2,3,4):
        cases=[s['cases'][METHODS[mi]+'_'+str(t)+'ps'] for t in plan['sigmas_ps']]
        gain=np.array([c['paired_gain_vs_time_pp'] for c in cases])
        limits=np.array([c['paired_gain_vs_time_stat95'] for c in cases])
        ax.fill_between(plan['sigmas_ps'],limits[:,0],limits[:,1],color=COLORS[mi],alpha=.18)
        ax.plot(plan['sigmas_ps'],gain,'o-',color=COLORS[mi],label=LABELS[mi],ms=4)
    ax.set_ylabel('Change in within-10-cm fraction [percentage points]',fontsize=12)
    ax.legend(fontsize=10);save(fig,'added_information_vs_time_only')
    fig,ax = _plot_axes('100 GeV pion; 50 ps\nSame accepted photons; failures included',
                       'Allowed absolute Z error [cm]',(0,200))
    edges=h['residual_edges_mm']; centers=(edges[:-1]+edges[1:])/2
    for mi,name in enumerate(METHODS):
        errors=np.sort(np.unique(abs(centers)))
        counts=h['residual'][mi,si]; denom=h['denominator'][mi,si].sum()
        cumulative=np.array([counts[abs(centers)<=e].sum() for e in errors])
        ax.plot(errors/10,100*cumulative/denom,color=COLORS[mi],label=LABELS[mi],ls='--' if mi==5 else '-')
    ax.set(ylim=(0,100),ylabel='Fraction of all accepted photons [%]')
    ax.legend(loc='lower right',fontsize=9);save(fig,'inclusive_error_comparison_50ps')
    with (target/'comparison_metrics.csv').open('x',newline='') as f:
        writer=csv.writer(f);writer.writerow(['method','sigma_ps','accepted_photons','within10cm_pct','stat95_low','stat95_high',
            'failed_pct','gain_vs_basic_pp','gain_stat95_low','gain_stat95_high','gain_vs_time_pp','gain_vs_time95_low','gain_vs_time95_high'])
        for name in METHODS:
            for t in plan['sigmas_ps']:
                c=s['cases'][name+'_'+str(t)+'ps'];writer.writerow([name,t,c['input'],c['within10cm_pct'],*c['stat95'][2],
                    c['failure_pct'],c['paired_gain_vs_axial_pp'],*c['paired_gain_stat95'],c['paired_gain_vs_time_pp'],*c['paired_gain_vs_time_stat95']])
    (target/'plot_manifest.json').write_text(json.dumps(dict(source_summary_sha256=file_sha256(directory/'summary.json'),
        plotter_sha256=file_sha256(__file__),
        plan_sha256=plan['_sha256'],files_sha256={p.name:file_sha256(p) for p in target.iterdir()}),indent=2)+'\n')


# Primary method study. The older direct-Z lookup above is an empirical control,
# NOT corrections 1-3. These entries call the actual dsipm.depth reconstruction.
TIMING_METHODS=('axial','global_factor','radius','radius_delay','neighbour','neighbour_guided',
                'empirical_time','empirical_radius','empirical_delay','empirical_neighbour','depth_prior')
TIMING_LABELS=('Basic axial inversion','Global timing-factor calibration','Correction 1: radius',
    'Correction 2: radius + delay','Correction 3: neighbour prior','Correction 3 + geometric heuristic',
    'Empirical Z lookup: time','Empirical Z lookup: + radius','Empirical Z lookup: + delay',
    'Empirical Z lookup: + neighbour time','No-timing depth prior')


def method_histograms(nm,ns,nt,nr):
    return dict(response=np.zeros((nm,ns,nt,nt),np.int64),failures=np.zeros((nm,ns,nt),np.int64),
        denominator=np.zeros((nm,ns,nt),np.int64),within=np.zeros((nm,ns,nt,3),np.int64),
        residual=np.zeros((nm,ns,nr),np.int64),fallback=np.zeros((nm,ns,5),np.int64),
        status_counts=np.zeros((nm,ns,3),np.int64))


def method_task(rows,plan,source_plan,target,*,fold,sigma):
    """One training-file fold and MATCHED timing sigma. Cached observations only."""
    if target.exists():raise FileExistsError('Immutable method task already exists')
    contract=DepthContract(**plan['contract']);spec=plan['empirical_features']
    directory=Path(plan['cache_directory']);train_rows=[r for r in rows if r['fold']!=fold]
    test_rows=[r for r in rows if r['fold']==fold]
    if not train_rows or not test_rows:raise ValueError('Empty train/evaluation fold')
    train_hist=np.zeros(dimensions(spec)+(len(spec['target_edges_mm'])-1,),np.int32)
    support=[np.zeros(dimensions(spec)[:j],np.int32) for j in range(1,5)]
    inputs={}
    def read(row):
        obs,truth,audit=load_observations(directory,row,source_plan)
        inputs[row['root_sha256']]=dict(observations_sha256=audit['observations_sha256'],
            audit_sha256=file_sha256(directory/'results'/Path(row['root_file']).stem/'audit.json'))
        return smear_observations(obs,sigma_t_ns=sigma/1000,seed=plan['noise_seed']),truth
    def training_batches():
        for i,row in enumerate(train_rows):
            obs,truth=read(row);feature=observable_features(obs,contract,spec)
            zbin=np.minimum(np.searchsorted(spec['target_edges_mm'],truth,side='right')-1,len(spec['target_edges_mm'])-2)
            np.add.at(train_hist,tuple(feature.T)+(zbin,),1)
            for j in range(1,5):
                unique=np.unique(feature[:,:j],axis=0);np.add.at(support[j-1],tuple(unique.T),1)
            if i%40==0:print('TRAIN',fold,sigma,i+1,len(train_rows),flush=True)
            yield obs,truth
    calibration=fit_timing_calibration_batches(training_batches(),contract,
        source_sha256=digest(dict(cache_plan_sha256=source_plan['_sha256'],sigma_ps=sigma,
            noise_seed=plan['noise_seed'],training_roots=[r['root_sha256'] for r in train_rows])),
        truth_z_validation=plan['truth_z_validation'],training_sigma_ns=sigma/1000,
        min_count=plan['factor_min_count'])
    tables,supported=fit_tables(train_hist,support,spec);del train_hist,support
    ze,re=np.array(plan['z_edges_mm']),np.array(plan['residual_edges_mm'])
    h=method_histograms(len(TIMING_METHODS),1,len(ze)-1,len(re)-1)
    blocks=np.zeros((len(test_rows),len(TIMING_METHODS),1,5),np.int64);options={}
    for fi,row in enumerate(test_rows):
        obs,truth=read(row)
        outputs=dict(axial=contract.inverse(obs.time_ns,1.,obs.sensor_x_mm,obs.sensor_y_mm),
            global_factor=contract.inverse(obs.time_ns,calibration.payload['global_median'],obs.sensor_x_mm,obs.sensor_y_mm))
        for name in ('radius','radius_delay','neighbour','neighbour_guided'):
            result=reconstruct(obs,calibration,contract,method='neighbour' if name=='neighbour_guided' else name,
                sigma_t_ns=sigma/1000,window=plan['window_rods'],clean_radius=plan['clean_radius'],
                geometric_guidance=name=='neighbour_guided')
            outputs[name]=result.z_mm;options[name]=result.options
            h['status_counts'][TIMING_METHODS.index(name),0]+=np.bincount(result.status,minlength=3)
        empirical,levels=predict_tables(observable_features(obs,contract,spec),tables,supported)
        for j,name in enumerate(TIMING_METHODS[6:10]):
            outputs[name]=empirical[j];h['fallback'][TIMING_METHODS.index(name),0]+=np.bincount(levels[j],minlength=5)
        outputs['depth_prior']=np.full(len(truth),float(tables[0]));h['fallback'][-1,0,0]+=len(truth)
        for mi,name in enumerate(TIMING_METHODS):blocks[fi,mi,0]=score(h,mi,0,truth,outputs[name],ze,re)
        with np.load(directory/'results'/Path(row['root_file']).stem/'histograms.npz',allow_pickle=False) as old:
            si=source_plan['sigmas_ps'].index(sigma)
            if (blocks[fi,0,0,0]!=old['denominator'][0,si].sum() or
                    blocks[fi,0,0,1]!=old['failures'][0,si].sum() or
                    blocks[fi,0,0,3]!=old['within'][0,si,:,1].sum()):
                raise ValueError('Matched-jitter baseline replay mismatch')
    if (not np.array_equal(h['response'].sum(axis=-1)+h['failures'],h['denominator']) or
            np.any(h['denominator']!=h['denominator'][0,0])):raise ValueError('Method task photon closure failed')
    stage=target.with_name(target.name+'.staging_'+os.environ['SLURM_JOB_ID']);stage.mkdir(parents=True)
    calibration.save(stage/'calibration.json')
    np.savez_compressed(stage/'histograms.npz',**h,blocks=blocks)
    np.savez_compressed(stage/'empirical_control.npz',**{f'table{i}':v for i,v in enumerate(tables)},
                        **{f'support{i}':v for i,v in enumerate(supported)})
    report=dict(status='complete_method_task_physics_provisional',fold=fold,sigma_ps=sigma,
        plan_sha256=plan['_sha256'],code_revision=plan['code_revision'],methods=list(TIMING_METHODS),
        evaluation_rows=test_rows,training_roots=[r['root_sha256'] for r in train_rows],
        reconstruction_options=options,calibration_sha256=calibration.sha256,
        training_sigma_ns=calibration.payload['training_sigma_ns'],source_inputs=inputs,
        factor_kind='timing_effective',optical_transport_validated=False,
        files_sha256={p.name:file_sha256(p) for p in stage.iterdir()})
    (stage/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');stage.rename(target)
    print('METHOD TASK COMPLETE',fold,sigma,len(test_rows),int(blocks[:,0,0,0].sum()),flush=True)
    return report


def inclusive_histogram_quantile(counts,edges,total,probability):
    """Failure-inclusive tolerance; upper bin edge, <= bin-width approximation."""
    rank=int(np.ceil(probability*total))
    if not total or rank>int(counts.sum()):return None
    absolute_upper=np.maximum(abs(edges[:-1]),abs(edges[1:]));order=np.argsort(absolute_upper)
    at=int(np.searchsorted(np.cumsum(counts[order]),rank,side='left'))
    return float(absolute_upper[order[at]])


def merge_methods(rows,plan,directory,target):
    if target.exists():raise FileExistsError('Immutable method comparison already exists')
    nm,ns,nt,nr=len(TIMING_METHODS),len(plan['sigmas_ps']),len(plan['z_edges_mm'])-1,len(plan['residual_edges_mm'])-1
    h=method_histograms(nm,ns,nt,nr);blocks=np.zeros((len(rows),nm,ns,5),np.int64)
    done=np.zeros((ns,len(rows)),bool);lookup={r['root_sha256']:i for i,r in enumerate(rows)}
    task_records=[];source_inputs=None
    for fold in range(plan['crossfit_folds']):
        for si,sigma in enumerate(plan['sigmas_ps']):
            d=directory/'method_results'/f'fold{fold}_{sigma}ps';s=json.loads((d/'summary.json').read_text())
            if (s['plan_sha256']!=plan['_sha256'] or s['code_revision']!=plan['code_revision'] or
                    s['fold']!=fold or s['sigma_ps']!=sigma or s['training_sigma_ns']!=sigma/1000 or
                    s['methods']!=list(TIMING_METHODS) or s['status']!='complete_method_task_physics_provisional'):
                raise ValueError('Method task provenance mismatch')
            for name,sha in s['files_sha256'].items():
                if Path(name).name!=name or file_sha256(d/name)!=sha:raise ValueError('Method task checksum mismatch')
            cal=Calibration.load(d/'calibration.json')
            if (cal.sha256!=s['calibration_sha256'] or cal.payload['training_sigma_ns']!=sigma/1000 or
                    cal.payload.get('training_entrance_sigma_ns',0)!=0 or
                    cal.payload.get('factor_kind')!='timing_effective'):
                raise ValueError('Calibration jitter/hash mismatch')
            for name in ('radius','radius_delay','neighbour','neighbour_guided'):
                option=s['reconstruction_options'][name]
                if (option['method']!=('neighbour' if name=='neighbour_guided' else name) or
                        option['sigma_t_ns']!=sigma/1000 or option['entrance_sigma_ns']!=0 or
                        option['factor_kind']!='timing_effective' or option['allow_training'] or
                        option['calibration_mode']!='jitter_matched' or
                        option['geometric_guidance']!=(name=='neighbour_guided')):
                    raise ValueError('Actual-method/jitter/guidance mismatch')
            expected=[r for r in rows if r['fold']==fold]
            if s['evaluation_rows']!=expected or set(s['training_roots'])!={r['root_sha256'] for r in rows if r['fold']!=fold}:
                raise ValueError('Source-file train/test overlap or missing source')
            training_events=set(cal.payload['training_events'])
            if any(row['root_sha256']+':'+str(i) in training_events for row in expected for i in range(row['events'])):
                raise ValueError('Calibration/evaluation event overlap')
            if source_inputs is None:source_inputs=s['source_inputs']
            if s['source_inputs']!=source_inputs or set(source_inputs)!=set(lookup):raise ValueError('Cached input drift between tasks')
            with np.load(d/'histograms.npz',allow_pickle=False) as f:
                for key in h:h[key][:,si]+=f[key][:,0]
                for j,row in enumerate(expected):
                    at=lookup[row['root_sha256']]
                    if done[si,at]:raise ValueError('Duplicate evaluation source')
                    blocks[at,:,si]=f['blocks'][j,:,0];done[si,at]=True
            task_records.append(dict(fold=fold,sigma_ps=sigma,summary_sha256=file_sha256(d/'summary.json'),
                files_sha256=s['files_sha256'],reconstruction_options=s['reconstruction_options'],
                undefined_calibration_factors=cal.payload['undefined_factor_photons'],
                factors_below_one=cal.payload['effective_factors_below_one']))
    if (not done.all() or np.any(h['denominator']!=h['denominator'][0,0]) or
            not np.array_equal(h['response'].sum(axis=-1)+h['failures'],h['denominator']) or
            not np.array_equal(h['residual'].sum(axis=-1),h['response'].sum(axis=(-1,-2))) or
            not np.array_equal(blocks.sum(axis=0)[:,:,0],h['denominator'].sum(axis=-1)) or
            not np.array_equal(blocks.sum(axis=0)[:,:,1],h['failures'].sum(axis=-1)) or
            not np.array_equal(blocks.sum(axis=0)[:,:,2:],h['within'].sum(axis=2)) or
            not np.array_equal(h['status_counts'][2:6].sum(axis=-1),h['denominator'][2:6].sum(axis=-1)) or
            not np.array_equal(h['status_counts'][2:6,:,-1],h['failures'][2:6].sum(axis=-1))):
        raise ValueError('Full method comparison coverage/count closure failed')
    weights=np.random.default_rng(plan['noise_seed']).multinomial(len(rows),np.ones(len(rows))/len(rows),size=500)
    boot=np.einsum('bf,fmsp->bmsp',weights.astype(float),blocks.astype(float));rates=100*boot[:,:,:,1:]/boot[:,:,:,0,None]
    totals=blocks.sum(axis=0);observed=100*totals[:,:,1:]/totals[:,:,0,None];cases={}
    for mi,name in enumerate(TIMING_METHODS):
        for si,sigma in enumerate(plan['sigmas_ps']):
            cases[f'{name}_{sigma}ps']=dict(input=int(totals[mi,si,0]),failure_pct=float(observed[mi,si,0]),
                within5cm_pct=float(observed[mi,si,1]),within10cm_pct=float(observed[mi,si,2]),within20cm_pct=float(observed[mi,si,3]),
                stat95=np.percentile(rates[:,mi,si],[2.5,97.5],axis=0).T.tolist(),
                inclusive_abs_error_q68_mm=inclusive_histogram_quantile(h['residual'][mi,si],np.array(plan['residual_edges_mm']),int(totals[mi,si,0]),.68),
                comparisons={ref:dict(gain_pp=float(observed[mi,si,2]-observed[r,si,2]),
                    stat95=np.percentile(rates[:,mi,si,2]-rates[:,r,si,2],[2.5,97.5]).tolist())
                    for r,ref in enumerate(TIMING_METHODS)})
    stage=target.with_name(target.name+'.staging_'+os.environ['SLURM_JOB_ID']);stage.mkdir()
    np.savez_compressed(stage/'comparison.npz',**h,blocks=blocks,z_edges_mm=plan['z_edges_mm'],residual_edges_mm=plan['residual_edges_mm'])
    report=dict(schema='timing-factor-method-comparison-v1',status='complete_methods_physics_provisional',
        events=sum(r['events'] for r in rows),files=len(rows),accepted_photons=int(totals[0,0,0]),closure=True,
        plan_sha256=plan['_sha256'],code_revision=plan['code_revision'],methods=list(TIMING_METHODS),sigmas_ps=plan['sigmas_ps'],
        cases=cases,tasks=task_records,source_inputs=source_inputs,optical_transport_validated=False,
        factor_definition='Timing-effective calibration parameter; NOT geometric path truth; includes emission/front, group-speed and jitter effects',
        bootstrap='500 paired source-file draws, conditional on fitted crossfit maps; excludes refit and transport/model uncertainty',
        inclusive_q68='Absolute-error tolerance containing 68% of ALL photons, failures at infinity; null when unreachable; <=2 mm histogram approximation',
        files_sha256={'comparison.npz':file_sha256(stage/'comparison.npz')})
    (stage/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');stage.rename(target)
    print('METHOD COMPARISON COMPLETE',report['events'],report['accepted_photons'],target,flush=True)
    return report


def plot_methods(directory,target,plan):
    """Separate actual inverse methods from direct-Z lookup controls."""
    from studyDepthReconstruction import _plot_axes
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm,LinearSegmentedColormap
    s=json.loads((directory/'summary.json').read_text())
    if (s['plan_sha256']!=plan['_sha256'] or s['events']!=5000 or s['files']!=200 or not s['closure'] or
            s['status']!='complete_methods_physics_provisional' or s['methods']!=list(TIMING_METHODS)):
        raise ValueError('Need the complete 5000-event jitter-matched method comparison')
    for name,sha in s['files_sha256'].items():
        if Path(name).name!=name or file_sha256(directory/name)!=sha:raise ValueError('Method plot data checksum mismatch')
    with np.load(directory/'comparison.npz',allow_pickle=False) as f:h={k:f[k] for k in f.files}
    if (not np.array_equal(h['response'].sum(axis=-1)+h['failures'],h['denominator']) or
            np.any(h['denominator']!=h['denominator'][0,0])):raise ValueError('Method plot closure mismatch')
    target.mkdir(parents=True,exist_ok=False);temp,_=_plot_axes('','',(0,200));plt.close(temp)
    colors=('#555555','#888888','#2171b5','#dd851b','#148060','#a51c7b','#399ac8','#17a398','#c7893b','#a51c7b','#9b895e')
    short=('Basic axial','Global factor','1: Radius','2: Radius + delay','3: Neighbour prior','3: Geometric-guided variant',
           'Z lookup: time','Z lookup: radius','Z lookup: delay','Z lookup: neighbour','No-timing prior')
    cmap=LinearSegmentedColormap.from_list('bird',['#352a87','#1062dc','#108cd2','#17acb9','#60bc8b','#babc62','#f7c638','#f9f90f']);cmap.set_bad('white')
    ze=h['z_edges_mm']/10;si=plan['sigmas_ps'].index(50)
    def save(fig,name):
        for ext in ('pdf','png'):fig.savefig(target/(name+'.'+ext),dpi=180,facecolor='white')
        plt.close(fig)
    for ids,title,name in [(range(6),'Actual timing inversion and corrections 1-3','actual_methods_Z_truth_vs_reconstructed_50ps'),
                           (range(6,11),'Empirical Z-lookup controls - NOT corrections 1-3','empirical_controls_Z_truth_vs_reconstructed_50ps')]:
        fig,axes=plt.subplots(2,3,figsize=(14.5,10));fig.subplots_adjust(left=.065,right=.90,bottom=.08,top=.815,wspace=.27,hspace=.36)
        fig.text(.065,.962,'FLASH',weight='bold',fontsize=20);fig.text(.145,.962,'Simulation Preliminary',style='italic',fontsize=16)
        fig.text(.94,.962,'100 GeV pion; 50 ps/photon; exact entrance clock',ha='right',fontsize=12)
        fig.text(.065,.917,title+'; same 5,000 events',fontsize=12)
        fig.text(.065,.888,'Jitter-matched calibration; whole-file holdout; timing-effective factor is not geometric path truth.',fontsize=11)
        for ax,mi in zip(axes.flat,ids):
            p=np.divide(100*h['response'][mi,si],h['denominator'][mi,si,:,None],out=np.zeros_like(h['response'][mi,si],float),where=h['denominator'][mi,si,:,None]>0)
            im=ax.pcolormesh(ze,ze,np.ma.masked_equal(p.T,0),cmap=cmap,norm=LogNorm(.01,100),rasterized=True)
            ax.plot([0,200],[0,200],'--',color='#ce007a',lw=1.1)
            ax.set(xlim=(0,200),ylim=(0,200),xlabel='True photon-production Z [cm]',ylabel='Reconstructed Z [cm]')
            ax.set_box_aspect(1);ax.minorticks_on();ax.tick_params(labelsize=10)
            c=s['cases'][TIMING_METHODS[mi]+'_50ps'];note=f"failed: {c['failure_pct']:.2f}%" if mi<6 else 'in-range median'
            ax.set_title(short[mi]+'\n'+f"Within 10 cm: {c['within10cm_pct']:.1f}%; "+note,fontsize=11,pad=8)
        if len(ids)<6:axes.flat[-1].set_visible(False)
        cax=fig.add_axes([.925,.18,.014,.61]);fig.colorbar(im,cax=cax).set_label('All inputs per true-Z bin [% / 1 cm reco bin]',weight='bold')
        save(fig,name)
    fig,ax=_plot_axes('100 GeV pion\nMatched training at every timing sigma\nExact entrance time','Assumed single-photon timing sigma [ps]',(-10,510))
    for mi in (0,1,2,3,4,5,10):
        cs=[s['cases'][TIMING_METHODS[mi]+f'_{t}ps'] for t in plan['sigmas_ps']];limits=np.array([c['stat95'][2] for c in cs])
        ax.fill_between(plan['sigmas_ps'],limits[:,0],limits[:,1],color=colors[mi],alpha=.13)
        ax.plot(plan['sigmas_ps'],[c['within10cm_pct'] for c in cs],marker='o',ls='--' if mi in (1,5,10) else '-',color=colors[mi],label=short[mi],ms=4)
    ax.set_xticks(plan['sigmas_ps'])
    ax.set(ylim=(0,100),ylabel='All accepted photons within 10 cm of true Z [%]');ax.legend(fontsize=8.5,loc='upper right')
    save(fig,'actual_methods_accuracy_vs_timing')
    fig,ax=_plot_axes('100 GeV pion\nEach added step versus the previous step\nPaired 95% file-bootstrap intervals','Assumed single-photon timing sigma [ps]',(-10,510))
    ax.axhline(0,color='#666666',ls='--',lw=1.2)
    for mi,ref,label in [(2,'global_factor','Add radius'),(3,'radius','Add same-rod delay'),(4,'radius_delay','Add neighbour prior'),(5,'neighbour','Add geometric heuristic')]:
        cs=[s['cases'][TIMING_METHODS[mi]+f'_{t}ps']['comparisons'][ref] for t in plan['sigmas_ps']];limits=np.array([c['stat95'] for c in cs])
        ax.fill_between(plan['sigmas_ps'],limits[:,0],limits[:,1],color=colors[mi],alpha=.18)
        ax.plot(plan['sigmas_ps'],[c['gain_pp'] for c in cs],'o-',color=colors[mi],label=label,ms=4)
    ax.set_xticks(plan['sigmas_ps'])
    ax.set_ylabel('Change in within-10-cm fraction [percentage points]',fontsize=12);ax.legend(fontsize=9)
    save(fig,'actual_methods_incremental_gain')
    fig,ax=_plot_axes('100 GeV pion; 50 ps\nFailures stay in the denominator','Allowed absolute Z error [cm]',(0,200))
    upper=np.maximum(abs(h['residual_edges_mm'][:-1]),abs(h['residual_edges_mm'][1:]));order=np.argsort(upper)
    for mi in (0,1,2,3,4,5,10):
        ax.plot(upper[order]/10,100*np.cumsum(h['residual'][mi,si][order])/h['denominator'][mi,si].sum(),color=colors[mi],label=short[mi],ls='--' if mi in (1,5,10) else '-')
    ax.set(ylim=(0,100),ylabel='Fraction of all accepted photons [%]');ax.legend(fontsize=8.5,loc='lower right')
    save(fig,'actual_methods_failure_inclusive_error_50ps')
    with (target/'all_method_metrics.csv').open('x',newline='') as f:
        w=csv.writer(f);w.writerow(['method','sigma_ps','accepted_photons','failure_pct','within5cm_pct','within10cm_pct','within10cm_stat95_low','within10cm_stat95_high','within20cm_pct','all_photon_q68_abs_error_cm'])
        for name in TIMING_METHODS:
            for t in plan['sigmas_ps']:
                c=s['cases'][name+f'_{t}ps'];q=c['inclusive_abs_error_q68_mm']
                w.writerow([name,t,c['input'],c['failure_pct'],c['within5cm_pct'],c['within10cm_pct'],*c['stat95'][2],c['within20cm_pct'],None if q is None else q/10])
    (target/'plot_manifest.json').write_text(json.dumps(dict(source_summary_sha256=file_sha256(directory/'summary.json'),
        plotter_sha256=file_sha256(__file__),plan_sha256=plan['_sha256'],files_sha256={p.name:file_sha256(p) for p in target.iterdir()}),indent=2)+'\n')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan',required=True,type=Path);ap.add_argument('--directory',required=True,type=Path)
    op=ap.add_mutually_exclusive_group(required=True)
    op.add_argument('--extract',action='store_true');op.add_argument('--compare',type=Path);op.add_argument('--plot',type=Path)
    op.add_argument('--method-task',type=int,nargs='?',const=-1)
    op.add_argument('--method-merge',type=Path)
    op.add_argument('--method-plot',type=Path)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());plan['_sha256']=file_sha256(args.plan)
    if args.method_plot:plot_methods(args.directory,args.method_plot,plan);return
    if args.method_task is not None or args.method_merge:
        if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Method fitting/scoring requires Slurm')
        if plan['schema']!='timing-factor-method-plan-v1':raise ValueError('Wrong method comparison plan')
        for name,sha in plan['code_sha256'].items():
            if file_sha256(Path(__file__).parent/name)!=sha:raise ValueError('Method code checksum mismatch: '+name)
        if file_sha256(plan['cache_plan'])!=plan['cache_plan_sha256']:raise ValueError('Cached extraction plan mismatch')
        source_plan=json.loads(Path(plan['cache_plan']).read_text());source_plan['_sha256']=plan['cache_plan_sha256']
        rows=plan['sources']
        if (rows!=source_plan['sources'] or sum(r['events'] for r in rows)!=5000 or
                len(rows)!=200 or len({r['root_sha256'] for r in rows})!=200 or
                len({r['offset']+e for r in rows for e in range(r['events'])})!=5000 or
                set(r['fold'] for r in rows)!=set(range(plan['crossfit_folds'])) or
                plan['contract']!=source_plan['contract'] or
                DepthContract(**plan['contract']).sha256!=plan['contract_sha256'] or
                plan['noise_seed']!=source_plan['noise_seed'] or
                plan['entrance_sigma_ps']!=0 or plan['sigmas_ps']!=source_plan['sigmas_ps']):
            raise ValueError('Method population/geometry/timing mismatch')
        if args.method_merge:merge_methods(rows,plan,args.directory,args.method_merge)
        else:
            task=args.method_task if args.method_task>=0 else int(os.environ.get('SLURM_ARRAY_TASK_ID','-1'))
            ns=len(plan['sigmas_ps'])
            if task<0 or task>=plan['crossfit_folds']*ns:raise ValueError('Invalid method task index')
            fold,sigma=task//ns,plan['sigmas_ps'][task%ns]
            method_task(rows,plan,source_plan,args.directory/'method_results'/f'fold{fold}_{sigma}ps',fold=fold,sigma=sigma)
        return
    if args.plot:
        plot(args.directory,args.plot,plan);return
    if not os.environ.get('SLURM_JOB_ID'): raise RuntimeError('Existing-photon analysis requires Slurm')
    for name,sha in plan['code_sha256'].items():
        if file_sha256((args.plan.parent if args.compare else Path(__file__).parent)/name)!=sha:
            raise ValueError('Extraction code hash mismatch: '+name)
    if args.compare and Path(__file__).parent.resolve()!=args.plan.parent.resolve():
        path=Path(__file__).parent/'aggregate_source.json';source=json.loads(path.read_text())
        if source['plan_sha256']!=plan['_sha256']: raise ValueError('Aggregate source plan mismatch')
        for name,sha in source['code_sha256'].items():
            if file_sha256(Path(__file__).parent/name)!=sha: raise ValueError('Aggregate source hash mismatch: '+name)
        plan['analysis_code_revision']=source['code_revision']
        plan['aggregate_source_manifest_sha256']=file_sha256(path)
    if file_sha256(plan['input_manifest'])!=plan['input_manifest_sha256']: raise ValueError('Input manifest hash mismatch')
    all_rows=json.loads(Path(plan['input_manifest']).read_text())
    rows=plan['sources']
    expected=[r for r in all_rows if r['sample']=='piplus_100gev']
    if ([{k:v for k,v in r.items() if k!='fold'} for r in rows] != expected or sum(r['events'] for r in rows)!=5000 or
            len({r['root_sha256'] for r in rows})!=len(rows) or
            len({r['offset']+e for r in rows for e in range(r['events'])})!=5000 or
            set(r['fold'] for r in rows)!=set(range(plan['crossfit_folds']))):
        raise ValueError('Incomplete/duplicate or altered pion manifest/folds')
    if file_sha256(Path(plan['reference_release'])/'plan.json')!=plan['reference_plan_sha256']:
        raise ValueError('Corrected reference-plan checksum mismatch')
    if args.compare:
        compare(rows,plan,args.directory,args.compare)
    else:
        task=int(os.environ['SLURM_ARRAY_TASK_ID']);chunk=plan['files_per_task'];start=task*chunk
        if start<0 or start>=len(rows): raise ValueError('Invalid extraction array index')
        for row in rows[start:start+chunk]:
            # The common reference table uses the original manifest, not its fold tag.
            reduce_one({k:v for k,v in row.items() if k!='fold'},plan,args.directory)


if __name__ == '__main__':
    main()
