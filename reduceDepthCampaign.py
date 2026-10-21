#!/usr/bin/env python3
"""Slurm-only all-event/all-photon baseline audit of existing production ROOTs.

No channel/event selection, fitting, raw export, simulation or production writes.
The historical fixed-path baseline excludes its 20 old calibration events from
its response scores, but processes and accounts for every input event/photon.
The unfitted axial baseline scores every event. Learned path models stay gated
on independently validated optical-path truth.
"""
import argparse
import csv
from dataclasses import replace
import gzip
import json
import os
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np

from dsipm.depth import DepthContract
from dsipm.geometry import dreamsim_front_contract
from dsipm.io import file_sha256, iter_corrected_root
from dsipm.sensors import (FullDetectorScenario, digitize_full_detector_event,
                           make_observations, smear_observations)


def make_histograms(plan):
    ze, re = np.asarray(plan['z_edges_mm'], float), np.asarray(plan['residual_edges_mm'], float)
    if (any(a.ndim != 1 or len(a)<3 or not np.isfinite(a).all() or np.any(np.diff(a)<=0) for a in (ze,re))
            or tuple(ze[[0,-1]]) != (0.,2000.) or re[0]>-2000. or re[-1]<2000.):
        raise ValueError('Histogram edges must cover the complete physical detector and residuals')
    sigmas = np.asarray(plan['sigmas_ps'], float)
    models = plan['baselines']
    if (sigmas.ndim != 1 or not len(sigmas) or not np.isfinite(sigmas).all() or np.any(sigmas<0)
            or len(set(sigmas))!=len(sigmas) or not models
            or len({m['name'] for m in models})!=len(models)
            or any(not np.isfinite(m['n_eff']) or m['n_eff']<=1. for m in models)):
        raise ValueError('Invalid timing scan or fixed-path baselines')
    nt, nr, nm, ns = len(ze)-1, len(re)-1, len(models), len(sigmas)
    return dict(response=np.zeros((nm,ns,nt,nt), np.int64),
        failures=np.zeros((nm,ns,nt), np.int64), denominator=np.zeros((nm,ns,nt), np.int64),
        residual=np.zeros((nm,ns,nr), np.int64), depth_residual=np.zeros((nm,ns,nt,nr), np.int64),
        within=np.zeros((nm,ns,nt,3), np.int64), stage_truth=np.zeros((6,nt), np.int64))


def accumulate(arrays, model, sigma, truth, reconstructed, zedges, redges):
    """Physical outputs and failures sum to ALL accepted inputs in each truth bin."""
    truth, reco = np.asarray(truth), np.asarray(reconstructed)
    if truth.ndim!=1 or reco.shape!=truth.shape:
        raise ValueError('Misaligned photon truth and reconstructed depth')
    valid = np.isfinite(reco) & (reco >= zedges[0]) & (reco <= zedges[-1])
    if not np.isfinite(truth).all() or np.any((truth < zedges[0]) | (truth > zedges[-1])):
        raise ValueError("Invalid true-Z target")
    row = (model, sigma)
    arrays['response'][row] += np.histogram2d(truth[valid], reco[valid], (zedges, zedges))[0].astype(np.int64)
    arrays['failures'][row] += np.histogram(truth[~valid], zedges)[0]
    arrays['denominator'][row] += np.histogram(truth, zedges)[0]
    error = reco[valid]-truth[valid]
    arrays['residual'][row] += np.histogram(error, redges)[0]
    arrays['depth_residual'][row] += np.histogram2d(truth[valid], error, (zedges, redges))[0].astype(np.int64)
    for i, tol in enumerate((50., 100., 200.)):
        arrays['within'][model, sigma, :, i] += np.histogram(truth[valid][abs(error) <= tol], zedges)[0]
    if not np.array_equal(arrays['response'][row].sum(axis=1)+arrays['failures'][row], arrays['denominator'][row]):
        raise ValueError("All-photon response closure failed")
    if (arrays['residual'][row].sum()!=arrays['response'][row].sum()
            or not np.array_equal(arrays['depth_residual'][row].sum(axis=1),arrays['response'][row].sum(axis=1))):
        raise ValueError('Residual histogram coverage failed')
    return dict(input=len(truth), reconstructed=int(valid.sum()), failed=int((~valid).sum()),
                within100mm=int(np.sum(abs(error) <= 100.)))


def reduce_one(row, plan, output):
    source = Path(row['root_file'])
    target = output/'results'/source.stem
    if target.exists():
        report = json.loads((target/'audit.json').read_text())
        if report['plan_sha256'] != plan['_sha256'] or report['input'] != row:
            raise ValueError("Existing output has different provenance")
        if file_sha256(target/'histograms.npz') != report['histograms_sha256']:
            raise ValueError("Existing output checksum mismatch")
        return report
    side = Path(str(source)+'.validation.json')
    if file_sha256(side) != row['sidecar_sha256']:
        raise ValueError("Sidecar checksum mismatch")
    metadata = json.loads(side.read_text())
    if metadata['status'] != 'valid' or metadata['entries'] != row['events']:
        raise ValueError("Invalid production sidecar")
    for k, v in plan['production_pins'].items():
        if metadata[k] != v:
            raise ValueError('Production provenance mismatch: '+k)
    table = Path(plan['reference_release'])/'reductions'/source.stem/(source.stem+'.events.csv.gz')
    validation = json.loads((table.parent/'release_validation.json').read_text())
    if validation['input'] != row or file_sha256(table) != validation['events_sha256']:
        raise ValueError("Reference event table/provenance mismatch")
    with gzip.open(table, 'rt') as f:
        expected_rows = list(csv.DictReader(f))
    expected = {int(r['source_entry']): r for r in expected_rows}
    if len(expected_rows)!=row['events'] or len(expected)!=row['events'] or set(expected)!=set(range(row['events'])):
        raise ValueError('Missing/duplicate reference event entries')
    ref = DepthContract.from_workflow(Path(__file__).parent/'sensor_workflow.json')
    contract, entrance, _ = dreamsim_front_contract(pde=ref.pde, pitch_um=ref.pitch_um,
        recovery_ns=ref.recovery_ns, workflow_sha256=ref.workflow_sha256)
    # Read every finite arrival time. The canonical sensor window remains an
    # explicit, separately accounted acquisition assumption, not a prompt cut.
    reader = replace(contract, gate_ns=(-1e300, 1e300))
    ze, re = np.asarray(plan['z_edges_mm']), np.asarray(plan['residual_edges_mm'])
    arrays = make_histograms(plan)
    records = []
    for key, (all_event, original_ids, truth, counts) in iter_corrected_root(source, row['root_sha256'], reader):
        entry = int(key.rsplit(':',1)[1])
        in_gate = (all_event.time_ns > contract.gate_ns[0]) & (all_event.time_ns < contract.gate_ns[1])
        order = np.flatnonzero(in_gate)
        order = order[np.argsort(all_event.time_ns[order], kind='stable')]
        event = SimpleNamespace(**{k:v[order] for k,v in vars(all_event).items()}, event_index=row['offset']+entry)
        z = truth['truth_z_mm'][order]
        digit = digitize_full_detector_event(event, FullDetectorScenario(ref.pde,ref.pitch_um,ref.recovery_ns), seed=4307)
        old = expected[entry]
        if len(order) != int(old['selected_cherenkov_photons']) or digit.n_detected != int(old['detected_reference']) or digit.n_accepted != int(old['accepted_reference']):
            raise ValueError('Reference photon-count mismatch at entry '+str(entry))
        batch = make_observations(event,key,contract,entrance_time_ns=entrance,
            photon_ids=original_ids[order],accepted_mask=digit.accepted_mask,external_detected_mask=digit.detected_mask)
        accepted_truth = z[batch.accepted_indices]
        stages = [truth['truth_z_mm'],truth['truth_z_mm'][~in_gate],z,z[~digit.detected_mask],
                  z[digit.detected_mask & ~digit.accepted_mask],accepted_truth]
        for i, values in enumerate(stages):
            hist = np.histogram(values, ze)[0]
            if hist.sum()!=len(values):
                raise ValueError('Photon-stage truth histogram lost entries')
            arrays['stage_truth'][i] += hist
        record = dict(entry=entry,event=row['offset']+entry,raw=counts['raw'],
            geometric_downstream=len(all_event.time_ns),nonincident_records=counts['rejected'],
            outside_sensor_window=int((~in_gate).sum()),incident_in_window=len(order),
            not_detected=len(order)-digit.n_detected,lost_holdoff=digit.n_detected-digit.n_accepted,
            accepted=digit.n_accepted,exit_outside_core=int(np.sum(np.hypot(event.local_x_mm,event.local_y_mm)>ref.core_radius_mm+1e-5)),
            max_arrival_ns=float(all_event.time_ns.max()) if len(all_event.time_ns) else None,baselines={})
        for si,sigma in enumerate(plan['sigmas_ps']):
            obs = smear_observations(batch.photons,sigma_t_ns=sigma/1000.,
                                     entrance_sigma_ns=plan['entrance_sigma_ps']/1000.,seed=plan['noise_seed'])
            np.testing.assert_array_equal(obs.photon_id,batch.photons.photon_id)
            for mi,model in enumerate(plan['baselines']):
                # Process, count and report the old calibration entries, while
                # keeping them out of that fitted baseline's evaluation scores.
                is_training = (model['name']=='historical_fixed' and row['root_sha256']==plan['old_calibration_root_sha256']
                               and row['offset']+entry in plan['old_calibration_event_indices'])
                reco = contract.inverse(obs.time_ns,model['n_eff']/contract.n_group,obs.sensor_x_mm,obs.sensor_y_mm)
                metric = accumulate(arrays,mi,si,accepted_truth[:0] if is_training else accepted_truth,
                                    reco[:0] if is_training else reco,ze,re)
                record['baselines'][model['name']+'_'+str(sigma)+'ps'] = dict(metric,calibration_entry=is_training)
        records.append(record)
    if len(records)!=row['events'] or {r['entry'] for r in records} != set(expected):
        raise ValueError("Missing/duplicate input events")
    stage = output/'staging'/(source.stem+'_'+os.environ['SLURM_JOB_ID'])
    stage.mkdir(parents=True,exist_ok=False)
    np.savez_compressed(stage/'histograms.npz',**arrays,z_edges_mm=ze,residual_edges_mm=re,sigmas_ps=plan['sigmas_ps'])
    report=dict(status='complete_accounting_physics_provisional',input=row,records=records,
        plan_sha256=plan['_sha256'],code_revision=plan['code_revision'],contract_sha256=contract.sha256,
        reader_contract_sha256=reader.sha256,histograms_sha256=file_sha256(stage/'histograms.npz'),
        slurm_job_id=os.environ['SLURM_JOB_ID'],source_table_sha256=file_sha256(table),
        target='all accepted individual photon production depths, not first hits or energy deposition',
        stage_names=['all_geometric_downstream','outside_sensor_window','incident_in_window','not_detected','lost_holdoff','accepted'],
        optical_transport_validated=False,trained_path_models=False)
    if not np.array_equal(arrays['stage_truth'][0],arrays['stage_truth'][1]+arrays['stage_truth'][2]):
        raise ValueError('Acquisition-window photon closure failed')
    if not np.array_equal(arrays['stage_truth'][2],arrays['stage_truth'][3]+arrays['stage_truth'][4]+arrays['stage_truth'][5]):
        raise ValueError('Sensor photon closure failed')
    with (stage/'audit.json').open('x') as f:
        json.dump(report,f,indent=2,allow_nan=False)
    target.parent.mkdir(parents=True,exist_ok=True)
    stage.rename(target)
    print('PASS',row['sample'],row['offset'],len(records),sum(r['accepted'] for r in records),flush=True)
    return report


def residual_quantiles(hist, edges):
    """Binned conditional quantiles; uncertainty from binning is <= bin width."""
    hist = np.asarray(hist)
    if hist.sum()==0:
        return [None]*3
    cdf = np.cumsum(hist)/hist.sum()
    return [float(np.interp(q,np.r_[0.,cdf],edges)) for q in (.16,.5,.84)]


def merge_campaign(rows, plan, output, checkpoint, *, bootstrap_draws=500):
    """Strict aggregate-only merge; resample complete independent source files."""
    if checkpoint.exists():
        raise FileExistsError('Immutable checkpoint already exists')
    totals, files, blocks = {}, [], {}
    ze, re = np.asarray(plan['z_edges_mm']), np.asarray(plan['residual_edges_mm'])
    for row in rows:
        directory = output/'results'/Path(row['root_file']).stem
        audit_file, hist_file = directory/'audit.json',directory/'histograms.npz'
        report = json.loads(audit_file.read_text())
        if (report['input']!=row or report['plan_sha256']!=plan['_sha256']
                or report['contract_sha256']!=plan['contract_sha256']
                or report['code_revision']!=plan['code_revision']
                or report['status']!='complete_accounting_physics_provisional'
                or report['optical_transport_validated'] or report['trained_path_models']
                or file_sha256(hist_file)!=report['histograms_sha256']):
            raise ValueError('Merge provenance/hash mismatch: '+str(directory))
        with np.load(hist_file,allow_pickle=False) as src:
            if (not np.array_equal(src['z_edges_mm'],ze) or not np.array_equal(src['residual_edges_mm'],re)
                    or not np.array_equal(src['sigmas_ps'],plan['sigmas_ps'])):
                raise ValueError('Incompatible histogram axes')
            a = {k:src[k] for k in make_histograms(plan)}
        shapes = make_histograms(plan)
        if any(v.shape!=shapes[k].shape or v.dtype.kind not in 'iu' or np.any(v<0) for k,v in a.items()):
            raise ValueError('Invalid histogram shape/counts')
        if (not np.array_equal(a['response'].sum(axis=-1)+a['failures'],a['denominator'])
                or not np.array_equal(a['depth_residual'].sum(axis=-1),a['response'].sum(axis=-1))
                or not np.array_equal(a['residual'].sum(axis=-1),a['response'].sum(axis=(-1,-2)))
                or not np.array_equal(a['stage_truth'][0],a['stage_truth'][1]+a['stage_truth'][2])
                or not np.array_equal(a['stage_truth'][2],a['stage_truth'][3:6].sum(axis=0))):
            raise ValueError('Merged response/stage/coverage closure failed')
        if np.any(np.diff(a['within'],axis=-1)<0) or np.any(a['within'][...,-1]>a['response'].sum(axis=-1)):
            raise ValueError('Nonphysical inclusive success counts')
        records = report['records']
        if (len(records)!=row['events'] or {r['entry'] for r in records}!=set(range(row['events']))
                or any(r['event']!=row['offset']+r['entry'] for r in records)
                or any(r['raw']!=r['geometric_downstream']+r['nonincident_records'] for r in records)):
            raise ValueError('Missing/duplicate event identity at merge')
        stages = ['geometric_downstream','outside_sensor_window','incident_in_window','not_detected','lost_holdoff','accepted']
        if not np.array_equal(a['stage_truth'].sum(axis=1),[sum(r[k] for r in records) for k in stages]):
            raise ValueError('Event/photon-stage totals disagree')
        metrics = np.zeros((len(plan['baselines']),len(plan['sigmas_ps']),4),np.int64)
        for mi,model in enumerate(plan['baselines']):
            for si,sigma in enumerate(plan['sigmas_ps']):
                tag = model['name']+'_'+str(sigma)+'ps'
                training = [r for r in records if model['name']=='historical_fixed'
                    and row['root_sha256']==plan['old_calibration_root_sha256']
                    and r['event'] in plan['old_calibration_event_indices']]
                expected = sum(r['accepted'] for r in records)-sum(r['accepted'] for r in training)
                metrics[mi,si] = [sum(r['baselines'][tag][k] for r in records)
                                  for k in ('input','reconstructed','failed','within100mm')]
                hmetrics = [a['denominator'][mi,si].sum(),a['response'][mi,si].sum(),
                            a['failures'][mi,si].sum(),a['within'][mi,si,:,1].sum()]
                if metrics[mi,si,0]!=expected or not np.array_equal(metrics[mi,si],hmetrics):
                    raise ValueError('All-photon scored denominator mismatch')
        sample = row['sample']
        if sample not in totals:
            totals[sample] = make_histograms(plan); blocks[sample] = []
        for k in a: totals[sample][k] += a[k]
        blocks[sample].append((metrics,a['residual']))
        files.append(dict(input=row,audit_sha256=file_sha256(audit_file),histograms_sha256=report['histograms_sha256']))
        if file_sha256(hist_file)!=report['histograms_sha256']:
            raise ValueError('Histogram changed during merge')
    stage = checkpoint.with_name(checkpoint.name+'.staging')
    stage.mkdir(parents=True,exist_ok=False)
    samples={}
    for sample in sorted(totals):
        a=totals[sample]; parts=blocks[sample]
        metrics=np.stack([b[0] for b in parts]); residuals=np.stack([b[1] for b in parts])
        rng=np.random.default_rng(plan['noise_seed'])
        weights=rng.multinomial(len(parts),np.ones(len(parts))/len(parts),size=bootstrap_draws)
        cases={}
        for mi,model in enumerate(plan['baselines']):
            for si,sigma in enumerate(plan['sigmas_ps']):
                m=metrics[:,mi,si].sum(axis=0); bm=weights @ metrics[:,mi,si]
                bh=weights.astype(float) @ residuals[:,mi,si].astype(float)
                quantiles=residual_quantiles(a['residual'][mi,si],re)
                bq=np.asarray([residual_quantiles(h,re) for h in bh],float)
                failure=np.divide(100*bm[:,2],bm[:,0],out=np.full(len(bm),np.nan),where=bm[:,0]>0)
                good=np.divide(100*bm[:,3],bm[:,0],out=np.full(len(bm),np.nan),where=bm[:,0]>0)
                interval=lambda v: np.nanpercentile(v,[2.5,97.5]).tolist() if np.isfinite(v).any() else [None,None]
                cases[model['name']+'_'+str(sigma)+'ps']=dict(input=int(m[0]),reconstructed=int(m[1]),failed=int(m[2]),
                    within100mm=int(m[3]),failure_pct=float(100*m[2]/m[0]) if m[0] else None,
                    within100mm_pct=float(100*m[3]/m[0]) if m[0] else None,
                    conditional_residual_q16_q50_q84_mm=quantiles,
                    conditional_sigma68_mm=None if quantiles[0] is None else (quantiles[2]-quantiles[0])/2,
                    file_bootstrap_95=dict(failure_pct=interval(failure),within100mm_pct=interval(good),
                        median_residual_mm=interval(bq[:,1]),sigma68_mm=interval((bq[:,2]-bq[:,0])/2)))
        name=sample+'.histograms.npz'
        np.savez_compressed(stage/name,**a,z_edges_mm=ze,residual_edges_mm=re,sigmas_ps=plan['sigmas_ps'])
        samples[sample]=dict(events=sum(r['input']['events'] for r in files if r['input']['sample']==sample),
            files=len(parts),photon_stages=dict(zip(['all_geometric_downstream','outside_sensor_window','incident_in_window',
            'not_detected','lost_holdoff','accepted'],map(int,a['stage_truth'].sum(axis=1)))),
            cases=cases,histograms_file=name,histograms_sha256=file_sha256(stage/name))
    report=dict(status='complete_accounting_physics_provisional',plan_sha256=plan['_sha256'],
        reduction_code_revision=plan['code_revision'],merge_code_sha256=file_sha256(__file__),
        files=len(files),events=sum(r['input']['events'] for r in files),samples=samples,inputs=files,
        optical_transport_validated=False,trained_path_models=False,closure=True,
        bootstrap=dict(unit='independent source ROOT file, all photons kept together',draws=bootstrap_draws,
            seed=plan['noise_seed'],interval='percentile 95%; statistical only, not transport/model uncertainty'),
        residual_quantile_bin_width_mm=float(np.diff(re).max()),
        selection=plan['selection'],scope=plan['scope'],target=plan['target'])
    with (stage/'summary.json').open('x') as f: json.dump(report,f,indent=2,allow_nan=False)
    stage.rename(checkpoint)
    print('MERGED',len(files),report['events'],checkpoint,flush=True)
    return report


def main():
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Raw campaign analysis requires Slurm')
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan',required=True); ap.add_argument('--output',required=True)
    ap.add_argument('--task',type=int,default=int(os.environ.get('SLURM_ARRAY_TASK_ID','0')))
    ap.add_argument('--merge',type=Path,help='New immutable aggregate checkpoint; no raw ROOT reads')
    args=ap.parse_args()
    plan_path=Path(args.plan); plan=json.loads(plan_path.read_text()); plan['_sha256']=file_sha256(plan_path)
    for name,digest in plan['code_sha256'].items():
        if file_sha256((plan_path.parent if args.merge else Path(__file__).parent)/name)!=digest:
            raise ValueError('Analysis code checksum mismatch: '+name)
    manifest=Path(plan['input_manifest'])
    if file_sha256(manifest)!=plan['input_manifest_sha256']:
        raise ValueError('Input manifest checksum mismatch')
    if file_sha256(Path(plan['reference_release'])/'plan.json')!=plan['reference_plan_sha256']:
        raise ValueError('Corrected reference-plan checksum mismatch')
    rows=json.loads(manifest.read_text())
    if len(rows)!=1648 or sum(r['events'] for r in rows)!=40000:
        raise ValueError('Campaign event/file count mismatch')
    if (len({r['root_file'] for r in rows})!=len(rows) or len({tuple(r['seeds']) for r in rows})!=len(rows)
            or len({r['root_sha256'] for r in rows})!=len(rows)
            or len({Path(r['root_file']).stem for r in rows})!=len(rows)):
        raise ValueError('Duplicate source/seed identity')
    samples = {s:sum(r['events'] for r in rows if r['sample']==s) for s in {r['sample'] for r in rows}}
    if len(samples)!=8 or set(samples.values())!={5000}:
        raise ValueError('Sample event count mismatch')
    identities = [(r['sample'],i) for r in rows for i in range(r['offset'],r['offset']+r['events'])]
    if len(set(identities))!=len(identities):
        raise ValueError('Duplicate sample event index')
    make_histograms(plan)
    if args.merge:
        merge_campaign(rows,plan,Path(args.output),args.merge)
        return
    chunk=plan['files_per_task']; start=args.task*chunk
    if start<0 or start>=len(rows):
        raise ValueError('Invalid array task')
    output=Path(args.output); errors=[]
    for row in rows[start:start+chunk]:
        try:
            reduce_one(row,plan,output)
        except Exception as exc:
            errors.append(dict(root=row['root_file'],exception=type(exc).__name__,message=str(exc)))
            print('FAIL',row['root_file'],repr(exc),flush=True)
    if errors:
        dest=output/'failures'; dest.mkdir(parents=True,exist_ok=True)
        with (dest/('task_'+str(args.task)+'_'+os.environ['SLURM_JOB_ID']+'.json')).open('x') as f:
            json.dump(errors,f,indent=2)
        sys.exit(1)


if __name__=='__main__':
    main()
