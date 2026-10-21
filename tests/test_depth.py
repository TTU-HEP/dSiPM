from dataclasses import replace
import copy
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
import numpy as np
from dsipm.depth import (C_MM_NS, Calibration, DepthContract, DepthStatus, PhotonObservations,
                        channel_delay, fit_calibration, reconstruct, evaluate, residual_summary, response_matrix)
from dsipm.sensors import make_observations, smear_observations
from dsipm.geometry import (FIBER_OFFSETS_MM, dreamsim_front_contract, localize_endpoints,
                           rotation, world_to_calorimeter, RayStatus, trace_cylinder)
from dsipm.io import save_cache, load_cache, iter_corrected_root, corrected_event_from_arrays

def observations(t, r=None, events=None, rx=None, ry=None):
    t = np.asarray(t, float); n = len(t)
    return PhotonObservations(np.array(["test"]*n) if events is None else np.asarray(events),
        np.arange(n, dtype=np.int64), np.zeros(n, np.int64) if rx is None else np.asarray(rx, np.int64),
        np.zeros(n, np.int64) if ry is None else np.asarray(ry, np.int64), t,
        np.full(n, .1) if r is None else np.asarray(r, float), np.zeros(n), np.zeros(n))


def calibration(contract=None, k=1.02, quantiles=None):
    contract = contract or DepthContract()
    o = observations(np.linspace(7, 9, 30), events=["train"]*30)
    c = fit_calibration(o, np.full(30, k), contract, source_sha256="a"*64)
    if quantiles is not None:
        p = copy.deepcopy(c.payload)
        p["joint_quantiles"] = np.broadcast_to(quantiles, np.asarray(p["joint_quantiles"]).shape).tolist()
        c = Calibration(p)
    return c


def event(n=700):
    rng = np.random.default_rng(91)
    phi = rng.uniform(0, 2*np.pi, n)
    radius = .39*np.sqrt(rng.random(n))
    return SimpleNamespace(time_ns=rng.uniform(7, 17, n), local_x_mm=radius*np.cos(phi),
        local_y_mm=radius*np.sin(phi), sensor_id=rng.integers(0, 4, n),
        sensor_x_mm=np.zeros(n), sensor_y_mm=np.zeros(n), rod_ix=np.zeros(n, np.int64), rod_iy=np.zeros(n, np.int64))


class DepthReconstructionTests(unittest.TestCase):
    def test_forward_inverse_and_group_index(self):
        c = DepthContract(n_group=1.55, n_phase=1.49)
        z = np.array([0., 200., 1000., 2000.]); k = np.array([1., 1.02, 1.5, 3.])
        np.testing.assert_allclose(c.inverse(c.forward_time(z, k), k), z, atol=3e-12)
        self.assertEqual(c.n_phase, 1.49)

    def test_explicit_fixed_path_baseline_matches_historical_inverse(self):
        c, _, _ = dreamsim_front_contract()
        o = observations([8., 9., 10.])
        neff = 1.553012152751959
        r = reconstruct(o, calibration(c), c, method="fixed", fixed_n_eff=neff)
        expected = c.inverse(o.time_ns, neff/c.n_group)
        expected[(expected < 0) | (expected > c.length_mm)] = np.nan
        np.testing.assert_allclose(r.z_mm, expected)
        with self.assertRaisesRegex(ValueError, "explicit"):
            reconstruct(o, calibration(c), c, method="fixed")

    def test_emission_delay_bias_is_not_clock_resolution(self):
        c = DepthContract(n_group=1.504)
        t = c.forward_time(1000, 1.02, emission_delay_ns=.1)
        z = c.inverse(t, 1.02)
        self.assertAlmostEqual(z-1000, -.1*C_MM_NS/(1.504*1.02-1))

    def test_canonical_workflow_is_pinned_and_incompatible_with_changes(self):
        from dsipm.io import file_sha256
        path = Path(__file__).resolve().parents[1]/"sensor_workflow.json"
        reference = json.loads(path.read_text())["reference_scenario"]
        c = DepthContract.from_workflow(path)
        self.assertEqual(c.workflow_sha256, file_sha256(path))
        self.assertEqual(c.pitch_um, reference["spad_pitch_um"])
        self.assertEqual(c.pde, reference["effective_pde"])
        changed = replace(c, workflow_sha256="f"*64)
        with self.assertRaisesRegex(ValueError, "contract"):
            reconstruct(observations([8.]), calibration(c), changed)

    def test_local_timing_uncertainties_are_separate_not_total_errors(self):
        c = DepthContract(); o = observations([c.forward_time(1000, 1.02)])
        r = reconstruct(o, calibration(c), c, sigma_t_ns=.05, entrance_sigma_ns=.01)
        self.assertAlmostEqual(r.photon_timing_sigma_z_mm[0], .05*C_MM_NS/(c.n_group*1.02-c.front_z))
        self.assertAlmostEqual(r.entrance_timing_sigma_z_mm[0], r.photon_timing_sigma_z_mm[0]/5)

    def test_exact_detector_endpoints_and_direct_exit_candidates(self):
        c = DepthContract()
        o = observations(c.forward_time(np.array([0., 2000.]), 1.02))
        r = reconstruct(o, calibration(c), c, method="radius_delay")
        np.testing.assert_allclose(r.z_mm, [0, 2000], atol=1e-9)
        # TIR does not apply to a high-angle ray generated microns from end.
        high = observations([c.forward_time(1999.999, 2.)])
        r = reconstruct(high, calibration(c, k=2.), c, method="radius")
        self.assertAlmostEqual(r.z_mm[0], 1999.999)

    def test_isolated_photon_has_exactly_zero_support(self):
        c = DepthContract(); o = observations([c.forward_time(1000, 1.02)])
        r = reconstruct(o, calibration(c), c, sigma_t_ns=.05)
        self.assertTrue(r.prior_seed[0])
        self.assertEqual(float(r.candidate_density.max()), 0)
        self.assertEqual(r.status[0], DepthStatus.CALIBRATED_FALLBACK)

    def test_invalid_seed_cannot_vote_at_boundary(self):
        c = DepthContract()
        o = observations([c.forward_time(-1000, 1.02), c.forward_time(1000, 1.5)], [.1, .95], rx=[0, 1])
        r = reconstruct(o, calibration(c, quantiles=np.linspace(1.1, 2., 10)), c)
        self.assertEqual(r.invalid_clean_seed_count, 1)
        self.assertFalse(r.prior_seed[0])
        self.assertEqual(float(r.candidate_density.max()), 0)

    def test_late_photon_is_unreconstructed_not_boundary(self):
        c = DepthContract(); o = observations([100.])
        r = reconstruct(o, calibration(c), c)
        self.assertTrue(np.isnan(r.z_mm[0]))
        self.assertEqual(r.status[0], DepthStatus.NO_PHYSICAL_CANDIDATE)
        self.assertEqual(residual_summary([1000], r.z_mm)["inclusive_within"]["100.0"], 0)

    def test_real_median_fallback_not_55_percent(self):
        c = DepthContract(); o = observations([c.forward_time(1000, 1.02)])
        r = reconstruct(o, calibration(c, quantiles=np.linspace(1.005, 1.035, 10)), c)
        self.assertAlmostEqual(r.z_mm[0], 1000)

    def test_neighbour_support_and_leave_rod_out(self):
        c = DepthContract(); o = observations(np.full(3, c.forward_time(1000, 1.02)), rx=[0, 0, 1])
        r = reconstruct(o, calibration(c), c, exclude_same_rod=True)
        self.assertTrue(np.all(r.status == DepthStatus.SUPPORTED))
        same = reconstruct(o.subset([0, 1]), calibration(c), c, exclude_same_rod=True)
        self.assertTrue(np.all(same.candidate_density == 0))

    def test_events_cannot_share_prior_and_order_does_not_matter(self):
        c = DepthContract(); t = c.forward_time(1000, 1.02)
        o = observations([t, t, t], events=["a", "a", "b"], rx=[0, 1, 0])
        a = reconstruct(o, calibration(c), c)
        self.assertEqual(a.candidate_density[2].max(), 0)
        order = np.array([2, 0, 1])
        b = reconstruct(o.subset(order), calibration(c), c)
        np.testing.assert_equal(a.z_mm[order], b.z_mm)

    def test_delay_no_packed_id_collisions_and_exact_zero_bin(self):
        o = observations([9, 8, 7, 7.005, 7.02], events=["e", "other", "e", "e", "e"], rx=[1000000, 0, 0, 0, 0])
        np.testing.assert_allclose(channel_delay(o), [0, 0, 0, .005, .02])
        c = calibration(); _, tb = c.bins(o)
        self.assertEqual(tb[2], 0); self.assertEqual(tb[3], 1)
        self.assertGreater(tb[4], 1)

    def test_train_test_overlap_and_model_mismatch_rejected(self):
        c = DepthContract(); cal = calibration(c)
        with self.assertRaisesRegex(ValueError, "overlap"):
            reconstruct(observations([8.], events=["train"]), cal, c)
        for changed in [replace(c, n_group=1.51), replace(c, gate_ns=(0, 20)), replace(c, pitch_um=50)]:
            with self.assertRaisesRegex(ValueError, "contract"):
                reconstruct(observations([8.]), cal, changed)

    def test_empty_and_input_validation(self):
        c = DepthContract()
        r = reconstruct(observations([]), calibration(c), c)
        self.assertEqual(r.z_mm.size, 0)
        self.assertEqual(evaluate(observations([]), [], r, c)["input"], 0)
        for t, rad in [([np.nan], [.1]), ([8.], [-.1])]:
            with self.assertRaises(ValueError): observations(t, rad)
        with self.assertRaises(ValueError): reconstruct(observations([8.]), calibration(c), c, window=2)
        with self.assertRaises(ValueError): DepthContract(n_group=.9)

    def test_artifact_hash_and_immutable_cache(self):
        c = DepthContract(); cal = calibration(c)
        with TemporaryDirectory() as d:
            p = Path(d)/"cal.json"; cal.save(p)
            self.assertEqual(Calibration.load(p).sha256, cal.sha256)
            with self.assertRaises(FileExistsError): cal.save(p)
            data = json.loads(p.read_text()); data["calibration"]["training_sigma_ns"] = 1
            p.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError, "checksum"): Calibration.load(p)
            target = Path(d)/"cache"
            save_cache(target, observations([8.]), c, provenance={"case": "unit-test"}, truth_z_mm=[500.])
            obs, cc, truth, _ = load_cache(target)
            self.assertEqual(cc.sha256, c.sha256); self.assertEqual(truth["truth_z_mm"][0], 500)
            with self.assertRaises(FileExistsError): save_cache(target, obs, c, provenance={})

    def test_inclusive_depth_closure_and_event_bootstrap(self):
        c = DepthContract(); truth = np.array([0., 500., 1000., 2000.])
        # Avoid floating-point inversion at z=0 for this bookkeeping test.
        o = observations(c.forward_time([10, 500, 1000, 1990], 1.02), events=["a", "a", "b", "b"])
        r = reconstruct(o, calibration(c), c)
        a = evaluate(o, truth, r, c, bootstrap_draws=25, truth_edges_mm=[100, 1000, 1900])
        b = evaluate(o, truth, r, c, bootstrap_draws=25, truth_edges_mm=[100, 1000, 1900])
        self.assertEqual(a, b); self.assertTrue(a["closure"])
        self.assertEqual(a["truth_underflow"], 1); self.assertEqual(a["truth_overflow"], 1)
        self.assertEqual(a["event_block_bootstrap"]["draws"], 25)


class DepthResponseTests(unittest.TestCase):
    def test_orientation_endpoints_and_failure_denominators(self):
        r = response_matrix([0., 100., 500., 2000.], [2000., np.nan, 500., 0.], bins=4)
        h = np.array(r["counts"])
        self.assertEqual(h[0, 3], 1)  # true on x, reco on y; never transposed
        self.assertEqual(h[3, 0], 1)
        self.assertEqual(r["input_by_truth"], [2, 1, 0, 1])
        self.assertEqual(r["unreconstructed_by_truth"], [1, 0, 0, 0])
        self.assertEqual(r["within100mm_by_truth"], [0, 1, 0, 0])
        np.testing.assert_equal(h.sum(axis=1)+r["unreconstructed_by_truth"], r["input_by_truth"])

    def test_empty_and_conditional_quantiles(self):
        self.assertEqual(response_matrix([], [], bins=2)["input"], 0)
        r = response_matrix([100., 200., 300.], [110., 230., np.nan], bins=2)
        self.assertEqual(r["conditional_residual_q16_q50_q84_mm"][0][1], 20.)
        self.assertEqual(r["input_by_truth"][0], 3)
        self.assertEqual(r["unreconstructed_by_truth"][0], 1)

    def test_out_of_detector_response_cannot_be_silently_dropped(self):
        for truth, z in [([500.], [-1.]), ([500.], [2001.]), ([-1.], [np.nan]), ([2001.], [5.])]:
            with self.assertRaisesRegex(ValueError, "physical"):
                response_matrix(truth, z)
        for bins in [1, 1.5, True]:
            with self.assertRaises(ValueError): response_matrix([10.], [20.], bins=bins)


class AllPhotonCampaignTests(unittest.TestCase):
    def test_all_inputs_include_nonphysical_and_nonfinite_failures(self):
        from reduceDepthCampaign import make_histograms, accumulate
        plan = dict(z_edges_mm=[0.,500.,1000.,2000.], residual_edges_mm=[-2000.,0.,2000.],
                    baselines=[dict(name='axial',n_eff=1.49)], sigmas_ps=[0.,50.])
        a = make_histograms(plan)
        r = accumulate(a,0,0,[0.,500.,2000.,1000.],[np.nan,500.,3000.,np.inf],
                       plan['z_edges_mm'],plan['residual_edges_mm'])
        self.assertEqual(r,dict(input=4,reconstructed=1,failed=3,within100mm=1))
        np.testing.assert_array_equal(a['response'][0,0].sum(axis=1)+a['failures'][0,0],a['denominator'][0,0])
        self.assertEqual(a['response'][0,0,1,1],1)
        accumulate(a,0,0,[0.,2000.],[2000.,0.],plan['z_edges_mm'],plan['residual_edges_mm'])
        self.assertEqual(a['residual'].sum(),3)  # both residual endpoints retained
        self.assertEqual(a['depth_residual'].sum(),3)
        self.assertEqual(a['denominator'].sum(),6)
        self.assertEqual(a['within'][0,0,:,1].sum(),1)
        self.assertEqual(a['denominator'][0,1].sum(),0)
        with self.assertRaisesRegex(ValueError,'Misaligned'):
            accumulate(a,0,0,[0.,100.],[0.],plan['z_edges_mm'],plan['residual_edges_mm'])
        with self.assertRaisesRegex(ValueError,'Invalid true'):
            accumulate(a,0,0,[-1.],[0.],plan['z_edges_mm'],plan['residual_edges_mm'])

    def test_plan_cannot_hide_histogram_overflows(self):
        from reduceDepthCampaign import make_histograms
        good = dict(z_edges_mm=[0.,1000.,2000.],residual_edges_mm=[-2000.,0.,2000.],
                    baselines=[dict(name='axial',n_eff=1.49)],sigmas_ps=[0.,50.])
        for bad in [dict(z_edges_mm=[0.,1000.,1999.]),dict(residual_edges_mm=[-1000.,0.,2000.]),
                    dict(sigmas_ps=[0.,0.]),dict(sigmas_ps=[-1.]),dict(z_edges_mm=[0.,2000.,1000.])]:
            with self.assertRaises(ValueError): make_histograms(dict(good,**bad))

    def test_campaign_requires_slurm_before_reading(self):
        from unittest.mock import patch
        import reduceDepthCampaign
        with patch.dict('os.environ',{},clear=True):
            with self.assertRaisesRegex(RuntimeError,'Slurm'): reduceDepthCampaign.main()

    def test_complete_synthetic_event_accounts_for_gate_and_calibration(self):
        from unittest.mock import patch
        import csv
        import gzip
        from dsipm.io import file_sha256
        from dsipm.sensors import digitize_full_detector_event, FullDetectorScenario
        from reduceDepthCampaign import reduce_one, merge_campaign
        raw = event()
        raw.time_ns[:2] = [-1.,45.]  # reader retains these; sensor window does not
        order = np.flatnonzero((raw.time_ns>0)&(raw.time_ns<40))
        order = order[np.argsort(raw.time_ns[order],kind='stable')]
        e = SimpleNamespace(**{k:v[order] for k,v in vars(raw).items()},event_index=0)
        ref = DepthContract.from_workflow()
        digit = digitize_full_detector_event(e,FullDetectorScenario(ref.pde,ref.pitch_um,ref.recovery_ns))
        with TemporaryDirectory() as tmp:
            root = Path(tmp); source = root/'fixture.root'
            side = Path(str(source)+'.validation.json')
            side.write_text(json.dumps(dict(status='valid',entries=1)))
            row = dict(root_file=str(source),sidecar_sha256=file_sha256(side),root_sha256='a'*64,
                       sample='fixture',events=1,offset=0)
            table = root/'reference/reductions/fixture/fixture.events.csv.gz'
            table.parent.mkdir(parents=True)
            with gzip.open(table,'wt') as f:
                writer = csv.DictWriter(f,fieldnames=['source_entry','selected_cherenkov_photons','detected_reference','accepted_reference'])
                writer.writeheader(); writer.writerow(dict(source_entry=0,selected_cherenkov_photons=len(order),
                    detected_reference=digit.n_detected,accepted_reference=digit.n_accepted))
            (table.parent/'release_validation.json').write_text(json.dumps(dict(input=row,events_sha256=file_sha256(table))))
            plan = dict(_sha256='b'*64,code_revision='fixture',production_pins={},reference_release=str(root/'reference'),
                z_edges_mm=[0.,1000.,2000.],residual_edges_mm=[-2000.,0.,2000.],sigmas_ps=[0,50],
                baselines=[dict(name='axial',n_eff=1.49),dict(name='historical_fixed',n_eff=1.553)],
                entrance_sigma_ps=0,noise_seed=20261021,old_calibration_root_sha256='a'*64,old_calibration_event_indices=[0],
                contract_sha256=dreamsim_front_contract(pde=ref.pde,pitch_um=ref.pitch_um,recovery_ns=ref.recovery_ns,
                    workflow_sha256=ref.workflow_sha256)[0].sha256,selection='synthetic',scope='test',target='test')
            truth = dict(truth_z_mm=np.linspace(0.,2000.,len(raw.time_ns)))
            yielded = [('a'*64+':0',(raw,np.arange(len(raw.time_ns)),truth,dict(raw=703,rejected=3)))]
            with patch.dict('os.environ',SLURM_JOB_ID='fixture'),patch('reduceDepthCampaign.iter_corrected_root',return_value=iter(yielded)) as read:
                report = reduce_one(row,plan,root/'output')
            self.assertEqual(read.call_args.args[2].gate_ns,(-1e300,1e300))
            r = report['records'][0]
            self.assertEqual(r['raw'],703); self.assertEqual(r['geometric_downstream'],700)
            self.assertEqual(r['outside_sensor_window'],2)
            self.assertEqual(r['incident_in_window'],r['not_detected']+r['lost_holdoff']+r['accepted'])
            self.assertTrue(r['baselines']['historical_fixed_50ps']['calibration_entry'])
            with np.load(root/'output/results/fixture/histograms.npz') as h:
                self.assertEqual(h['denominator'][0,0].sum(),digit.n_accepted)
                self.assertEqual(h['denominator'][0,1].sum(),digit.n_accepted)
                self.assertEqual(h['denominator'][1].sum(),0)
                np.testing.assert_array_equal(h['stage_truth'].sum(axis=1),[700,2,698,698-digit.n_detected,digit.n_detected-digit.n_accepted,digit.n_accepted])
            self.assertEqual(reduce_one(row,plan,root/'output'),report)
            merged = merge_campaign([row],plan,root/'output',root/'checkpoint',bootstrap_draws=5)
            self.assertEqual(merged['events'],1)
            self.assertEqual(merged['samples']['fixture']['photon_stages']['accepted'],digit.n_accepted)
            self.assertEqual(merged['samples']['fixture']['cases']['axial_50ps']['input'],digit.n_accepted)
            self.assertFalse(merged['optical_transport_validated'])
            with self.assertRaises(FileExistsError): merge_campaign([row],plan,root/'output',root/'checkpoint')
            (root/'output/results/fixture/histograms.npz').write_bytes(b'corrupt')
            with self.assertRaisesRegex(ValueError,'checksum'): reduce_one(row,plan,root/'output')
            with self.assertRaisesRegex(ValueError,'hash mismatch'): merge_campaign([row],plan,root/'output',root/'bad_checkpoint')

    def test_residual_quantiles_keep_empty_and_exact_endpoints(self):
        from reduceDepthCampaign import residual_quantiles
        self.assertEqual(residual_quantiles([0,0],[-2,0,2]),[None]*3)
        np.testing.assert_allclose(residual_quantiles([1,1],[-2,0,2]),[-1.36,0.,1.36])

    def test_campaign_plot_requires_full_checkpoint_plan_and_hashes(self):
        from unittest.mock import patch
        from dsipm.io import file_sha256
        from reduceDepthCampaign import make_histograms
        from studyDepthReconstruction import plot_campaign
        with TemporaryDirectory() as tmp:
            root=Path(tmp)
            plan=dict(baselines=[dict(name='axial',n_eff=1.49)],z_edges_mm=[0,1000,2000],
                residual_edges_mm=[-2000,0,2000],sigmas_ps=[50],contract={},entrance_sigma_ps=0,scope='test',selection='test')
            path=root/'plan.json';path.write_text(json.dumps(plan))
            a=make_histograms(plan);a['response'][0,0,0,1]=1;a['response'][0,0,1,0]=1
            a['failures'][0,0,0]=1;a['denominator'][0,0]=[2,1]
            a['depth_residual'][0,0,0,1]=1;a['depth_residual'][0,0,1,0]=1
            samples={}
            for p in ('eplus','piplus'):
                for e in (10,30,100,300):
                    sample=p+'_'+str(e)+'gev';name=sample+'.npz'
                    np.savez_compressed(root/name,**a,z_edges_mm=plan['z_edges_mm'],
                        residual_edges_mm=plan['residual_edges_mm'],sigmas_ps=[50])
                    samples[sample]=dict(events=5000,histograms_file=name,histograms_sha256=file_sha256(root/name),
                        cases={'axial_50ps':dict(input=3,reconstructed=2,failed=1,within100mm=0,within100mm_pct=0.,
                            failure_pct=100/3,conditional_residual_q16_q50_q84_mm=[-1000.,0.,1000.],conditional_sigma68_mm=1000.,
                            file_bootstrap_95=dict(within100mm_pct=[0.,0.],failure_pct=[100/3,100/3]))})
            summary=dict(status='complete_accounting_physics_provisional',closure=True,files=1648,events=40000,
                samples=samples,plan_sha256=file_sha256(path),optical_transport_validated=False,residual_quantile_bin_width_mm=2000)
            (root/'summary.json').write_text(json.dumps(summary))
            args=SimpleNamespace(input=str(root),plan=str(path),output=str(root/'plots'),sigma_ps=50,baseline='axial')
            with patch('studyDepthReconstruction.render_response') as render:
                plot_campaign(args)
            self.assertEqual(render.call_count,8)
            self.assertEqual(render.call_args.args[0]['input_by_truth'],[2,1])
            self.assertEqual(len((root/'plots/all_resolution_results.csv').read_text().splitlines()),9)
            args.output=str(root/'must_not_exist')
            (root/'piplus_300gev.npz').write_bytes(b'corrupt')
            with self.assertRaisesRegex(ValueError,'hash mismatch'): plot_campaign(args)
            self.assertFalse(Path(args.output).exists())
            path.write_text('{}')
            with self.assertRaisesRegex(ValueError,'plan checksum'): plot_campaign(args)


class DepthSensorTests(unittest.TestCase):
    def test_matched_zero_and_nonzero_noise(self):
        c = DepthContract(); e = event()
        batches = [make_observations(e, "rootsha:7", c, sigma_t_ns=s) for s in [0, .05, .2, .5]]
        for b in batches[1:]:
            np.testing.assert_equal(b.accepted_indices, batches[0].accepted_indices)
            np.testing.assert_equal(b.detected_mask, batches[0].detected_mask)
        a, b, d = [x.photons.time_ns for x in batches[:3]]
        np.testing.assert_allclose((b-a)*4, d-a, atol=5e-15)

    def test_order_and_partition_invariant_variates(self):
        c = DepthContract(); e = event(); order = np.random.default_rng(4).permutation(len(e.time_ns))
        a = make_observations(e, "source:3", c, sigma_t_ns=.05)
        p = SimpleNamespace(**{k: v[order] for k, v in vars(e).items()})
        b = make_observations(p, "source:3", c, sigma_t_ns=.05, photon_ids=order)
        ia, ib = np.argsort(a.photons.photon_id), np.argsort(b.photons.photon_id)
        np.testing.assert_equal(a.photons.photon_id[ia], b.photons.photon_id[ib])
        np.testing.assert_equal(a.photons.time_ns[ia], b.photons.time_ns[ib])
        sub = a.photons.subset(np.arange(0, len(ia), 2))
        s = smear_observations(a.photons, sigma_t_ns=.1)
        ss = smear_observations(sub, sigma_t_ns=.1)
        np.testing.assert_equal(s.time_ns[::2], ss.time_ns)

    def test_clock_error_is_event_common_and_independent_of_sigma(self):
        o = observations([8., 9., 10.], events=["a", "a", "b"])
        r = smear_observations(o, entrance_sigma_ns=.05)
        self.assertAlmostEqual(r.time_ns[0]-8, r.time_ns[1]-9)
        self.assertNotAlmostEqual(r.time_ns[0]-8, r.time_ns[2]-10)

    def test_holdoff_nonparalyzable_and_fibers_independent(self):
        e = event(4)
        e.time_ns = np.array([7., 8., 11.3, 7.1]); e.sensor_id = np.array([1, 1, 1, 2])
        e.local_x_mm[:] = 0; e.local_y_mm[:] = 0
        b = make_observations(e, "a", DepthContract(pde=1.))
        np.testing.assert_equal(b.accepted_mask, [True, False, True, True])
        self.assertEqual(b.counts["detected"], b.counts["accepted"]+b.counts["lost_holdoff"])

    def test_truth_never_read_by_observation_adapter(self):
        e = event(10); e.production_z_mm = "poisoned unused truth"
        b = make_observations(e, "a", DepthContract(pde=1.))
        self.assertGreater(len(b.photons.time_ns), 0)

    def test_external_mask_does_not_invent_detection_counts(self):
        e = event(3)
        b = make_observations(e, "a", DepthContract(), accepted_mask=np.array([True, False, True]))
        self.assertIsNone(b.counts["detected"])
        self.assertIsNone(b.counts["lost_holdoff"])
        with self.assertRaisesRegex(ValueError, "closure"):
            make_observations(e, "a", DepthContract(), accepted_mask=np.array([True, False, True]),
                              external_detected_mask=np.array([False, False, True]))


class DepthGeometryTests(unittest.TestCase):
    def test_production_core_and_outer_radii_are_distinct(self):
        # DREAMSim 54b2447, B4DetectorConstruction.cc: core_C_rMax=.39;
        # clad_C_rMax=.40 mm. These are independent source constants.
        c, _, _ = dreamsim_front_contract()
        self.assertEqual(c.core_radius_mm, .39)
        self.assertEqual(c.fiber_outer_radius_mm, .40)

    def test_production_contract_rejects_outer_radius_as_core(self):
        with self.assertRaisesRegex(ValueError, "core"):
            dreamsim_front_contract(core_radius_mm=.4)

    def test_source_geometry_roundtrip_all_fibers_and_edges(self):
        for rx, ry in [(-45, -40), (-25, -25), (0, 0), (15, 25), (44, 39)]:
            for f in range(4):
                for dx, dy in [(0, 0), (.39, 0), (0, -.39), (.2, .2)]:
                    calo = np.array([[2+4*rx+FIBER_OFFSETS_MM[f, 0]+dx, 2+4*ry+FIBER_OFFSETS_MM[f, 1]+dy, 1000.]])
                    world = calo @ rotation()
                    local, _, index, _, back = localize_endpoints(world, np.array([f]))
                    np.testing.assert_allclose(local, [[dx, dy]], atol=1e-10)
                    np.testing.assert_equal(index, [[rx, ry]])
                    np.testing.assert_allclose(back, calo, atol=1e-10)

    def test_rotated_front_inverse_matches_world_flight(self):
        c, tface, face = dreamsim_front_contract()
        r = rotation(); self.assertGreater(tface, 0)
        for z in [0, 100, 1000, 2000]:
            xyz = np.array([20., -12., z-1000])
            world = r.T @ xyz
            front_time = (world[2]-face[2])/C_MM_NS
            t = front_time+c.n_group*1.02*(2000-z)/C_MM_NS
            self.assertAlmostEqual(float(c.inverse(t, 1.02, xyz[0], xyz[1])), z, places=9)
        np.testing.assert_allclose(world_to_calorimeter(np.array([face]))[0, 2], -1000., atol=1e-9)

    def test_raw_adapter_requires_slurm(self):
        import os
        from unittest.mock import patch
        with patch.dict(os.environ, {}, clear=True):
            with self.assertRaisesRegex(RuntimeError, "Slurm"):
                next(iter_corrected_root("nonexistent", "a"*64, DepthContract()))


class IdealTransportTests(unittest.TestCase):
    @staticmethod
    def reference(pos, u, distance, radius):
        pos, u = np.array(pos, float), np.array(u, float)
        for count in range(10000):
            dot = float(pos @ u)
            wall = -dot+np.sqrt(max(0., dot*dot+radius*radius-float(pos @ pos)))
            if distance <= wall+1e-12:
                return pos+distance*u, count
            pos += wall*u; distance -= wall
            u -= 2*float(u @ pos)*pos/radius**2
        raise AssertionError("Independent reflection reference did not finish")

    def test_zero_impact_two_reflections(self):
        d = np.array([[.1, 0, np.sqrt(.99)]])
        length = 1.2*d[0, 2]/d[0, 0]
        r = trace_cylinder([[.1, 0, 0]], d, length_mm=length, radius_mm=.39, n_phase=1.504)
        np.testing.assert_allclose(r.exit_xy_mm[0], [-.26, 0], atol=1e-12)
        ref, nr = self.reference([.1, 0], [1, 0], 1.2, .39)
        np.testing.assert_allclose(r.exit_xy_mm[0], ref, atol=1e-10)
        self.assertEqual(r.reflections[0], nr)

    def test_nonzero_impact_matches_independent_reflections(self):
        rng = np.random.default_rng(2)
        for _ in range(64):
            radius, phi = np.sqrt(rng.random())*.35, rng.uniform(0, 2*np.pi)
            pos = [radius*np.cos(phi), radius*np.sin(phi)]
            theta = rng.uniform(0, 2*np.pi); u = np.array([np.cos(theta), np.sin(theta)])
            pt = rng.uniform(.03, .15); dz = np.sqrt(1-pt*pt); length = rng.uniform(1, 25)
            r = trace_cylinder([pos+[0]], [list(pt*u)+[dz]], length_mm=length)
            ref, count = self.reference(pos, u, length*pt/dz, .4)
            np.testing.assert_allclose(r.exit_xy_mm[0], ref, atol=1e-8)
            self.assertEqual(r.reflections[0], count)

    def test_face_before_wall_and_axial_and_backward(self):
        r = trace_cylinder([[0, 0, 1999.999], [0, 0, 0], [0, 0, 100]],
            [[np.sqrt(.75), 0, .5], [0, 0, 1], [0, 0, -1]])
        np.testing.assert_equal(r.status, [RayStatus.DIRECT_FACE, RayStatus.DIRECT_FACE, RayStatus.BACKWARD])
        self.assertTrue(np.isnan(r.arrival_ns[2])); self.assertEqual(r.reflections[0], 0)

    def test_non_guided_wall_and_invalid_source(self):
        r = trace_cylinder([[0, 0, 0]], [[np.sqrt(.75), 0, .5]])
        self.assertEqual(r.status[0], RayStatus.ESCAPED_WALL)
        self.assertTrue(np.isnan(r.arrival_ns[0]))
        with self.assertRaises(ValueError): trace_cylinder([[.5, 0, 0]], [[0, 0, 1]])


REPO = Path(__file__).resolve().parents[1]


def fixture(prefix, events=4, n=120):
    """Synthetic unit-test arrays, NOT generated Geant4/shower events."""
    rng = np.random.default_rng(8)
    c = DepthContract()
    keys = np.repeat([f"{prefix}:{i}" for i in range(events)], n)
    ids = np.tile(np.arange(n, dtype=np.int64), events)
    r = rng.uniform(.05, .8, len(keys))
    k = 1.01+.015*r+.005*rng.random(len(keys))
    x = rng.integers(-2, 3, len(keys)); y = rng.integers(-2, 3, len(keys))
    # Two separated populations test that the method returns a candidate, not
    # the unoccupied midpoint between depth modes.
    z = np.where(np.arange(len(keys)) % 2, 450., 1450.)+rng.normal(0, 15, len(keys))
    obs = PhotonObservations(keys, ids, x, y, c.forward_time(z, k), r, 4*x.astype(float), 4*y.astype(float))
    return obs, c, z, k


class DepthPipelineTests(unittest.TestCase):
    def test_cli_calibration_and_all_resolution_methods(self):
        with TemporaryDirectory() as d:
            d = Path(d)
            a, c, z, k = fixture("train")
            b, _, zb, _ = fixture("test")
            save_cache(d/"train", a, c, provenance={"truth_validation": "synthetic fixture exact kappa", "timestamp_sigma_ns": 0},
                       truth_z_mm=z, kappa_truth=k, calibration_eligible=True)
            save_cache(d/"test", b, c, provenance={"timestamp_sigma_ns": 0}, truth_z_mm=zb)
            subprocess.run([sys.executable, str(REPO/"studyDepthReconstruction.py"), "calibrate", "--input", str(d/"train"),
                            "--output", str(d/"cal.json")], check=True, capture_output=True, text=True)
            cmd = [sys.executable, str(REPO/"studyDepthReconstruction.py"), "scan", "--input", str(d/"test"),
                   "--calibration", str(d/"cal.json"), "--output", str(d/"scan"), "--bootstrap-draws", "8"]
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            report = json.loads((d/"scan/summary.json").read_text())
            self.assertEqual(len(report["results"]), 12)
            for r in report["results"]:
                self.assertTrue(r["evaluation"]["closure"])
                self.assertEqual(r["evaluation"]["input"], len(zb))
                matrix = r["evaluation"]["response"]
                self.assertTrue(matrix["closure"])
                self.assertEqual(matrix["input"], len(zb))
                with np.load(d/"scan"/(r["name"]+".npz"), allow_pickle=False) as f:
                    np.testing.assert_equal(f["photon_id"], b.photon_id)
                    np.testing.assert_equal(f["event_key"], b.event_key)
            self.assertNotEqual(subprocess.run(cmd, capture_output=True).returncode, 0)

    def test_plot_rejects_cache_and_result_mismatches(self):
        from dsipm.io import file_sha256
        with TemporaryDirectory() as d:
            d = Path(d)
            a, c, z, _ = fixture("eval", events=2, n=8)
            save_cache(d/"cache", a, c, provenance={"timestamp_sigma_ns": 0}, truth_z_mm=z)
            _, _, _, meta = load_cache(d/"cache")
            (d/"run").mkdir()
            record = d/"run/radius_50ps_draw0.npz"
            np.savez_compressed(record, event_key=a.event_key, photon_id=a.photon_id, z_mm=z, status=np.zeros(len(z),int))
            row = dict(name=record.stem, sha256=file_sha256(record), options=dict(method="radius",sigma_t_ns=.05))
            summary = dict(source_cache=meta, results=[row])
            cmd = [sys.executable,str(REPO/"studyDepthReconstruction.py"),"plot","--input",str(d/"cache"),
                   "--run",str(d/"run"),"--output",str(d/"plots"),"--label","synthetic test"]
            summary["source_cache"] = dict(meta, events=["wrong"])
            (d/"run/summary.json").write_text(json.dumps(summary))
            self.assertIn("does not match",subprocess.run(cmd,capture_output=True,text=True).stderr)
            summary["source_cache"] = meta
            row["sha256"] = "0"*64
            (d/"run/summary.json").write_text(json.dumps(summary))
            self.assertIn("hash mismatch",subprocess.run(cmd,capture_output=True,text=True).stderr)
            np.savez_compressed(record, event_key=a.event_key, photon_id=a.photon_id[::-1], z_mm=z, status=np.zeros(len(z),int))
            row["sha256"] = file_sha256(record)
            (d/"run/summary.json").write_text(json.dumps(summary))
            self.assertIn("alignment",subprocess.run(cmd,capture_output=True,text=True).stderr)
            self.assertFalse((d/"plots").exists())

    def test_plot_cli_writes_three_individual_figures(self):
        import importlib.util
        if importlib.util.find_spec("matplotlib") is None:
            self.skipTest("Optional matplotlib not installed")
        from dsipm.io import file_sha256
        with TemporaryDirectory() as d:
            d = Path(d)
            a, c, z, _ = fixture("eval",events=2,n=8)
            save_cache(d/"cache",a,c,provenance={"timestamp_sigma_ns":0},truth_z_mm=z)
            _, _, _, meta = load_cache(d/"cache")
            (d/"run").mkdir()
            result=d/"run/radius_50ps_draw0.npz"
            np.savez_compressed(result,event_key=a.event_key,photon_id=a.photon_id,z_mm=z,status=np.zeros(len(z),int))
            (d/"run/summary.json").write_text(json.dumps(dict(source_cache=meta,results=[dict(name=result.stem,
                sha256=file_sha256(result),options=dict(method="radius",sigma_t_ns=.05))])))
            cmd=[sys.executable,str(REPO/"studyDepthReconstruction.py"),"plot","--input",str(d/"cache"),
                 "--run",str(d/"run"),"--output",str(d/"plots"),"--label","synthetic test","--sigma-ps","50"]
            subprocess.run(cmd,check=True,capture_output=True,text=True)
            self.assertEqual(len(list((d/"plots").rglob("*.pdf"))),3)
            self.assertEqual(len(list((d/"plots").rglob("*.png"))),3)
            self.assertNotEqual(subprocess.run(cmd,capture_output=True).returncode,0)

    def test_unvalidated_truth_cannot_be_trained_by_cli(self):
        with TemporaryDirectory() as d:
            d = Path(d); a, c, z, k = fixture("x")
            save_cache(d/"input", a, c, provenance={"timestamp_sigma_ns": 0}, truth_z_mm=z, kappa_truth=k)
            p = subprocess.run([sys.executable, str(REPO/"studyDepthReconstruction.py"), "calibrate", "--input", str(d/"input"),
                                "--output", str(d/"cal.json")], capture_output=True, text=True)
            self.assertNotEqual(p.returncode, 0); self.assertIn("not validated", p.stderr)
            self.assertFalse((d/"cal.json").exists())

    def test_cache_corruption_rejected(self):
        with TemporaryDirectory() as d:
            d = Path(d); a, c, z, k = fixture("a")
            save_cache(d/"input", a, c, provenance={})
            with (d/"input/photons.npz").open("ab") as f:
                f.write(b"changed")
            with self.assertRaisesRegex(ValueError, "hash mismatch"): load_cache(d/"input")

    def test_two_depth_populations_do_not_force_a_mean_between_modes(self):
        a, c, z, k = fixture("cal")
        b, _, zb, _ = fixture("eval")
        cal = fit_calibration(a, k, c, source_sha256="b"*64)
        r = reconstruct(b, cal, c)
        # Loose structural regression, NOT a detector-performance threshold.
        self.assertLess(np.mean(np.isfinite(r.z_mm) & (r.z_mm > 800) & (r.z_mm < 1100)), .05)

    def test_raw_corrected_adapter_geometry_and_observable_separation(self):
        c, tface, _ = dreamsim_front_contract(pde=1.)
        f = np.arange(4, dtype=np.int32)
        calo = np.column_stack([2+FIBER_OFFSETS_MM[:, 0], 2+FIBER_OFFSETS_MM[:, 1], np.full(4, 1000.)])
        world = calo @ rotation()/10
        birth = calo.copy(); birth[:, 2] = -500.
        birth_world = birth @ rotation()/10
        mom = np.tile([0., 0., 1.], (4, 1)) @ rotation()
        data = dict(OP_time_final=c.forward_time(np.full(4, 500.), 1., calo[:, 0], calo[:, 1])+tface,
            OP_time_produced=np.zeros(4), OP_productionFiber=f, OP_finalFiber=f.copy(), OP_isCoreC=np.ones(4, bool))
        for j, axis in enumerate("xyz"):
            data["OP_pos_final_"+axis] = world[:, j]
            data["OP_pos_produced_"+axis] = birth_world[:, j]
            data["OP_mom_final_"+axis] = mom[:, j]
        e, ids, diagnostic, counts = corrected_event_from_arrays(data, c)
        np.testing.assert_allclose(e.local_x_mm, 0, atol=1e-10)
        np.testing.assert_allclose(diagnostic["truth_z_mm"], 500, atol=1e-10)
        np.testing.assert_allclose(diagnostic["exit_direction_kappa_proxy"], 1, atol=1e-10)
        batch = make_observations(e, "hash:0", c, entrance_time_ns=tface, photon_ids=ids)
        self.assertEqual(counts["selected"], 4); self.assertEqual(len(batch.photons.time_ns), 4)
        data["OP_finalFiber"][:] = -99
        _, _, _, unknown = corrected_event_from_arrays(data, c)
        self.assertEqual(unknown["unavailable_final_fiber"], 4)
        data["OP_finalFiber"][0] = 3
        with self.assertRaisesRegex(ValueError, "fiber mismatch"): corrected_event_from_arrays(data, c)


if __name__ == '__main__':
    unittest.main()
