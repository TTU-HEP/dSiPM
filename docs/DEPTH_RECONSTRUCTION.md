# Fiber-depth reconstruction: development integration

This is a new implementation, not a replacement of historical reductions.
It uses existing photons; there is no Geant4 generator or production submission
in the workflow. Requires Python >=3.11 and NumPy. ROOT reading additionally
requires uproot and awkward, and is Slurm-only. No SciPy is required.

## Forward model and what is actually inferred

Depth z is in mm from the entrance face; timestamps are ns relative to an
entrance timing layer. For a prompt beta=1 production front,

```
c (t_arrival - t_face) = a_z z + a_x x + a_y y + b + n_group kappa (L-z)
z = [n_group kappa L + a_x x + a_y y + b - c t] / (n_group kappa - a_z)
```

kappa is optical path / axial distance, not measured by the sensor. Radius is
the exit pixel centre's distance from the known physical fiber centre,
normalized by the core radius. It constrains possible paths, not a unique path.
The ideal TIR bound uses the **phase** indices, while flight time uses the
**group** index. Calibration candidates incompatible with this bound are
excluded, not projected onto it. Near-face direct exits are allowed without
TIR when their transverse distance can fit an unreflected chord. This is a
necessary geometric feasibility test, not knowledge of the exit azimuth.

The default n_phase=n_group=1.49 and n_clad=1.42 are an explicitly constant-index
development assumption matching the nominal pinned clear-fiber material, NOT
the supplied simplified study's n=1.504. Calibrations cannot be mixed between
them. Default core radius is 0.4 mm, length 2000 mm. These are not parameters
automatically inferred from a ROOT filename.

The full emission model has an additional delta_emit:

```
t_emit - t_face = t_prompt_front(z,x,y) + delta_emit
t_arrival - t_face = t_prompt_front + delta_emit + n_group kappa (L-z)/c
```

The current inverse assumes delta_emit=0. The saved first-step production time
is **not** guaranteed to be the true photon-birth timestamp: pinned DREAMSim
records a pre-step position with track global time in the first stepping action.
Keep that diagnostic; do not use it as an exact truth correction. End-face exit
momentum is only a direction proxy (including possible refraction), NOT an
automatically validated in-core optical path. The unresolved transport problem
cannot be fixed or certified by training lookup tables on its own delayed output.

## Estimators and failure accounting

1. `radius`: median calibration kappa conditional on exit radius.
2. `radius_delay`: median conditional on radius and measured delay to the first
   observed photon in the same rod/event. Exact zero has a distinct delay bin.
3. `neighbour`: ten equal-mass midpoint kappa quantiles conditional on those
   observables; choose the physical candidate with largest neighbouring clean-
   light density. No additional calibration density is multiplied into those
   equal-probability candidates. Default 5x5 rods, 20 mm depth bins and Gaussian
   smoothing at least 30 mm; these are **untuned development defaults**.

Only finite, physical clean-light inverses seed the density. A negative depth
cannot vote at the front face. The identical finite Gaussian kernel is used
for smoothing and leave-one-photon-out subtraction. `exclude_same_rod=True`
is an additional leave-one-rod-out ablation to expose within-rod correlations.
All context comes from the same event; rods use signed tuple identities, not
collision-prone packed multiplication. No smoothing in truth space is used.

No support: return the **actual calibrated median**, not the 55th-percentile
candidate. If that median is invalid but a physical quantile exists, choose
the nearest quantile in kappa. No physical estimate: z=NaN, status=2. A fallback
is status=1; positive neighbour support is status=0. Positive support does not
imply a correct reconstruction or calibrated confidence. Sub-nanometre floating-
point endpoint roundoff is normalized; no unphysical depth is clipped into the
detector. Multimodal candidates are retained; no forced mean between deposits.

For comparison with the older fixed-path inverse, add `--fixed-n-eff VALUE` to
the scan (API: `method='fixed', fixed_n_eff=VALUE`). The value must be frozen
outside the test sample; no implicit universal effective index is supplied.
This baseline uses the same entrance/rotation model and physical depth check,
but not the radius/TIR constraint. It is still a per-photon baseline, not the
older earliest-photon-per-rod observable; do not compare unlike truth targets.

`photon_timing_sigma_z_mm` and `entrance_timing_sigma_z_mm` report local linear
propagation c*sigma_t/(n_group*kappa-a_z), **conditional on the chosen kappa**.
They omit path/emission/calibration ambiguity and candidate switching, and are
not total uncertainties or calibrated confidence intervals. Entrance-clock
error is common to all photons of one event and must not average as 1/sqrt(N).

## Geometry and sensor integration

Use `depth_geometry.dreamsim_front_contract()` for DREAMSim 54b2447, two-degree
x/y passive placement rotations, world-z beam at (25,-25,-1050) mm. It returns
the contract, true entrance time relative to the generation clock, and entrance
point. World-to-calorimeter is Ry @ Rx. Depth is calo_z+1000, without the
historical factor-of-two error. Physical fiber offsets, 4 mm rod pitch, replica
limits and the downstream plane are checked. Sensor-centre x/y approximates
sub-fiber front-time position (<~0.1 ps at this angle); it is not truth x/y.

`make_observations` accepts a corrected FullDetectorPhotonEvent-compatible
object. It never reads its truth-depth field. Four fibers have independent
SPAD holdoff domains; neighbours are rods, not shared HG-DREAM optical channels.
PDE, strict stored-time gate, pixel grid and nonparalyzable holdoff are applied
before jitter. The default pixel origin matches the historical 3 mm sensor
grid. Event keys must contain source identity and original entry, e.g.
`<ROOT_SHA256>:<entry>`; original photon indices are preserved.

The SHA-keyed counter streams make PDE, photon jitter and event-clock jitter
independent and invariant to input order/subsets/batching. 0, 50, 200 and 500 ps
scans have exactly the same accepted identities at fixed pitch. Increasing
pitch reruns holdoff but shares the per-photon PDE draws. This new keyed response
is not claimed bit-identical to the historical NumPy RNG; an existing validated
mask can instead be supplied. Its detected count remains unknown unless the
detected mask is supplied too. Do not silently relabel old response results.

For canonical sensor provenance, construct `DepthContract.from_workflow()`;
this reads PDE/pitch/holdoff from `sensor_workflow.json` and pins its SHA-256.
For rotated production geometry, pass those fields and that workflow hash to
`dreamsim_front_contract`. A changed workflow hash invalidates compatibility
even if a particular sensor parameter stayed the same. Bare `DepthContract()`
is a standalone assumption, not proof that the canonical workflow was loaded.

The observations in this first integration are **accepted-SPAD timestamps and
pixel positions**, before shared-channel timestamp export. They assume that
location information is available at that scale. Actual timestamp buffers,
80-micrometre grouping, resolving time, count records and optical coupling are
not simulated here. A workflow checksum does not mean those electronics were
applied. Evaluate the exported-record model before describing the result as
achievable FLASH hardware performance; coarser measured position must be an
explicitly different observation/sensor contract, not substituted silently.

The optional raw reader `iter_corrected_root` verifies the complete ROOT SHA-256
before and after reading and requires SLURM_JOB_ID. It uses a new reader, not
a monkeypatch of any historical module. It preserves the known exitFiber=-99
sentinel as unavailable. Actual non-sentinel fiber mismatches, production-fiber
containment failures and wrong geometry raise errors. Exhaust the iterator to
complete its after-read hash check. The reader exposes recorded first-step time
and exit-direction kappa ONLY as separate diagnostics; neither is automatically
eligible calibration truth. Keep the validated sidecar/production manifest
check upstream of this reader and archive it in the cache provenance.

## Calibration, caching and reproducible scans

`fit_calibration(observations, kappa_truth, contract, source_sha256=...)` is the
ONLY estimator API taking truth. Its tables retain training event IDs, bin
counts, global/radius fallback levels, training timestamp sigma and checksums.
Calibration saves are immutable and hash-checked. Inference rejects overlapping
training/evaluation events and any optical, geometry, gate, sensor or holdoff
contract mismatch. `allow_training=True` is a labelled training-only diagnostic,
never held-out evaluation. Hold out tuning events too; do not choose windows,
cuts or examples using the final test score.

`save_cache` writes a new directory with `photons.npz` and `metadata.json`.
`PhotonObservations` contains NO truth. Optional truth columns remain separate
in the cache. A calibration cache must be explicitly `calibration_eligible=True`
with documented `provenance['truth_validation']`; the CLI rejects a direction
proxy or unvalidated cache. This marker is a scientific audit decision, not a
test that the software can independently prove.

Example after making separately validated `calibration_cache/` and
`heldout_cache/` (timestamp_sigma_ns=0 in provenance):

```bash
python studyDepthReconstruction.py calibrate --input calibration_cache --output calibration.json
python studyDepthReconstruction.py scan --input heldout_cache --calibration calibration.json \
  --output results/depth-fixed-cal --sigmas-ps 0 50 200 500 --clock-ps 0
```

This is a **fixed zero-jitter calibration** scan; do not describe it as separately
optimized at every resolution. Use `--draws` for matched noise repetitions and
`--leave-rod-out` for the correlation ablation. Repeated draws are not independent
showers and their bootstrap intervals must not be combined as independent data.
Outputs retain per-photon IDs, candidates, statuses, conditional timing scales,
code/input/calibration hashes and complete truth-depth coverage counts. NumPy
arrays use no pickle. Existing output directories are never overwritten.

Evaluation reports bias/68% width conditional on reconstruction, the inclusive
fraction within 5/10/20 cm of truth, explicit failures, depth-bin losses and
event-block bootstrap intervals. Inclusive here means **every accepted photon**,
not every generated photon or all deposited energy. Report gate/PDE/holdoff
counts alongside it. Channel medians require >=5 accepted AND reconstructed
photons and compare to accepted-photon truth medians; they do not mean the
pre-holdoff population or a resolved energy-deposition coordinate.

## Checks and release gates

```bash
python -m unittest discover -s tests -v
```

The suite covers the reviewed defects, all 80 source-based geometry cases,
world-front timing, event-common clock error, matched response, sparse and
two-depth inputs, artifact compatibility, end-to-end scans and histogram/failure
closure. An independent reflection-by-reflection reference checks 64 off-axis
rays plus the zero-impact two-reflection case. Direct face arrival before first
wall, axial rays and nonguided losses are explicit. This ideal-cylinder library
does not repair/replace stored Geant4 photons.

Required before a detector-performance claim: resolve the production transport
and exact creation-time semantics, define trustworthy path-factor calibration,
freeze disjoint calibration/tuning/test manifests, run a small existing-data
Slurm pilot, compare each ablation and geometry/model transfer by species/energy,
and measure coverage/tails/failures with event/file-block uncertainties. No
broad new event generation is needed to test this implementation first.

## Method provenance

Concepts adapted from the user-supplied `fiber_sim_study` and companion
*Timing-Based Depth Reconstruction in Fiber Calorimeters: Method and Performance*
(received 2026-10-21). The source package and its summary numbers were preserved
and reviewed separately. This is a new independently tested implementation, not
an exact replay of its absent photon caches or a claim to reproduce its 10.8 cm.
Reviewed source SHA-256 values:

- `depth_reco.py`: `599f9b182a762f3b702895d10bb7fb2837b974f780986020aeaa1874e6ccceab`
- `fibers.py`: `6722987846af3c6ab390aae98da03a2ded944519bc0edb7245f71893ad1234f8`

Historical f12de72/schema-v3 checkpoints, production commits, final40000
correction release and paper inputs remain unchanged. The production scope
remains central-window Cherenkov, not all-fiber optical response. No timing
resolution, explicit scintillation, installed HG-DREAM coupling or prototype
performance is established by the software tests.
