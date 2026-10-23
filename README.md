# DREAM / FLASH photon analysis

A small, NumPy-based core for sensor response and timing-based fiber-depth
reconstruction. It reads existing photons; it does not generate new showers.
The original main-branch scripts are left untouched as legacy references.
For new upstream DREAMSim ROOT files, the detector definition comes from the
producer's committed source and actual run macro, not the historical CPAD
constants or an independently tuned dSiPM geometry.

## Core

```text
dsipm/
  geometry.py   corrected coordinates, DREAM/HG-DREAM layouts, ideal ray reference
  sensors.py    FLASH holdoff, analog recharge, matched clocks, readout/pulse proxies
  depth.py      model contract, calibration, reconstruction, truth-based evaluation
  io.py         immutable caches and hash-checked, Slurm-only existing-ROOT reader
studyDepthReconstruction.py   compact-cache entry point: calibrate / scan / plot
reduceDepthCampaign.py       Slurm-only complete-campaign baseline reduction / merge
compareDepthCorrections.py   pion-only, held-out depth methods and empirical controls
sensor_workflow.json         one canonical sensor configuration
tests/                      focused physics and end-to-end regressions
```

Calibration/evaluation truth stays outside the measured-observation interface,
even though related functions now live together. Imports do not load ROOT,
generate events, submit jobs, or produce plots.

```bash
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python studyDepthReconstruction.py --help
```

The optional ROOT reader additionally needs `uproot` and `awkward`. Python
3.11 or newer is required. NumPy is the only required dependency for the core.

## Reconstruction

z is depth in mm from the entrance face. t is arrival time in ns relative to
the entrance timing layer. The prompt beta=1 production-front approximation is

```text
c t = a_z z + a_x x + a_y y + b + n_group kappa (L-z)
z = (n_group kappa L + a_x x + a_y y + b - c t) / (n_group kappa - a_z)
```

kappa is optical path / axial distance, not a measured quantity. Phase indices
set the ideal guidance bound; the independent group index sets flight time.
Standalone/historical defaults are L=2000 mm, core/outer radii=0.39/0.40 mm and
indices 1.49/1.42. They are **not** the current upstream detector. There is no implicit universal
effective index for the optional historical fixed-path comparison.

The core and cladding boundary must remain distinct. The source audit corrected
the earlier 0.40 mm core assumption; the v2 algorithm/contract rejects old
calibration tables rather than silently reusing their radius normalization.

`DepthContract.from_workflow()` loads canonical PDE/pitch/holdoff and pins the
complete configuration hash. It also accepts a workflow JSON **path**. For the
rotated production geometry, use `geometry.dreamsim_front_contract`, passing
the workflow's sensor values/hash; this returns the contract, entrance time
and entrance point. World-to-calorimeter is Ry @ Rx for the pinned two-degree
rotations. Depth is calo_z+1000, not the historical doubled scale. Coordinates,
fiber containment, replica limits and the downstream plane are checked.

## Current upstream detector: DREAMSim is the authority

Use `geometry.DreamsimDetector.from_checkout(checkout, actual_run_macro)` for
new upstream files. It reads the committed detector/material/placement source
and run parameters, pins the full source revision and file hashes, and rejects
uncommitted source changes, overlapping fibers, unsupported source forms,
nonconstant refractive-index tables and unsupported beam configurations.
It does not execute or modify C++ or silently fit missing group indices.

The upstream model currently has a 2500 mm length, 90 rods by 80 layers at
4 mm pitch, 0.39/0.40 mm core/outer radii, two plastic Cherenkov fibers with
constant core index 1.504, two quartz fibers with index 1.468, and three
scintillation fibers. The Cherenkov cladding index is 1.42. Group index equals
phase index **only because these producer tables are verified constant**.
Scintillation reconstruction is not implemented or claimed.

`detector.contract('plastic')` and `detector.contract('quartz')` return separate
contracts, entrance clocks and entrance points. Length, indices, fiber identity,
replica dimensions, rotation, beam front and explicit sample-rod scope are
source-derived. Sensor PDE/pitch/recovery may come from `sensor_workflow.json`;
its historical geometry block never overrides this producer definition. Do
not change the locked historical workflow or reuse its calibrations as an
upstream calibration. Source profile/material/gate/workflow changes invalidate
calibration compatibility automatically.

Upstream `020d84c` contains overlapping quartz/scintillation placements. The
reader rejects that source rather than compensating in Python. The isolated
DREAMSim correction restores the quartz copy-2 y sign to the slot indicated
by the adjacent commented placement and previous geometry: **one line in one
C++ file**. No detector redesign, tracking-window enlargement or transport
rewrite is part of this migration.

The upstream default explicitly tracks one configured rod/layer. With the
whole-detector-capable producer, set **both** `opSampleRod` and `opSampleLayer`
to **-1** in the actual run macro to track optical photons in all 90 x 80 rods.
The source parser verifies that this producer implements the opt-in; the old
single-rod source cannot be relabelled as whole-detector data. Mixed negative
indices are invalid. Scope changes invalidate calibration compatibility.
`detector.tracks_all_rods` records this source-derived distinction.

Generate the complete event with all rods tracked together for a whole-detector
Z study. The neighbour estimator then receives observed light from adjacent
rods in the **same event**; combining unrelated single-rod simulations is not an
equivalent test. Plastic and quartz remain separate calibration/evaluation
populations. Sensor holdoff remains per fiber, not one shared detector-wide SPAD.

In single-rod runs, other optical records carry birth/analytical information,
not measured sensor arrivals. The reader uses **only recorded downstream Cherenkov-core exits**, never
`OP_analyticalArrivalTime_*`. It separates plastic and quartz, cross-checks
endpoint localization against exported replica IDs, preserves original photon
IDs, and accounts separately for unrecorded, non-downstream, outside-gate and
selected records. Single-rod products must not be described as all-detector
optical performance. Whole-detector transport is a simulation scope, not proof
of accurate depth inference or validated end-face/sensor coupling. Benchmark a
small run before scaling: all-rod optical transport can be substantially slower.

Run extraction under Slurm with a producer ledger binding `root_sha256`,
`producer_revision`, `detector_sha256` and `run_macro_sha256` to the actual
source and run configuration. Any command-line geometry/beam/sampling overrides
must be reflected in that run macro/profile. Current upstream does not export
rotation branches, so the rotation is pinned by this producer ledger and run
macro, not falsely reported as independently measured in each ROOT event.
Exported beam position/direction are checked; rotation branches are checked too
if a future schema supplies them.

```bash
python studyDepthReconstruction.py extract \
  --input fresh.root --sha256 ROOT_SHA256 \
  --dreamsim /path/to/clean/DREAMSim --macro actual_run.mac \
  --provenance producer.json --material plastic --output plastic_cache
# Repeat with --material quartz and a different output directory.
```

The resulting immutable caches feed the existing `calibrate` and `scan` tools
without a new estimator implementation. Truth remains separate from measured
observations. Optical-path truth is not declared validated: a source migration
and a geometry fix do not establish that late-photon transport is correct.
The old `reduceDepthCampaign.py` CPAD reduction and old fixed-2000-mm comparison
plans remain historical entry points; do not feed new-model ROOTs into them.

## Depth estimators

Three estimators share exactly the same accepted photons:

- `radius`: calibrated median kappa conditional on pixel-centre exit radius.
- `radius_delay`: also conditions on measured delay to the first photon of that
  rod/event; exact zero has its own bin.
- `neighbour`: selects among equal-probability kappa quantiles using same-event
  neighbouring clean light; retains candidates for ambiguity diagnostics.

Only physical inverses seed the neighbour density. Smoothing and self-photon
subtraction use the same kernel. Leave-one-rod-out is an explicit ablation.
No support uses the actual calibrated 50% median, not the 55% candidate. If
there is no physical solution, z is NaN with an unreconstructed status; no
unphysical distance is clipped into the detector. A tiny endpoint-roundoff
tolerance is not a boundary-repair procedure. Near-face direct exits do not
require TIR before a wall they never encounter.

`fit_calibration(observations, kappa_truth, contract, source_sha256=...)` fits
validated geometric path factors. Separately, `fit_timing_calibration` (or its
streaming batch variant) fits an **effective timing factor** from measured
arrival time and geometry-checked production Z, using the exact inverse model:

```text
kappa_eff = (c*t - a_z*z_true - a_x*x - a_y*y - b) / (n_group*(L-z_true))
```

This needs neither photon direction nor exact creation time. It absorbs emission
front mismatch, mixed propagation speeds and timestamp noise; it is **not**
geometric path truth or a repair of transport physics. Finite factors below one
(even negative factors under noise) are retained. Undefined Z=L calibration
labels are counted, never used to remove evaluation photons. The geometric
TIR guidance remains the default for geometric calibrations; timing-effective
models default to physical-Z-range checks only. An explicit geometric-guidance
variant tests that upper-bound heuristic without asserting it is a physical
constraint on an effective coefficient.

Inference rejects overlapping calibration
and evaluation events and mismatched geometry, optics, workflow, gate or sensor
models. Photon jitter and entrance-clock jitter must match the calibration's
recorded training values. Only an explicit `allow_jitter_mismatch` override
permits a fixed-calibration robustness study. Hold out tuning events as well;
never choose cuts/windows on final test
scores. Defaults are development assumptions, not an optimized performance.

The original `compareDepthCorrections.py --compare` study directly calibrated conditional medians of geometry-checked
production Z, not kappa: measured time, then pixel exit radius, same-rod arrival
delay, and other-rod timing in a 5x5 window. Five source-file folds keep every
evaluation event out of its model's training. Feature bins/support thresholds
are fixed before evaluation; all methods retain identical accepted photons and
matched 0/50/200/500 ps jitter, with fixed zero-jitter training. Sparse cells
explicitly fall back to less detailed models. A no-timing depth-prior control
is mandatory: direct conditional estimates are physical by construction, not
because they solve a validated optical inverse. These are sample-specific,
empirical development comparisons, not proof of optical-path correction or
energy-deposition resolution. Paired file-bootstrap errors condition on the
fitted crossfit maps; they exclude refit and transport/model uncertainties.
It is an empirical control, **not** a test of corrections 1-3. The new
`--method-task` / `--method-merge` comparison invokes the actual `reconstruct`
radius, radius-delay and neighbour methods on the same 5,000 pion events. It
refits both the effective-factor calibration and empirical controls separately
at each timing sigma, with five whole-source-file holdouts. It also reports a
global-factor reference, range-only versus geometric-guidance ablations, and a
no-timing prior. Radius bins resolve 0.9, 0.95 and 0.975 and retain outer-fiber
pixel observations. Grouped sorting replaces repeated full-array bin masks.

Headline metrics are failure-inclusive fractions within 5/10/20 cm and the
absolute-error tolerance containing 68% of all accepted photons. Failed returns
have infinite error in this CDF; that tolerance is undefined/unreachable if
too many photons fail. A failed estimate has no signed residual or finite RMS
without an explicitly declared penalty. Conditional residual widths remain
labeled diagnostics, never substitutes for all-photon performance.

## Data and matched scans

`sensors.make_observations` accepts a corrected FullDetectorPhotonEvent-compatible
object. Only measured time, local endpoint, sensor/rod identity and sensor-centre
coordinates are read. Use `<source_SHA256>:<original_entry>` event keys and keep
original photon indices. No inference function reads true z or emission time.

PDE, stored-time gate and per-fiber SPAD holdoff precede timestamp smearing.
PDE, photon jitter and the event-common entrance-clock error use separate keyed
streams. The same photons and standard-normal draws are used at 0/50/200/500 ps;
input order and batching cannot change selection. Existing validated masks may
be reused, with unknown detection/loss counts left unknown unless supplied.

`io.save_cache` creates a new directory with `photons.npz` and `metadata.json`.
Truth columns are separate. A calibration cache must explicitly be marked
eligible with documented `provenance['truth_validation']`. That is a scientific
audit decision, not proof the software can establish automatically. Do not mark
an exit-direction proxy or transport-bug-derived path factor as validated truth.
The raw reader requires Slurm and checks complete ROOT hashes before/after an
exhausted read. Check the validated sidecar/production manifest upstream too.

After preparing disjoint caches with `timestamp_sigma_ns=0` and documented
`truth_z_validation` in the calibration cache provenance:

```bash
for timing_sigma_ps in 0 50 200 500; do
  python studyDepthReconstruction.py calibrate --input calibration_cache \
    --factor-kind timing_effective --sigma-ps "$timing_sigma_ps" --output "cal_${timing_sigma_ps}.json"
done
python studyDepthReconstruction.py scan --input heldout_cache --calibration cal_0.json cal_50.json cal_200.json cal_500.json \
  --output results/depth-scan --sigmas-ps 0 50 200 500 --clock-ps 0
```

This uses one fixed zero-jitter calibration, not retraining at each resolution.
`--draws` adds matched noise repetitions; they are not independent showers.
`--leave-rod-out` tests within-rod correlations. `--fixed-n-eff VALUE` adds an
explicitly frozen fixed-path baseline. No existing output is overwritten.

Scoring retains failures in the accepted-photon denominator, reports conditional
bias/68% width and inclusive fractions within 5/10/20 cm, and uses event-block
bootstrap intervals with depth-bin underflow/overflow accounting. Gate/PDE/
holdoff losses are separate. Channel medians are accepted-photon truth medians,
not energy-deposition truth. Local timing-only z-error contributions exclude
path/emission ambiguity and are not total calibrated uncertainty intervals.

Every evaluated scan also saves a true-Z/reconstructed-Z response matrix with
explicit failure counts in each truth bin. True Z is horizontal, reconstructed
Z vertical. Response-column normalization retains unreconstructed inputs.
Optional plotting checks the cache/result hashes and original photon identities
before writing individual square plots for response, residual versus truth,
and inclusive success/failure versus truth (one PDF and PNG per plot):

```bash
python -m pip install matplotlib  # optional, not needed by the NumPy core
python studyDepthReconstruction.py plot --input heldout_cache --run results/depth-scan \
  --output results/depth-plots --method neighbour --sigma-ps 50 --label "pi+ 100 GeV"
```

Keep particle and energy samples separate. A plot of reconstructed versus true
photon-production depth is the primary depth deliverable; arrival time versus
truth is only a transport diagnostic. Neither establishes energy-deposition
resolution or cures the unresolved optical-transport issue.

`reduceDepthCampaign.py` is the Slurm-only existing-ROOT baseline driver for the
complete pinned 40,000-event campaign. It evaluates **every accepted individual
photon**, not selected channels or first hits. It preserves the same photons and
noise draws across timing assumptions, records out-of-window photons separately,
and closes incident = undetected + holdoff loss + accepted. Physical inversions
and failures sum to all scored photons; residual histograms must also close.
Per-source immutable products permit file-block uncertainties without exporting
raw photons. The unfitted axial baseline scores all events; the historical fixed
index baseline excludes its original 20 calibration events from scores only.
Neither baseline is a trained path model or a certified final detector result.
The plan pins the combined manifest, code hashes, reference table, timing scan,
optics and calibration exclusions; use a new output directory for every new plan.
Its `--merge NEW_CHECKPOINT` mode verifies every per-source hash, event identity,
sensor-stage count and response/residual closure before publishing a checkpoint.
It reports 500-draw source-file bootstrap intervals, separately by particle and
energy. Those intervals describe statistical uncertainty only; they do not cover
the unresolved transport model, path ambiguity or electronics assumptions.
After this strict merge, `studyDepthReconstruction.py campaign-plot --input
CHECKPOINT --plan all_photon_plan.json --output NEW_DIRECTORY --baseline axial
--sigma-ps 50` writes three individual figures for every sample. It requires all
40,000 events, checks the pinned plan and each histogram hash/closure before
writing any plot, and records model and acquisition assumptions in each manifest.
An additional accuracy-versus-timestamp-resolution figure compares both frozen
baselines, with source-file statistical intervals and failed inversions retained
in the denominator. The exported CSV contains all samples and resolution points.
`reduceDepthCampaign.py --depth-support NEW_DIRECTORY --source-checkpoint CHECKPOINT
--summary-sha256 SHA256` (with the same `--plan` and `--output`) uses only the
hash-checked per-source aggregates to compute depth-bin file-bootstrap intervals.
It records independent-file coverage and qualifies the within-10-cm rate only
when at least 30 source files contribute and relative bootstrap error is <=5%.
This is a display-precision flag, not an event/photon cut or a physics-validity
certificate. Undefined or sparsely supported bins remain explicitly recorded.

## Comparison models and scope

The core retains 19-DREAM-PMT / 448-HG-DREAM-C-channel transverse projections,
independent-fiber FLASH response, the conventional-SiPM partial-recharge proxy,
shared-channel resource counting and explicit electronics-display convolution.
The unchanged layout-projection APIs use their documented rod-centred reference
frame. They are not a separate DREAM detector simulation or a model of actual
HG-DREAM shared optical transport. Corrected raw localization uses the physical
calorimeter frame separately; do not mix the two coordinate conventions.

Accepted-SPAD timestamps/positions are an optimistic observation model, NOT
actual FLASH exported records. Shared-channel count/timestamp resource proxies
are separate explicit calls, not automatically applied to reconstruction. A
convolved waveform proxy is not measured voltage or a dSiPM analog waveform.
Historical seeded comparison APIs retain their RNG convention; matched depth
observations use the newer keyed convention. Neither is silently relabelled.

The unresolved production optical-transport issue and exact photon-birth timing
remain physics-validation gates. Saved first-step time is not exact birth time;
exit momentum may include end refraction and is not in-core path truth. Learning
from delayed output cannot certify the bug away. No detector timing resolution,
explicit scintillation, all-fiber optical response or prototype performance is
claimed. The production scope remains central-window Cherenkov.

## Preservation and checks

The full pre-cleanup tree, including all paper/talk sources, study drivers,
plotters, provenance notes and production tooling, is preserved at
[`cpad2026-complete-before-cleanup-20261021`](https://github.com/TTU-HEP/dSiPM/tree/cpad2026-complete-before-cleanup-20261021)
(commit `bfac679`). Nothing was force-pushed or removed from Git history.
HPCC checkouts/checkpoints and local paper/figure products are untouched.
Those historical rebuilds should run from that frozen version, not this core.

This cleanup moves computational definitions rather than rederiving the physics.
All current depth-regression cases and the relevant comparison/pulse/readout
tests remain. A frozen before/after numerical comparison checks sensor masks,
charges, layouts, noise, calibration hashes, all reconstruction methods and their
evaluation metrics. Software closure is not a new physics-performance result.
The cleanup passes 48 focused tests, retaining all 35 current depth regressions;
165 frozen arrays/serialized metric records match the pre-cleanup version exactly.

Import migration: `depth_reconstruction` and `depth_evaluation` -> `depth`;
`depth_geometry` and `fiber_transport` -> `geometry`; `depth_observations`,
sensor response and `pulse_shape` -> `sensors`; `depth_io` -> `io`.
Device catalog/legacy plotting/rebuild/production workflows remain in the archive,
not hidden behind compatibility wrappers in the new core.
