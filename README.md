# DREAM / FLASH photon analysis

A small, NumPy-based core for sensor response and timing-based fiber-depth
reconstruction. It reads existing photons; it does not generate new showers.
The original main-branch scripts are left untouched as legacy references.

## Core

```text
dsipm/
  geometry.py   corrected coordinates, DREAM/HG-DREAM layouts, ideal ray reference
  sensors.py    FLASH holdoff, analog recharge, matched clocks, readout/pulse proxies
  depth.py      model contract, calibration, reconstruction, truth-based evaluation
  io.py         immutable caches and hash-checked, Slurm-only existing-ROOT reader
studyDepthReconstruction.py   compact-cache entry point: calibrate / scan / plot
reduceDepthCampaign.py       Slurm-only complete-campaign baseline reduction / merge
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
Default assumptions are L=2000 mm, core/outer radii=0.39/0.40 mm and indices 1.49/1.42,
not the supplied simplified study's 1.504 model. There is no implicit universal
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

`fit_calibration(observations, kappa_truth, contract, source_sha256=...)` is the
only estimator input taking truth. Inference rejects overlapping calibration
and evaluation events and mismatched geometry, optics, workflow, gate or sensor
models. Hold out tuning events as well; never choose cuts/windows on final test
scores. Defaults are development assumptions, not an optimized performance.

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

After preparing disjoint, validated caches with `timestamp_sigma_ns=0`:

```bash
python studyDepthReconstruction.py calibrate --input calibration_cache --output calibration.json
python studyDepthReconstruction.py scan --input heldout_cache --calibration calibration.json \
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
