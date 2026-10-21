# Depth-integration verification, 21 October 2026

Status: software regression and adapter closure; **not certified physics**.

Command: `python -m unittest discover -s tests -v`

Result: **74 tests passed** (39 pre-existing + 35 new). Verified with Python
3.12 and NumPy 2.3.5. The CLI tests build only small synthetic array fixtures in
temporary directories, not Geant4 events or new production samples.

New coverage includes:

- Forward/inverse timing, phase/group separation and the emission-delay bias.
- Full detector rotations, all 80 fiber/rod endpoint cases and world-front time.
- Explicit unreconstructed status, endpoint coverage and inclusive denominators.
- Invalid-prior-seed exclusion and exact self-photon removal.
- Same-event neighbourhoods, sparse cases, leave-one-rod-out and two-depth input.
- Actual 50% fallback, zero-delay bin and collision-free event/channel grouping.
- PDE/holdoff independence across fibers and matched 0/50/200/500 ps inputs.
- Input-order/partition invariance, photon noise and common entrance-clock error.
- Separately reported local timing contributions, not total depth errors.
- Strict calibration/event split, workflow compatibility, artifact checksums,
  immutable output and rejection of unvalidated calibration truth.
- End-to-end calibration/scan/closure, failure accounting and event bootstrap.
- Optional explicitly frozen fixed-path baseline, not tuned to the test sample.
- Ideal-cylinder transport checked with an independent reflection-by-reflection
  implementation: 64 off-axis rays, meridional two-reflection, axial/direct face
  arrivals, nonguided escape and backward rays.
- Corrected raw-array reader with unavailable exitFiber sentinel and deliberate
  mismatch rejection; enforcement of the Slurm requirement before ROOT access.

A separate read-only adapter check used a previously extracted, corrected
9,217-photon event (not a training sample). Reusing its original masks retained
1,844 detected / 1,842 accepted photons, identical accepted identities and
identical scaled jitter draws at 0/50/200/500 ps. It did **not** fit a calibration,
measure a reconstruction resolution or establish a representative performance.
No production data are included in this branch.

The pinned historical HPCC f12de72 checkout was checked read-only and remains
clean. No generation/reduction jobs were submitted for this integration; no
historical outputs or paper inputs changed.

Not yet verified: a live Slurm ROOT-reader pilot, a trustworthy production path-
factor calibration after the optical-transport/creation-time gates, independent
test-sample physics performance, cross-sample transfer, uncertainty coverage or
the actual exported FLASH timestamp/channel electronics. Those remain release
gates, not conclusions inferred from these passing software tests.
