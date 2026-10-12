# Studies after the October 2026 audit

## What is complete

- The historical 2025 ROOT products reproduce exactly, including the original
  one-hit-per-event efficiency.
- The coordinate audit explains why the old labelled pitch acted as
  `pitch / 3` in physical DREAMSim coordinates.
- The corrected study uses a centered 3 x 3 mm sensor in physical millimeters
  and scans 5--100 um pitch, 10--100% constant PDE, and 1, 4.3, and 40 ns
  recovery on both 5,000-event samples.
- Input cutflows, event-level survival, pulse-shape distortion, and a first
  train/test timing-depth proxy are now machine readable and plotted.
- The full-detector study maps photons across the rod lattice, preserves fiber
  identity, and propagates pitch/recovery response to event-level Cherenkov
  observables.
- A representative conventional C-Series SiPM proxy is evaluated on the same
  photons, separating PDE selection from partial microcell recharge.
- Physical pitch, nominal active-area fill factor, spectral junction response,
  SPAD recovery, and shared-channel timestamp retention are now reported as
  separate response stages.
- A shared-channel stress test scans grouping, 50 ps resolving time, an
  eight-word timestamp buffer, and an eight-bit counter.
- The detector and HPCC-reduction logic has 35 passing unit tests.

## Current result and limits

For the 1 TeV pion stress sample at 20% PDE and 4.3 ns recovery:

| Physical pitch | Usable / incident | Accepted / detected | 10th-percentile event survival |
|---:|---:|---:|---:|
| 10 um | 19.64% | 97.93% | 96.7% |
| 20 um | 18.59% | 92.73% | 90.3% |
| 25 um | 17.95% | 89.50% | 86.2% |

These values support studying a 10--20 um physical-pitch region under the
stated PDE and recovery assumptions. They do not define a universal pitch.
The one-fiber 100 GeV positron and 1 TeV pion files remain unmatched-energy
stress samples.

The available full-detector samples contain 1,000 100 GeV positron events and
100 90 GeV pion events. At 20 um pitch, 20% PDE, and 4.3 ns recovery, electron
calibration gives a +1.67% pion Cherenkov shift relative to the matched
PDE-only baseline. The apparent narrowing of the pion response is saturation
compression. This establishes differential linearity as the first detector
requirement.

The representative Onsemi MicroFC-30035 proxy uses 31% PDE at 420 nm, 4,774
data-sheet microcells, and an 82 ns recharge constant. It retains 88.54% of
the electron charge and 96.93% of the pion charge, producing a +9.47% pion C
shift after electron calibration relative to a matched 31% PDE-only baseline.
The corresponding FLASH point retains 97.79% and 99.43% and shifts by +1.67%.
This analog result is a sensor-only proxy; measured installed-device and
waveform-chain parameters are still required.

The explicit photon branch is Cherenkov-only. Stored scintillation sums exist,
but explicit scintillation photons do not. Present HG-DREAM and FLASH response
must therefore be digitized from a new common C/S photon sample before making
a final dual-readout resolution claim.

The first timing proxy uses the 10th percentile of accepted photon arrival
times to predict mean selected-photon production depth. The ideal test RMSE is
236 mm, with a 59 mm median absolute error and a long residual tail. Adding a
50 ps response and 50 ps TDC changes the RMSE by less than 0.1 mm. A 1.145 ns
response plus 37 ps TDC adds about 12 mm. This proxy does not yet justify a
10--50 ps system requirement.

At 20 um pitch and 70% nominal square active area, the new spectral proxy gives
18.74% detected photons per incident photon, 93.06% recovery survival, and
16.79 exported timestamps per 100 incident photons for an approximately 80 um
timing channel. The 30% peak junction response and active square are study
assumptions, not prototype measurements.

## Priority 1: close the full-detector C/S loop

Generate one full-detector configuration with explicit Cherenkov and
scintillation photons and sufficient event statistics for electrons and pions.
Digitize two readout chains from the same photons:

1. the present HG-DREAM SiPM plus waveform chain;
2. the FLASH PDE, recovery, time-stamp, masking, and output chain.

Reconstruct $C$, $S$, $C/S$, and dual-readout energy with one calibration
procedure. The key pass/fail quantities are electron and pion linearity,
stochastic and constant terms, transverse moments, timing observables, and
response tails. This is the study that directly answers whether HG-DREAM can
be replaced by dSiPMs.

## Priority 2: matched-energy simulation production (running)

The current production generates positron and pion samples with one pinned
DREAMSim commit and configuration at 10, 30, 100, and 300 GeV, with 5,000
events per point. Every ROOT file has deterministic seeds, a validation
sidecar, and commit/configuration checksums. Muons and explicit scintillation
photons remain a later production step.

First outputs:

- selected-photon yield and event-tail distributions by energy and particle;
- physical PDE-by-pitch-by-recovery surfaces at every energy;
- shower-depth truth distributions and simple timing baselines;
- full-tile occupancy maps rather than one collapsed fiber.

The first matched-energy CPAD/JINST figure set is fixed before inspecting the
results:

1. selected Cherenkov photons versus beam energy, with event-bootstrap bands;
2. recovery survival versus energy for positrons and pions at each pitch;
3. the 10th, 50th, and 90th percentiles of event occupancy and recovery loss;
4. electron-calibrated Cherenkov linearity for ideal photons, the FLASH
   reference point, and the conventional-SiPM proxy;
5. RMS/mean and central-68% resolution versus energy, with saturation
   compression shown separately from genuine resolution improvement;
6. incident/detected/accepted/lost photon-arrival profiles versus energy;
7. full-detector and fiber-local x-y occupancy versus energy;
8. matched representative event displays in truth, DREAM, HG-DREAM, and FLASH
   views.

Raw ROOT files remain on HPCC. The checked-in reduction pipeline writes
sample-global event tables and additive histograms, verifies them against the
production manifest, and transfers only compact plot inputs to the local
machine.

## Priority 3: measured optical response

Replace constant PDE and Gaussian timing assumptions with measurements or
traceable reference curves:

- PDE versus wavelength, overvoltage, and incidence angle;
- instrument-response function versus wavelength and avalanche amplitude;
- recovery and partial-recharge behavior from a two-pulse delay scan;
- dark counts, crosstalk, afterpulsing, hot-cell fraction, and masking cost;
- fixed channel skew and TDC differential/integral nonlinearity.

The prototype measurement plan should use the same parameterization consumed
by the simulation so bench results can update the requirement plots directly.

## Priority 4: information versus output volume (started)

Compare identical reconstruction targets under these exported records:

1. all accepted timestamps;
2. first timestamp per timing group;
3. global first/last timestamps plus counts;
4. short per-group time histograms;
5. count only.

Scan timing-group size, timestamp bits, buffer depth, and arbitration. Report
shower-depth bias/resolution and particle-identification performance versus
bits per event and dropped-hit probability. This is more useful for architecture
decisions than a standalone TDC-bin comparison.

The first grouping/buffer scan is complete for a common all-timestamps output
model. The next iteration must vary buffer depth and output policy and add the
installed waveform bit count as the comparison denominator.

## Priority 5: waveform and beam-test closure

Convolve the simulated photon times with a measured SiPM impulse response and
the current HG-DREAM DRS chain. Compare pulse shape, charge, leading-edge time,
and late-light fractions with representative electron, pion, and muon beam
data. Only after this closure should the same simulated events be reduced to
FLASH-like sparse hits.

## Machine-learning work after baselines

Before extending the student neural-network study:

- establish count-only and timing-only linear/tree baselines;
- freeze preprocessing and sensor random seeds;
- keep independent Geant4 training and prediction batches;
- use cross-validation and event-bootstrap uncertainties;
- test generalization across energy and sensor assumptions;
- compare against ideal-photon and present-waveform upper/lower baselines.

The first targets should be shower depth, particle identification, and the
incremental improvement to a dual-readout energy correction.
