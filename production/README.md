# DREAMSim event production for FLASH

This directory contains the event-production layer for the FLASH/dSiPM study.
It is adapted from Samuel McKinley's public `dSiPM_SPAD_NN` workflow at commit
`80810a7e18eb7e1946436d47c49007310708f5e5`.

## What this production represents

The initial campaigns intentionally reproduce Chris Madrid's June 2025 dSiPM
production configuration:

- DREAMSim commit `54b2447c323dd7eb5093fda5baf9ee433bad1bf6`;
- Geant4 11.2.2 in the TTU HPCC Singularity/Apptainer environment;
- the full HG-DREAM shower geometry;
- optical-photon tracking restricted in `B4bSteppingAction.cc` to rod indices
  20--60 and layer indices 15--65.

These samples simulate the complete calorimeter shower but retain explicit
optical photons only in the central 41-by-51 rod/layer window. They are the
correct basis for the current full-detector event-display and response-transfer
studies, but they must not be labelled as all-fiber optical readout.

The pinned DREAMSim configuration also has no scintillation emission spectrum
or yield in the plastic-fiber material. `OP_isScintillation` is consequently
expected to be empty. A later, separately validated DREAMSim production is
required for explicit Cherenkov and scintillation photon readout across the
full detector.

## Campaigns

`campaigns/cpad2026_pilot_100gev.tsv`

- 20 positron and 20 pion events at 100 GeV;
- split into two 10-event jobs per species;
- validates the build, random seeds, output branches, particle/energy labels,
  and file-size scaling before substantial cluster use.

`campaigns/cpad2026_one_fiber_5k.tsv`

- positrons and pions at 10, 30, 100, and 300 GeV;
- 5,000 events for every particle/energy point;
- an energy-dependent split derived from the validated 100 GeV pilot: 100
  events/job at 10 and 30 GeV, 25 events/job at 100 GeV, and 10 events/job at
  300 GeV. This keeps the projected wall time below 75% of the 12-hour request
  for the slower positron sample while preserving exactly 5,000 events at each
  point.

Do not submit the 5,000-event campaign until the pilot ROOT files pass
`validate_dreamsim_output.py` and the HPCC storage estimate has been checked.

## HPCC layout

The submission host needs:

```text
<parent>/dSiPM/
<parent>/DREAMSim/
```

DREAMSim must be built at the pinned commit. The expected executable is
`DREAMSim/sim/build/exampleB4b`. The default TTU container is:

```text
/lustre/research/hep/yofeng/SimulationEnv/alma9forgeant4_sbox
```

## Dry run

From the dSiPM repository:

```bash
SUBMIT=0 \
DREAMSIM_DIR=/path/to/DREAMSim \
OUTPUT_BASE=/tmp/flash-cpad2026-pilot \
./production/submit_dreamsim_campaign.sh \
  production/campaigns/cpad2026_pilot_100gev.tsv
```

This verifies the campaign and writes the SLURM scripts without submitting
jobs. Generated job scripts and logs live under `OUTPUT_BASE` rather than in
the repository.

## Submit the pilot on TTU HPCC

```bash
SUBMIT=1 \
DREAMSIM_DIR=/home/$USER/DREAMSim \
VALIDATOR_PYTHON=/lustre/work/$USER/flash_dreamsim/venv/bin/python3 \
OUTPUT_BASE=/lustre/work/$USER/flash_dreamsim/cpad2026_pilot_100gev \
./production/submit_dreamsim_campaign.sh \
  production/campaigns/cpad2026_pilot_100gev.tsv
```

Every ROOT file is validated inside its SLURM job. A JSON sidecar records the
entry count, required branches, beam label, photon counts, DREAMSim commit,
the exact dSiPM analysis commit, campaign/launcher/validator/executable
checksums, macro checksum, SLURM job ID, execution host, and seeds. The
validator is copied into each job directory when the campaign is prepared, so
a later checkout change cannot alter an already-submitted job. Failed
validation makes the SLURM job fail.

Submit only from a clean, committed dSiPM checkout. For CPAD 2026 production,
the local development branch and the clean HPCC checkout should resolve to the
same commit. The submission script snapshots the campaign, launcher, and
validator under `OUTPUT_BASE/provenance/` and refuses to reuse that output
directory if any snapshot differs.

## After the pilot

Inspect the queue and results with:

```bash
squeue -u $USER
find /lustre/work/$USER/flash_dreamsim/cpad2026_pilot_100gev \
  -name '*.validation.json' -print
```

Only after all four pilot jobs pass should the matched 5,000-event campaign be
submitted. Large ROOT files should stay on HPCC; reduce them there to the
event-level tables and analysis histograms needed locally.

## HPCC reduction environment

The simulation validator uses the existing Python 3.9 environment because it
is deliberately standalone. The shared `dsipm` analysis package requires
Python 3.11 or newer. Create its separate, pinned HPCC environment with:

```bash
/usr/bin/python3.11 -m venv /lustre/work/$USER/flash_dreamsim/analysis-venv-py311
/lustre/work/$USER/flash_dreamsim/analysis-venv-py311/bin/pip install \
  -r production/requirements-hpcc-analysis.txt
```

Once a ROOT file has a valid sidecar, reduce it without copying its photon
vectors off the cluster:

```bash
PYTHONPATH=$PWD \
/lustre/work/$USER/flash_dreamsim/analysis-venv-py311/bin/python \
  production/reduce_dreamsim_output.py /path/to/file.root \
  --output-dir /lustre/work/$USER/flash_dreamsim/cpad2026_reduced \
  --sample eplus_100gev \
  --event-index-offset 0
```

Each input produces a compressed event table, additive histogram archive, and
JSON provenance summary. The event table contains the full pitch/recovery scan,
the FLASH reference response, the conventional-SiPM proxy, timing quantiles,
spatial moments, active-fiber occupancy quantiles, and high-occupancy fiber
counts. Histograms retain the incident/detected/accepted/lost arrival-time and
x-y views, production-depth timing, photon-weighted fiber occupancy, and
detected-SPAD multiplicity for every scanned pitch. The occupancy products
close exactly to the photon totals while keeping the raw photon vectors on
HPCC. The archive also retains the accepted-photon time profile for every
pitch/recovery point and the detected/partial-recharge charge profile of the
conventional-SiPM proxy, allowing matched pulse-shape plots without copying
raw vectors. The reducer uses the same full-detector coordinate and digitization
functions as the local plotting workflow; only the ROOT reader changes from
PyROOT to uproot.
