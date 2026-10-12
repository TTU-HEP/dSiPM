# CPAD 2026 narrative and evidence map

## The continuity statement

The 2025 CPAD talk introduced digital SiPM readout for HG-DREAM and used a
one-fiber simulation to study photon occupancy, SPAD segmentation, dead time,
and an early timing-based reconstruction. The 2026 talk reports the next
stage of the same program. A funded SLAC LDRD now supports FLASH development,
and the simulation has advanced from a local design study to an audited,
detector-level comparison.

## What CPAD 2025 established

1. HG-DREAM provides a demanding optical-calorimeter benchmark.
2. Per-SPAD digital readout can preserve photon timing while reducing waveform
   data volume.
3. Local photon density, SPAD size, and recovery determine count loss.
4. The original analysis proposed a maximum 25 x 25 um2 SPAD scale and showed
   first timing/ML studies.

## What changed for CPAD 2026

1. The historical calculation was reproduced exactly.
2. The coordinate convention was audited. The original pitch label acted as
   one third of the physical pitch in DREAMSim coordinates.
3. The corrected study now scans physical pitch, PDE, and recovery with
   deterministic random numbers and event bootstrap intervals.
4. The response model now runs over every physical Cherenkov fiber in the
   available full-detector samples.
5. A representative conventional analog SiPM and FLASH are evaluated on the
   same photons with matched-PDE linear baselines.
6. The next LDRD measurements are now connected directly to simulation
   parameters and calorimeter observables.

## Main claims and evidence

| Claim | Evidence | Limitation that stays visible |
|---|---|---|
| The 2025 result is reproducible | historical closure plot and unit tests | exact reproduction does not validate the physical coordinate interpretation |
| One physical pitch does not define the sensor | physical PDE-pitch-recovery surface | one-fiber positron and pion samples have different energies |
| Recovery can bias time structure before total count becomes dramatic | pulse-shape comparison | photon arrival profiles are not analog voltage waveforms |
| Full-detector loss is concentrated in bright fibers | occupancy and event-display plots | current samples are single-energy and unmatched |
| FLASH reduces recovery-driven response transfer in the present comparison | calibrated performance scan | the conventional chain is a sensor proxy, not the installed waveform system |
| A final dual-readout claim requires new production | validation and response plots | explicit scintillation photons are absent in the present full-detector files |

## Funding-progress language

Use concrete progress statements:

- The LDRD now has a validated simulation chain from Geant4 photons to
  detector-level response.
- The study has converted a qualitative pixel-size idea into measurable
  requirements for PDE, recovery, timing response, and output policy.
- The next prototype measurements have named simulation inputs and pass/fail
  observables.

Avoid claiming that the prototype performance is already measured or that the
present samples establish a final calorimeter resolution.
