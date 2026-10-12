# FLASH/HG-DREAM JINST paper draft

This folder contains the audit-first JINST draft based on the corrected
physical HG-DREAM occupancy, event-tail, pulse-shape, and timing studies.

## Build

```bash
pdflatex flash_hgdream_jinst.tex
bibtex flash_hgdream_jinst
pdflatex flash_hgdream_jinst.tex
pdflatex flash_hgdream_jinst.tex
```

The official `jinstpub.sty` and `JHEP.bst` files are included so that the
source can be compiled and later packaged for submission without relying on a
local TeX installation of the JINST style.

## Scientific scope of this draft

The central result is that the historical labelled pitch acted three times
finer than physical pitch because a 3 mm segmentation was applied to a 1 mm
display coordinate. The paper retains the historical curve as software
closure and reports the corrected PDE-by-pitch-by-recovery surface, event-level
tails, prompt pulse distortion, and a first longitudinal timing proxy. It does
not claim a complete FLASH or calorimeter performance prediction.

## Items requiring author decisions before circulation

1. Confirm title, author list, order, affiliations, and corresponding email.
2. Confirm the complete SLAC LDRD funding language and identifier.
3. Treat matched-energy production as mandatory before a particle-dependent
   performance claim; decide whether it must precede the first submission.
4. Obtain and document the exact DREAMSim commit, geometry configuration, and
   production provenance for both large ROOT files.
5. Decide whether the timing-depth proof of concept is mature enough for the
   first submission or belongs in a follow-up.
6. Replace presentation citations with archival references wherever a paper
   becomes available.
7. Confirm the wording of the AI-assistance disclosure against the final
   journal and institutional policy.
8. Tag and archive the exact software and machine-readable outputs, then add
   the persistent identifier to the data-availability statement.
