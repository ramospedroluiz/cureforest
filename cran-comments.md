## R CMD check results

- Local environment: Windows 11 x64, R 4.6.0 (2026-04-24 ucrt)
- `R CMD check --as-cran`: 0 errors, 0 warnings, 1 note
- Downstream dependencies: none (new submission)

The single note is the expected incoming-feasibility note for a first
submission:

```text
Maintainer: 'Pedro Luiz Ramos <pedro.ramos@uc.cl>'
New submission
```

## New submission

This is the first submission of `cureforest`. The package implements honest
cure-fraction-oriented survival trees and random survival forests for
right-censored long-term survival data.

The package has no network access, does not write outside the R temporary
directory, and defaults to one core. Examples and vignette computations use one
core and small data sets.

## Notes for the reviewer

The word "cure" denotes a population-level long-term survival estimand. The
documentation explicitly states that it is not an observed individual label
and that adequate tail follow-up is required.
