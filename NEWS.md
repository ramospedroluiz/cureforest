# cureforest 0.1.1

- Added the built-in `cure_latency_demo` data set for a self-contained,
  train/test illustration of cure-latency separation.
- Added `compare_cureforest()` to fit `cureforest` and `randomForestSRC` on the
  same split and report held-out Harrell C-index and marginal-IPCW Brier scores.
- Added `cureforest_demo()` as a one-command reproducible demonstration.

# cureforest 0.1.0

- Initial CRAN candidate release.
- Added `survcure()` for honest cure-fraction survival trees.
- Added `randomforestcure()` with OOB prediction and Windows parallelism.
- Added KM tail averages with linear-time Greenwood variance calculations.
- Added marginal, cross-fitted Cox, and user-supplied IPCW censoring weights.
- Added coherent predictions for population survival, susceptible survival,
  and dynamic cure probabilities.
- Added split-count, total-gain, and root-frequency variable importance.
- Added the fixed-size, without-replacement calibration factor to
  infinitesimal-jackknife variances and influence contributions. Simulations
  support variance calibration around the forest expectation; intervals for a
  cure target still require negligible forest and tail bias.
- Added a reproducible introductory vignette, package citation, CRAN submission
  notes, and portable API tests.
