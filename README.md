# cureforest

[![R-CMD-check](https://github.com/ramospedroluiz/cureforest/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ramospedroluiz/cureforest/actions/workflows/R-CMD-check.yaml)

`cureforest` fits survival trees and random survival forests whose split rule
targets differences in a finite-tail survival probability. With sufficiently
mature follow-up and negligible residual susceptible survival in the chosen
window, this target can be interpreted as a population-level cure fraction.

The package does not classify an observed individual as biologically cured.
Its estimand is a tail-supported probability, and the maturity controls are an
essential part of the analysis.

## Installation

Install the development version from GitHub with:

```r
install.packages("remotes")
remotes::install_github("ramospedroluiz/cureforest")
```

After a CRAN release, installation will be available through:

```r
install.packages("cureforest")
```

## Two-minute comparison with randomForestSRC

The package includes a prespecified train/test data set in which one marker
changes cure probability and another changes failure timing among susceptible
subjects. The two forests are fitted explicitly, just as in a standard
random-forest example:

```r
install.packages("randomForestSRC")
library(cureforest)
library(survival)

data(cure_latency_demo)
train <- subset(cure_latency_demo, cohort == "train")
test <- subset(cure_latency_demo, cohort == "test")
form <- Surv(time, event) ~
  cure_marker + latency_marker + age + stage + biomarker + noise

fit.cure <- cureforest(
  form, train, ntree = 200, mtry = 6, nodesize = 80, maxdepth = 1,
  tail.start = 7, tail.end = 9.5, n.cores = 6, seed = 20260725
)

fit.rsf <- randomForestSRC::rfsrc(
  form, train, ntree = 200, mtry = 6, nodesize = 80, nodedepth = 1,
  splitrule = "logrank", nthread = 6, seed = 20260725
)

pred.cure <- predict(fit.cure, test, type = "cure")
head(pred.cure)

comparison <- compare_forest_fits(fit.cure, fit.rsf, test, n.cores = 6)
comparison$performance
```

The returned table reports the held-out C-index for the proposed cure forest
and for the ordinary log-rank random survival forest, together with Brier
scores and IBS. The complete script is installed with the package:

```r
system.file("examples", "cureforest-vs-rsf.R", package = "cureforest")
```

The example is a controlled illustration of the difference between cure and
log-rank splitting, not an empirical superiority claim.

## Cure-fraction survival tree

```r
library(survival)
library(cureforest)

fit_tree <- curetree(
  Surv(time, status) ~ age + stage + biomarker,
  data = cancer_data,
  tail.start = 8, tail.end = 10,
  nodesize = 50, maxdepth = 3, nsplit = 10
)

print(fit_tree)
summary(fit_tree)
plot(fit_tree)
predict(fit_tree, newdata = cancer_data, type = "cure")
```

## Cure-fraction random survival forest

```r
fit_forest <- cureforest(
  Surv(time, status) ~ age + stage + biomarker,
  data = cancer_data,
  ntree = 500,
  mtry = 2,
  nodesize = 50,
  maxdepth = 4,
  nsplit = 10,
  sampsize = 0.8,
  honesty = TRUE,
  n.cores = 4,
  seed = 2026
)

fit_forest$predicted                  # out-of-bag cure predictions
cure_importance(fit_forest)
predict(fit_forest, cancer_data, type = "cure")
predict(fit_forest, cancer_data, type = "survival", times = c(5, 8, 10))
predict(fit_forest, cancer_data, type = "susceptible", times = c(5, 8, 10))
predict(fit_forest, cancer_data, type = "dynamic", times = c(5, 8, 10))
```

## Covariate-dependent censoring

For censoring that depends on observed covariates, a cross-fitted Cox model can
estimate the censoring distribution:

```r
fit_ipcw <- randomforestcure(
  Surv(time, status) ~ age + stage + biomarker,
  data = cancer_data,
  estimator = "ipcw",
  ipcw.model = "cox",
  ipcw.folds = 5,
  ntree = 500
)
```

Alternatively, `ipcw.model = "user"` accepts a function `g.fun(t, newdata)`
that returns `P(C > t | X)` for every row of `newdata`.

## Predictions

`predict()` supports:

- `"cure"`: long-term survival or cure probability;
- `"survival"`: population survival;
- `"susceptible"`: survival among susceptible individuals;
- `"dynamic"`: probability of cure conditional on survival to time `t`;
- `"terminal"`: terminal-node identifiers;
- `"all"`: all available functionals.

Within every tree, population survival is constrained to be no smaller than
that tree's estimated tail probability; the coherent tree curves are then
averaged. This order matters in finite forests. The returned population and
susceptible curves satisfy the mixture identity numerically.

## Practical guidance

- Choose the tail window from scientific considerations and follow-up support,
  not by maximizing test-set performance.
- Inspect `summary()` for terminal-node fallback and maturity diagnostics.
- Use `splitrule = "logrank"` with the same controls for a matched comparison.
- Treat IJ standard errors as experimental diagnostics. Requesting them
  requires a fit made with `inference = TRUE`, which enforces one-shot
  fixed-size subsampling without replacement, a sublinear integer subsample
  schedule, prespecified tail times, positive score stabilization, and
  external user-supplied IPCW weights when IPCW is used. Calibration around
  the forest expectation does not remove localization or residual-tail bias,
  so target-centered intervals require additional bias rates;
  `se.fit = FALSE` remains the default.
- When IJ output is requested, inspect `ij.unstable` and
  `finite.B.correction.fraction`; a warning means that the finite-tree
  subtraction exhausted the raw IJ variance and more trees are needed before
  using the corrected standard error.
- `nsplit` values between 8 and 20 are often a useful starting range.
- `n.cores > 1` uses PSOCK workers and is supported on Windows, macOS, and
  Unix-like systems.

The original functions `survcure()` and `randomforestcure()` remain available
for backward compatibility.

See `vignette("cureforest-introduction", package = "cureforest")` for a complete
reproducible example.
