if (file.exists(file.path("R", "internals.R"))) {
  source(file.path("R", "internals.R"))
  source(file.path("R", "fit.R"))
  source(file.path("R", "predict.R"))
} else {
  library(cureforest)
}

simulate_cure_data <- function(n = 700, seed = 1) {
  set.seed(seed)
  x1 <- runif(n)
  x2 <- factor(rbinom(n, 1, 0.5), labels = c("low", "high"))
  x3 <- rnorm(n)
  pi0 <- plogis(-0.8 + 3.2 * (x1 - 0.5))
  cured <- rbinom(n, 1, pi0)
  rate <- ifelse(x2 == "high", 0.36, 0.15)
  susceptible_time <- pmin(rexp(n, rate), 9)
  event_time <- ifelse(cured == 1, Inf, susceptible_time)
  censor_k <- 1.2 + 1.0 * x1
  censor_time <- 12 * runif(n)^(1 / censor_k)
  data.frame(
    time = pmin(event_time, censor_time),
    status = as.integer(event_time <= censor_time),
    x1 = x1,
    x2 = x2,
    x3 = x3,
    censor_k = censor_k,
    pi0 = pi0,
    cured = cured
  )
}

data <- simulate_cure_data()
formula <- survival::Surv(time, status) ~ x1 + x2 + x3

tree <- survcure(
  formula,
  data,
  tail.start = 7,
  tail.end = 9,
  tail.times = seq(7, 9, length.out = 5),
  nodesize = 30,
  min.tail.at.risk = 8,
  min.tail.end.at.risk = 3,
  maxdepth = 2,
  nsplit = 8,
  seed = 2026,
  keep.data = TRUE
)
stopifnot(inherits(tree, "survcure"), nrow(tree$tree$nodes) >= 1)
stopped_tree <- survcure(
  formula,
  data,
  tail.start = 7,
  tail.end = 9,
  nodesize = 30,
  min.tail.at.risk = 8,
  min.tail.end.at.risk = 3,
  maxdepth = 2,
  min.gain = 1e12,
  seed = 2026
)
stopifnot(isTRUE(stopped_tree$tree$nodes$terminal[1L]))

tree_cure <- predict(tree, data[1:25, ], type = "cure")
tree_response <- predict(tree, data[1:25, ], type = "response")
tree_survival <- predict(tree, data[1:25, ], type = "survival", times = c(0, 4, 8))
stopifnot(
  length(tree_cure) == 25,
  isTRUE(all.equal(tree_cure, tree_response)),
  all(is.finite(tree_cure)), all(tree_cure >= 0 & tree_cure <= 1),
  identical(dim(tree_survival), c(25L, 3L)),
  all(tree_survival >= 0 & tree_survival <= 1)
)

new_missing <- data[1:5, ]
new_missing$x1[1] <- NA_real_
missing_prediction <- predict(tree, new_missing, type = "cure")
stopifnot(all(is.finite(missing_prediction)))

forest <- randomforestcure(
  formula,
  data,
  ntree = 24,
  mtry = 2,
  nodesize = 25,
  maxdepth = 2,
  nsplit = 8,
  sampsize = 0.8,
  tail.start = 7,
  tail.end = 9,
  tail.times = seq(7, 9, length.out = 5),
  min.tail.at.risk = 7,
  min.tail.end.at.risk = 3,
  seed = 77,
  oob = TRUE,
  importance = TRUE,
  inference = TRUE,
  keep.data = TRUE,
  n.cores = 1
)
stopifnot(
  inherits(forest, "randomforestcure"),
  length(forest$trees) == 24,
  mean(forest$oob.count > 0) > 0.95,
  nrow(cure_importance(forest, aggregate = FALSE)) == forest$p,
  nrow(cure_importance(forest, aggregate = TRUE)) == 3,
  isTRUE(all.equal(forest$tail_grid, seq(7, 9, length.out = 5)))
)

all_prediction <- predict(
  forest, data[1:40, ], type = "all", times = c(0, 3, 6, 9)
)
reconstructed <- matrix(all_prediction$cure, 40, 4) +
  (1 - matrix(all_prediction$cure, 40, 4)) * all_prediction$susceptible
stopifnot(
  max(abs(reconstructed - all_prediction$survival)) < 1e-10,
  all(diff(all_prediction$survival[1, ]) <= 1e-10),
  all(all_prediction$dynamic >= 0 & all_prediction$dynamic <= 1)
)

ij_prediction <- predict(
  forest, data[1:4, ], type = "cure", se.fit = TRUE, se.type = "ij",
  return.influence = TRUE
)
both_prediction <- predict(
  forest, data[1:4, ], type = "cure", se.fit = TRUE, se.type = "both"
)
expected_fpc <- ((forest$n - 1) / forest$n) *
  (forest$n / (forest$n - forest$sampsize))^2
stopifnot(
  length(ij_prediction$fit) == 4L,
  all(is.finite(ij_prediction$se)),
  all(ij_prediction$se >= 0),
  identical(dim(ij_prediction$conf.int), c(4L, 2L)),
  all(ij_prediction$conf.int >= 0 & ij_prediction$conf.int <= 1),
  identical(dim(ij_prediction$influence), c(4L, nrow(data))),
  all(is.finite(ij_prediction$ij.variance.raw)),
  all(is.finite(ij_prediction$finite.B.correction.fraction)),
  is.logical(ij_prediction$ij.unstable),
  length(ij_prediction$ij.unstable) == length(ij_prediction$fit),
  isTRUE(all.equal(
    ij_prediction$finite.population.multiplier, expected_fpc
  )),
  isTRUE(all.equal(
    unname(rowSums(ij_prediction$influence^2)),
    unname(ij_prediction$ij.variance.raw),
    tolerance = 1e-12
  )),
  all(both_prediction$se + 1e-12 >= both_prediction$sampling.se)
)

forest_parallel <- randomforestcure(
  formula,
  data,
  ntree = 24,
  mtry = 2,
  nodesize = 25,
  maxdepth = 2,
  nsplit = 8,
  sampsize = 0.8,
  tail.start = 7,
  tail.end = 9,
  tail.times = seq(7, 9, length.out = 5),
  min.tail.at.risk = 7,
  min.tail.end.at.risk = 3,
  seed = 77,
  oob = TRUE,
  importance = TRUE,
  n.cores = 2
)
serial_prediction <- predict(forest, data[1:30, ], type = "cure")
parallel_prediction <- predict(forest_parallel, data[1:30, ], type = "cure")
stopifnot(isTRUE(all.equal(serial_prediction, parallel_prediction, tolerance = 0)))

g_fun <- function(t, newdata) {
  ratio <- pmin(pmax(t / 12, 0), 1)
  pmax(1 - ratio^newdata$censor_k, 1e-6)
}
tree_ipcw <- survcure(
  formula,
  data,
  estimator = "ipcw",
  ipcw.model = "user",
  g.fun = g_fun,
  tail.start = 7,
  tail.end = 9,
  nodesize = 30,
  min.tail.at.risk = 8,
  min.tail.end.at.risk = 3,
  maxdepth = 1,
  seed = 90
)
stopifnot(
  tree_ipcw$estimator == "ipcw",
  all(is.finite(predict(tree_ipcw, data[1:20, ], type = "cure")))
)

tree_ipcw_cox <- survcure(
  formula,
  data,
  estimator = "ipcw",
  ipcw.model = "cox",
  ipcw.folds = 3,
  tail.start = 7,
  tail.end = 9,
  nodesize = 30,
  min.tail.at.risk = 8,
  min.tail.end.at.risk = 3,
  maxdepth = 1,
  seed = 91
)
stopifnot(
  tree_ipcw_cox$ipcw_model == "cox",
  all(is.finite(predict(tree_ipcw_cox, data[1:20, ], type = "cure")))
)

temporary <- tempfile(fileext = ".rds")
saveRDS(forest, temporary)
restored <- readRDS(temporary)
stopifnot(isTRUE(all.equal(
  predict(forest, data[1:10, ], type = "cure"),
  predict(restored, data[1:10, ], type = "cure")
)))
unlink(temporary)

cat("survcure API tests passed\n")
