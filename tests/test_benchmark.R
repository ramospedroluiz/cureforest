if (file.exists(file.path("R", "internals.R"))) {
  source(file.path("R", "internals.R"))
  source(file.path("R", "fit.R"))
  source(file.path("R", "predict.R"))
  source(file.path("R", "aliases.R"))
  source(file.path("R", "benchmark.R"))
  load(file.path("data", "cure_latency_demo.rda"))
} else {
  library(cureforest)
  data(cure_latency_demo)
}

if (requireNamespace("randomForestSRC", quietly = TRUE)) {
  train_rows <- which(cure_latency_demo$cohort == "train")[1:300]
  test_rows <- which(cure_latency_demo$cohort == "test")[1:250]
  small_demo <- cure_latency_demo[c(train_rows, test_rows), , drop = FALSE]

  benchmark <- compare_cureforest(
    survival::Surv(time, event) ~
      cure_marker + latency_marker + age + stage + biomarker + noise,
    data = small_demo,
    ntree = 4L,
    mtry = 6L,
    nodesize = 20L,
    maxdepth = 1L,
    nsplit = 2L,
    n.cores = 1L,
    seed = 71L
  )

  stopifnot(
    inherits(benchmark, "cureforest_benchmark"),
    identical(benchmark$performance$Method,
              c("cureforest", "randomForestSRC")),
    all(is.finite(benchmark$performance$C_index)),
    all(is.finite(benchmark$performance$Brier_8y)),
    nrow(benchmark$contrast) == 3L
  )

  training <- small_demo[small_demo$cohort == "train", , drop = FALSE]
  testing <- small_demo[small_demo$cohort == "test", , drop = FALSE]
  direct_formula <- stats::as.formula(
    paste(
      "Surv(time, event) ~ cure_marker + latency_marker +",
      "age + stage + biomarker + noise"
    ),
    env = asNamespace("survival")
  )
  fit_cure <- cureforest(
    direct_formula,
    data = training,
    ntree = 4L,
    mtry = 6L,
    nodesize = 20L,
    maxdepth = 1L,
    nsplit = 2L,
    tail.start = 7,
    tail.end = 9.5,
    n.cores = 1L,
    seed = 72L
  )
  fit_rsf <- randomForestSRC::rfsrc(
    direct_formula,
    data = training,
    ntree = 4L,
    mtry = 6L,
    nodesize = 20L,
    nodedepth = 1L,
    nsplit = 2L,
    splitrule = "logrank",
    nthread = 1L,
    seed = 72L
  )
  direct_comparison <- compare_forest_fits(
    fit_cure,
    fit_rsf,
    newdata = testing,
    n.cores = 1L
  )
  stopifnot(
    inherits(direct_comparison, "cureforest_benchmark"),
    identical(
      direct_comparison$performance$Method,
      c("cureforest", "randomForestSRC")
    ),
    all(is.finite(direct_comparison$performance$C_index)),
    all(is.finite(direct_comparison$performance$Brier_8y))
  )
}
