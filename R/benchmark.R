.cfr_benchmark_step_survival <- function(survival, source_times, target_times) {
  index <- findInterval(target_times, source_times)
  output <- matrix(1, nrow(survival), length(target_times))
  use <- index > 0L
  if (any(use)) {
    output[, use] <- survival[, index[use], drop = FALSE]
  }
  colnames(output) <- format(target_times, trim = TRUE, scientific = FALSE)
  output
}

.cfr_benchmark_outcome <- function(formula, data) {
  frame <- stats::model.frame(formula, data = data, na.action = stats::na.fail)
  outcome <- stats::model.response(frame)
  if (!inherits(outcome, "Surv") || ncol(outcome) != 2L) {
    stop(
      "formula must have a right-censored survival::Surv(time, status) response.",
      call. = FALSE
    )
  }
  list(
    time = as.numeric(outcome[, 1L]),
    status = as.integer(outcome[, 2L])
  )
}

.cfr_benchmark_concordance <- function(time, status, risk) {
  outcome <- survival::Surv(time, status)
  as.numeric(
    survival::concordance(outcome ~ risk, reverse = TRUE)$concordance
  )
}

.cfr_benchmark_censoring_fit <- function(time, status) {
  survival::survfit(survival::Surv(time, 1L - status) ~ 1)
}

.cfr_benchmark_censoring_survival <- function(fit, times, before = FALSE,
                                               lower = 0.02) {
  query <- as.numeric(times)
  if (before) {
    query <- query -
      sqrt(.Machine$double.eps) * pmax(1, abs(query))
  }
  index <- findInterval(query, fit$time)
  probability <- rep(1, length(query))
  use <- index > 0L
  if (any(use)) probability[use] <- fit$surv[index[use]]
  pmax(probability, lower)
}

.cfr_benchmark_brier <- function(time, status, predicted_survival, horizon,
                                 censoring_fit) {
  after_horizon <- time > horizon
  event_by_horizon <- status == 1L & time <= horizon
  weights <- numeric(length(time))
  weights[after_horizon] <- 1 / .cfr_benchmark_censoring_survival(
    censoring_fit, horizon
  )
  if (any(event_by_horizon)) {
    weights[event_by_horizon] <- 1 / .cfr_benchmark_censoring_survival(
      censoring_fit,
      time[event_by_horizon],
      before = TRUE
    )
  }
  target_survival <- as.numeric(after_horizon)
  mean(weights * (target_survival - predicted_survival)^2)
}

.cfr_benchmark_integral <- function(times, values) {
  if (length(times) < 2L || diff(range(times)) <= 0) return(NA_real_)
  sum(diff(times) * (utils::head(values, -1L) + utils::tail(values, -1L)) / 2) /
    diff(range(times))
}

.cfr_benchmark_integer <- function(x, name, lower = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < lower || x != as.integer(x)) {
    stop(name, " must be a single integer not smaller than ", lower, ".",
         call. = FALSE)
  }
  as.integer(x)
}

#' Compare a cure forest with a standard random survival forest
#'
#' Fits `cureforest` and `randomForestSRC` on the same training cohort and
#' reports held-out Harrell C-index and marginal-IPCW Brier scores. The
#' C-index uses each model's late-survival risk score. No files are written.
#'
#' @param formula A formula with a `survival::Surv(time, status)` response.
#' @param data A data frame containing both training and test observations.
#' @param split Name of the column identifying the two cohorts.
#' @param train,test Values of `split` identifying training and test rows.
#' @param horizons Evaluation horizons for Brier scores.
#' @param tail Two numbers giving the late-survival averaging window.
#' @param ntree Number of trees in each forest.
#' @param mtry Number of candidate variables at each split. By default all
#'   formula terms are considered.
#' @param nodesize Minimum terminal-node size used by both forests.
#' @param maxdepth Maximum tree depth used by both forests.
#' @param nsplit Number of candidate cutpoints for continuous variables.
#' @param n.cores Number of parallel workers.
#' @param seed Random seed shared by the two fits.
#' @return An object of class `cureforest_benchmark`. Its `performance`
#'   component contains the held-out metrics, and its `fits` component contains
#'   both fitted forests.
#' @export
compare_cureforest <- function(formula,
                               data,
                               split = "cohort",
                               train = "train",
                               test = "test",
                               horizons = c(3, 5, 8),
                               tail = c(7, 9.5),
                               ntree = 200L,
                               mtry = NULL,
                               nodesize = 80L,
                               maxdepth = 1L,
                               nsplit = 8L,
                               n.cores = 1L,
                               seed = 20260725L) {
  if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
    stop(
      "Package 'randomForestSRC' is required. Install it with ",
      "install.packages(\"randomForestSRC\").",
      call. = FALSE
    )
  }
  if (!inherits(formula, "formula")) {
    stop("formula must be a formula.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("data must be a data frame.", call. = FALSE)
  }
  if (!is.character(split) || length(split) != 1L || !split %in% names(data)) {
    stop("split must name a column in data.", call. = FALSE)
  }
  ntree <- .cfr_benchmark_integer(ntree, "ntree", 2L)
  nodesize <- .cfr_benchmark_integer(nodesize, "nodesize", 2L)
  maxdepth <- .cfr_benchmark_integer(maxdepth, "maxdepth", 1L)
  nsplit <- .cfr_benchmark_integer(nsplit, "nsplit", 1L)
  n.cores <- .cfr_benchmark_integer(n.cores, "n.cores", 1L)
  seed <- .cfr_benchmark_integer(seed, "seed", 1L)
  if (!is.numeric(horizons) || length(horizons) < 1L ||
      any(!is.finite(horizons)) || any(horizons <= 0)) {
    stop("horizons must contain positive finite values.", call. = FALSE)
  }
  horizons <- sort(unique(as.numeric(horizons)))
  if (!is.numeric(tail) || length(tail) != 2L ||
      any(!is.finite(tail)) || tail[1L] <= 0 || tail[2L] <= tail[1L]) {
    stop("tail must contain two increasing positive values.", call. = FALSE)
  }
  if (max(horizons) > tail[2L]) {
    stop("Evaluation horizons cannot exceed the end of the tail window.",
         call. = FALSE)
  }

  training <- data[data[[split]] == train, , drop = FALSE]
  testing <- data[data[[split]] == test, , drop = FALSE]
  if (nrow(training) < 2L * nodesize || nrow(testing) < 20L) {
    stop("The requested train/test split is too small.", call. = FALSE)
  }
  training_outcome <- .cfr_benchmark_outcome(formula, training)
  testing_outcome <- .cfr_benchmark_outcome(formula, testing)
  if (sum(training_outcome$status) < 10L ||
      sum(testing_outcome$status) < 5L) {
    stop("Both cohorts must contain observed events.", call. = FALSE)
  }

  term_count <- length(attr(stats::terms(formula), "term.labels"))
  if (is.null(mtry)) mtry <- max(1L, term_count)
  mtry <- .cfr_benchmark_integer(mtry, "mtry", 1L)
  rsf_formula <- formula
  rsf_response <- rsf_formula[[2L]]
  rsf_response[[1L]] <- as.name("Surv")
  rsf_formula[[2L]] <- rsf_response
  environment(rsf_formula) <- asNamespace("survival")
  min_tail_risk <- max(8L, min(20L, floor(0.01 * nrow(training))))
  min_tail_end_risk <- max(3L, min(8L, floor(0.004 * nrow(training))))
  min_events <- max(6L, min(8L, floor(0.004 * nrow(training))))
  tail_grid <- seq(tail[1L], tail[2L], length.out = 6L)

  start_cure <- proc.time()[["elapsed"]]
  fit_cure <- cureforest(
    formula = formula,
    data = training,
    ntree = ntree,
    mtry = mtry,
    nodesize = nodesize,
    maxdepth = maxdepth,
    nsplit = nsplit,
    sampsize = 0.8,
    replace = FALSE,
    honesty = TRUE,
    honesty.fraction = 0.5,
    splitrule = "cure",
    estimator = "km",
    tail.start = tail[1L],
    tail.end = tail[2L],
    tail.times = tail_grid,
    tail.grid.size = 24L,
    curve.grid.size = 120L,
    min.tail.at.risk = min_tail_risk,
    min.tail.end.at.risk = min_tail_end_risk,
    min.events.before.tail = min_events,
    oob = FALSE,
    importance = FALSE,
    n.cores = n.cores,
    seed = seed,
    keep.data = FALSE
  )
  cure_seconds <- proc.time()[["elapsed"]] - start_cure

  start_rsf <- proc.time()[["elapsed"]]
  fit_rsf <- randomForestSRC::rfsrc(
    formula = rsf_formula,
    data = training,
    ntree = ntree,
    mtry = mtry,
    nodesize = nodesize,
    nodedepth = maxdepth,
    nsplit = nsplit,
    splitrule = "logrank",
    importance = "none",
    forest = TRUE,
    nthread = n.cores,
    seed = seed
  )
  rsf_seconds <- proc.time()[["elapsed"]] - start_rsf

  cure_survival <- predict(
    fit_cure,
    newdata = testing,
    type = "survival",
    times = horizons
  )
  cure_probability <- as.numeric(
    predict(fit_cure, newdata = testing, type = "cure")
  )
  rsf_prediction <- stats::predict(
    fit_rsf,
    newdata = testing,
    outcome = "test",
    nthread = n.cores
  )
  rsf_survival <- .cfr_benchmark_step_survival(
    rsf_prediction$survival,
    rsf_prediction$time.interest,
    horizons
  )
  rsf_tail_survival <- .cfr_benchmark_step_survival(
    rsf_prediction$survival,
    rsf_prediction$time.interest,
    tail_grid
  )
  rsf_cure_proxy <- rowMeans(rsf_tail_survival)

  censoring_fit <- .cfr_benchmark_censoring_fit(
    testing_outcome$time, testing_outcome$status
  )
  cure_brier <- vapply(
    seq_along(horizons),
    function(j) {
      .cfr_benchmark_brier(
        testing_outcome$time,
        testing_outcome$status,
        cure_survival[, j],
        horizons[j],
        censoring_fit
      )
    },
    numeric(1L)
  )
  rsf_brier <- vapply(
    seq_along(horizons),
    function(j) {
      .cfr_benchmark_brier(
        testing_outcome$time,
        testing_outcome$status,
        rsf_survival[, j],
        horizons[j],
        censoring_fit
      )
    },
    numeric(1L)
  )

  performance <- data.frame(
    Method = c("cureforest", "randomForestSRC"),
    C_index = c(
      .cfr_benchmark_concordance(
        testing_outcome$time,
        testing_outcome$status,
        1 - cure_probability
      ),
      .cfr_benchmark_concordance(
        testing_outcome$time,
        testing_outcome$status,
        1 - rsf_cure_proxy
      )
    ),
    stringsAsFactors = FALSE
  )
  for (j in seq_along(horizons)) {
    name <- paste0("Brier_", format(horizons[j], trim = TRUE), "y")
    performance[[name]] <- c(cure_brier[j], rsf_brier[j])
  }
  performance$IBS <- c(
    .cfr_benchmark_integral(horizons, cure_brier),
    .cfr_benchmark_integral(horizons, rsf_brier)
  )
  performance$Runtime_seconds <- c(cure_seconds, rsf_seconds)

  primary_brier <- ncol(performance) - 2L
  brier_reduction <- 100 * (
    performance[2L, primary_brier] - performance[1L, primary_brier]
  ) / performance[2L, primary_brier]
  contrast <- data.frame(
    Metric = c(
      "C-index difference",
      paste0("Brier reduction at ", max(horizons), " years (%)"),
      "IBS reduction (%)"
    ),
    Value = c(
      performance$C_index[1L] - performance$C_index[2L],
      brier_reduction,
      100 * (performance$IBS[2L] - performance$IBS[1L]) /
        performance$IBS[2L]
    ),
    stringsAsFactors = FALSE
  )

  output <- list(
    call = match.call(),
    performance = performance,
    contrast = contrast,
    fits = list(cureforest = fit_cure, randomForestSRC = fit_rsf),
    predictions = list(
      cureforest_cure = cure_probability,
      randomForestSRC_tail = rsf_cure_proxy,
      horizons = horizons
    ),
    cohorts = data.frame(
      Cohort = c("Training", "Test"),
      N = c(nrow(training), nrow(testing)),
      Events = c(sum(training_outcome$status), sum(testing_outcome$status))
    ),
    settings = list(
      tail = tail,
      ntree = ntree,
      mtry = mtry,
      nodesize = nodesize,
      maxdepth = maxdepth,
      nsplit = nsplit,
      n.cores = n.cores,
      seed = seed
    )
  )
  class(output) <- "cureforest_benchmark"
  output
}

#' Run the self-contained cure-latency demonstration
#'
#' Loads `cure_latency_demo`, fits the two forests, and prints held-out C-index
#' and Brier scores. This controlled benchmark illustrates the consequence of
#' changing the split target; it is not an empirical superiority claim.
#'
#' @inheritParams compare_cureforest
#' @return An object of class `cureforest_benchmark`, invisibly.
#' @export
cureforest_demo <- function(ntree = 200L,
                            n.cores = 1L,
                            seed = 20260725L) {
  data_environment <- new.env(parent = emptyenv())
  utils::data(
    "cure_latency_demo",
    package = "cureforest",
    envir = data_environment
  )
  formula <- stats::as.formula(
    paste(
      "Surv(time, event) ~ cure_marker + latency_marker +",
      "age + stage + biomarker + noise"
    ),
    env = asNamespace("survival")
  )
  result <- compare_cureforest(
    formula = formula,
    data = data_environment$cure_latency_demo,
    split = "cohort",
    train = "train",
    test = "test",
    horizons = c(3, 5, 8),
    tail = c(7, 9.5),
    ntree = ntree,
    mtry = 6L,
    nodesize = 80L,
    maxdepth = 1L,
    nsplit = 8L,
    n.cores = n.cores,
    seed = seed
  )
  print(result)
  invisible(result)
}

#' @export
print.cureforest_benchmark <- function(x, digits = 4L, ...) {
  cat("Held-out cureforest benchmark\n")
  print(x$cohorts, row.names = FALSE)
  cat("\nHigher C-index and lower Brier/IBS are better.\n\n")
  performance <- x$performance
  numeric_columns <- vapply(performance, is.numeric, logical(1L))
  performance[numeric_columns] <- lapply(
    performance[numeric_columns], round, digits = digits
  )
  print(performance, row.names = FALSE)
  cat("\nPositive contrasts favor cureforest.\n\n")
  contrast <- x$contrast
  contrast$Value <- round(contrast$Value, digits = digits)
  print(contrast, row.names = FALSE)
  invisible(x)
}
