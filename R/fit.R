.cfr_make_context <- function(prepared,
                              splitrule,
                              estimator,
                              tail_start,
                              tail_end,
                              tail_grid_size,
                              tail_times,
                              curve_grid_size,
                              nodesize,
                              min_tail_at_risk,
                              min_tail_end_at_risk,
                              min_events_before_tail,
                              maxdepth,
                              mtry,
                              nsplit,
                              epsilon,
                              lambda,
                              min_gain,
                              ipcw_model,
                              g_fun,
                              ipcw_folds,
                              g_lower,
                              seed) {
  event_time <- prepared$time[prepared$status == 1L]
  tail_grid <- if (is.null(tail_times)) {
    tail_values <- event_time[event_time >= tail_start & event_time <= tail_end]
    .cfr_make_grid(
      tail_values, tail_grid_size, include = c(tail_start, tail_end)
    )
  } else {
    tail_times
  }
  curve_grid <- .cfr_make_grid(
    event_time, curve_grid_size,
    include = c(0, tail_start, tail_end, max(prepared$time))
  )

  context <- list(
    x = prepared$x,
    time = prepared$time,
    status = prepared$status,
    variable_names = prepared$variable_names,
    splitrule = splitrule,
    estimator = estimator,
    tail_start = tail_start,
    tail_end = tail_end,
    tail_grid = tail_grid,
    curve_grid = curve_grid,
    nodesize = nodesize,
    min_tail_at_risk = min_tail_at_risk,
    min_tail_end_at_risk = min_tail_end_at_risk,
    min_events_before_tail = min_events_before_tail,
    maxdepth = maxdepth,
    mtry = mtry,
    nsplit = nsplit,
    epsilon = epsilon,
    lambda = lambda,
    min_gain = min_gain,
    ipcw_model = if (estimator == "ipcw") ipcw_model else NA_character_,
    ipcw_clipped_fraction = 0
  )

  if (estimator == "ipcw") {
    combined_grid <- sort(unique(c(curve_grid, tail_grid)))
    ipcw <- .cfr_prepare_ipcw(
      prepared, combined_grid, ipcw_model, g_fun, ipcw_folds, g_lower, seed
    )
    context$ipcw_curve <- ipcw$contribution[, match(curve_grid, combined_grid),
                                             drop = FALSE]
    context$ipcw_tail <- ipcw$contribution[, match(tail_grid, combined_grid),
                                            drop = FALSE]
    context$ipcw_tail_mean <- rowMeans(context$ipcw_tail)
    context$ipcw_clipped_fraction <- ipcw$clipped_fraction
  }
  context
}

.cfr_prepare_fit <- function(formula,
                             data,
                             na.action,
                             splitrule,
                             estimator,
                             tail.start,
                             tail.quantile,
                             tail.end,
                             tail.times,
                             tail.grid.size,
                             curve.grid.size,
                             nodesize,
                             min.tail.at.risk,
                             min.tail.end.at.risk,
                             min.events.before.tail,
                             maxdepth,
                             mtry,
                             nsplit,
                             epsilon,
                             lambda,
                             min.gain,
                             ipcw.model,
                             g.fun,
                             ipcw.folds,
                             g.lower,
                             seed) {
  splitrule <- match.arg(splitrule, c("cure", "logrank"))
  estimator <- match.arg(estimator, c("km", "ipcw"))
  ipcw.model <- match.arg(ipcw.model, c("cox", "marginal", "user"))
  prepared <- .cfr_prepare_data(formula, data, na.action)

  .cfr_assert_scalar(tail.quantile, "tail.quantile", 0, 1, inclusive = FALSE)
  event_time <- prepared$time[prepared$status == 1L]
  if (!is.null(tail.times)) {
    if (!is.numeric(tail.times) || !length(tail.times) ||
        any(!is.finite(tail.times)) || any(tail.times < 0)) {
      stop("tail.times must contain finite nonnegative values.", call. = FALSE)
    }
    tail.times <- sort(unique(as.numeric(tail.times)))
  }
  if (is.null(tail.start)) {
    if (!is.null(tail.times)) {
      tail.start <- min(tail.times)
    } else {
    tail.start <- as.numeric(stats::quantile(
      event_time, tail.quantile, names = FALSE, type = 8
    ))
    }
  }
  .cfr_assert_scalar(tail.start, "tail.start", 0, max(prepared$time))
  if (is.null(tail.end)) {
    if (!is.null(tail.times)) {
      tail.end <- max(tail.times)
    } else {
      tail.end <- as.numeric(stats::quantile(
        prepared$time, 0.90, names = FALSE, type = 8
      ))
      tail.end <- max(tail.start, tail.end)
    }
  }
  .cfr_assert_scalar(tail.end, "tail.end", tail.start, max(prepared$time))
  if (!is.null(tail.times) &&
      any(tail.times < tail.start | tail.times > tail.end)) {
    stop("tail.times must lie between tail.start and tail.end.", call. = FALSE)
  }

  if (is.null(nodesize)) nodesize <- max(20L, floor(sqrt(prepared$n) / 2))
  if (is.null(min.tail.at.risk)) min.tail.at.risk <- max(10L, ceiling(nodesize / 3))
  if (is.null(min.tail.end.at.risk)) {
    min.tail.end.at.risk <- max(3L, ceiling(min.tail.at.risk / 5))
  }
  if (is.null(mtry)) mtry <- max(1L, floor(sqrt(prepared$p)))

  .cfr_assert_scalar(nodesize, "nodesize", 2, prepared$n, integer = TRUE)
  .cfr_assert_scalar(min.tail.at.risk, "min.tail.at.risk", 1, prepared$n,
                     integer = TRUE)
  .cfr_assert_scalar(min.tail.end.at.risk, "min.tail.end.at.risk", 1,
                     prepared$n, integer = TRUE)
  .cfr_assert_scalar(min.events.before.tail, "min.events.before.tail", 1,
                     prepared$n, integer = TRUE)
  .cfr_assert_scalar(maxdepth, "maxdepth", 0, 20, integer = TRUE)
  .cfr_assert_scalar(mtry, "mtry", 1, prepared$p, integer = TRUE)
  .cfr_assert_scalar(nsplit, "nsplit", 1, 1000, integer = TRUE)
  .cfr_assert_scalar(tail.grid.size, "tail.grid.size", 1, 1000, integer = TRUE)
  .cfr_assert_scalar(curve.grid.size, "curve.grid.size", 2, 5000, integer = TRUE)
  .cfr_assert_scalar(epsilon, "epsilon", 0, 0.1, inclusive = FALSE)
  .cfr_assert_scalar(lambda, "lambda", 0, Inf)
  .cfr_assert_scalar(min.gain, "min.gain", 0, Inf)
  .cfr_assert_scalar(ipcw.folds, "ipcw.folds", 1, 50, integer = TRUE)
  .cfr_assert_scalar(g.lower, "g.lower", 0, 1, inclusive = FALSE)
  .cfr_assert_scalar(seed, "seed", 0, .Machine$integer.max, integer = TRUE)

  context <- .cfr_make_context(
    prepared = prepared,
    splitrule = splitrule,
    estimator = estimator,
    tail_start = tail.start,
    tail_end = tail.end,
    tail_grid_size = tail.grid.size,
    tail_times = tail.times,
    curve_grid_size = curve.grid.size,
    nodesize = as.integer(nodesize),
    min_tail_at_risk = as.integer(min.tail.at.risk),
    min_tail_end_at_risk = as.integer(min.tail.end.at.risk),
    min_events_before_tail = as.integer(min.events.before.tail),
    maxdepth = as.integer(maxdepth),
    mtry = as.integer(mtry),
    nsplit = as.integer(nsplit),
    epsilon = epsilon,
    lambda = lambda,
    min_gain = min.gain,
    ipcw_model = ipcw.model,
    g_fun = g.fun,
    ipcw_folds = as.integer(ipcw.folds),
    g_lower = g.lower,
    seed = as.integer(seed)
  )
  root <- .cfr_node_estimate(seq_len(prepared$n), context)
  if (!isTRUE(root$valid)) {
    stop(
      "The complete sample fails the tail-maturity rules. Reduce tail.start/tail.end or the maturity minima.",
      call. = FALSE
    )
  }
  list(prepared = prepared, context = context, root = root)
}

.cfr_fit_metadata <- function(prepared, context, call, keep.data) {
  list(
    call = call,
    formula = prepared$formula,
    terms = prepared$terms,
    xterms = prepared$xterms,
    xlevels = prepared$xlevels,
    contrasts = prepared$contrasts,
    variable_names = prepared$variable_names,
    variable_map = prepared$variable_map,
    n = prepared$n,
    p = prepared$p,
    analysis_rows = prepared$row_index,
    omitted = prepared$omitted,
    tail_start = context$tail_start,
    tail_end = context$tail_end,
    tail_grid = context$tail_grid,
    curve_grid = context$curve_grid,
    splitrule = context$splitrule,
    estimator = context$estimator,
    ipcw_model = context$ipcw_model,
    ipcw_clipped_fraction = context$ipcw_clipped_fraction,
    control = context[c(
      "nodesize", "min_tail_at_risk", "min_tail_end_at_risk",
      "min_events_before_tail", "maxdepth", "mtry", "nsplit",
      "epsilon", "lambda", "min_gain"
    )],
    training_data = if (keep.data) prepared$data else NULL
  )
}

#' Fit a cure-fraction survival tree
#'
#' Fits one honest survival tree whose split criterion targets differences in
#' daughter-node cure fractions.
#'
#' @param formula A formula with a `survival::Surv(time, status)` response.
#' @param data A data frame.
#' @param splitrule Either `"cure"` or `"logrank"`.
#' @param estimator Either `"km"` or `"ipcw"`.
#' @return An object of class `survcure`.
#' @export
survcure <- function(formula,
                     data,
                     splitrule = c("cure", "logrank"),
                     estimator = c("km", "ipcw"),
                     tail.start = NULL,
                     tail.quantile = 0.80,
                     tail.end = NULL,
                     tail.times = NULL,
                     tail.grid.size = 20L,
                     curve.grid.size = 100L,
                     nodesize = NULL,
                     min.tail.at.risk = NULL,
                     min.tail.end.at.risk = NULL,
                     min.events.before.tail = 5L,
                     maxdepth = 3L,
                     mtry = NULL,
                     nsplit = 10L,
                     honesty = TRUE,
                     honesty.fraction = 0.5,
                     epsilon = 1e-4,
                     lambda = 0.05,
                     min.gain = 0,
                     ipcw.model = c("cox", "marginal", "user"),
                     g.fun = NULL,
                     ipcw.folds = 5L,
                     g.lower = 0.02,
                     seed = 1L,
                     na.action = stats::na.omit,
                     keep.data = FALSE) {
  call <- match.call()
  .cfr_assert_scalar(honesty.fraction, "honesty.fraction", 0, 1,
                     inclusive = FALSE)
  fit <- .cfr_prepare_fit(
    formula, data, na.action, splitrule, estimator, tail.start, tail.quantile,
    tail.end, tail.times, tail.grid.size, curve.grid.size, nodesize, min.tail.at.risk,
    min.tail.end.at.risk, min.events.before.tail, maxdepth, mtry, nsplit,
    epsilon, lambda, min.gain, ipcw.model, g.fun, ipcw.folds, g.lower, seed
  )
  task <- list(
    seed = as.integer(seed),
    sampsize = fit$prepared$n,
    replace = FALSE,
    honesty = isTRUE(honesty),
    honesty_fraction = honesty.fraction,
    oob = FALSE,
    keep_inbag = FALSE,
    max_attempts = 20L
  )
  start <- proc.time()[[3L]]
  tree <- .cfr_tree_task(task, fit$context)
  tree$oob <- NULL
  metadata <- .cfr_fit_metadata(fit$prepared, fit$context, call, keep.data)
  object <- c(
    list(tree = tree, elapsed = proc.time()[[3L]] - start,
         honesty = isTRUE(honesty), honesty_fraction = honesty.fraction),
    metadata
  )
  class(object) <- "survcure"
  object
}

#' Fit a cure-fraction random survival forest
#'
#' Fits an ensemble of honest cure-fraction survival trees with randomized
#' covariate selection and optional parallel construction.
#'
#' @param ntree Number of trees.
#' @param sampsize Sampling fraction in `(0,1]` or an integer sample size.
#' @param n.cores Number of parallel R workers.
#' @return An object of class `randomforestcure`.
#' @export
randomforestcure <- function(formula,
                             data,
                             ntree = 500L,
                             mtry = NULL,
                             nodesize = NULL,
                             maxdepth = 4L,
                             nsplit = 10L,
                             sampsize = 0.80,
                             replace = FALSE,
                             honesty = TRUE,
                             honesty.fraction = 0.5,
                             splitrule = c("cure", "logrank"),
                             estimator = c("km", "ipcw"),
                             tail.start = NULL,
                             tail.quantile = 0.80,
                             tail.end = NULL,
                             tail.times = NULL,
                             tail.grid.size = 20L,
                             curve.grid.size = 100L,
                             min.tail.at.risk = NULL,
                             min.tail.end.at.risk = NULL,
                             min.events.before.tail = 5L,
                             epsilon = 1e-4,
                             lambda = 0.05,
                             min.gain = 0,
                             ipcw.model = c("cox", "marginal", "user"),
                             g.fun = NULL,
                             ipcw.folds = 5L,
                             g.lower = 0.02,
                             oob = TRUE,
                             importance = TRUE,
                             n.cores = 1L,
                             seed = 1L,
                             na.action = stats::na.omit,
                             keep.data = FALSE,
                             inference = FALSE,
                             keep.inbag = inference) {
  call <- match.call()
  .cfr_assert_scalar(ntree, "ntree", 1, 100000, integer = TRUE)
  .cfr_assert_scalar(n.cores, "n.cores", 1, 256, integer = TRUE)
  if (tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = "false")) %in%
      c("true", "yes", "1")) {
    n.cores <- min(as.integer(n.cores), 2L)
  }
  .cfr_assert_scalar(honesty.fraction, "honesty.fraction", 0, 1,
                     inclusive = FALSE)
  if (isTRUE(inference) && (!isTRUE(honesty) || isTRUE(replace))) {
    stop("inference = TRUE requires honesty = TRUE and replace = FALSE.",
         call. = FALSE)
  }
  if (isTRUE(replace) && isTRUE(honesty)) {
    warning("Sampling with replacement can place repeated subjects in both honest subsamples.",
            call. = FALSE)
  }

  fit <- .cfr_prepare_fit(
    formula, data, na.action, splitrule, estimator, tail.start, tail.quantile,
    tail.end, tail.times, tail.grid.size, curve.grid.size, nodesize, min.tail.at.risk,
    min.tail.end.at.risk, min.events.before.tail, maxdepth, mtry, nsplit,
    epsilon, lambda, min.gain, ipcw.model, g.fun, ipcw.folds, g.lower, seed
  )
  n <- fit$prepared$n
  sample_n <- if (length(sampsize) == 1L && is.finite(sampsize) &&
                  sampsize > 0 && sampsize <= 1) {
    max(2L, floor(n * sampsize))
  } else {
    .cfr_assert_scalar(sampsize, "sampsize", 2, if (replace) Inf else n,
                       integer = TRUE)
    as.integer(sampsize)
  }
  if (isTRUE(inference) && sample_n >= n) {
    stop("inference = TRUE requires sampsize to be smaller than n.",
         call. = FALSE)
  }
  if (honesty && floor(sample_n * min(honesty.fraction, 1 - honesty.fraction)) <
      2L * fit$context$nodesize) {
    warning("The honest subsamples are small relative to nodesize; many trees may remain shallow.",
            call. = FALSE)
  }

  set.seed(seed)
  tree_seeds <- sample.int(.Machine$integer.max, ntree)
  tasks <- lapply(seq_len(ntree), function(b) list(
    seed = tree_seeds[b],
    sampsize = sample_n,
    replace = isTRUE(replace),
    honesty = isTRUE(honesty),
    honesty_fraction = honesty.fraction,
    oob = isTRUE(oob),
    keep_inbag = isTRUE(keep.inbag),
    max_attempts = 20L
  ))

  start <- proc.time()[[3L]]
  n.cores <- min(as.integer(n.cores), as.integer(ntree))
  if (n.cores == 1L) {
    trees <- lapply(tasks, .cfr_tree_task, context = fit$context)
  } else {
    cluster <- parallel::makePSOCKcluster(
      n.cores,
      rscript_args = "--vanilla"
    )
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    api_environment <- environment(randomforestcure)
    helper_names <- grep(
      "^\\.cfr_", ls(api_environment, all.names = TRUE), value = TRUE
    )
    parallel::clusterExport(cluster, helper_names, envir = api_environment)
    local_context <- fit$context
    parallel::clusterExport(cluster, "local_context", envir = environment())
    trees <- parallel::parLapplyLB(
      cluster, tasks, function(task) .cfr_tree_task(task, local_context)
    )
    parallel::stopCluster(cluster)
    on.exit(NULL, add = FALSE)
  }

  oob_prediction <- rep(NA_real_, n)
  oob_count <- integer(n)
  if (isTRUE(oob)) {
    oob_sum <- numeric(n)
    for (tree in trees) {
      index <- tree$oob
      if (!length(index)) next
      prediction <- .cfr_predict_tree_raw(
        tree, fit$context$x[index, , drop = FALSE], survival = FALSE
      )$cure
      oob_sum[index] <- oob_sum[index] + prediction
      oob_count[index] <- oob_count[index] + 1L
    }
    use <- oob_count > 0L
    oob_prediction[use] <- oob_sum[use] / oob_count[use]
  }
  if (!keep.inbag) {
    for (b in seq_along(trees)) trees[[b]]$inbag <- NULL
  }
  for (b in seq_along(trees)) trees[[b]]$oob <- NULL

  metadata <- .cfr_fit_metadata(fit$prepared, fit$context, call, keep.data)
  object <- c(
    list(
      trees = trees,
      ntree = as.integer(ntree),
      sampsize = sample_n,
      replace = isTRUE(replace),
      honesty = isTRUE(honesty),
      honesty_fraction = honesty.fraction,
      inference = isTRUE(inference),
      n.cores = n.cores,
      oob = isTRUE(oob),
      predicted = oob_prediction,
      oob.count = oob_count,
      importance = if (importance) {
        .cfr_build_importance(trees, fit$prepared$variable_names)
      } else NULL,
      elapsed = proc.time()[[3L]] - start
    ),
    metadata
  )
  class(object) <- "randomforestcure"
  object
}

#' Variable importance for a cure forest
#' @export
cure_importance <- function(object, aggregate = TRUE) {
  if (!inherits(object, "randomforestcure")) {
    stop("object must inherit from randomforestcure.", call. = FALSE)
  }
  importance <- if (is.null(object$importance)) {
    .cfr_build_importance(object$trees, object$variable_names)
  } else object$importance
  importance <- importance[order(importance$total_gain, decreasing = TRUE), ,
                           drop = FALSE]
  if (!isTRUE(aggregate)) return(importance)
  .cfr_aggregate_importance(importance, object$variable_map)
}
