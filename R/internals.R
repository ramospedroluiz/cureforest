.cfr_clip <- function(x, lower = 0, upper = 1) {
  pmin(pmax(x, lower), upper)
}

.cfr_assert_scalar <- function(x, name, lower = -Inf, upper = Inf,
                               integer = FALSE, inclusive = TRUE) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop(name, " must be one finite value.", call. = FALSE)
  }
  valid <- if (inclusive) x >= lower && x <= upper else x > lower && x < upper
  if (!valid || (integer && x != as.integer(x))) {
    stop(name, " is outside its admissible range.", call. = FALSE)
  }
  invisible(TRUE)
}

.cfr_make_grid <- function(values, size, include = numeric(0)) {
  values <- sort(unique(values[is.finite(values)]))
  if (length(values) > size) {
    values <- unique(as.numeric(stats::quantile(
      values, probs = seq(0, 1, length.out = size), names = FALSE, type = 8
    )))
  }
  sort(unique(c(include, values)))
}

.cfr_prepare_data <- function(formula, data, na.action) {
  if (!inherits(formula, "formula")) {
    stop("formula must be a two-sided formula with a Surv response.", call. = FALSE)
  }
  if (!is.data.frame(data)) data <- as.data.frame(data)
  if (nrow(data) < 2L) stop("data must contain at least two rows.", call. = FALSE)

  tt <- stats::terms(formula, data = data)
  mf <- stats::model.frame(
    tt, data = data, na.action = na.action, drop.unused.levels = TRUE
  )
  y <- stats::model.response(mf)
  if (!inherits(y, "Surv") || ncol(y) != 2L || attr(y, "type") != "right") {
    stop("The response must be survival::Surv(time, status) with right censoring.",
         call. = FALSE)
  }

  time <- as.numeric(y[, 1L])
  status <- as.integer(y[, 2L])
  if (any(!is.finite(time)) || any(time < 0)) {
    stop("Observed times must be finite and nonnegative.", call. = FALSE)
  }
  if (any(!status %in% c(0L, 1L))) {
    stop("The event indicator must use 0 for censoring and 1 for an event.",
         call. = FALSE)
  }
  if (sum(status) < 2L) stop("At least two observed events are required.", call. = FALSE)

  xterms <- stats::delete.response(tt)
  x <- stats::model.matrix(xterms, mf)
  assignment <- attr(x, "assign")
  contrast_specification <- attr(x, "contrasts")
  if ("(Intercept)" %in% colnames(x)) {
    keep <- colnames(x) != "(Intercept)"
    x <- x[, keep, drop = FALSE]
    assignment <- assignment[keep]
  }
  if (ncol(x) == 0L) stop("The formula must contain at least one predictor.", call. = FALSE)
  storage.mode(x) <- "double"

  variable <- vapply(seq_len(ncol(x)), function(j) {
    z <- x[, j]
    any(is.finite(z)) && length(unique(z[is.finite(z)])) > 1L
  }, logical(1))
  if (!all(variable)) {
    x <- x[, variable, drop = FALSE]
    assignment <- assignment[variable]
  }
  if (ncol(x) == 0L) stop("All encoded predictors are constant.", call. = FALSE)

  row_index <- match(rownames(mf), rownames(data))
  if (anyNA(row_index)) {
    stop("Could not map the model frame back to the supplied data.", call. = FALSE)
  }
  omitted <- attr(mf, "na.action")
  if (!is.null(omitted)) {
    warning(length(omitted), " incomplete row(s) were omitted from model fitting.",
            call. = FALSE)
  }

  list(
    formula = formula,
    terms = tt,
    xterms = xterms,
    xlevels = stats::.getXlevels(tt, mf),
    contrasts = contrast_specification,
    x = x,
    time = time,
    status = status,
    data = data[row_index, , drop = FALSE],
    row_index = row_index,
    omitted = omitted,
    variable_names = colnames(x),
    variable_map = stats::setNames(
      attr(tt, "term.labels")[assignment], colnames(x)
    ),
    n = nrow(x),
    p = ncol(x)
  )
}

.cfr_new_matrix <- function(object, newdata) {
  if (!is.data.frame(newdata)) newdata <- as.data.frame(newdata)
  mf <- stats::model.frame(
    object$xterms,
    data = newdata,
    na.action = stats::na.pass,
    xlev = object$xlevels
  )
  x <- stats::model.matrix(object$xterms, mf, contrasts.arg = object$contrasts)
  if ("(Intercept)" %in% colnames(x)) {
    x <- x[, colnames(x) != "(Intercept)", drop = FALSE]
  }

  missing_columns <- setdiff(object$variable_names, colnames(x))
  if (length(missing_columns)) {
    add <- matrix(0, nrow(x), length(missing_columns),
                  dimnames = list(NULL, missing_columns))
    x <- cbind(x, add)
  }
  x <- x[, object$variable_names, drop = FALSE]
  storage.mode(x) <- "double"
  x
}

.cfr_components_sorted <- function(time, status) {
  n <- length(time)
  if (n == 0L) {
    return(list(time = numeric(0), survival = numeric(0), greenwood = numeric(0)))
  }
  runs <- rle(time)
  group <- rep.int(seq_along(runs$lengths), runs$lengths)
  events_all <- as.numeric(rowsum(status, group, reorder = FALSE))
  start <- c(1L, head(cumsum(runs$lengths), -1L) + 1L)
  risk_all <- n - start + 1L
  event <- events_all > 0
  if (!any(event)) {
    return(list(time = numeric(0), survival = numeric(0), greenwood = numeric(0)))
  }

  event_time <- runs$values[event]
  events <- events_all[event]
  risk <- risk_all[event]
  survival <- cumprod(1 - events / risk)
  increment <- ifelse(risk > events, events / (risk * (risk - events)), 0)
  greenwood <- cumsum(increment)
  list(time = event_time, survival = survival, greenwood = greenwood)
}

.cfr_eval_step <- function(components, times, field = "survival", initial = 1) {
  if (length(times) == 0L) return(numeric(0))
  if (length(components$time) == 0L) return(rep(initial, length(times)))
  index <- findInterval(times, components$time)
  values <- components[[field]]
  ifelse(index == 0L, initial, values[pmax(index, 1L)])
}

.cfr_km_estimate_sorted <- function(time, status, context) {
  n <- length(time)
  n_tail <- sum(time >= context$tail_start)
  n_tail_end <- sum(time >= context$tail_end)
  n_events <- sum(status == 1L & time <= context$tail_start)
  if (n == 0L || n_tail < context$min_tail_at_risk ||
      n_tail_end < context$min_tail_end_at_risk ||
      n_events < context$min_events_before_tail) {
    return(list(valid = FALSE, reason = "maturity", pi = NA_real_,
                var_pi = NA_real_, var_logit = NA_real_, n_tail = n_tail,
                n_tail_end = n_tail_end, n_events = n_events,
                tail_slope = NA_real_))
  }

  components <- .cfr_components_sorted(time, status)
  survival <- .cfr_eval_step(components, context$tail_grid, "survival", 1)
  greenwood <- .cfr_eval_step(components, context$tail_grid, "greenwood", 0)
  if (any(!is.finite(survival)) || any(!is.finite(greenwood))) {
    return(list(valid = FALSE, reason = "greenwood", pi = NA_real_,
                var_pi = NA_real_, var_logit = NA_real_, n_tail = n_tail,
                n_tail_end = n_tail_end, n_events = n_events,
                tail_slope = NA_real_))
  }

  weights <- rep.int(1 / length(survival), length(survival))
  weighted_survival <- weights * survival
  greenwood_increment <- c(greenwood[1L], diff(greenwood))
  reverse_sum <- rev(cumsum(rev(weighted_survival)))
  var_pi <- sum(pmax(greenwood_increment, 0) * reverse_sum^2)
  pi_hat <- mean(survival)
  p <- .cfr_clip(pi_hat, context$epsilon, 1 - context$epsilon)
  var_logit <- var_pi / (p^2 * (1 - p)^2)

  list(
    valid = is.finite(pi_hat) && is.finite(var_logit),
    reason = "ok",
    pi = .cfr_clip(pi_hat),
    var_pi = var_pi,
    var_logit = var_logit,
    n_tail = n_tail,
    n_tail_end = n_tail_end,
    n_events = n_events,
    tail_slope = survival[1L] - survival[length(survival)]
  )
}

.cfr_ipcw_estimate <- function(index, context) {
  n <- length(index)
  time <- context$time[index]
  status <- context$status[index]
  n_tail <- sum(time >= context$tail_start)
  n_tail_end <- sum(time >= context$tail_end)
  n_events <- sum(status == 1L & time <= context$tail_start)
  if (n == 0L || n_tail < context$min_tail_at_risk ||
      n_tail_end < context$min_tail_end_at_risk ||
      n_events < context$min_events_before_tail) {
    return(list(valid = FALSE, reason = "maturity", pi = NA_real_,
                var_pi = NA_real_, var_logit = NA_real_, n_tail = n_tail,
                n_tail_end = n_tail_end, n_events = n_events,
                tail_slope = NA_real_))
  }

  h <- context$ipcw_tail_mean[index]
  pi_hat <- mean(h)
  var_pi <- if (n > 1L) stats::var(h) / n else 0
  p <- .cfr_clip(pi_hat, context$epsilon, 1 - context$epsilon)
  var_logit <- var_pi / (p^2 * (1 - p)^2)
  tail_survival <- colMeans(context$ipcw_tail[index, , drop = FALSE])
  list(
    valid = is.finite(pi_hat) && is.finite(var_logit),
    reason = "ok",
    pi = .cfr_clip(pi_hat),
    var_pi = var_pi,
    var_logit = var_logit,
    n_tail = n_tail,
    n_tail_end = n_tail_end,
    n_events = n_events,
    tail_slope = tail_survival[1L] - tail_survival[length(tail_survival)]
  )
}

.cfr_node_estimate <- function(index, context, sorted = FALSE) {
  if (context$estimator == "ipcw") return(.cfr_ipcw_estimate(index, context))
  time <- context$time[index]
  status <- context$status[index]
  if (!sorted) {
    order_time <- order(time)
    time <- time[order_time]
    status <- status[order_time]
  }
  .cfr_km_estimate_sorted(time, status, context)
}

.cfr_terminal_survival <- function(index, context) {
  if (context$estimator == "ipcw") {
    return(.cfr_clip(colMeans(context$ipcw_curve[index, , drop = FALSE])))
  }
  order_time <- order(context$time[index])
  components <- .cfr_components_sorted(
    context$time[index][order_time], context$status[index][order_time]
  )
  .cfr_eval_step(components, context$curve_grid, "survival", 1)
}

.cfr_unrestricted_root_estimate <- function(index, context) {
  if (!length(index)) {
    return(list(
      pi = 0.5, var_pi = 0, n_tail = 0L, n_tail_end = 0L,
      n_events = 0L, tail_slope = 0
    ))
  }
  time <- context$time[index]
  status <- context$status[index]
  if (context$estimator == "ipcw") {
    tail_survival <- colMeans(
      context$ipcw_tail[index, , drop = FALSE]
    )
    contribution <- context$ipcw_tail_mean[index]
    pi_hat <- mean(contribution)
    var_pi <- if (length(contribution) > 1L) {
      stats::var(contribution) / length(contribution)
    } else 0
  } else {
    order_time <- order(time)
    components <- .cfr_components_sorted(
      time[order_time], status[order_time]
    )
    tail_survival <- .cfr_eval_step(
      components, context$tail_grid, "survival", 1
    )
    greenwood <- .cfr_eval_step(
      components, context$tail_grid, "greenwood", 0
    )
    weights <- rep.int(1 / length(tail_survival), length(tail_survival))
    weighted_survival <- weights * tail_survival
    greenwood_increment <- c(greenwood[1L], diff(greenwood))
    reverse_sum <- rev(cumsum(rev(weighted_survival)))
    var_pi <- sum(pmax(greenwood_increment, 0) * reverse_sum^2)
    pi_hat <- mean(tail_survival)
  }
  if (!is.finite(pi_hat)) pi_hat <- 0.5
  if (!is.finite(var_pi)) var_pi <- 0
  list(
    pi = .cfr_clip(pi_hat),
    var_pi = max(var_pi, 0),
    n_tail = sum(time >= context$tail_start),
    n_tail_end = sum(time >= context$tail_end),
    n_events = sum(status == 1L & time <= context$tail_start),
    tail_slope = if (length(tail_survival)) {
      tail_survival[1L] - tail_survival[length(tail_survival)]
    } else 0
  )
}

.cfr_one_shot_fallback_tree <- function(structure, estimation, context,
                                        reason = "tree construction failed") {
  estimate <- .cfr_unrestricted_root_estimate(estimation, context)
  nodes <- data.frame(
    node = 1L,
    depth = 0L,
    split_variable = NA_integer_,
    split_value = NA_real_,
    split_score = NA_real_,
    left = NA_integer_,
    right = NA_integer_,
    default_left = TRUE,
    cure = estimate$pi,
    variance = estimate$var_pi,
    estimate_valid = FALSE,
    estimate_reason = paste0("one_shot_fallback: ", reason),
    used_parent_fallback = FALSE,
    n_structure = length(structure),
    n_estimation = length(estimation),
    n_tail = estimate$n_tail,
    n_tail_end = estimate$n_tail_end,
    tail_slope = estimate$tail_slope,
    terminal = TRUE
  )
  survival <- matrix(
    .cfr_terminal_survival(estimation, context),
    nrow = 1L,
    dimnames = list(
      NULL,
      format(context$curve_grid, trim = TRUE, scientific = FALSE)
    )
  )
  list(
    nodes = nodes,
    survival = survival,
    curve_grid = context$curve_grid,
    variable_names = context$variable_names,
    root_cure = estimate$pi,
    one_shot_fallback = TRUE
  )
}

.cfr_logrank_sorted <- function(time, status, left) {
  n <- length(time)
  runs <- rle(time)
  group <- rep.int(seq_along(runs$lengths), runs$lengths)
  start <- c(1L, head(cumsum(runs$lengths), -1L) + 1L)
  risk <- n - start + 1L
  risk_left <- rev(cumsum(rev(as.integer(left))))[start]
  events <- as.numeric(rowsum(status, group, reorder = FALSE))
  events_left <- as.numeric(rowsum(status * as.integer(left), group, reorder = FALSE))
  use <- events > 0L & risk > 1L
  if (!any(use)) return(NA_real_)
  expected <- events[use] * risk_left[use] / risk[use]
  variance <- events[use] * (risk_left[use] / risk[use]) *
    (1 - risk_left[use] / risk[use]) *
    (risk[use] - events[use]) / (risk[use] - 1)
  v <- sum(variance)
  if (!is.finite(v) || v <= 0) return(NA_real_)
  u <- sum(events_left[use] - expected)
  as.numeric(u^2 / v)
}

.cfr_cutpoints <- function(x, nsplit) {
  values <- sort(unique(x[is.finite(x)]))
  if (length(values) <= 1L) return(numeric(0))
  if (length(values) == 2L) return(mean(values))
  probs <- seq(0, 1, length.out = nsplit + 2L)[-c(1L, nsplit + 2L)]
  cuts <- unique(as.numeric(stats::quantile(x, probs, na.rm = TRUE,
                                             names = FALSE, type = 8)))
  cuts[cuts > values[1L] & cuts < values[length(values)]]
}

.cfr_find_best_split <- function(structure_index, estimation_index, context) {
  p <- ncol(context$x)
  variables <- if (context$mtry >= p) seq_len(p) else sample.int(p, context$mtry)
  x_structure <- context$x[structure_index, , drop = FALSE]
  x_estimation <- context$x[estimation_index, , drop = FALSE]
  time_structure <- context$time[structure_index]
  status_structure <- context$status[structure_index]
  time_order <- order(time_structure)
  sorted_time <- time_structure[time_order]
  sorted_status <- status_structure[time_order]
  best <- NULL
  best_score <- -Inf

  for (variable in variables) {
    cutpoints <- .cfr_cutpoints(x_structure[, variable], context$nsplit)
    if (!length(cutpoints)) next
    for (cutpoint in cutpoints) {
      left_structure <- x_structure[, variable] <= cutpoint
      left_estimation <- x_estimation[, variable] <= cutpoint
      if (anyNA(left_structure) || anyNA(left_estimation)) next
      n_left_structure <- sum(left_structure)
      n_left_estimation <- sum(left_estimation)
      if (n_left_structure < context$nodesize ||
          length(left_structure) - n_left_structure < context$nodesize ||
          n_left_estimation < context$nodesize ||
          length(left_estimation) - n_left_estimation < context$nodesize) next

      balance <- mean(left_structure) * mean(!left_structure)
      if (context$splitrule == "logrank") {
        score <- .cfr_logrank_sorted(
          sorted_time, sorted_status, left_structure[time_order]
        )
        pi_left <- pi_right <- NA_real_
      } else if (context$estimator == "km") {
        left_sorted <- left_structure[time_order]
        left_estimate <- .cfr_km_estimate_sorted(
          sorted_time[left_sorted], sorted_status[left_sorted], context
        )
        right_estimate <- .cfr_km_estimate_sorted(
          sorted_time[!left_sorted], sorted_status[!left_sorted], context
        )
        if (!isTRUE(left_estimate$valid) || !isTRUE(right_estimate$valid)) next
        pi_left <- left_estimate$pi
        pi_right <- right_estimate$pi
        difference <- stats::qlogis(.cfr_clip(pi_left, context$epsilon,
                                              1 - context$epsilon)) -
          stats::qlogis(.cfr_clip(pi_right, context$epsilon,
                                  1 - context$epsilon))
        score <- balance * difference^2 /
          (left_estimate$var_logit + right_estimate$var_logit + context$lambda)
      } else {
        left_estimate <- .cfr_ipcw_estimate(
          structure_index[left_structure], context
        )
        right_estimate <- .cfr_ipcw_estimate(
          structure_index[!left_structure], context
        )
        if (!isTRUE(left_estimate$valid) || !isTRUE(right_estimate$valid)) next
        pi_left <- left_estimate$pi
        pi_right <- right_estimate$pi
        difference <- stats::qlogis(.cfr_clip(pi_left, context$epsilon,
                                              1 - context$epsilon)) -
          stats::qlogis(.cfr_clip(pi_right, context$epsilon,
                                  1 - context$epsilon))
        score <- balance * difference^2 /
          (left_estimate$var_logit + right_estimate$var_logit + context$lambda)
      }

      if (is.finite(score) && score > best_score) {
        best_score <- score
        best <- list(
          variable = variable,
          cutpoint = cutpoint,
          score = score,
          balance = balance,
          pi_left = pi_left,
          pi_right = pi_right,
          left_structure = left_structure,
          left_estimation = left_estimation,
          default_left = n_left_estimation >= length(left_estimation) / 2
        )
      }
    }
  }
  best
}

.cfr_grow_tree <- function(structure_index, estimation_index, context) {
  max_nodes <- 2^(context$maxdepth + 1L) - 1L
  node_store <- vector("list", max_nodes)
  survival_store <- vector("list", max_nodes)
  node_count <- 0L

  root_estimate <- .cfr_node_estimate(estimation_index, context)
  if (!isTRUE(root_estimate$valid)) {
    stop("The root node does not satisfy the requested tail-maturity rules.",
         call. = FALSE)
  }

  grow <- function(structure, estimation, depth, parent_cure) {
    node_count <<- node_count + 1L
    node_id <- node_count
    estimate <- .cfr_node_estimate(estimation, context)
    cure <- if (isTRUE(estimate$valid)) estimate$pi else parent_cure
    node_store[[node_id]] <<- data.frame(
      node = node_id,
      depth = depth,
      split_variable = NA_integer_,
      split_value = NA_real_,
      split_score = NA_real_,
      left = NA_integer_,
      right = NA_integer_,
      default_left = TRUE,
      cure = cure,
      variance = if (isTRUE(estimate$valid)) estimate$var_pi else NA_real_,
      estimate_valid = isTRUE(estimate$valid),
      estimate_reason = as.character(estimate$reason),
      used_parent_fallback = !isTRUE(estimate$valid),
      n_structure = length(structure),
      n_estimation = length(estimation),
      n_tail = estimate$n_tail,
      n_tail_end = estimate$n_tail_end,
      tail_slope = estimate$tail_slope,
      terminal = TRUE
    )

    can_split <- depth < context$maxdepth &&
      length(structure) >= 2L * context$nodesize &&
      length(estimation) >= 2L * context$nodesize
    split <- if (can_split) {
      .cfr_find_best_split(structure, estimation, context)
    } else NULL

    if (is.null(split) || split$score <= context$min_gain) {
      survival_store[[node_id]] <<- .cfr_terminal_survival(estimation, context)
      return(node_id)
    }

    left_id <- grow(
      structure[split$left_structure],
      estimation[split$left_estimation],
      depth + 1L,
      cure
    )
    right_id <- grow(
      structure[!split$left_structure],
      estimation[!split$left_estimation],
      depth + 1L,
      cure
    )
    node_store[[node_id]]$split_variable <<- split$variable
    node_store[[node_id]]$split_value <<- split$cutpoint
    node_store[[node_id]]$split_score <<- split$score
    node_store[[node_id]]$left <<- left_id
    node_store[[node_id]]$right <<- right_id
    node_store[[node_id]]$default_left <<- split$default_left
    node_store[[node_id]]$terminal <<- FALSE
    node_store[[node_id]]$cure <<- cure
    left_id <- left_id
    right_id <- right_id
    survival_store[[node_id]] <<- rep(NA_real_, length(context$curve_grid))
    node_id
  }

  grow(structure_index, estimation_index, 0L, root_estimate$pi)
  nodes <- do.call(rbind, node_store[seq_len(node_count)])
  rownames(nodes) <- NULL
  survival <- do.call(rbind, survival_store[seq_len(node_count)])
  colnames(survival) <- format(context$curve_grid, trim = TRUE, scientific = FALSE)
  list(
    nodes = nodes,
    survival = survival,
    curve_grid = context$curve_grid,
    variable_names = context$variable_names,
    root_cure = root_estimate$pi,
    one_shot_fallback = FALSE
  )
}

.cfr_tree_task <- function(task, context) {
  set.seed(task$seed)
  n <- nrow(context$x)
  tree <- NULL
  last_error <- NULL
  attempts <- if (isTRUE(task$one_shot)) 1L else task$max_attempts
  for (attempt in seq_len(attempts)) {
    sampled <- sample.int(n, task$sampsize, replace = task$replace)
    sampled_unique <- unique(sampled)
    if (task$honesty) {
      sampled <- sample(sampled, length(sampled), replace = FALSE)
      n_structure <- max(2L, floor(task$honesty_fraction * length(sampled)))
      structure <- sampled[seq_len(n_structure)]
      estimation <- sampled[-seq_len(n_structure)]
    } else {
      structure <- estimation <- sampled
    }
    candidate <- try(.cfr_grow_tree(structure, estimation, context), silent = TRUE)
    if (!inherits(candidate, "try-error")) {
      tree <- candidate
      break
    }
    last_error <- candidate
  }
  if (is.null(tree) && isTRUE(task$one_shot)) {
    tree <- .cfr_one_shot_fallback_tree(
      structure, estimation, context, as.character(last_error)
    )
  }
  if (is.null(tree)) {
    stop("Tree construction failed after ", attempts,
         " resampling attempts: ", as.character(last_error), call. = FALSE)
  }
  tree$subsample_attempts <- as.integer(attempt)
  oob <- if (task$oob) setdiff(seq_len(n), sampled_unique) else integer(0)
  tree$oob <- oob
  if (task$keep_inbag) tree$inbag <- sampled_unique
  tree
}

.cfr_route_tree <- function(tree, x) {
  n <- nrow(x)
  node <- rep.int(1L, n)
  repeat {
    active_nodes <- unique(node[!tree$nodes$terminal[node]])
    if (!length(active_nodes)) break
    for (node_id in active_nodes) {
      rows <- which(node == node_id)
      variable <- tree$nodes$split_variable[node_id]
      values <- x[rows, variable]
      go_left <- values <= tree$nodes$split_value[node_id]
      go_left[is.na(go_left)] <- tree$nodes$default_left[node_id]
      node[rows[go_left]] <- tree$nodes$left[node_id]
      node[rows[!go_left]] <- tree$nodes$right[node_id]
    }
  }
  node
}

.cfr_predict_tree_raw <- function(tree, x, survival = FALSE) {
  node <- .cfr_route_tree(tree, x)
  out <- list(node = node, cure = tree$nodes$cure[node])
  if (survival) out$survival <- tree$survival[node, , drop = FALSE]
  out
}

.cfr_build_importance <- function(trees, variable_names) {
  count <- numeric(length(variable_names))
  gain <- numeric(length(variable_names))
  root <- numeric(length(variable_names))
  for (tree in trees) {
    internal <- which(!tree$nodes$terminal)
    if (!length(internal)) next
    variable <- tree$nodes$split_variable[internal]
    count <- count + tabulate(variable, nbins = length(variable_names))
    for (j in seq_along(variable)) {
      gain[variable[j]] <- gain[variable[j]] + tree$nodes$split_score[internal[j]]
    }
    root_variable <- tree$nodes$split_variable[1L]
    if (is.finite(root_variable)) root[root_variable] <- root[root_variable] + 1
  }
  data.frame(
    variable = variable_names,
    split_count = count,
    total_gain = gain,
    mean_gain = ifelse(count > 0, gain / count, 0),
    root_frequency = root / length(trees),
    row.names = NULL
  )
}

.cfr_aggregate_importance <- function(importance, variable_map) {
  mapped <- unname(variable_map[importance$variable])
  mapped[is.na(mapped) | !nzchar(mapped)] <- importance$variable[
    is.na(mapped) | !nzchar(mapped)
  ]
  groups <- split(seq_len(nrow(importance)), mapped)
  rows <- lapply(names(groups), function(variable) {
    index <- groups[[variable]]
    split_count <- sum(importance$split_count[index])
    total_gain <- sum(importance$total_gain[index])
    data.frame(
      variable = variable,
      split_count = split_count,
      total_gain = total_gain,
      mean_gain = if (split_count > 0) total_gain / split_count else 0,
      root_frequency = sum(importance$root_frequency[index]),
      row.names = NULL
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$total_gain, decreasing = TRUE), , drop = FALSE]
}

.cfr_marginal_censoring <- function(time, status, grid) {
  order_time <- order(time)
  components <- .cfr_components_sorted(time[order_time], 1L - status[order_time])
  .cfr_eval_step(components, grid, "survival", 1)
}

.cfr_cox_g <- function(time, status, x, train, predict, grid) {
  if (sum(status[train] == 0L) < 5L) stop("Too few censoring events for Cox IPCW.")
  x_names <- paste0("x", seq_len(ncol(x)))
  training <- data.frame(
    .time = time[train], .censor = 1L - status[train],
    x[train, , drop = FALSE], check.names = FALSE
  )
  names(training)[-(1:2)] <- x_names
  prediction <- data.frame(x[predict, , drop = FALSE], check.names = FALSE)
  names(prediction) <- x_names
  fit <- suppressWarnings(survival::coxph(
    survival::Surv(.time, .censor) ~ .,
    data = training,
    ties = "breslow",
    singular.ok = TRUE,
    control = survival::coxph.control(iter.max = 30),
    model = FALSE,
    x = FALSE,
    y = FALSE
  ))
  baseline <- survival::basehaz(fit, centered = FALSE)
  if (!nrow(baseline)) stop("The censoring Cox model has no baseline hazard.")
  index <- findInterval(grid, baseline$time)
  h0 <- ifelse(index == 0L, 0, baseline$hazard[pmax(index, 1L)])
  lp <- as.numeric(stats::predict(fit, newdata = prediction, type = "lp",
                                  reference = "zero"))
  if (any(!is.finite(lp))) stop("Non-finite censoring-model linear predictors.")
  exp(-outer(exp(.cfr_clip(lp, -20, 20)), h0))
}

.cfr_prepare_ipcw <- function(prepared, grid, model, g_fun, folds, g_lower, seed) {
  n <- prepared$n
  if (model == "user" && !is.function(g_fun)) {
    stop("g.fun must be supplied when ipcw.model = 'user'.", call. = FALSE)
  }

  if (model == "user") {
    g <- matrix(NA_real_, n, length(grid))
    for (j in seq_along(grid)) {
      value <- g_fun(grid[j], prepared$data)
      if (length(value) == 1L) value <- rep(value, n)
      if (length(value) != n || any(!is.finite(value))) {
        stop("g.fun must return one finite censoring-survival value per row.",
             call. = FALSE)
      }
      g[, j] <- value
    }
  } else if (model == "marginal") {
    marginal <- .cfr_marginal_censoring(prepared$time, prepared$status, grid)
    g <- matrix(rep(marginal, each = n), nrow = n)
  } else {
    folds <- min(as.integer(folds), n)
    if (folds < 2L) folds <- 1L
    set.seed(seed + 7919L)
    fold_id <- if (folds == 1L) rep.int(1L, n) else
      sample(rep(seq_len(folds), length.out = n))
    g <- matrix(NA_real_, n, length(grid))
    for (fold in seq_len(folds)) {
      prediction <- which(fold_id == fold)
      training <- if (folds == 1L) seq_len(n) else which(fold_id != fold)
      value <- try(.cfr_cox_g(
        prepared$time, prepared$status, prepared$x,
        training, prediction, grid
      ), silent = TRUE)
      if (inherits(value, "try-error") || any(!is.finite(value))) {
        marginal <- .cfr_marginal_censoring(
          prepared$time[training], prepared$status[training], grid
        )
        value <- matrix(rep(marginal, each = length(prediction)),
                        nrow = length(prediction))
        warning("A censoring Cox fold failed; marginal censoring weights were used for that fold.",
                call. = FALSE)
      }
      g[prediction, ] <- value
    }
  }

  clipped <- mean(g < g_lower)
  g <- .cfr_clip(g, g_lower, 1)
  contribution <- matrix(0, n, length(grid))
  for (j in seq_along(grid)) {
    contribution[, j] <- as.numeric(prepared$time > grid[j]) / g[, j]
  }
  list(contribution = contribution, clipped_fraction = clipped)
}
