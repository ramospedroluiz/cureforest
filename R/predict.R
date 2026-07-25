.cfr_select_times <- function(survival, grid, times) {
  if (is.null(times)) {
    attr(survival, "times") <- grid
    return(survival)
  }
  if (!is.numeric(times) || any(!is.finite(times)) || any(times < 0)) {
    stop("times must contain finite nonnegative values.", call. = FALSE)
  }
  index <- findInterval(times, grid)
  out <- matrix(1, nrow(survival), length(times))
  use <- index > 0L
  if (any(use)) out[, use] <- survival[, index[use], drop = FALSE]
  colnames(out) <- format(times, trim = TRUE, scientific = FALSE)
  attr(out, "times") <- times
  out
}

.cfr_derived_prediction <- function(cure, survival, type, epsilon = 1e-8) {
  cure_matrix <- matrix(cure, nrow(survival), ncol(survival))
  if (type == "susceptible") {
    denominator <- pmax(1 - cure_matrix, epsilon)
    return(.cfr_clip((survival - cure_matrix) / denominator))
  }
  if (type == "dynamic") {
    return(.cfr_clip(cure_matrix / pmax(survival, epsilon)))
  }
  survival
}

#' Predict from a cure-fraction survival tree
#' @export
predict.survcure <- function(object,
                             newdata = NULL,
                             type = c("cure", "response", "survival", "susceptible",
                                      "dynamic", "terminal", "all"),
                             times = NULL,
                             ...) {
  type <- match.arg(type)
  if (type == "response") type <- "cure"
  if (is.null(newdata)) {
    if (is.null(object$training_data)) {
      stop("newdata is required because the training data were not retained.",
           call. = FALSE)
    }
    newdata <- object$training_data
  }
  x <- .cfr_new_matrix(object, newdata)
  need_survival <- type %in% c("survival", "susceptible", "dynamic", "all")
  raw <- .cfr_predict_tree_raw(object$tree, x, survival = need_survival)
  names(raw$cure) <- rownames(newdata)
  if (type == "cure") return(raw$cure)
  if (type == "terminal") {
    names(raw$node) <- rownames(newdata)
    return(raw$node)
  }

  survival <- .cfr_select_times(raw$survival, object$curve_grid, times)
  survival_dimension <- dim(survival)
  survival <- matrix(
    .cfr_clip(pmax(as.numeric(survival), rep(raw$cure, survival_dimension[2L]))),
    nrow = survival_dimension[1L], ncol = survival_dimension[2L]
  )
  attr(survival, "times") <- if (is.null(times)) object$curve_grid else times
  rownames(survival) <- rownames(newdata)
  if (type == "survival") return(survival)
  susceptible <- .cfr_derived_prediction(raw$cure, survival, "susceptible")
  dynamic <- .cfr_derived_prediction(raw$cure, survival, "dynamic")
  if (type == "susceptible") return(susceptible)
  if (type == "dynamic") return(dynamic)
  list(
    cure = raw$cure,
    survival = survival,
    susceptible = susceptible,
    dynamic = dynamic,
    terminal = raw$node,
    times = attr(survival, "times")
  )
}

#' Predict from a cure-fraction random survival forest
#' @export
predict.randomforestcure <- function(object,
                                     newdata = NULL,
                                     type = c("cure", "response", "survival", "susceptible",
                                              "dynamic", "terminal", "all"),
                                     times = NULL,
                                     se.fit = FALSE,
                                     se.type = c("monte_carlo", "ij", "both"),
                                     conf.level = 0.95,
                                     finite.B = TRUE,
                                     return.influence = FALSE,
                                     ...) {
  type <- match.arg(type)
  se.type <- match.arg(se.type)
  if (type == "response") type <- "cure"
  if (isTRUE(se.fit) && type != "cure") {
    stop("se.fit is currently available only for type = \"cure\".",
         call. = FALSE)
  }
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("conf.level must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  use_ij <- isTRUE(se.fit) && se.type %in% c("ij", "both")
  if (use_ij) {
    if (!isTRUE(object$honesty) || isTRUE(object$replace)) {
      stop("IJ inference requires an honest forest sampled without replacement.",
           call. = FALSE)
    }
    has_inbag <- vapply(
      object$trees, function(tree) !is.null(tree$inbag), logical(1L)
    )
    if (!all(has_inbag)) {
      stop("IJ inference requires a fit with inference = TRUE or keep.inbag = TRUE.",
           call. = FALSE)
    }
  }
  if (is.null(newdata)) {
    if (type == "cure" && isTRUE(object$oob)) return(object$predicted)
    if (is.null(object$training_data)) {
      stop("newdata is required for this prediction type because the training data were not retained.",
           call. = FALSE)
    }
    newdata <- object$training_data
  }
  x <- .cfr_new_matrix(object, newdata)
  n <- nrow(x)
  need_survival <- type %in% c("survival", "susceptible", "dynamic", "all")
  cure_sum <- numeric(n)
  tree_cure <- if (se.fit) matrix(NA_real_, n, object$ntree) else NULL
  survival_sum <- if (need_survival) {
    matrix(0, n, length(object$curve_grid))
  } else NULL
  terminal <- if (type %in% c("terminal", "all")) {
    matrix(NA_integer_, n, object$ntree)
  } else NULL

  for (b in seq_along(object$trees)) {
    raw <- .cfr_predict_tree_raw(object$trees[[b]], x, survival = need_survival)
    cure_sum <- cure_sum + raw$cure
    if (se.fit) tree_cure[, b] <- raw$cure
    if (need_survival) survival_sum <- survival_sum + raw$survival
    if (!is.null(terminal)) terminal[, b] <- raw$node
  }
  cure <- cure_sum / object$ntree
  names(cure) <- rownames(newdata)
  if (type == "cure") {
    if (!se.fit) return(cure)
    tree_centered <- tree_cure - cure
    tree_variance <- rowMeans(tree_centered^2)
    monte_carlo_variance <- tree_variance / object$ntree
    sampling_variance <- rep(NA_real_, n)
    ij_variance_raw <- rep(NA_real_, n)
    finite_B_correction <- rep(NA_real_, n)
    finite_B_correction_fraction <- rep(NA_real_, n)
    ij_unstable <- rep(FALSE, n)
    influence <- NULL
    inbag_variance_sum <- NA_real_
    finite_population_multiplier <- NA_real_

    if (use_ij) {
      if (object$sampsize >= object$n) {
        stop("IJ inference requires sampsize to be smaller than n.",
             call. = FALSE)
      }
      inbag <- matrix(0, nrow = object$n, ncol = object$ntree)
      for (b in seq_along(object$trees)) {
        inbag[object$trees[[b]]$inbag, b] <- 1
      }
      inbag_centered <- inbag - rowMeans(inbag)
      inclusion_covariance <- inbag_centered %*% t(tree_centered) / object$ntree
      finite_population_multiplier <-
        ((object$n - 1) / object$n) *
        (object$n / (object$n - object$sampsize))^2
      ij_variance_raw <-
        finite_population_multiplier * colSums(inclusion_covariance^2)

      # Leading finite-forest inflation of squared empirical covariances.
      inbag_variance_sum <- sum(rowMeans(inbag_centered^2))
      finite_B_correction <- finite_population_multiplier *
        inbag_variance_sum * tree_variance / object$ntree
      finite_B_correction_fraction <- finite_B_correction /
        pmax(ij_variance_raw, .Machine$double.eps)
      ij_unstable <- isTRUE(finite.B) &
        finite_B_correction >= ij_variance_raw
      sampling_variance <- if (isTRUE(finite.B)) {
        pmax(ij_variance_raw - finite_B_correction, 0)
      } else {
        ij_variance_raw
      }
      if (any(ij_unstable)) {
        warning(
          "The finite-B correction exhausted the raw IJ variance for ",
          sum(ij_unstable), " of ", length(ij_unstable),
          " prediction(s). Increase ntree or treat the corrected IJ standard ",
          "errors as unstable; see $ij.unstable.",
          call. = FALSE
        )
      }
      if (isTRUE(return.influence)) {
        influence <- sqrt(finite_population_multiplier) *
          t(inclusion_covariance)
        rownames(influence) <- rownames(newdata)
        colnames(influence) <- paste0("observation", seq_len(object$n))
      }
    }

    variance <- switch(
      se.type,
      monte_carlo = monte_carlo_variance,
      ij = sampling_variance,
      both = sampling_variance + monte_carlo_variance
    )
    standard_error <- sqrt(pmax(variance, 0))
    probability <- pmin(pmax(cure, 1e-8), 1 - 1e-8)
    eta <- stats::qlogis(probability)
    eta_se <- standard_error / (probability * (1 - probability))
    critical <- stats::qnorm(1 - (1 - conf.level) / 2)
    conf_int <- cbind(
      lower = stats::plogis(eta - critical * eta_se),
      upper = stats::plogis(eta + critical * eta_se)
    )
    rownames(conf_int) <- rownames(newdata)
    return(list(
      fit = cure,
      se = standard_error,
      variance = variance,
      se.type = se.type,
      sampling.se = sqrt(pmax(sampling_variance, 0)),
      monte.carlo.se = sqrt(pmax(monte_carlo_variance, 0)),
      ij.variance.raw = ij_variance_raw,
      finite.B.correction = finite_B_correction,
      finite.B.correction.fraction = finite_B_correction_fraction,
      ij.unstable = ij_unstable,
      inbag.variance.sum = inbag_variance_sum,
      finite.population.multiplier = finite_population_multiplier,
      conf.int = conf_int,
      conf.level = conf.level,
      influence = influence,
      tree.predictions = if (isTRUE(return.influence)) tree_cure else NULL,
      note = if (se.type == "monte_carlo") {
        "The standard error measures finite-forest Monte Carlo variation only."
      } else if (se.type == "ij") {
        paste(
          "Finite-B and without-replacement corrected infinitesimal-jackknife",
          "sampling standard error; cure-target coverage also requires negligible bias.",
          "Inspect ij.unstable before using corrected intervals."
        )
      } else {
        paste(
          "Experimental IJ sampling variation combined with finite-forest Monte Carlo variation;",
          "finite-sample calibration is not guaranteed."
        )
      }
    ))
  }
  if (type == "terminal") {
    rownames(terminal) <- rownames(newdata)
    colnames(terminal) <- paste0("tree", seq_len(object$ntree))
    return(terminal)
  }

  survival <- .cfr_select_times(
    survival_sum / object$ntree, object$curve_grid, times
  )
  survival_dimension <- dim(survival)
  survival <- matrix(
    .cfr_clip(pmax(as.numeric(survival), rep(cure, survival_dimension[2L]))),
    nrow = survival_dimension[1L], ncol = survival_dimension[2L]
  )
  attr(survival, "times") <- if (is.null(times)) object$curve_grid else times
  rownames(survival) <- rownames(newdata)
  if (type == "survival") return(survival)
  susceptible <- .cfr_derived_prediction(cure, survival, "susceptible")
  dynamic <- .cfr_derived_prediction(cure, survival, "dynamic")
  if (type == "susceptible") return(susceptible)
  if (type == "dynamic") return(dynamic)
  list(
    cure = cure,
    survival = survival,
    susceptible = susceptible,
    dynamic = dynamic,
    terminal = terminal,
    times = attr(survival, "times")
  )
}

#' @export
print.survcure <- function(x, ...) {
  terminal <- sum(x$tree$nodes$terminal)
  cat("Cure-fraction survival tree\n")
  cat("  Observations:", x$n, "\n")
  cat("  Encoded predictors:", x$p, "\n")
  cat("  Split rule / estimator:", x$splitrule, "/", x$estimator, "\n")
  cat("  Tail window:", format(x$tail_start, digits = 5), "to",
      format(x$tail_end, digits = 5), "\n")
  cat("  Nodes / terminal nodes:", nrow(x$tree$nodes), "/", terminal, "\n")
  cat("  Root cure estimate:", format(x$tree$root_cure, digits = 4), "\n")
  cat("  Elapsed seconds:", format(x$elapsed, digits = 4), "\n")
  invisible(x)
}

#' @export
summary.survcure <- function(object, ...) {
  nodes <- object$tree$nodes
  internal <- nodes[!nodes$terminal, , drop = FALSE]
  split_table <- if (nrow(internal)) {
    data.frame(
      variable = object$variable_names[internal$split_variable],
      split_value = internal$split_value,
      gain = internal$split_score,
      depth = internal$depth,
      row.names = NULL
    )
  } else data.frame()
  out <- list(
    call = object$call,
    n = object$n,
    p = object$p,
    splitrule = object$splitrule,
    estimator = object$estimator,
    tail_window = c(object$tail_start, object$tail_end),
    nodes = nrow(nodes),
    terminal_nodes = sum(nodes$terminal),
    maximum_depth = max(nodes$depth),
    root_cure = object$tree$root_cure,
    minimum_terminal_tail = min(nodes$n_tail[nodes$terminal], na.rm = TRUE),
    fallback_terminal_fraction = mean(
      nodes$used_parent_fallback[nodes$terminal], na.rm = TRUE
    ),
    splits = split_table,
    elapsed = object$elapsed
  )
  class(out) <- "summary.survcure"
  out
}

#' @export
print.summary.survcure <- function(x, ...) {
  cat("Cure-fraction survival tree summary\n")
  cat("  n =", x$n, "; p =", x$p, "; nodes =", x$nodes,
      "; terminal =", x$terminal_nodes, "\n")
  cat("  maximum depth =", x$maximum_depth,
      "; root cure =", format(x$root_cure, digits = 4), "\n")
  cat("  terminal fallback fraction =",
      format(x$fallback_terminal_fraction, digits = 4), "\n")
  if (nrow(x$splits)) {
    cat("\nSplits:\n")
    print(x$splits, row.names = FALSE)
  } else cat("  No split was accepted.\n")
  invisible(x)
}

#' @export
print.randomforestcure <- function(x, ...) {
  terminal <- vapply(x$trees, function(tree) sum(tree$nodes$terminal), numeric(1))
  cat("Cure-fraction random survival forest\n")
  cat("  Observations / encoded predictors:", x$n, "/", x$p, "\n")
  cat("  Trees:", x$ntree, "\n")
  cat("  Split rule / estimator:", x$splitrule, "/", x$estimator, "\n")
  cat("  Tail window:", format(x$tail_start, digits = 5), "to",
      format(x$tail_end, digits = 5), "\n")
  cat("  Mean terminal nodes:", format(mean(terminal), digits = 4), "\n")
  if (isTRUE(x$oob)) {
    cat("  OOB coverage:", format(mean(x$oob.count > 0), digits = 4), "\n")
  }
  cat("  Workers / elapsed seconds:", x$n.cores, "/",
      format(x$elapsed, digits = 5), "\n")
  invisible(x)
}

#' @export
summary.randomforestcure <- function(object, ...) {
  node_count <- vapply(object$trees, function(tree) nrow(tree$nodes), numeric(1))
  terminal_count <- vapply(
    object$trees, function(tree) sum(tree$nodes$terminal), numeric(1)
  )
  fallback_fraction <- vapply(object$trees, function(tree) {
    terminal <- tree$nodes$terminal
    mean(tree$nodes$used_parent_fallback[terminal], na.rm = TRUE)
  }, numeric(1))
  root_variable <- vapply(object$trees, function(tree) {
    variable <- tree$nodes$split_variable[1L]
    if (is.na(variable)) NA_character_ else object$variable_names[variable]
  }, character(1))
  out <- list(
    call = object$call,
    n = object$n,
    p = object$p,
    ntree = object$ntree,
    splitrule = object$splitrule,
    estimator = object$estimator,
    tail_window = c(object$tail_start, object$tail_end),
    mean_nodes = mean(node_count),
    mean_terminal_nodes = mean(terminal_count),
    mean_terminal_fallback = mean(fallback_fraction),
    no_root_split = mean(is.na(root_variable)),
    root_variables = sort(table(root_variable), decreasing = TRUE),
    oob_coverage = if (object$oob) mean(object$oob.count > 0) else NA_real_,
    importance = cure_importance(object),
    elapsed = object$elapsed
  )
  class(out) <- "summary.randomforestcure"
  out
}

#' @export
print.summary.randomforestcure <- function(x, ...) {
  cat("Cure-fraction random survival forest summary\n")
  cat("  n =", x$n, "; p =", x$p, "; trees =", x$ntree, "\n")
  cat("  mean nodes =", format(x$mean_nodes, digits = 4),
      "; mean terminal nodes =", format(x$mean_terminal_nodes, digits = 4),
      "\n")
  cat("  mean terminal fallback fraction =",
      format(x$mean_terminal_fallback, digits = 4), "\n")
  cat("  root no-split frequency =", format(x$no_root_split, digits = 4),
      "; OOB coverage =", format(x$oob_coverage, digits = 4), "\n")
  cat("\nLeading variables by total gain:\n")
  print(utils::head(x$importance, 10L), row.names = FALSE)
  invisible(x)
}

#' @export
plot.survcure <- function(x, digits = 3L, cex = 0.75, ...) {
  nodes <- x$tree$nodes
  position <- matrix(NA_real_, nrow(nodes), 2L,
                     dimnames = list(NULL, c("x", "y")))
  leaf_counter <- 0L
  place <- function(node) {
    position[node, "y"] <<- -nodes$depth[node]
    if (nodes$terminal[node]) {
      leaf_counter <<- leaf_counter + 1L
      position[node, "x"] <<- leaf_counter
    } else {
      place(nodes$left[node])
      place(nodes$right[node])
      position[node, "x"] <<- mean(position[c(nodes$left[node], nodes$right[node]), "x"])
    }
  }
  place(1L)
  graphics::plot(position, type = "n", axes = FALSE, xlab = "", ylab = "", ...)
  internal <- which(!nodes$terminal)
  for (node in internal) {
    for (child in c(nodes$left[node], nodes$right[node])) {
      graphics::segments(position[node, 1L], position[node, 2L],
                         position[child, 1L], position[child, 2L])
    }
  }
  label <- vapply(seq_len(nrow(nodes)), function(node) {
    if (nodes$terminal[node]) {
      paste0("Node ", node, "\npi=", round(nodes$cure[node], digits),
             "\nn=", nodes$n_estimation[node])
    } else {
      paste0(x$variable_names[nodes$split_variable[node]], " <= ",
             round(nodes$split_value[node], digits))
    }
  }, character(1))
  graphics::text(position[, 1L], position[, 2L], labels = label, cex = cex)
  graphics::title(main = "Cure-fraction survival tree")
  invisible(position)
}
