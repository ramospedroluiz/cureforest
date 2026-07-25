library(cureforest)

expect_error <- function(expr, pattern = NULL) {
  error <- tryCatch({
    force(expr)
    NULL
  }, error = identity)
  if (is.null(error)) stop("Expected an error, but the expression succeeded.")
  if (!is.null(pattern) && !grepl(pattern, conditionMessage(error), fixed = TRUE)) {
    stop("Unexpected error message: ", conditionMessage(error))
  }
  invisible(error)
}

set.seed(501)
n <- 220
x <- runif(n)
group <- factor(rbinom(n, 1, 0.5), labels = c("a", "b"))
cured <- rbinom(n, 1, plogis(-0.6 + 2.5 * (x - 0.5)))
event_time <- ifelse(cured == 1L, Inf, pmin(rexp(n, 0.22), 7))
censor_time <- runif(n, 0, 10)
dat <- data.frame(
  time = pmin(event_time, censor_time),
  status = as.integer(event_time <= censor_time),
  x = x,
  group = group
)

form <- survival::Surv(time, status) ~ x + group
fit <- curetree(
  form, dat,
  tail.start = 5,
  tail.end = 7,
  nodesize = 20,
  min.tail.at.risk = 5,
  min.tail.end.at.risk = 2,
  maxdepth = 1,
  keep.data = FALSE,
  seed = 4
)

expect_error(predict(fit, type = "cure"), "newdata is required")
expect_error(survcure(status ~ x, dat), "Surv(time, status)")
expect_error(randomforestcure(form, dat, ntree = 0), "ntree")
expect_error(
  randomforestcure(form, dat, ntree = 2, inference = TRUE, honesty = FALSE),
  "requires honesty"
)
expect_error(survcure(form, dat, maxdepth = 99), "maxdepth")
expect_error(
  survcure(form, dat, tail.start = 5, tail.end = 7, tail.times = c(4, 6)),
  "between tail.start and tail.end"
)
expect_error(survcure(form, transform(dat, time = -abs(time))),
             "nonnegative")

all_censored <- transform(dat, status = 0L)
expect_error(survcure(form, all_censored), "At least two observed events")

immature <- dat
expect_error(
  survcure(
    form, immature, tail.start = 8, tail.end = 9,
    min.tail.at.risk = n, min.tail.end.at.risk = n
  ),
  "tail-maturity"
)

new_level <- dat[1:3, ]
new_level$group <- factor("new", levels = "new")
expect_error(predict(fit, new_level, type = "cure"), "new level")

forest <- cureforest(
  form, dat,
  ntree = 8,
  mtry = 2,
  nodesize = 18,
  maxdepth = 1,
  tail.start = 5,
  tail.end = 7,
  min.tail.at.risk = 5,
  min.tail.end.at.risk = 2,
  n.cores = 1,
  seed = 8,
  keep.data = TRUE
)

mc <- predict(forest, dat[1:6, ], type = "cure", se.fit = TRUE)
terminal <- predict(forest, dat[1:6, ], type = "terminal")
stopifnot(
  all(c("fit", "se", "monte.carlo.se", "conf.int", "note") %in% names(mc)),
  length(mc$fit) == 6L,
  all(is.finite(mc$se)),
  identical(dim(terminal), c(6L, 8L)),
  inherits(summary(forest), "summary.randomforestcure"),
  inherits(summary(fit), "summary.survcure")
)

cat("cureforest validation tests passed\n")
