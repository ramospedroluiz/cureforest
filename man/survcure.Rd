\name{survcure}
\alias{survcure}
\title{Fit a Cure-Fraction Survival Tree}
\usage{
survcure(formula, data, splitrule = c("cure", "logrank"),
  estimator = c("km", "ipcw"), tail.start = NULL,
  tail.quantile = 0.8, tail.end = NULL, tail.times = NULL,
  tail.grid.size = 20L,
  curve.grid.size = 100L, nodesize = NULL, min.tail.at.risk = NULL,
  min.tail.end.at.risk = NULL, min.events.before.tail = 5L,
  maxdepth = 3L, mtry = NULL, nsplit = 10L, honesty = TRUE,
  honesty.fraction = 0.5, epsilon = 1e-4, lambda = 0.05,
  min.gain = 0, ipcw.model = c("cox", "marginal", "user"),
  g.fun = NULL, ipcw.folds = 5L, g.lower = 0.02, seed = 1L,
  na.action = stats::na.omit, keep.data = FALSE)
}
\arguments{
  \item{formula}{A two-sided formula with a \code{Surv(time, status)} response.}
  \item{data}{A data frame containing the response and predictors.}
  \item{splitrule}{The cure-oriented or log-rank split rule.}
  \item{estimator}{The KM or IPCW node survival estimator.}
  \item{tail.start, tail.end}{Start and end of the tail window.}
  \item{tail.quantile}{Event-time quantile used when \code{tail.start} is omitted.}
  \item{tail.times}{Optional prespecified deterministic times used for the
  tail average. When supplied, these replace the event-adaptive tail grid.}
  \item{tail.grid.size, curve.grid.size}{Maximum tail and prediction grid sizes.}
  \item{nodesize}{Minimum structure and estimation observations per daughter node.}
  \item{min.tail.at.risk, min.tail.end.at.risk}{Tail maturity minima.}
  \item{min.events.before.tail}{Minimum failures observed by the tail start.}
  \item{maxdepth}{Maximum tree depth.}
  \item{mtry}{Number of encoded predictors considered at a node.}
  \item{nsplit}{Number of quantile cutpoints considered per predictor.}
  \item{honesty}{Whether independent subsamples select splits and estimate nodes.}
  \item{honesty.fraction}{Fraction assigned to split selection.}
  \item{epsilon, lambda}{Probability truncation and variance stabilization constants.}
  \item{min.gain}{Minimum split score required to divide a node.}
  \item{ipcw.model}{Censoring model: cross-fitted Cox, marginal KM, or user supplied.}
  \item{g.fun}{Function of time and data returning conditional censoring survival.}
  \item{ipcw.folds}{Number of folds for the censoring Cox model.}
  \item{g.lower}{Lower truncation bound for censoring survival probabilities.}
  \item{seed}{Random seed.}
  \item{na.action}{Training-data missing-value action.}
  \item{keep.data}{Whether to retain the analysis data in the fitted object.}
}
\value{An object of class \code{survcure}.}
\description{
Fits one survival tree using stabilized contrasts between daughter-node
tail-averaged cure-fraction estimates.
}
\details{
With \code{splitrule = "cure"}, candidate daughter nodes are compared through
a balanced, variance-stabilized contrast between their tail-averaged long-term
survival estimates. With \code{splitrule = "logrank"}, the same tree controls
are retained but the split score is the usual log-rank statistic.

When \code{honesty = TRUE}, different observations select the partition and
estimate terminal-node quantities. A terminal estimate that does not satisfy
the maturity constraints falls back to the nearest mature ancestor. The
summary method reports the frequency of this fallback.

The estimated cure fraction is a population-level tail estimand and is not an
observed individual cure label.

For formal inference, \code{tail.times} can be used to keep the late-survival
estimand fixed across samples and trees.
}
\examples{
library(survival)
lung2 <- transform(lung, status = as.integer(status == 2))
fit <- survcure(
  Surv(time, status) ~ age + sex,
  data = lung2,
  tail.start = 300,
  tail.end = 600,
  maxdepth = 1,
  nodesize = 20,
  min.tail.at.risk = 8,
  min.tail.end.at.risk = 3,
  seed = 7
)
predict(fit, lung2[1:5, ], type = "cure")
}
\seealso{\code{\link{randomforestcure}}, \code{\link{predict.survcure}}}
\keyword{survival}
\keyword{tree}
