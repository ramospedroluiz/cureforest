\name{randomforestcure}
\alias{randomforestcure}
\title{Fit a Cure-Fraction Random Survival Forest}
\usage{
randomforestcure(formula, data, ntree = 500L, mtry = NULL,
  nodesize = NULL, maxdepth = 4L, nsplit = 10L, sampsize = 0.8,
  replace = FALSE, honesty = TRUE, honesty.fraction = 0.5,
  splitrule = c("cure", "logrank"), estimator = c("km", "ipcw"),
  tail.start = NULL, tail.quantile = 0.8, tail.end = NULL,
  tail.times = NULL, tail.grid.size = 20L, curve.grid.size = 100L,
  min.tail.at.risk = NULL, min.tail.end.at.risk = NULL,
  min.events.before.tail = 5L, epsilon = 1e-4, lambda = 0.05,
  min.gain = 0, ipcw.model = c("cox", "marginal", "user"),
  g.fun = NULL, ipcw.folds = 5L, g.lower = 0.02, oob = TRUE,
  importance = TRUE, n.cores = 1L, seed = 1L,
  na.action = stats::na.omit, keep.data = FALSE, inference = FALSE,
  inference.exponent = 0.9, keep.inbag = inference)
}
\arguments{
  \item{formula, data}{A right-censored survival formula and data frame.}
  \item{ntree}{Number of trees.}
  \item{mtry}{Encoded predictors sampled at each node.}
  \item{nodesize, maxdepth, nsplit}{Tree size and split-search controls.}
  \item{sampsize}{A sampling fraction or integer sample size. For
  \code{inference = TRUE}, this must be an integer satisfying the sublinear
  schedule described below.}
  \item{replace}{Whether trees use sampling with replacement.}
  \item{honesty, honesty.fraction}{Honest sample-splitting controls.}
  \item{splitrule, estimator}{Split criterion and node estimator.}
  \item{tail.start, tail.quantile, tail.end}{Tail-window controls.}
  \item{tail.times}{Optional prespecified deterministic times for the tail
  average.}
  \item{tail.grid.size, curve.grid.size}{Stored time-grid controls.}
  \item{min.tail.at.risk, min.tail.end.at.risk, min.events.before.tail}{Maturity controls.}
  \item{epsilon, lambda, min.gain}{Split-score stabilization and stopping controls.}
  \item{ipcw.model, g.fun, ipcw.folds, g.lower}{IPCW controls.}
  \item{oob}{Whether to compute out-of-bag cure predictions.}
  \item{importance}{Whether to compute split-based importance.}
  \item{n.cores}{Number of PSOCK workers.}
  \item{seed}{Random seed.}
  \item{na.action}{Training-data missing-value action.}
  \item{keep.data}{Whether to retain the analysis data.}
  \item{inference}{Enable the honest subsampling architecture and retain
  inclusion indicators required for infinitesimal-jackknife inference.}
  \item{inference.exponent}{A number strictly between zero and one. In
  inference mode, the integer subsample size cannot exceed
  \code{floor(n^inference.exponent)}.}
  \item{keep.inbag}{Whether to retain tree-specific subject-inclusion
  indicators. Defaults to the value of \code{inference}.}
}
\value{An object of class \code{randomforestcure}.}
\description{
Fits a randomized ensemble of cure-fraction survival trees and returns OOB
cure predictions, split diagnostics, and coherent survival functionals.
}
\details{
Each tree is grown from a randomized sample and a randomized set of encoded
predictors at each node. The default sampling is without replacement. When
\code{honesty = TRUE}, the sampled observations are divided into independent
structure and estimation subsamples.

Out-of-bag predictions in \code{predicted} estimate the cure or long-term
survival probability. They are available only for observations with positive
\code{oob.count}. The \code{n.cores} argument uses PSOCK workers. During CRAN
checks the package respects the two-core check-farm limit.

Setting \code{inference = TRUE} activates the restricted kernel covered by
the accompanying asymptotic argument. It requires honest sampling without
replacement, one fixed-size integer subsample per tree, a sublinear subsample
schedule, prespecified \code{tail.times}, positive \code{lambda}, and retained
subject-inclusion indicators. A failed tree is replaced by a bounded root
estimate from the same subsample; the implementation does not reject that
subsample and draw another. If IPCW is requested, only external
user-supplied censoring weights are permitted in this mode.

The tail average is a finite-tail survival target. Calling it an exact cure
fraction additionally requires residual susceptible survival over the chosen
window to be negligible.
}
\examples{
library(survival)
lung2 <- transform(lung, status = as.integer(status == 2))
fit <- randomforestcure(
  Surv(time, status) ~ age + sex,
  data = lung2,
  ntree = 5,
  mtry = 2,
  nodesize = 18,
  maxdepth = 1,
  tail.start = 300,
  tail.end = 600,
  min.tail.at.risk = 6,
  min.tail.end.at.risk = 2,
  n.cores = 1,
  seed = 19
)
predict(fit, lung2[1:4, ], type = "survival", times = c(180, 365, 600))
}
\seealso{\code{\link{survcure}}, \code{\link{cure_importance}},
\code{\link{predict.randomforestcure}}}
\keyword{survival}
\keyword{tree}
