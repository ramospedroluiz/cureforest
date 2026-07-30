\name{predict.survcure}
\alias{predict.survcure}
\alias{predict.randomforestcure}
\title{Predict Cure and Survival Functionals}
\usage{
\method{predict}{survcure}(object, newdata = NULL,
  type = c("cure", "response", "survival", "susceptible", "dynamic", "terminal", "all"),
  times = NULL, ...)

\method{predict}{randomforestcure}(object, newdata = NULL,
  type = c("cure", "response", "survival", "susceptible", "dynamic", "terminal", "all"),
  times = NULL, se.fit = FALSE,
  se.type = c("monte_carlo", "ij", "both"), conf.level = 0.95,
  finite.B = TRUE, return.influence = FALSE, ...)
}
\arguments{
  \item{object}{A fitted tree or forest.}
  \item{newdata}{A data frame containing the predictors.}
  \item{type}{The requested cure or survival functional. \code{"response"} is an alias for \code{"cure"}.}
  \item{times}{Nonnegative prediction times.}
  \item{se.fit}{Return Monte Carlo or experimental IJ standard errors for cure predictions.}
  \item{se.type}{Use finite-forest Monte Carlo variation, infinitesimal-jackknife
  sampling variation, or their sum.}
  \item{conf.level}{Confidence level for pointwise logit-scale intervals.}
  \item{finite.B}{Apply the leading finite-tree correction to the IJ variance.
  The fixed-size, without-replacement factor is applied independently of this
  option.}
  \item{return.influence}{Return the subject-by-prediction IJ influence
  contributions and tree predictions for covariance-aware contrasts.}
  \item{...}{Reserved for future use.}
}
\value{
A numeric vector for \code{"cure"}; a survival matrix for
\code{"survival"}, \code{"susceptible"}, or \code{"dynamic"}; terminal-node
identifiers for \code{"terminal"}; or a named list for \code{"all"}. For a
forest with \code{se.fit = TRUE} and \code{type = "cure"}, a list containing
\code{fit}, \code{se}, the sampling and Monte Carlo components, and a
pointwise confidence interval is returned.
}
\description{Predicts cure probabilities and related long-term survival quantities.}
\details{
For a forest, each tree survival curve is first constrained to be no smaller
than that tree's tail-supported probability, and these coherent tree curves
are then averaged. The susceptible curve is derived from the returned forest
survival and average tail probability. The dynamic quantity is the probability
of belonging to the long-term-survivor component conditional on survival to
each requested time. These transformations are truncated to the unit interval
for numerical stability.

New factor levels not present during model fitting are not supported. Missing
numeric or encoded predictor values follow the default branch stored in each
tree.

IJ inference requires a forest fitted explicitly with
\code{inference = TRUE}. That mode enforces honest, one-shot, fixed-size
subsampling without replacement, a sublinear integer subsample schedule,
prespecified tail times, and the remaining restrictions documented in
\code{\link{randomforestcure}}. The IJ component
estimates repeated-sample variation; the Monte Carlo component measures
variation from using finitely many randomized trees. Setting
\code{return.influence = TRUE} returns the inclusion covariances needed to
estimate variances of paired prediction contrasts. The IJ output is
calibrated for fixed-size subsampling without replacement and is centered on
the forest expectation. Centering an interval on the finite-tail target, and
then on exact cure, requires additional forest-bias and residual-tail-bias
rates described in the accompanying manuscript.

The returned list includes \code{finite.B.correction.fraction}, the ratio of
the finite-\eqn{B} correction to the raw IJ variance, and
\code{ij.unstable}. The latter is \code{TRUE} when the correction exhausts the
raw variance. In that case the function warns and the zero-truncated corrected
standard error should not be used as evidence of negligible uncertainty.
}
\examples{
library(survival)
lung2 <- transform(lung, status = as.integer(status == 2))
fit <- survcure(
  Surv(time, status) ~ age + sex, data = lung2,
  tail.start = 300, tail.end = 600, maxdepth = 1,
  nodesize = 20, min.tail.at.risk = 8, min.tail.end.at.risk = 3
)
predict(fit, lung2[1:3, ], type = "all", times = c(180, 365, 600))
}
\keyword{survival}
\keyword{predict}
