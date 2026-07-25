\name{compare_cureforest}
\alias{compare_cureforest}
\alias{cureforest_demo}
\alias{print.cureforest_benchmark}
\title{Compare Cure-Oriented and Standard Random Survival Forests}
\description{
Fits \code{cureforest} and \code{randomForestSRC} on the same training cohort
and reports held-out Harrell C-index and marginal-IPCW Brier scores.
\code{cureforest_demo()} runs the complete built-in cure-latency example.
}
\usage{
compare_cureforest(
  formula,
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
  seed = 20260725L
)

cureforest_demo(ntree = 200L, n.cores = 1L, seed = 20260725L)

\method{print}{cureforest_benchmark}(x, digits = 4L, ...)
}
\arguments{
  \item{formula}{A formula with a
  \code{survival::Surv(time, status)} response.}
  \item{data}{A data frame containing training and test observations.}
  \item{split}{Name of the column identifying the two cohorts.}
  \item{train, test}{Values of \code{split} identifying training and test rows.}
  \item{horizons}{Evaluation horizons for the Brier scores.}
  \item{tail}{Two increasing numbers defining the late-survival averaging
  window.}
  \item{ntree}{Number of trees in each forest.}
  \item{mtry}{Number of candidate variables at each split. If \code{NULL}, all
  formula terms are considered.}
  \item{nodesize}{Minimum terminal-node size used by both forests.}
  \item{maxdepth}{Maximum tree depth used by both forests.}
  \item{nsplit}{Number of candidate cutpoints for continuous variables.}
  \item{n.cores}{Number of parallel workers.}
  \item{seed}{Random seed shared by the two fits.}
  \item{x}{A \code{cureforest_benchmark} object.}
  \item{digits}{Number of digits shown in printed metrics.}
  \item{...}{Additional arguments, currently unused.}
}
\details{
The C-index uses each method's late-survival risk score: one minus the estimated
cure probability for \code{cureforest}, and one minus the average RSF survival
over \code{tail} for \code{randomForestSRC}. Brier scores use the predicted
population survival at each horizon and marginal Kaplan--Meier weights for the
test-cohort censoring distribution.

The built-in example is a controlled benchmark designed to make cure and
latency signals distinct. It illustrates the consequence of changing the
split target and is not an empirical superiority claim.
}
\value{
An object of class \code{cureforest_benchmark}. The \code{performance}
component contains the held-out metrics; \code{contrast} contains direct
comparisons; and \code{fits} contains both fitted forests.
}
\examples{
if (requireNamespace("randomForestSRC", quietly = TRUE)) {
  result <- cureforest_demo(ntree = 40, n.cores = 1)
  result$performance
}
}
