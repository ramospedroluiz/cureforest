\name{compare_forest_fits}
\alias{compare_forest_fits}
\title{Compare Two Already Fitted Survival Forests}
\description{
Evaluates an already fitted \code{cureforest} object and an ordinary
\code{randomForestSRC} random survival forest on the same held-out data.
}
\usage{
compare_forest_fits(
  cure_fit,
  rsf_fit,
  newdata,
  time = "time",
  status = "event",
  horizons = c(3, 5, 8),
  tail = c(7, 9.5),
  n.cores = 1L
)
}
\arguments{
  \item{cure_fit}{An object returned by \code{cureforest()}.}
  \item{rsf_fit}{An ordinary random survival forest returned by
  \code{randomForestSRC::rfsrc()}.}
  \item{newdata}{Held-out test data containing the outcome and predictors.}
  \item{time}{Name of the observed-time column.}
  \item{status}{Name of the event-indicator column.}
  \item{horizons}{Evaluation horizons for Brier scores.}
  \item{tail}{Two increasing numbers defining the late-survival window used
  for the ordinary RSF risk score.}
  \item{n.cores}{Number of threads used when predicting from the ordinary
  RSF.}
}
\details{
The C-index uses one minus the predicted cure probability for
\code{cureforest}. For the ordinary RSF, it uses one minus average predicted
survival over the late-survival window. The function also calculates
marginal-IPCW Brier scores and the integrated Brier score.
}
\value{
An object of class \code{cureforest_benchmark}. Its \code{performance}
component contains the held-out C-index, Brier scores, and IBS.
}
\examples{
\dontrun{
library(survival)
data(cure_latency_demo)
train <- subset(cure_latency_demo, cohort == "train")
test <- subset(cure_latency_demo, cohort == "test")
form <- Surv(time, event) ~ cure_marker + latency_marker + age + stage

fit.cure <- cureforest(form, train, ntree = 200)
fit.rsf <- randomForestSRC::rfsrc(form, train, ntree = 200)
compare_forest_fits(fit.cure, fit.rsf, test)
}
}
