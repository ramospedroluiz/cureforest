\name{cure_latency_demo}
\alias{cure_latency_demo}
\docType{data}
\title{Controlled Cure-Latency Survival Example}
\description{
A simulated train/test survival data set in which one marker controls the cure
probability and another controls failure timing among susceptible subjects.
The data are included for a self-contained comparison of cure-oriented and
log-rank random survival forests.
}
\usage{cure_latency_demo}
\format{
A data frame with 4000 observations and 9 variables:
\describe{
  \item{\code{time}}{Observed event or censoring time.}
  \item{\code{event}}{Event indicator.}
  \item{\code{cure_marker}}{Binary marker governing cure probability.}
  \item{\code{latency_marker}}{Continuous marker governing susceptible-event
  timing.}
  \item{\code{age}}{Simulated age in years.}
  \item{\code{stage}}{Simulated disease stage.}
  \item{\code{biomarker}}{Continuous noise biomarker.}
  \item{\code{noise}}{Additional continuous noise variable.}
  \item{\code{cohort}}{Prespecified \code{"train"} or \code{"test"} cohort.}
}
}
\details{
The data-generating mechanism is recorded in
\code{data-raw/cure_latency_demo.R}. The cure marker changes long-term event-free
probability, whereas the latency marker strongly changes early failure times
without changing cure status. Administrative censoring occurs at time 10.
}
\source{Simulated by the cureforest authors for package illustration.}
\examples{
data(cure_latency_demo)
table(cure_latency_demo$cohort)
}
