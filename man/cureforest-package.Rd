\name{cureforest-package}
\alias{cureforest-package}
\docType{package}
\title{Cure-Fraction Survival Trees and Random Survival Forests}
\description{
Fits honest survival trees and randomized ensembles whose split criterion can
target differences in tail-supported long-term survival probabilities.
}
\details{
The main functions are \code{curetree()} for a single tree and
\code{cureforest()} for an ensemble. The original functions
\code{survcure()} and \code{randomforestcure()} remain available for backward
compatibility. Node estimates use either a
tail-averaged Kaplan--Meier estimator or an inverse-probability-of-censoring
weighted estimator. Explicit maturity constraints control the number at risk
at the beginning and end of the tail window.

The cure probability is a population-level statistical estimand. It should not
be interpreted as an observed biological-cure label for an individual.
}
\author{
Pedro Luiz Ramos and Lucas Osses.
}
\seealso{
\code{\link{survcure}}, \code{\link{randomforestcure}},
\code{\link{predict.survcure}}, and \code{\link{cure_importance}}.
}
\keyword{survival}
\keyword{tree}
