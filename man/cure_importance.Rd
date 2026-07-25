\name{cure_importance}
\alias{cure_importance}
\title{Split-Based Importance for a Cure Forest}
\usage{cure_importance(object, aggregate = TRUE)}
\arguments{
  \item{object}{A fitted \code{randomforestcure} object.}
  \item{aggregate}{Aggregate encoded factor levels and basis columns to their original formula terms.}
}
\value{A data frame with split counts, total and mean gain, and root frequency.}
\description{Summarizes how often and how strongly encoded predictors split cure forests.}
\details{
This is a split-based diagnostic, not permutation importance. With
\code{aggregate = TRUE}, dummy variables created from one formula term are
combined before ranking.
}
\seealso{\code{\link{randomforestcure}}}
\keyword{survival}
\keyword{tree}
