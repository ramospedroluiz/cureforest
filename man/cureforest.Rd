\name{cureforest}
\alias{cureforest}
\alias{curetree}
\title{Concise Interfaces for Cure-Fraction Trees and Forests}
\usage{
curetree(...)

cureforest(...)
}
\arguments{
  \item{...}{Arguments passed to \code{\link{survcure}} for
  \code{curetree()}, or to \code{\link{randomforestcure}} for
  \code{cureforest()}.}
}
\value{
\code{curetree()} returns an object of class \code{survcure}.
\code{cureforest()} returns an object of class \code{randomforestcure}.
}
\description{
User-facing aliases for fitting cure-fraction survival trees and random
survival forests. The original fitting functions remain available for
backward compatibility.
}
\seealso{
\code{\link{survcure}}, \code{\link{randomforestcure}},
\code{\link{cure_importance}}
}
