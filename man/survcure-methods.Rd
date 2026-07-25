\name{survcure-methods}
\alias{print.survcure}
\alias{summary.survcure}
\alias{print.summary.survcure}
\alias{plot.survcure}
\alias{print.randomforestcure}
\alias{summary.randomforestcure}
\alias{print.summary.randomforestcure}
\title{Print, Summarize, and Plot survcure Models}
\usage{
\method{print}{survcure}(x, ...)
\method{summary}{survcure}(object, ...)
\method{print}{summary.survcure}(x, ...)
\method{plot}{survcure}(x, digits = 3L, cex = 0.75, ...)
\method{print}{randomforestcure}(x, ...)
\method{summary}{randomforestcure}(object, ...)
\method{print}{summary.randomforestcure}(x, ...)
}
\arguments{
  \item{x}{A fitted model or summary object.}
  \item{object}{A fitted \code{survcure} or \code{randomforestcure} object.}
  \item{digits}{Number of digits used in terminal-node cure labels.}
  \item{cex}{Character expansion factor for tree labels.}
  \item{...}{Additional arguments passed to the underlying print or plot method.}
}
\value{
The print methods return their input invisibly. The summary methods return a
list of model dimensions, maturity diagnostics, split information, and timing.
The plot method invisibly returns the node coordinates.
}
\description{
Compact methods for inspecting tree structure, forest maturity, root-variable
stability, and split-based importance.
}
\keyword{survival}
\keyword{tree}
