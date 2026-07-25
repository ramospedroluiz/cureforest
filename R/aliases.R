#' Fit a cure-fraction survival tree
#'
#' `curetree()` is the concise front end for [survcure()]. All arguments are
#' passed unchanged.
#'
#' @param ... Arguments passed to [survcure()].
#' @return An object of class `survcure`.
#' @export
curetree <- function(...) {
  survcure(...)
}

#' Fit a cure-fraction random survival forest
#'
#' `cureforest()` is the concise front end for [randomforestcure()]. All
#' arguments are passed unchanged.
#'
#' @param ... Arguments passed to [randomforestcure()].
#' @return An object of class `randomforestcure`.
#' @export
cureforest <- function(...) {
  randomforestcure(...)
}
