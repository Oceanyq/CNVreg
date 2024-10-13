#' Wrapper around plogis
#' @noRd
#' @param X A numeric vector.
#' 
#' @returns A numeric vector
#' 
#' @importFrom stats plogis
#' @keywords internal
.probPred <- function(X) {
  stopifnot(
    "`X` must be a numeric vector" = !missing(X) && .isNumericVector(X)
  )

  p <- stats::plogis(X)

  if (isTRUE(all.equal(p, rep(1.0, length(X)))) || 
      isTRUE(all.equal(p, rep(0.0, length(X))))) {
    warning("probabilities are all 0 or all 1", call. = FALSE)
  }

  p
}



