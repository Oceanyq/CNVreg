#' Calculate log likelihood for binomial logistic type of regression data
#' Take X, beta, and binary Y,  calculate log likelihood
#' 
#' @noRd
#' @param X a matrix n \times p, numeric
#' @param beta a vector of coefficients, p \times 1, numeric
#' @param Y a vector of (0,1) combination
#' 
#' @returns A scalar numeric
#' 
#' @include utils.R linearPred.R
#' @keywords internal
.logLH <- function(X,  beta,  Y) {
  
  stopifnot(
    "`X` must be provided" = !missing(X),
    "`beta` must be provided" = !missing(beta),
    "`Y` must be a vector of number (0,1)" = !missing(Y) && 
      .isNumericVector(Y, nrow(X))
  )
  
  Y <- .confirmBinary(Y)
  
  if (all(Y == 0L) || all(Y == 1L)) {
    warning("All Y are all equal", call. = FALSE)
  }
  
  linear_pred <- .linearPred(X = X, beta = beta)
  sum(linear_pred * Y - log(1.0 + exp(linear_pred)))
}