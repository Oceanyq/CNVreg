#' Take X and beta, calculate linear combination of X%*%beta,
#' 
#' @noRd
#' @param X A Matrix object with dimension n \times p
#' @param beta A numeric vector of coefficients, {p + 1}
#' 
#' @returns A numeric vector.
.linearPred <- function(X, beta) {
  stopifnot(
    "`X` must be a Matrix with at least 1 column" =
      !missing(X) && inherits(X, "Matrix"),
    "`beta` must be a vector with numeric values, beta has intercept"=
      !missing(beta) && .isNumericVector(beta, ncol(X) + 1L)
  )
  
  {beta[1L] + X %*% beta[-1L]} |> as.matrix() |> drop()
}