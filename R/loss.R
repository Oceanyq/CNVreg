#' Calculate loss based on beta coefficients, x and y (binary/continuous).
#' takes X(CNV+Cov),  Y(outcome), A(weight matrix), lmd1(candidate lambda1), lmd2 as input
#' options family="gaussian"/"binomial"
#' output coefficients
#'
#' @noRd
#' @param X A numeric matrix. The CNV fragments and covariates. Must be of
#'   dimension n x p, i.e., should not include an intercept.
#' @param Y A numeric vector. The binary or continuous outcome for each
#'   participant.
#' @param beta A numeric vector. The current estimated parameters. Must be of
#'   length p + 1.
#' @param family A character. The family of the outcome. Must be one of 
#'   "gaussian", "binomial"
#'   
#' @returns A scalar numeric
#'   
#' @include linearPred.R logLH.R
#'
#' @keywords internal
.loss <- function(X, Y, beta,  family = c("gaussian", "binomial")) {

  # take the first value as default
  family <- match.arg(family)

  stopifnot(
    "`X` must be provided" = !missing(X),
    "`Y` must be a numeric vector" = !missing(Y) && 
      is.vector(Y, mode = "numeric") && nrow(X) == length(Y),
    "`beta` vector must be provided" = !missing(beta)
  )

  if (family == "gaussian") {
    Y_hat <- .linearPred(X = X, beta = beta)
    loss <- {crossprod(Y_hat - Y) / length(Y)} |> drop()
  } else {
    loss <- .logLH(X = X, beta = beta, Y = Y)
  }

  loss
}
