#' Inner Function `.loss()`
#' 
#' Calculate loss based on X, Y and beta coefficients for a binary/continuous outcome Y.
#' 
#'
#' @noRd
#' @param X A numeric matrix of the CNV fragments and covariates. Must be of
#'   dimension n x p (n is the number of sample size, p is the total number of CNV fragments 
#'   + the total number of covariates). Intercept would be considered during calculation, no
#'   need to include here.
#' @param Y A numeric vector. The binary/continuous outcome for each participant.
#' @param beta A numeric vector. The current estimated coefficients of length 1+p, 
#'   which should include the intercept and
#'   p is the total number of CNV fragments + the total number of covariates.
#' @param family A character. The family of the outcome. Must be one of "gaussian"/"binomial". 
#'   Default value is the "gaussian" for continuous outcomes.
#'   
#' @returns A scalar numeric
#'   
#' @include linearPred.R logLH.R
#'
#' @keywords internal
#' 
.loss <- function(X, Y, beta,  family = c("gaussian", "binomial")) {
  
  # take the first value as default
  family <- match.arg(family)
  
  stopifnot(
    "`X` must be provided" = !missing(X),
    "`Y` must be a numeric vector with matching sample size" = !missing(Y) && 
      is.vector(Y, mode = "numeric") && nrow(X) == length(Y),
    "`beta` vector must be provided" = !missing(beta)
  )
  
  if (family == "gaussian") {
    Y_hat <- .linearPred(X = X, beta = beta)
    loss <- {crossprod(Y_hat - Y)/length(Y)} |> drop()
  } else {
    loss <- (-1/nrow(X)) * .logLH(X = X, Y = Y, beta = beta)
  }
  
  loss
}