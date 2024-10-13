#' @noRd
#' @param Y A numeric vector.
#' @param cv.control A list object. Allows user to control cross-validation
#'   procedure. Allowed elements are "n.fold", the number of cross-validation
#'   folds; "n.core", the number of cores to use in procedure; and
#'   "stratified", if TRUE and family = "binomial", the folds will be
#'   stratified (this option should be used if either category is
#'   "rare.")
#'
#' @returns An integer vector containing fold membership assignments
#' 
#' @include helpful_tests.R utils.R
.nfoldSplit <- function(Y, cv.control) {
  
  # take the first value as default
  family <- match.arg(family)

  stopifnot(
    "`Y` must be a vector" = !missing(Y) && .isNumericVector(Y),
    "`cv.control` must be a list; allowed elements are n.fold, n.core, and stratified" = 
      .isNamedList(cv.control, c("n.fold", "n.core", "stratified")),
  )
  
  if (cv.control$stratified) {
    tr <- rep(0L, length(Y))
    tr[Y == 0L] <- sample(seq_len(cv.control$n.fold), sum(Y == 0L), replace = TRUE)
    tr[Y != 0L] <- sample(seq_len(cv.control$n.fold), sum(Y != 0L), replace = TRUE)
  } else{
    tr <- sample(seq_len(cv.control$n.fold), length(Y), replace = TRUE)
  }

  tr
}
