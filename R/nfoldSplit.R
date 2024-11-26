#' @noRd
#' @param Y A numeric vector.
#' @param X.id A vector.
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
.nfoldSplit <- function(Y, X.id, cv.control) {
  
  stopifnot(
    "`Y` must be a vector" = !missing(Y) && .isNumericVector(Y),
    "`X.id` must be a vector" = !missing(X.id) && is.vector(X.id) && !is.list(X.id),
    "`cv.control` must be a list; allowed elements are n.fold, n.core, and stratified" = 
      !missing(cv.control) && .isNamedList(cv.control, c("n.fold", "n.core", "stratified"))
  )
  
  Y_not_in_x <- !{names(Y) %in% X.id}
  
  if (cv.control$stratified) {
    if (sum(!Y_not_in_x & Y == 0L) < cv.control$n.fold ||
        sum(!Y_not_in_x & Y != 0L) < cv.control$n.fold) {
      warning("cannot stratify -- too few IDs in X")
      cv.control$stratified <- FALSE
    }
  }
  
  if (cv.control$stratified) {
    Y <- .confirmBinary(Y)
    tr <- rep(0L, length(Y))
    tr[Y == 0L & Y_not_in_x] <- sample(seq_len(cv.control$n.fold), sum(Y == 0L & Y_not_in_x), replace = TRUE)
    tr[Y != 0L & Y_not_in_x] <- sample(seq_len(cv.control$n.fold), sum(Y != 0L & Y_not_in_x), replace = TRUE)
    if (any(!Y_not_in_x)) {
      for (i in 1L:5L){
        tr[Y == 0L & !Y_not_in_x] <- sample(seq_len(cv.control$n.fold), sum(Y == 0L & !Y_not_in_x), 
                                            replace = sum(Y == 0L & !Y_not_in_x) > cv.control$n.fold)
        tr[Y != 0L & !Y_not_in_x] <- sample(seq_len(cv.control$n.fold), sum(Y != 0L & !Y_not_in_x), 
                                            replace = sum(Y != 0L & !Y_not_in_x) > cv.control$n.fold)
        
        if (all(seq_len(cv.control$n.fold) %in% tr[Y == 0L & !Y_not_in_x]) &
            all(seq_len(cv.control$n.fold) %in% tr[Y != 0L & !Y_not_in_x])) break
        if (i == 5L) warning("unable to stratify IDs in x across folds")
      }
    }
  } else{
    tr <- sample(seq_len(cv.control$n.fold), length(Y), replace = TRUE)
  }
  
  tr
}