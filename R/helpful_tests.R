#' Inner functions to test data format meets requirements.
#' 
#' Test if provided object is a numeric vector of specified length
#'
#' @noRd
#' @param x An R object.
#' @param length An integer or NULL. If NULL, length is tested to be non-zero;
#'   if positive, length of x must match.
#'
#' @returns A logical. TRUE if x is a numeric vector of appropriate length
#' @keywords internal
.isNumericVector <- function(x, length = NULL) {
  tst <- is.numeric(x) && is.vector(x, mode = "numeric") && length(x) > 0L
  if (!is.null(length)) tst <- tst && length(x) == length
  tst
}

#' Test if provided object is a character vector of specified length
#'
#' @noRd
#' @param x An R object.
#' @param length An integer or NULL. If NULL, length is tested to be non-zero;
#'   if positive, length of x must match.
#'
#' @returns A logical. TRUE if x is a character vector of appropriate length.
#' @keywords internal
.isCharacterVector <- function(x, length = NULL) {
  if (is.null(length)) {
    is.character(x) && is.vector(x, mode = "character") && length(x) > 0L
  } else {
    is.character(x) && is.vector(x, mode = "character") && length(x) == length
  }
}

#' Test if provided object is an integer vector of specified length
#'
#' @noRd
#' @param x An R object.
#' @param length An integer or NULL. If NULL, length is tested to be non-zero;
#'   if positive, length of x must match.
#'
#' @returns A logical. TRUE if x is an integer vector of appropriate length.
#' @keywords internal
.isIntegerVector <- function(x, length = NULL) {
  if (is.null(length)) {
    is.integer(x) && is.vector(x, mode = "integer") && length(x) > 0L
  } else {
    is.integer(x) && is.vector(x, mode = "integer") && length(x) == length
  }
}

#' Test if provided object is a logical vector of specified length
#'
#' @noRd
#' @param x An R object.
#' @param length An integer or NULL. If NULL, length is tested to be non-zero;
#'   if positive, length of x must match.
#'
#' @returns A logical. TRUE if x is a logical vector of appropriate length.
#' @keywords internal
.isLogicalVector <- function(x, length = NULL) {
  if (is.null(length)) {
    is.logical(x) && is.vector(x, mode = "logical") && length(x) > 0L
  } else {
    is.logical(x) && is.vector(x, mode = "logical") && length(x) == length
  }
}

#' Test if provided object is a list containing at a minimum elements nms
#'
#' @noRd
#' @param x An R object.
#' @param nms A character vector of element names
#'
#' @returns A logical. TRUE if x is a list containing the nms elements
#' @keywords internal
.isNamedList <- function(x, nms) {
  tst <- is.vector(x, mode = "list") && length(x) >= length(nms)
  tst <- tst && all(nms %in% names(x))
  tst
}

.isCNV <- function(CNV) {
  
  if (missing(CNV)) stop("`CNV` must be provided", call. = FALSE)
  
  if (!is.data.frame(CNV) || !all(c("ID", "CHR", "BP1", "BP2", "TYPE") %in% colnames(CNV))) {
    stop("`CNV` data.frame in PLINK format with columns: `ID`, `CHR`, `BP1`, `BP2`,`TYPE`",
         call. = FALSE)
  }
  
  if (!is.vector(CNV$CHR, mode = "integer") ||
      length(unique(CNV$CHR)) != 1L || !all(CNV$CHR %in% c(1L:22L))) {
    stop("`CNV$CHR` must take a single value in [1, 22]", call. = FALSE)
  }
  
  if (!is.vector(CNV$BP1, mode = "numeric") || !all(CNV$BP1 >= 0L)) {
      #||
      #!isTRUE(all.equal(CNV$BP1, round(CNV$BP1)))) {
    stop("`CNV$BP1 must contain non-negative numbers", call. = FALSE)
  }
  
  if (!is.vector(CNV$BP2, mode = "numeric") || !all(CNV$BP2 >= 0L)){
      #||
      #!isTRUE(all.equal(CNV$BP2, round(CNV$BP2)))) {
    stop("`CNV$BP2 must contain non-negative numbers", call. = FALSE)
  }
  
  if (any(CNV$BP1 > CNV$BP2)) {
    stop("`CNV$BP2 must be great than CNV$BP1 for all records", call. = FALSE)
  }
  
  if (!is.vector(CNV$TYPE, mode = "integer") ||
      !all(CNV$TYPE %in% c(0L, 1L, 3L:9999L)) ||
      sum(CNV$TYPE == 2L) != 0L) {
    stop("`CNV$Type` must be integer in [0, 1, 3, ..., 9999]", call. = FALSE)
  }

  if (!isTRUE(all.equal(CNV$BP1, round(CNV$BP1)))) {
    warning("Non-integer CNV break points exist in BP1")
  }
  
  if (!isTRUE(all.equal(CNV$BP2, round(CNV$BP2)))) {
    warning("Non-integer CNV break points exist in BP2")
  }
  if (any(CNV$BP2 - CNV$BP1 < 50)) {
    warning(" There are CNV events with length <50bp, double check the data or units of data")
  }
  
  TRUE
  
}

.testIterControl <- function(iter.control) {
  
  default_iter <- list(max.iter = 8L, tol.beta = 1e-3, tol.loss = 1e-6)
  
  if (missing(iter.control)) stop("`iter.control` must be provided")
  if (is.null(iter.control)) iter.control <- default_iter
  
  if (!is.vector(iter.control, "list")) stop("`iter.control` must be a list")
  
  tst <- names(iter.control) %in% c("max.iter", "tol.beta", "tol.loss")
  if (any(!{tst})) {
    warning("ignored iter.control elements ",
            paste(names(iter.control)[!{tst}], collapse = ", ")
    )
    iter.control <- iter.control[tst]
  }
  
  default_iter[names(iter.control)] <- iter.control
  iter.control <- default_iter
  
  if (!.isNumericVector(iter.control$max.iter, 1L) ||
      iter.control$max.iter < 2.0) {
    stop("iter.control$max.iter is not appropriately defined")
  } else {
    iter.control$max.iter <- as.integer(round(iter.control$max.iter))
  }
  
  if (!.isNumericVector(iter.control$tol.beta, 1L) ||
      iter.control$tol.beta <= 0.0 || iter.control$tol.beta > 1e-2) {
    stop("iter.control$tol.beta is not appropriately defined (0 < tol.beta <= 1e-2)")
  }
  
  if (!.isNumericVector(iter.control$tol.loss, 1L) ||
      iter.control$tol.loss <= 0.0 || iter.control$tol.loss > 1e-3) {
    stop("iter.control$tol.loss is not appropriately defined (0 < tol.loss <= 1e-3)")
  }
  
  iter.control
}

.testCVControl <- function(cv.control, family) {
  
  default_cv <- list(n.fold = 5L, n.core = 1L, stratified = FALSE)
  
  if (missing(cv.control)) stop("`cv.control` must be provided")
  if (is.null(cv.control)) cv.control <- default_cv
  
  if (!is.vector(cv.control, "list")) stop("`cv.control` must be a list")
  
  tst <- names(cv.control) %in% c("n.fold", "n.core", "stratified")
  if (any(!{tst})) {
    warning("ignored cv.control elements ",
            paste(names(cv.control)[!{tst}], collapse = ", ")
    )
    cv.control <- cv.control[tst]
  }
  
  default_cv[names(cv.control)] <- cv.control
  cv.control <- default_cv
  
  if (!.isNumericVector(cv.control$n.fold, 1L) ||
      cv.control$n.fold < 2L) {
    stop("cv.control$n.fold is not appropriately defined.")
  } else {
    cv.control$n.fold <- as.integer(cv.control$n.fold)
  }
  
  if (!.isNumericVector(cv.control$n.core, 1L)) {
    stop("cv.control$n.core is not appropriately defined.")
  } else {
    cv.control$n.core <- max(as.integer(cv.control$n.core), 1L)
  }
  
  if (!.isLogicalVector(cv.control$stratified, 1L)) {
    stop("cv.control$stratified is not appropriately defined.")
  }
  
  if (family == "gaussian") cv.control$stratified <- FALSE
  
  cv.control
}