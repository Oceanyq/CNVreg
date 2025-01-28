#' Inner Functions for regression related data format transformation
#' 
#' Combine CNV design with covariates
#' 
#' @noRd
#' @param data The data object as returned by prep()
#' @returns The XZ matrix as a Matrix object
#' @import Matrix
#' @keywords internal
.combineXZ <- function(data) {
  
  stopifnot(
    "`data` must be a WTsmth.data object" = !missing(data) && inherits(data, "WTsmth.data")
  )
  
  X <-  data$design
  
  mrg_dimnames <- list(rownames(data$Z), c(colnames(X), colnames(data$Z)))
  
  idx <- match(rownames(X), rownames(data$Z))
  if (any(is.na(idx))) { stop("IDs in X not found in Z", call. = FALSE) }
  
  XZ <- Matrix::Matrix(0.0, 
                       nrow = length(mrg_dimnames[[1L]]), 
                       ncol = length(mrg_dimnames[[2L]]), 
                       dimnames = mrg_dimnames)
  XZ[idx, seq_len(ncol(X))] <- X
  XZ[, seq_len(ncol(data$Z)) + ncol(X)] <- data$Z
  XZ[is.na(XZ)] <- 0.0
  XZ
}

#' @noRd
#' @param data The data object as returned by prep()
#' @param weight The specific weight requested
#' @returns The A matrix as a Matrix object
#' @import Matrix
#' @keywords internal
.createA <- function(data, weight) {
  
  stopifnot(
    "`data` must be a WTsmth.data object" = !missing(data) && 
      inherits(data, "WTsmth.data") && "XZ" %in% names(data),
    "`weight` must be a character" = !missing(weight) && 
      is.vector(weight, mode = "character") && length(weight) == 1L
  )
  
  if (!weight %in% rownames(data$weight.options)) {
    stop("weight not recognized", call. = FALSE)
  }
  A <- data$weight.options[weight, ] * data$weight.structure
  rownames(A) <- paste0("A", seq_len(nrow(A)))
  
  A_matrix = Matrix::Matrix(0.0, 
                            nrow = nrow(A), ncol = ncol(data$XZ), 
                            dimnames = list(rownames(A), colnames(data$XZ)),
                            sparse = TRUE)
  A_matrix[, 1L:ncol(A)] <- A
  A_matrix[is.na(A_matrix)] <- 0.0
  
  A_matrix
}

.expandWTsmth <- function(data, weight) {
  
  stopifnot(
    "`data` must be a WTsmth.data object" = !missing(data) && 
      inherits(data, "WTsmth.data"),
    "`weight` must be one of eql, keql, wcs, kwcs, wif, kwif" =
      !missing(weight) && .isCharacterVector(weight, 1L) &&
      weight %in% c("eql", "keql", "wcs", "kwcs", "wif", "kwif")
  )
  
  data$XZ <- .combineXZ(data)
  data$A <- .createA(data, weight)
  data
}

.confirmBinary <- function(Y) {
  dn <- names(Y)
  if (!.isIntegerVector(Y) || !all(Y %in% c(0L, 1L))) {
    if (is.numeric(Y) && !isTRUE(all.equal(Y, round(Y)))) {
      stop("Y is not integer-like", call. = FALSE)
    }
    Y <- factor(Y)
    if (nlevels(Y) != 2L) {
      warning("Y does not appear to be binary", call. = FALSE)
    }
    Y <- {unclass(Y) - 1L} |> as.integer()
  }

  names(Y) <- dn
  Y
  
}

.confirmContinuous <- function(Y) {
  if (length(unique(Y)) <= 2) {
    warning("Y does not appear to be continuous", call. = FALSE)
  }
  Y
}