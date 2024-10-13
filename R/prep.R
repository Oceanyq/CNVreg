.internalStep1 <- function(itv.long, n.samples, rare.out, type) {
  
  # A sparseMatrix of dimension n x n_fragments, where n is the number
  # of unique participant IDs of the del/dup type
  wide_raw_sum <- .wideDataRaw(itv.long = itv.long)
  
  freqs <- .wideFrequency(wide.raw = wide_raw_sum, 
                          sample.size = n.samples, 
                          rare.out = rare.out)
  
  # keep only those fragments that contain non-rare events
  wide_common_data <- wide_raw_sum[, freqs$not.rare.idx, drop = FALSE]
  colnames(wide_common_data) <- paste0(type, colnames(wide_common_data))
  
  weight_deldup <- .weightMatrix(wide.data = wide_raw_sum, 
                                 not.rare.idx = freqs$not.rare.idx, 
                                 freq = freqs$freq)
  
  weight_deldup$wide.data <- wide_common_data
  weight_deldup
}

.mergeMatrices <- function(old, new) {
  
  dimnames <- list(union(rownames(old), rownames(new)), 
                   c(colnames(old), colnames(new)))
  
  mrg_wide <- Matrix::Matrix(0.0, 
                             nrow = length(dimnames[[1L]]), 
                             ncol = length(dimnames[[2L]]), 
                             dimnames = dimnames)
  
  mrg_wide[rownames(old), colnames(old) ] <- old
  mrg_wide[rownames(new), colnames(new) ] <- new
  
  mrg_wide
}

#' Preprocess CNV, Covariate, and Outcome Data into Package Format
#'
#' Converts CNV in Plink format into fragments, filters out rare CNV events,
#'   aligns data objects, and creates weight matrices.
#'
#' @param CNV A data.frame in PLINK format. Specifically, must contain
#' columns: 
#' \itemize{
#'   \item "ID": character, unique identity for each sample
#'   \item "CHR": integer, allowed range 1-22 NOTE: only 1 CHR can be present
#'   \item "BP1": integer, CNV event starting position
#'   \item "BP2": integer, CNV event ending position, each record must have 
#'         BP1 <= BP2, CNV at least 1bp (or other unit length)
#'   \item "TYPE": integer, range 0, 1, 3, 4, and larger allowed, i.e.,
#'        2 is not allowed.
#'   }
#' @param Z A data.frame. Must include column "ID". All other columns are covariates,
#'   which can be continuous, binary, or categorical variables. At a minimum, 
#'   Z must contain all unique CNV$ID values.
#' @param Y A data.frame. Must include column "ID". Must have 2 columns. For binary,
#'   values must be 0 (control) or 1 (case). For continuous, values must be 
#'   real. Y$ID must contain all Z$ID.
#' @param  rare.out  A scalar numeric in (0, 0.5) event rates below which are
#'   filtered out of the data.
#'
#' @returns An S3 object of class "WTsmth.data" extending a list object containing
#' \itemize{
#'   \item \code{design} CNV data converted to design matrix.
#'   \item \code{Z} The processed covariate matrix.
#'   \item \code{Y} The processed response vector.
#'   \item \code{weight.structure} A Matrix object. The structure of the weight matrix.
#'   \item \code{weight.options} A matrix object. Each row is the multiplicative
#'     vector to obtain each available weight. Specifically, the A matrix is
#'     obtained as weight_option[i, ] * weight.structure where i = 1-6 with
#'     1="eql", 2="keql", 3="wcs", 4="kwcs", 5="wif", and 6="kwif".
#'  \item \code{CNVR.info} A data.frame containing details about the fragment
#'    structure.
#' }
#' 
#' @examples
#' # example code
#' 
#' @include breakCNV.R weightMatrix.R wideDataRaw.R wideFrequency.R
#' @import Matrix
#'
#' @keywords internal
prep <- function(CNV, Z, Y, rare.out = 0.05) {
  
  # QUESTION: can Z be empty -- i.e., no covariates used?
  
  stopifnot(
    "`CNV` must be provided" = !missing(CNV) && .isCNV(CNV),
    "`Z` must be a data.frame with a column named `ID`" =
      !missing(Z) && is.data.frame(Z) && {"ID" %in% colnames(Z)} &&
      all(CNV$ID %in% Z$ID),
    "`Y` must be a data.frame with 2 columns" =
      !missing(Y) && is.data.frame(Y) && {"ID" %in% colnames(Y)} &&
      all(Y$ID %in% Z$ID) && all(Z$ID %in% Y$ID) && ncol(Y) == 2L,
    "`rare.out` is a number in (0, 0.5)" = !missing(rare.out) && 
      is.vector(rare.out, "numeric") && length(rare.out) == 1L &&
      rare.out > 0.0 && rare.out < 0.5
  )
  
  # order all provided data using common logic
  CNV <- CNV[order(CNV$ID, CNV$BP1), ]
  Z <- Z[order(Z$ID), ]
  Y <- Y[order(Y$ID), ]
  
  # align Z and Y such that the top elements are in the same order as CNV
  IDs <- c(intersect(Z$ID, CNV$ID), setdiff(Z$ID, CNV$ID))
  idx <- match(IDs, Z$ID)
  Z <- Z[idx, ]
  Y <- Y[idx, ]
  
  # remove ID columns from Z and Y and convert to data matrices with IDs
  # as rownames
  Z$ID <- NULL
  Y$ID <- NULL
  Z <- data.matrix(Z)
  rownames(Z) <- IDs
  Y <- data.matrix(Y)
  rownames(Y) <- IDs
  
  n_samples <- length(IDs)
  
  out_wide_deldup <- Matrix::Matrix(seq_len(n_samples), 
                                    ncol  = 1L,
                                    sparse = TRUE,
                                    dimnames = list(IDs, "mid"))
  
  weight_any <- matrix(0.0, nrow = 0L, ncol = 0L)
  weight_options <- matrix(0.0, nrow = 6L, ncol = 0L,
                           dimnames = list(c("eql", "keql", 
                                             "wcs", "kwcs", 
                                             "wif", "kwif"), NULL))
  CNVRinfo = matrix(nrow = 0L, ncol = 0L)
  
  itv_data <- .breakCNV(CNV)
  
  itv_long <- itv_data$CNV_frag_l
  itv_info <- itv_data$ITV_info

  for (deldup_i in c("del", "dup")) {
    
    tst <- itv_long$deldup == deldup_i

    if (!any(tst)) next
    
    res <- .internalStep1(itv_long[tst, ], n.samples = n_samples, 
                          rare.out = rare.out, type = deldup_i)
    
    out_wide_deldup <- .mergeMatrices(out_wide_deldup, res$wide.data)
    weight_options <- cbind(weight_options, res$weight.options)
    weight_any <- bdiag(weight_any, res$weight.structure)
    
    CNVR_ITV <- merge(res$CNVR.summary, itv_info, by = "idx", all.x = TRUE)
    CNVR_ITV$deldup <- deldup_i
    CNVRinfo <- rbind(CNVRinfo, CNVR_ITV)
    
  }
  
  out_wide_deldup <- out_wide_deldup[ , colnames(out_wide_deldup) != "mid"]
  
  res <- list("design" = out_wide_deldup,
              "Z" = as(Z, "sparseMatrix"),
              "Y" = Y,
              "weight.structure" = weight_any,
              "weight.options" = weight_options,
              "CNVR.info" = CNVRinfo)
  
  class(res) <- "WTsmth.data"
  res
}