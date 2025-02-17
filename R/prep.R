#' Function `prep()` and inner dependent functions.


#dependent functions
#' Y and Z have the same samples
#' CNV can show on part of all samples
#' The function will order the sample according to id
#'  and matches the order in CNV, Z and Y

#' @noRd
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
#' 
#' @returns A list object containing 
#'   \item{CNV} {ordered by character(ID), BP1}
#'   \item{Z} {ordered by CNV$ID, order(setdiff(Z$ID, CNV$ID))}
#'   \item{Y} {ordered by CNV$ID, order(setdiff(Y$ID, CNV$ID))}
#'   
#' @importFrom stats model.matrix
#' @keywords internal
#' 

.orderData <- function(CNV, Z, Y) {
  
  stopifnot(
    "`CNV` must be provided" = !missing(CNV) && .isCNV(CNV),
    "`Z` must be a data.frame with a column named `ID`" =
      !missing(Z) && {is.null(Z) ||
          {is.data.frame(Z) && {"ID" %in% colnames(Z)} &&
              all(as.character(CNV$ID) %in% as.character(Z$ID))}},
    "`Y` must be a data.frame with 2 columns" =
      !missing(Y) && is.data.frame(Y) && {"ID" %in% colnames(Y)} && ncol(Y) >= 2L &&
      {is.null(Z) || 
          {all(as.character(Y$ID) %in% as.character(Z$ID)) && 
              all(as.character(Z$ID) %in% as.character(Y$ID))}}
  )
  
  if (is.null(Z)) Z <- data.frame(ID = Y$ID)
  
  # convert all IDs to common type
  CNV$ID <- as.character(CNV$ID)
  Z$ID <- as.character(Z$ID)
  Y$ID <- as.character(Y$ID)
  
  # order all provided data 
  CNV <- CNV[order(CNV$ID, CNV$BP1), ]
  Z <- Z[order(Z$ID), , drop = FALSE]
  Y <- Y[order(Y$ID), ]
  
  # align Z and Y such that the top elements are in the same order as CNV
  IDs <- c(intersect(Z$ID, CNV$ID), setdiff(Z$ID, CNV$ID))
  idx <- match(IDs, Z$ID)
  Z <- Z[idx, , drop = FALSE]
  Y <- Y[idx, ]
  
  # remove ID columns from Z and Y and convert to data matrices with IDs
  # as rownames
  Z$ID <- NULL
  Y$ID <- NULL
  if (ncol(Z) > 0L) {
    Z <- stats::model.matrix(~., Z)[, -1L, drop=FALSE]
    rownames(Z) <- IDs
  } else {
    Z <- matrix(NA, nrow = nrow(Z), ncol = 0L,
                dimnames = list(IDs, NULL))
  }
  Y <- data.matrix(Y)
  rownames(Y) <- IDs
  
  list("CNV" = CNV, "Z" = Z, "Y" = Y)
  
}

#' Process CNV region with continuous fragments, build wide data shape, remove rare, 
#' build weight matrix 
#' @noRd
#' @param cnv.long CNV data converted to long format
#' @param n.samples Number of samples in dataset
#' @param rare.out Allowed range for rare events
#' @param type The type of CNV event (del/dup)


.internalStep1 <- function(cnv.long, n.samples, rare.out, type) {
  
  # A sparseMatrix of dimension n x n_fragments, where n is the number
  # of unique participant IDs of the del/dup type
  wide_raw_sum <- .wideDataRaw(cnv.long = cnv.long)
  #wide.raw = wide_raw_sum; sample.size = n.samples; rare.out = rare.out
  freqs <- .wideFrequency(wide.raw = wide_raw_sum, 
                          sample.size = n.samples, 
                          rare.out = rare.out)
  
  # keep only those fragments that contain non-rare events
  wide_common_data <- wide_raw_sum[, freqs$not.rare.idx, drop = FALSE]
  colnames(wide_common_data) <- paste0(type, colnames(wide_common_data))
  
  #wide.data = wide_raw_sum; not.rare.idx = freqs$not.rare.idx; freq = freqs$freq
  weight_deldup <- .weightMatrix(wide.data = wide_raw_sum, 
                                 not.rare.idx = freqs$not.rare.idx, 
                                 freq = freqs$freq)
  
  weight_deldup$wide.data <- wide_common_data
  weight_deldup
}

# combine Matrices of different CNVR
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




#' Function `prep()`
#' 
#'  Function `prep()` converts an individual's CNV events within a genomic 
#'  region to a CNV profile curve, further processes it as CNV fragments, 
#'  and filter out rare fragments. It analyzes the adjacency relationship 
#'  between CNV fragments and prepares different options of weight matrices
#'  for the penalized regression analysis. In addition, it 
#'  formats data into `WTsmth.data` format for regression analysis.
#' 

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
#' @param  rare.out  A scalar numeric in (0, 0.5); event rates below this value are
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
#' @export
#' 
#' @examples
#' # Note that the example data set is smaller than actual CNV data to accommodate a fast example.
#' # Real data analysis would take a little bit longer time to run. 
#' 
#' # load provided illustrative toy dataset with a continuous outcome and a binary outcome
#' library("CNVreg")
#' data("CNVCOVY")
#' #prepare data format for regression analysis
#' 
#' ## first try with the continuous outcome Y_QT
#' frag_data <- prep(CNV = CNV, Y = Y_QT, Z = Cov, rare.out = 0.05)
#' 
#' ## Second, to prepare for the binary outcome Y_BT,
#' # We can directly replace frag_data$Y with Y_BT in the correct format.
#' 
#' rownames(Y_BT) <- Y_BT$ID
#' #
#' ##order the sample in Y_BT as in frag_data$Y and name it
#' frag_data$Y <- Y_BT[names(frag_data$Y), "Y"] |> drop()
#' names(frag_data$Y) <- rownames(frag_data$Z) 
#' 
#' # Or, we can also repeat the procedure using prep() function 
#' # frag_data <- prep(CNV = CNV, Y = Y_BT, Z = Cov, rare.out = 0.05)
#' 
#' 
#' @include breakCNV.R weightMatrix.R wideDataRaw.R wideFrequency.R
#' @import Matrix
#'
#' @keywords internal
prep <- function(CNV, Y, Z = NULL, rare.out = 0.05) {
  
  # QUESTION: can Z be empty -- i.e., no covariates used?
  
  stopifnot(
    "`CNV` must be provided" = !missing(CNV),
    "`Z` must be NULL or a data.frame" = is.null(Z) || is.data.frame(Z),
    "`Y` must be provided" = !missing(Y),
    "`rare.out` is a number in (0, 0.5)" = is.vector(rare.out, "numeric") && 
      length(rare.out) == 1L && rare.out > 0.0 && rare.out < 0.5
  )
  
  # order data object to common structure
  ordered_data <- .orderData(CNV = CNV, Z = Z, Y = Y)
  
  # number of samples in the data
  n_samples <- nrow(ordered_data$Z)
  
  wide_deldup <- Matrix::Matrix(seq_len(n_samples), 
                                ncol  = 1L,
                                sparse = TRUE,
                                dimnames = list(rownames(ordered_data$Y), "mid"))
  
  # weight definitions
  weight_any <- matrix(0.0, nrow = 0L, ncol = 0L)
  weight_options <- matrix(0.0, nrow = 6L, ncol = 0L,
                           dimnames = list(c("eql", "keql", 
                                             "wcs", "kwcs", 
                                             "wif", "kwif"), NULL))
  
  # CNV
  CNVRinfo = matrix(nrow = 0L, ncol = 0L)
  
  cnv_long_form <- .breakCNV(ordered_data$CNV)
  
  cnv_long <- cnv_long_form$long.cnv
  grid_info <- cnv_long_form$grid.info
  #deldup_i ="del"
  for (deldup_i in c("del", "dup")) {
    
    tst <- cnv_long$deldup == deldup_i
    
    if (!any(tst)) next
    #cnv.long =cnv_long[tst, ]; n.samples = n_samples; rare.out = rare.out; type = deldup_i
    res <- .internalStep1(cnv.long = cnv_long[tst, ], n.samples = n_samples, 
                          rare.out = rare.out, type = deldup_i)
    
    wide_deldup <- .mergeMatrices(wide_deldup, res$wide.data)
    weight_options <- cbind(weight_options, res$weight.options)
    weight_any <- bdiag(weight_any, res$weight.structure)
    CNVR_ITV <- merge(res$CNVR.summary, grid_info, by = "grid.id", all.x = TRUE)
    CNVR_ITV$deldup <- deldup_i
    CNVRinfo <- rbind(CNVRinfo, CNVR_ITV)
    
  }
  
  wide_deldup <- wide_deldup[ , colnames(wide_deldup) != "mid"]
  
  res <- list("design" = wide_deldup,
              "Z" = Matrix::Matrix(ordered_data$Z[ , , drop = FALSE ], sparse = TRUE),
              "Y" = ordered_data$Y |> drop(),
              "weight.structure" = weight_any,
              "weight.options" = weight_options,
              "CNVR.info" = CNVRinfo)
  
  class(res) <- "WTsmth.data"
  res
}