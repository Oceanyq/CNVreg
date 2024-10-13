#' Construct A Weight Matrix
#'
#' @noRd
#' @param wide.data A sparseMatrix of dimension n (unique participants) x n_fragments;
#'   as returned by .wideDataRaw().
#' @param not.rare.idx An integer vector. The fragments (columns of wide.data)
#'   that contain non-rare events; as returned by .wideFrequency()
#' @param freq A numeric vector. The frequency of each fragment (column of
#'   wide.data).
#'
#' @returns A list object containing
#'   \itemize{
#'     \item \code{weight.structure} A Matrix object defining the general
#'       structure of the weight matrix. (n.fragment - 1 x n.fragment)
#'     \item \code{weight.options} A matrix object 6 x n.fragment - 1
#'       The multiplier to convert weight.structure to the desired weight matrix.
#'     \item \code{CNVR.summary} A matrix. Summary information about the selected
#'       fragments. Contains an id, fragment index, and fragment frequency
#'       each "id" may have multiple rows
#'   }
#'   
#' @import Matrix
#' 
#' @keywords internal
.weightMatrix <- function(wide.data, not.rare.idx, freq) {

  stopifnot(
    "`wide.data` must be a Matrix" = !missing(wide.data) && 
      inherits(wide.data, "Matrix"),
    "`not.rare.idx` must be a vector" =
      !missing(not.rare.idx) && is.vector(not.rare.idx, "integer"),
    "`freq` must be a vector of length ncol(wide.data)" =
      !missing(freq) && is.vector(freq, "integer") && 
      length(freq) == ncol(wide.data)
  )
  
  CNV_summary <- matrix(nrow = 0L, ncol = 3L)

  weight_structure <- Matrix::Matrix(nrow = 0L, ncol = 0L, sparse = TRUE)
  
  weight_options <- matrix(nrow = 6L, ncol = 0L,
                           dimnames = list(c("eql", "keql",
                                             "wcs", "kwcs",
                                             "wif", "kwif"), NULL))
  
  seq_s <- 1L
  CNVR_id <- 1L
  while (length(not.rare.idx) > 0L) {
    
    diffs <- diff(not.rare.idx)
    
    # define a sequence of fragments with ending positions
    seq_e <- min(length(diffs) + 1L, 
                 which(diffs > 1L)[1L], na.rm = TRUE)
    
    seq_frag <- not.rare.idx[seq_s:seq_e]
    n_frag <- length(seq_frag)
    
    CNVR_block <- wide.data[, seq_frag, drop = FALSE]
    
    CNV_summary <- rbind(CNV_summary, 
                         cbind(CNVR_id, seq_frag, freq[seq_frag]))
    
    if (n_frag == 1L) {
      weight_structure <- Matrix::bdiag(weight_structure, 0.0)
      weight_options <- cbind(weight_options, rep(0.0, 6L))
    } else {
      CNVR_block <- CNVR_block[rowSums(CNVR_block) > 10e-8, ]
      K <- nrow(CNVR_block)
      
      wgt_matrix <- matrix(nrow = n_frag - 1L, ncol = n_frag)
      
      for (i in seq_len(n_frag - 1L)) {
        term_wcs <- {{CNVR_block[, i] %*% CNVR_block[, i + 1L]} / 
          sqrt(crossprod(CNVR_block[, i]) * crossprod(CNVR_block[, i + 1L]))} |> 
          drop()
        term_wif <- 1.0 / sum(CNVR_block[, i] & CNVR_block[, i + 1L])

        weight_options <- cbind(weight_options,
                                c(1.0, K * 1.0, 
                                  term_wcs, K * term_wcs, 
                                  term_wif, K * term_wif))
        
        wgt_matrix[i, i] <- -1.0
        wgt_matrix[i, i + 1L] <- 1.0
      }
      
      weight_structure <- Matrix::bdiag(weight_structure, wgt_matrix)
    }
    
    CNVR_id = CNVR_id + 1L
    
    not.rare.idx <- not.rare.idx[-(seq_s:seq_e)]
    
  }
  
  colnames(CNV_summary) <- c("CNV.id", "idx", "freq")
  
  # factor of 2 is because there are 2 elements in each row
  norm_nonzero <- rowSums(abs(weight_options) > 1e-12, na.rm = TRUE) * 2.0
  norm_sum <- rowSums(abs(weight_options)) * 2.0
  norm_sum[norm_sum < 1e-8] <- 1.0
  
  weight_options <- weight_options * norm_nonzero / norm_sum
  
  list("weight.structure" = weight_structure,
       "weight.options" = weight_options,
       "CNVR.summary" = CNV_summary)
}
