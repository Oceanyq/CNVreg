#' input data with interval tags for 1 type of CNV
#' returns CNV intervals (all possible intervals, including those
#'  with no CNV activity) as CNV fragments in wide format (as defined in CNVRular)
#'
#' @noRd
#' @param cnv.long A data.frame object containin
#'   \itemize{
#'     \item ID The CNV IDs replicated for each grid unit spanned by the record
#'     \item CHR The chromosome of the CNV data
#'     \item idx Grid unit id spanned by the record
#'     \item TYPE The TYPE value of the record
#'     \item deldup Character indicating if TYPE is < 2 (del) or a dup
#'     \item AUC The abs(2 - TYPE) * width_of_grid
#'   }
#'
#' @returns A sparseMatrix of dimension n (unique participants) x n_fragments
#' 
#' @import Matrix
#' 
#' @keywords internal
.wideDataRaw <- function(cnv.long) {
  
  stopifnot(
    "`cnv.long` data.frame is missing or not appropriately defined" =
      !missing(cnv.long) && is.data.frame(cnv.long) &&
      all(c("ID", "CHR", "grid.id", "TYPE", "deldup", "AUC") %in% colnames(cnv.long))
  )
  
  cnv.long <- cnv.long[order(cnv.long$ID, cnv.long$grid.id, -cnv.long$AUC), ]
  
  test <- duplicated(cnv.long[, c("ID", "grid.id")])
  
  if (any(test)) {
    warning("One or more samples have more than 1 record for a CNV del/dup event at the same location, kept the highest dosage")
    cnv.long <- cnv.long[!test, ]
  }
  
  #spread data as wide table for GLMNET analysis
  ID_mid <- data.frame("ID" = unique(cnv.long$ID))
  ID_mid$mid <- seq_len(nrow(ID_mid))
  cnv.long <- merge(cnv.long, ID_mid, by = "ID", all.x = TRUE)
  
  cnv.long <- cnv.long[order(cnv.long$grid.id, cnv.long$ID), ]
  
  Matrix::sparseMatrix(i = cnv.long$mid,
                       j = cnv.long$grid.id, 
                       x = cnv.long$AUC, 
                       dims = c(max(cnv.long$mid), 
                                max(cnv.long$grid.id)),
                       dimnames = list(ID_mid$ID,
                                       seq_len(max(cnv.long$grid.id))))
}