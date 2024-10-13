#' input data with interval tags for 1 type of CNV
#' returns CNV intervals (all possible intervals, including those
#'  with no CNV activity) as CNV fragments in wide format (as defined in CNVRular)
#'
#' @noRd
#' @param itv.long A data.frame object containin
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
.wideDataRaw <- function(itv.long) {
  
  stopifnot(
    "`itv.long` data.frame is missing or not appropriately defined" =
      !missing(itv.long) && is.data.frame(itv.long) &&
      all(c("ID", "CHR", "idx", "TYPE", "deldup", "AUC") %in% colnames(itv.long))
  )

  itv.long <- itv.long[order(itv.long$ID, itv.long$idx, -itv.long$AUC), ]

  test <- duplicated(itv.long[, c("ID", "idx")])

  if (any(test)) {
    warning("One or more samples have more than 1 record for a CNV del/dup event at the same location, kept the highest dosage")
    itv.long <- itv.long[!test, ]
  }

  #spread data as wide table for GLMNET analysis
  ID_mid <- data.frame("ID" = unique(itv.long$ID))
  ID_mid$mid <- seq_len(nrow(ID_mid))
  itv.long <- merge(itv.long, ID_mid, by = "ID", all.x = TRUE)

  itv.long <- itv.long[order(itv.long$idx, itv.long$ID), ]

  Matrix::sparseMatrix(i = itv.long$mid,
                       j = itv.long$idx, 
                       x = itv.long$AUC, 
                       dims = c(max(itv.long$mid), 
                                max(itv.long$idx)),
                       dimnames = list(ID_mid$ID,
                                       seq_len(max(itv.long$idx))))
}
