.createITVData <- function(CNV) {
  # Identify boundaries of all grid units
  BPs <- c(CNV$BP1, CNV$BP2) |> unique() |> sort()
  
  # Track boundaries for each grid unit
  n_units <- length(BPs) - 1L
  ITV_info <- data.frame("CHR" = rep(CNV$CHR[1L], n_units))
  ITV_info$idx <- seq_len(n_units)
  ITV_info$ITV1 <- BPs[-{n_units + 1L}]
  ITV_info$ITV2 <- BPs[-1L]
  
  ITV_info
}

.createLongData <- function(CNV, itv.data) {
  
  BPs <- c(itv.data$ITV1, max(itv.data$ITV2))
  
  # idx of Match starting BP of each record to its grid unit
  ITVs <- findInterval(CNV$BP1, BPs)
  
  # idx of Match ending BP of each record to its grid unit
  ITVe <- findInterval(CNV$BP2, BPs, left.open = TRUE)
  
  
  # create long data format, for each idx of interval from ITVs to ITVe
  long <- data.frame("ID" = rep(CNV$ID, times = ITVe - ITVs + 1L))
  long$CHR <- CNV$CHR[1L]
  long$idx <- mapply(seq, ITVs, ITVe, MoreArgs = list(by = 1L)) |> unlist()
  long$TYPE <- rep(CNV$TYPE, times = ITVe - ITVs + 1L)
  long$deldup <- ifelse(long$TYPE < 2L, "del", "dup")
  
  itv1 <- itv.data$ITV1[long$idx]
  itv2 <- itv.data$ITV2[long$idx]
  long$AUC <- abs(2.0 - long$TYPE) * abs(itv2 - itv1)
  
  long <- long[!duplicated(long), ]
  
  long[order(long$idx), ]
}

#' Convert CNV Records to Grid
#'
#' @noRd
#' @param CNV A data.frame in PLINK format. Specifically, must contain
#' columns: 
#' \itemize{
#'   \item "ID": character, unique identity for each sample
#'   \item "CHR": integer, range 1-22
#'   \item "BP1": integer, CNV event starting position
#'   \item "BP2": integer, CNV event ending position, each rewcsd must have BP1 <= BP2, CNV at least 1bp(or other unit length)
#'   \item "TYPE": integer, range 0,1, 3,4, and larger allowed, 2(is not allowed)
#'   }
#'   
#' @returns A list object with elements
#' \itemize{
#'   \item CNV_frag_l The input CNV data.frame augments with
#'   \itemize{
#'     \item ID The CNV IDs replicated for each grid unit spanned by the record
#'     \item CHR The chromosome of the CNV data
#'     \item idx Grid unit id spanned by the record
#'     \item TYPE The TYPE value of the record
#'     \item deldup Character indicating if TYPE is < 2 (del) or a dup
#'     \item AUC The abs(2 - TYPE) * width_of_grid
#'   }
#'   \item ITV_info A data.frame. Links the grid information to the CNV data
#'   \itemize{
#'     \item CHR The chromosome
#'     \item idx The grid unit id
#'     \item ITV1 The BP value of the left grid boundary
#'     \item ITV2 The BP value of the right grid boundary
#'.  }
#' }
#' 
#' @include helpful_tests.R
#' @keywords internal
.breakCNV <- function(CNV) {

  stopifnot(
    "`CNV` data.frame in PLINK format with columns: `ID`, `CHR`, `BP1`, `BP2`,`TYPE`" =
      !missing(CNV) && .isCNV(CNV)
  )
  
  CNV <- CNV[order(CNV$ID, CNV$BP1), ]
  
  ITV_info <- .createITVData(CNV)

  # create long data format, for each idx of interval from ITVs to ITVe
  CNV_frag_l <- .createLongData(CNV, ITV_info)

  list("CNV_frag_l" = CNV_frag_l,
       "ITV_info" = ITV_info)
}
