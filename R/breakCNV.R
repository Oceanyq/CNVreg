#' Inner Function `.breakCNV()` and the dependent functions
#' 
#' Use CNV breakpoints (BP1 and BP2) on the same chromosome to construct CNV fragments/grids.
#' Process CNV data in PLINK format 1 chromosome at a time using some dependent functions. 



#' Create CNV fragments (grids) and define the boundary of each fragment.
#'  
#' @noRd
#' @param CNV A data.frame in PLINK format. Specifically, must contain
#' columns: 
#' \itemize{
#'   \item "ID": character, unique identity for each sample
#'   \item "CHR": integer, range 1-22
#'   \item "BP1": integer, CNV event starting position
#'   \item "BP2": integer, CNV event ending position, each record must have BP1 <= BP2, i.e., data at least 1bp (or data can be other unit length)
#'   \item "TYPE": integer, range 0, 1, 3, 4, and larger allowed, 2 is not allowed.
#'   }
#'
#' @returns A list object containing
#'   \item{CHR}{The chromosome under analysis}
#'   \item{id}{An integer vector of ids for each grid unit}
#'   \item{lower.boundary}{An integer vector of the lower boundary of each grid}
#'   \item{upper.boundary}{An integer vector of the upper boundary of each grid}
#' @keywords internal
#' 
.createGrid <- function(CNV) {
  
  # Identify boundaries of all grid units
  BPs <- c(CNV$BP1, CNV$BP2) |> unique() |> sort()
  
  # Track boundaries for each grid unit
  n_units <- length(BPs) - 1L
  ITV_info <- data.frame("CHR" = rep(CNV$CHR[1L], n_units))
  ITV_info$grid.id <- seq_len(n_units)
  ITV_info$lower.boundary <- BPs[-{n_units + 1L}]
  ITV_info$upper.boundary <- BPs[-1L]
  
  ITV_info
}

#'
#'  Construct a long-shaped matrix linking CNV BP1/BP2 to units of the fragments(grid),
#'  One record of CNV can match to >=1 grid(s).
#'
#' @noRd
#' @param CNV A data.frame in PLINK format. Specifically, must contain
#' columns: 
#' \itemize{
#'   \item "ID": character, unique identity for each sample
#'   \item "CHR": integer, range 1-22
#'   \item "BP1": integer, CNV event starting position
#'   \item "BP2": integer, CNV event ending position, each record must have BP1 <= BP2, i.e., data at least 1bp (or data can be other unit length)
#'   \item "TYPE": integer, range 0, 1, 3, 4, and larger allowed, 2 is not allowed.
#'   }
#' @param grid A list object containing the grid details
#' \itemize{
#'   \item CHR The chromosome under analysis
#'   \item id An integer vector of ids for each grid unit
#'   \item lower.boundary An integer vector of the lower boundary of each grid
#'   \item upper.boundary An integer vector of the upper boundary of each grid
#' }
#' 
#' @returns A data.frame ordered by grid.id containing
#'   \item{ID}{The sample ID}
#'   \item{CHR}{The chromosome under analysis}
#'   \item{grid.id}{The id of a grid unit}
#'   \item{TYPE}{The CNV value for the grid unit}
#'   \item{deldup}{A character, del/dup, indicating if deletion or duplication}
#'   \item{AUC}{Numeric AUC value}
#'
#' @keywords internal
.createLongData <- function(CNV, grid) {
  
  BPs <- c(grid$lower.boundary, max(grid$upper.boundary))
  
  # idx of Match starting BP of each record to its grid unit
  starts <- findInterval(CNV$BP1, BPs)
  
  # idx of Match ending BP of each record to its grid unit
  ends <- findInterval(CNV$BP2, BPs, left.open = TRUE)
  
  # create long data format, for each idx of interval from ITVs to ITVe
  long <- data.frame("ID" = rep(CNV$ID, times = ends - starts + 1L))
  long$CHR <- CNV$CHR[1L]
  long$grid.id <- mapply(seq, starts, ends, MoreArgs = list(by = 1L)) |> unlist()
  long$TYPE <- rep(CNV$TYPE, times = ends - starts + 1L)
  long$deldup <- ifelse(long$TYPE < 2L, "del", "dup")
  
  lower_boundary <- grid$lower.boundary[long$grid.id]
  upper_boundary <- grid$upper.boundary[long$grid.id]
  long$AUC <- abs(2.0 - long$TYPE) * abs(upper_boundary - lower_boundary)
  
  long <- long[!duplicated(long), ]
  
  long[order(long$deldup, long$grid.id), ]
}

#' Convert CNV Records to Grids/fragments
#'
#' @noRd
#' @param CNV A data.frame in PLINK format. Specifically, must contain
#' columns: 
#' \itemize{
#'   \item "ID": character, unique identity for each sample
#'   \item "CHR": integer, range 1-22
#'   \item "BP1": integer, CNV event starting position
#'   \item "BP2": integer, CNV event ending position,  each record must have BP1 <= BP2, i.e., data at least 1bp (or data can be other unit length)
#'   \item "TYPE": integer, range 0, 1, 3, 4, and larger allowed, 2 is not allowed.
#'   }
#'   
#' @returns A list object with elements
#' \itemize{
#'   \item long.cnv CNV data converted to long format
#'   \itemize{
#'     \item{ID}{The sample ID}
#'     \item{CHR}{The chromosome under analysis}
#'     \item{grid.id}{The id of a grid unit}
#'     \item{TYPE}{The CNV value for the grid unit}
#'     \item{deldup}{A character, del/dup, indicating if deletion or duplication}
#'     \item{AUC}{Numeric AUC value}
#'   }
#'   \item grid.info A data.frame. Links the grid information to the CNV data
#'   \itemize{
#'     \item CHR The chromosome
#'     \item id The grid unit id
#'     \item lower.boundary The BP value of the left grid boundary
#'     \item upper.boundary The BP value of the right grid boundary
#'.  }
#' }
#' 
#' @include helpful_tests.R
#' @keywords internal
#' 
#' 
#' Construct CNV fragments for CNVs in PLINK format, spread CNV records to a long table,
#' and output the long-shaped data.frame for all CNV fragments and the fragment boundary information. 
.breakCNV <- function(CNV) {
  
  stopifnot(
    "`CNV` data.frame in PLINK format with columns: `ID`, `CHR`, `BP1`, `BP2`,`TYPE`" =
      !missing(CNV) && .isCNV(CNV)
  )
  
  # create grid
  grid_info <- .createGrid(CNV)
  
  # create long data format
  CNV_frag_l <- .createLongData(CNV = CNV, grid = grid_info)
  
  list("long.cnv" = CNV_frag_l,
       "grid.info" = grid_info)
}