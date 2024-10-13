#' @noRd
#' @param wide.raw a sparse matrix, result from function Wide_Type_Raw
#' @param sample.size an integer, from Z covariates and Y outcomes file, number of unique IDs of cov
#' @param  rare.out  is a number, to filter in common CNV"
#' 
#' @returns A list object containing "not.rare.idx", a vector of the column
#'   indices of wide.raw (corresponding to fragments) that are not considered 
#'   rare events; and "freq", the number of non-rare events in each fragment.
#'   
#' @include wideDataRaw.R
#' @keywords internal
.wideFrequency <- function(wide.raw, sample.size, rare.out){

  stopifnot(
    "`wide.raw` must be a Matrix object" = !missing(wide.raw) && 
      inherits(wide.raw, "Matrix"),
    "`sample.size` be a scalar" = !missing(sample.size) && 
      is.vector(sample.size, "numeric") && length(sample.size) == 1L,
    "`rare.out` is a number in (0, 0.5)" = !missing(rare.out) && 
      is.vector(rare.out, "numeric") && length(rare.out) == 1L &&
      rare.out > 0.0 && rare.out < 0.5
  )

  # the number of non-zero values for each fragment,  find interval separated by non-CNV regions
  wide_freq <- apply(wide.raw, 2L, function(x)  sum({x > 1e-8}) )
  
  tst <- which(wide_freq > {sample.size * rare.out})
  if (length(tst) == 0L) stop("no cnv regions satisfy rare.out condition", call. = FALSE)

  list("not.rare.idx" = tst,
       "freq" = wide_freq)
}

