#' Simulated data with copy number variants (CNV), Covariates (Cov), and outcomes traits (Y_QT for a continuous outcome
#' and Y_BT for a binary outcome)for the illustration of CNV association analysis with penalized regression in CNVreg.  
#' 
#' @docType data
#' @keywords datasets
#' @aliases CNV Cov Y_QT Y_BT
#' @name CNVCOVY
#' @usage data('CNVCOVY')
#' @format CNVCOVY.RData provides 4 datasets: CNV, Cov, Y_QT, and Y_BT.
#' 
#' \itemize{
#'   \item CNV. A data frame of 2680 CNV records, it has 5 variables:
#'        \itemize{
#'          \item ID. The sample ID of the CNV records in each row. There are 797 unique IDs.
#'          \item CHR. One integer that indicates the chromosome number of CNV records.
#'          \item BP1. A number, show the starting breakpoint of the CNV records. 
#'          \item BP2. A number, show the ending breakpoint of the CNV records. 
#'          \item TYPE. An integer variable, describe how many copies of the CNV present. 
#'        }
#'   \item Cov. A data frame contains covariats of 900 samples 
#'   (including 797 samples in the CNV data set).
#'   Cov has 3 variables.
#'        \itemize{
#'          \item ID. The sample ID of 900 individuals (900 unique IDs).
#'          \item Sex. An integer covariate, sex of each sample: 0 male, 1 female.
#'          \item Age. A numeric covariate, age of each sample. 
#'        }
#'   \item Y_QT and Y_BT are two data frames for outcomes traits. Y_QT contains a continuous trait. 
#'   Y_BT contains a binary trait. Both have 2 variables
#'        \itemize{
#'          \item ID. The sample ID of 900 individuals (900 unique IDs).
#'          \item Y. Y_QT has a numeric variable range (-4.89 -- 16.70)
#'                   Y_BT has an integer variable with controls 0 and cases 1.
#'        }
#' }
#'
NULL
