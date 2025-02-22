% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prep.R
\name{prep}
\alias{prep}
\title{Function `prep()`}
\usage{
prep(CNV, Y, Z = NULL, rare.out = 0.05)
}
\arguments{
\item{CNV}{A data.frame in PLINK format. Specifically, must contain
columns: 
\itemize{
  \item "ID": character, unique identity for each sample
  \item "CHR": integer, allowed range 1-22 NOTE: only 1 CHR can be present,
   which means this function processes one chromosome at a time.
  \item "BP1": integer, CNV event starting position, 
  \item "BP2": integer, CNV event ending position, each record must have 
        BP1 <= BP2, CNV at least 1bp (or other unit length)
  \item "TYPE": integer, range 0, 1, 3, 4, and larger allowed, i.e.,
       2 is not allowed.
  }}

\item{Y}{A data.frame. Must include column "ID". Must have 2 columns. For binary,
values must be 0 (control) or 1 (case). For continuous, values must be 
real number. Y$ID must contain all unique CNV$ID. Y and Z have the same IDs.}

\item{Z}{A data.frame. Must include column "ID". All other columns are covariates,
which can be continuous, binary, or categorical variables. At a minimum, 
Z must contain all unique CNV$ID values.}

\item{rare.out}{A scalar numeric in the range [0, 0.5); event rates below this value are
filtered out of the data.}
}
\value{
An S3 object of class "WTsmth.data" extending a list object containing
\itemize{
  \item \code{design} CNV data converted to design matrix.
  \item \code{Z} The processed covariate matrix.
  \item \code{Y} The processed response vector.
  \item \code{weight.structure} A Matrix object. The structure of the weight matrix.
  \item \code{weight.options} A matrix object. Each row is the multiplicative
    vector to obtain each available weight. Specifically, the A matrix is
    obtained as weight_option[i, ] * weight.structure where i = 1-6 with
    1="eql", 2="keql", 3="wcs", 4="kwcs", 5="wif", and 6="kwif".
 \item \code{CNVR.info} A data.frame containing details about the fragment
   structure.
}
}
\description{
Function `prep()` converts an individual's CNV events within a genomic 
 region (from one chromosome) to a CNV profile curve, further processes it as CNV fragments, 
 and filter out rare fragments. It analyzes the adjacency relationship 
 between CNV fragments and prepares different options of weight matrices
 for the penalized regression analysis. In addition, it 
 formats data into `WTsmth.data` format for regression analysis.
}
\examples{
# Note that the example data set is smaller than actual CNV data to accommodate 
# a fast example.
# Real data analysis would take a little bit longer time to run. 

# load provided illustrative toy dataset with a continuous outcome and a 
# binary outcome
library("CNVreg")
data("CNVCOVY")
#prepare data format for regression analysis

## first try with the continuous outcome Y_QT
frag_data <- prep(CNV = CNV, Y = Y_QT, Z = Cov, rare.out = 0.05)

## Second, to prepare for the binary outcome Y_BT,
# We can directly replace frag_data$Y with Y_BT in the correct format.

rownames(Y_BT) <- Y_BT$ID
#
##order the sample in Y_BT as in frag_data$Y and name it
frag_data$Y <- Y_BT[names(frag_data$Y), "Y"] |> drop()
names(frag_data$Y) <- rownames(frag_data$Z) 

# Or, we can repeat the procedure using prep() function with Y=Y_BT
# frag_data <- prep(CNV = CNV, Y = Y_BT, Z = Cov, rare.out = 0.05)


}
