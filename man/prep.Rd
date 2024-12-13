% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prep.R
\name{prep}
\alias{prep}
\title{Preprocess CNV, Covariate, and Outcome Data into Package Format}
\usage{
prep(CNV, Z, Y, rare.out = 0.05)
}
\arguments{
\item{CNV}{A data.frame in PLINK format. Specifically, must contain
columns: 
\itemize{
  \item "ID": character, unique identity for each sample
  \item "CHR": integer, allowed range 1-22 NOTE: only 1 CHR can be present
  \item "BP1": integer, CNV event starting position
  \item "BP2": integer, CNV event ending position, each record must have 
        BP1 <= BP2, CNV at least 1bp (or other unit length)
  \item "TYPE": integer, range 0, 1, 3, 4, and larger allowed, i.e.,
       2 is not allowed.
  }}

\item{Z}{A data.frame. Must include column "ID". All other columns are covariates,
which can be continuous, binary, or categorical variables. At a minimum, 
Z must contain all unique CNV$ID values.}

\item{Y}{A data.frame. Must include column "ID". Must have 2 columns. For binary,
values must be 0 (control) or 1 (case). For continuous, values must be 
real. Y$ID must contain all Z$ID.}

\item{rare.out}{A scalar numeric in (0, 0.5) event rates below which are
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
Converts CNV in Plink format into fragments, filters out rare CNV events,
  aligns data objects, and creates weight matrices.
}
\examples{
# example code

}
\keyword{internal}
