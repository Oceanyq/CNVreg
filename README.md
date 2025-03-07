
 `CNVreg` Package


## Introduction

The *CNVreg* package provides functions to perform copy number variants (CNV) association analysis with penalized regression model. 

This package convert CNVs over a genomic region as a piecewise constant curve to capture the dosage and length of CNVs. The association analysis is then evaluated by regressing outcome traits on all CNV fragments in the region while adjusting for covariates. The corresponding CNV effects are obtained at each genome position. The penalized regression model with Lasso and weighted fusion penalties would perform variable selection and encourage adjacent CNVs to share similar effect size.

This package has 3 main functions: 

* `prep()`: Data preprocessing and format conversion.

* `cvfit_WTSMTH()`: Model fitting and effect estimate with cross-validation(CV). The CV procedure is to tune an optimal model by selecting the best pair of candidate tuning parameters.

* `fit_WTSMTH()`: Model fitting and effect estimate with a given pair of tuning parameters. 


We have a more detailed tutorial for all functions using an example data included in the `CNVreg` package. 
Please see the **CNVreg** vignette for detailed usage of the CNVreg package.
