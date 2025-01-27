#' Inner Function `.ctnsSolution()`
#' 
#' Penalized regression with Lasso and weighted fusion penalties for a continuous 
#' outcome. Specifically, it performs linear regression with Lasso penalties for
#' augmented data (which incorporates weighted fusion penalty and a given lambda_2)
#' for a continuous outcome. Get coefficient estimates at the given lambda 1. 
#'
#' @noRd
#' @param X.aug A sparseMatrix. CNV fragments (X) and covariates (Z),
#'   ID as rownames, unique ID colnames are fragments and covariant names
#'   each row represents a unique sample , covariates and cnv are matched by ID
#'   CNV fragments columns: 0 or positive values, length * dosage
#'   covariates columns: all numeric values (convert to numeric vales )
#' @param Y.aug A numeric vector.
#' @param lambda1 A scalar numeric
#'
#' @include helpful_tests.R
#' @importFrom glmnet glmnet coef.glmnet
#'
#' @keywords internal
.ctnsSolution <- function(data, X.app, Y.app, lambda1) {
  
  stopifnot(
    "`data` must be a 'WTsmth.data' object" = !missing(data) && 
      inherits(data, "WTsmth.data") && !is.null(data$XZ),
    " `lambda1 must be a scalar numeric" =
      !missing(lambda1) && .isNumericVector(lambda1, 1L)
  )
  

  XZ_aug <- rbind(data$XZ_update, X.app)
  Y_aug <- c(data$Y_update, Y.app) 

  
  p_fac <- rep(1.0, ncol(XZ_aug))
  
  # set the first penalty factor as 0, for intercept
  p_fac[1L] <- 0
  
  fit_tr <- tryCatch(glmnet::glmnet(x = XZ_aug, y = Y_aug, intercept = FALSE, 
                                    penalty.factor = p_fac),
                     error = function(e){
                       stop("glmnet encountered errors\n\t", call. = FALSE, e$message)
                     })
  
  glmnet::coef.glmnet(fit_tr, s = 2.0^lambda1, exact = FALSE)[-1L, ]
}