#' Linear regression with Lasso and Weighted fusion penalties for
#' continuous outcomes to get coefficient estimates
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
    "`data` must be a 'WTsth.data' object" = !missing(data) && 
      inherits(data, "WTsth.data") && !is.null(data$XZ),
    " `lambda1 must be a scalar numeric" =
      !missing(lambda1) && .isNumericVector(lambda1, 1L),
  )
  
  intercept_name <- "(Intercept)"
  while (intercept_name %in% colnames(data$XZ)) {
    intercept_name <- sample(LETTERS, 10, TRUE)
  }
  XZ_colnames <- c(intercept_name, colnames(data$XZ))
  
  # continuous outcomes
  Xp1 <- cbind(1.0, data$XZ)
  colnames(Xp1) <- XZ_colnames
  
  X_aug <- rbind(Xp1, X.app)
  Y_aug <- c(data$Y, Y.app)

  p_fac <- rep(1.0, ncol(X_aug))
  
  # set the first penalty factor as 0, for intercept
  p_fac[1L] <- 0
  
  fit_tr <- tryCatch(glmnet::glmnet(x = X_aug, y = Y_aug, intercept = FALSE, 
                                    penalty.factor = p_fac),
                     error = function(e){
                       stop("glmnet encountered errors\n\t", call. = FALSE, e$message)
                     })
  
  glmnet::coef.glmnet(fit_tr, s = 2.0^lambda1, exact = FALSE)[-1L, ]
}