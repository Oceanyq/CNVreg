cnv <- data.frame("ID" = 11L:20L,
                  "CHR" = rep(22L, 10L),
                  "BP1" = seq(1L, 20L, 2L),
                  "BP2" = seq(3L, 22L, 2L),
                  "TYPE" = c(0L, 0L, 1L, 1L,
                             0L, 3L, 3L, 3L,
                             4L, 4L))
Z <- data.frame("ID" = seq(20, 1, by = -1),
                "X1" = withr::with_seed(1234, rnorm(20)),
                "X2" = withr::with_seed(2345, rnorm(20)),
                "X3" = factor(withr::with_seed(3456, sample(1:3, 20, replace = TRUE))))

Y <- data.frame("ID" = seq(20, 1, by = -1), 
                "Y" = withr::with_seed(3456, rnorm(20)))

data <- prep(cnv, Z, Y, rare.out = 0.03)

test_that("`fit_WTSMTH()` returns expected errors", {
  
  expect_error(fit_WTSMTH(),
               "`data` must be a 'WTsmth.data' object")
  expect_error(fit_WTSMTH(NULL),
               "`data` must be a 'WTsmth.data' object")
  expect_error(fit_WTSMTH(unclass(data)),
               "`data` must be a 'WTsmth.data' object")
  
  expect_error(fit_WTSMTH(data),
               "`lambda1 must be a numeric vector")
  expect_error(fit_WTSMTH(data, NULL),
               "`lambda1 must be a numeric vector")
  expect_error(fit_WTSMTH(data, numeric(0)),
               "`lambda1 must be a numeric vector")
  expect_error(fit_WTSMTH(data, character(3)),
               "`lambda1 must be a numeric vector")
  
  expect_error(fit_WTSMTH(data, lambda1 = 1),
               "`lambda2 must be a numeric vector")
  expect_error(fit_WTSMTH(data, lambda1 = 1, NULL),
               "`lambda2 must be a numeric vector")
  expect_error(fit_WTSMTH(data, lambda1 = 1, numeric(0)),
               "`lambda2 must be a numeric vector")
  expect_error(fit_WTSMTH(data, lambda1 = 1, character(3)),
               "`lambda2 must be a numeric vector")
  
  expect_error(fit_WTSMTH(data, 1, 1, family = "gaussIAN"),
               "'arg' should be one of “gaussian”, “binomial”")
  expect_error(fit_WTSMTH(data, 1, 1, family = c("gaussian", "Binomial")),
               "'arg' must be of length 1")
  
  expect_error(fit_WTSMTH(data, 1, 1, family = "gaussian"),
               "`iter.control` must be a list; allowed elements are max.iter, tol.beta, and tol.loss")
  expect_error(fit_WTSMTH(data, 1, 1, family = "gaussian",
                          iter.control = c("max.iter" = 10, "tol.beta" = 1e-4, "tol.loss" = 1e-3)),
               "`iter.control` must be a list; allowed elements are max.iter, tol.beta, and tol.loss")
  expect_error(fit_WTSMTH(data, 1, 1, family = "gaussian",
                          iter.control = list("maxiter" = 10, "tol.beta" = 1e-4, "tol.loss" = 1e-3)),
               "`iter.control` must be a list; allowed elements are max.iter, tol.beta, and tol.loss")
})

test_that("`fit_WTSMTH()` returns expected results", {
  
  expect_equal(fit_WTSMTH(data, 2, 2, iter.control = list("max.iter" = 10, "tol.beta" = 1e-4, "tol.loss" = 1e-3)), 
               fit_WTSMTH(.expandWTsmth(data), 2, 2), iter.control = list("max.iter" = 10, "tol.beta" = 1e-4, "tol.loss" = 1e-3))
  
})

fit_WTSMTH <- function(data, lambda1, lambda2,
                       family = c("gaussian", "binomial"), 
                       iter.control = list(max_iter = 8L, 
                                           tol_beta = 10^(-3), 
                                           tol_loss = 10^(-6)),
                       ...) {
  
  extras <- list(...)
  
  # take the first value as default
  family <- match.arg(family)
  
  stopifnot(
    "`data` must be a 'WTsth.data' object" = !missing(data) && inherits(data, "WTsth.data"),
    "`lambda1 must be a numeric vector" = !missing(lambda1) && .isNumericVector(lambda1),
    "`lambda2 must be a numeric vector" = !missing(lambda2) && .isNumericVector(lambda2),
    "`iter.control` must be a list; allowed elements are max.iter, tol.beta, and tol.loss" = 
      .isNamedList(iter.control, c("max.iter", "tol.beta", "tol.loss")),
  )
  
  iter.control <- .testIterControl(iter.control)  
  
  if (family == "binomial") data$Y <- .confirmBinary(data$Y)
  
  if (is.null(data$XZ)) data <- .expandWTsmth(data)
  
  if (!is.null(extras$subset)) {
    data$Y <- data$Y[subset]
    data$XZ <- data$XZ[subset, ]
  }
  
  X_app <- cbind(0.0, sqrt(2.0^lambda2) * data$A)
  rownames(X_app) <- rownames(data$A)
  Y_app <- rep.int(0L, nrow(data$A))
  names(Y_app) <- rownames(data$A)
  
  if (family == "gaussian") {
    .ctnsSolution(data = data, X.app = X_app, Y.app = Y_app, lambda1 = lambda1)
  } else {
    .rwlsSolution(data = data,
                  X.app = X_app, Y.app = Y_app, 
                  lambda1 = lambda1,
                  iter.control = iter.control)
  }
}
