cnv <- data.frame("ID" = 11L:20L,
                  "CHR" = rep(22L, 10L),
                  "BP1" = seq(100L, 2000L, 200L),
                  "BP2" = seq(300L, 2200L, 200L),
                  "TYPE" = c(0L, 0L, 1L, 1L,
                             0L, 3L, 3L, 3L,
                             4L, 4L))
Z <- data.frame("ID" = seq(20, 1, by = -1),
                "X1" = withr::with_seed(1234, rnorm(20)),
                "X2" = withr::with_seed(2345, rnorm(20)),
                "X3" = factor(withr::with_seed(3456, sample(1:3, 20, replace = TRUE))))

Y <- data.frame("ID" = seq(20, 1, by = -1), 
                "Y" = withr::with_seed(3456, rnorm(20)))

data <- prep(cnv, Y, Z, rare.out = 0.03)

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
  
  expect_error(fit_WTSMTH(data, 1, 1, weight = NA, family = "gaussian"),
               "`weight` must be one of eql, keql, wcs, kwcs, wif, kwif")
  expect_error(fit_WTSMTH(data, 1, 1, family = "gaussian", weight = "EQL"),
               "`weight` must be one of eql, keql, wcs, kwcs, wif, kwif")
  
  expect_error(fit_WTSMTH(data, 1, 1, family = "gaussian", weight = 'eql',
                          iter.control = c("max.iter" = 10, "tol.beta" = 1e-4, "tol.loss" = 1e-3)),
               "`iter.control` must be a list; allowed elements are max.iter, tol.beta, and tol.loss")
  expect_error(fit_WTSMTH(data, 1, 1, family = "gaussian", weight = 'eql',
                          iter.control = list("maxiter" = 10, "tol.beta" = 1e-4, "tol.loss" = 1e-3)),
               "`iter.control` must be a list; allowed elements are max.iter, tol.beta, and tol.loss")
})

test_that("`fit_WTSMTH()` returns expected results", {
  
  expect_warning(out <- fit_WTSMTH(.expandWTsmth(data, weight = 'eql'), 2, 2, 
                                   weight = 'eql', 
                                   iter.control = list("max.iter" = 10, 
                                                       "tol.beta" = 1e-4, 
                                                       "tol.loss" = 1e-3)),
                 "`weight` input ignored; data already expanded")
  
  expect_equal(fit_WTSMTH(data, 2, 2, weight = 'eql', 
                          iter.control = list("max.iter" = 10, 
                                              "tol.beta" = 1e-4, 
                                              "tol.loss" = 1e-3)), 
               out)
  
})

test_that("`fit_WTSMTH()` mismatched family", {
  expect_error(fit_WTSMTH(data, 2, 2, family = "binomial", weight = "eql",
                          iter.control = list("max.iter" = 10, "tol.beta" = 1e-4, "tol.loss" = 1e-3)),
               "Y is not integer-like")
})

test_that("`fitWTSMTH()` returns expected results", {

  lambda2 <- 2
  lambda1 <- (-2.5)
  iter.control <- .testIterControl(list("max.iter" = 10, 
                                        "tol.beta" = 1e-4, 
                                        "tol.loss" = 1e-3))  
  data2 <- .expandWTsmth(data, weight = 'eql')
  data2$XZ_update <- cbind(1.0, data2$XZ)
  data2$Y_update <- data2$Y
  X_app <- cbind(0.0, 2.0 * data2$A)
  rownames(X_app) <- rownames(data2$A)
  Y_app <- rep.int(0L, nrow(data2$A))
  names(Y_app) <- rownames(data2$A)
  
  expected <- .ctnsSolution(data = data2, X.app = X_app, Y.app = Y_app, lambda1 = lambda1) |> as.numeric()
 
  expect_equal(fit_WTSMTH(data, -2.5, 2, weight = "eql",
                          family = "gaussian",
                          iter.control = list("max.iter" = 10, 
                                              "tol.beta" = 1e-4, 
                                              "tol.loss" = 1e-3))[,"coef"]|> as.numeric(), 
               expected)
  
 
 expected <- data.frame(Vnames= c("(Intercept)", colnames(data$design), colnames(data$Z)),
                        coef=expected)
 row.names(expected) <- NULL
  
  expect_equal(fit_WTSMTH(data2, -2.5, 2, 
                          family = "gaussian",
                          iter.control = list("max.iter" = 10, 
                                              "tol.beta" = 1e-4, 
                                              "tol.loss" = 1e-3))[ ,c("Vnames", "coef")], 
               expected)
})

test_that("`fitWTSMTH()` returns expected results", {
  
  Y <- data.frame("ID" = seq(20, 1, by = -1), 
                  "Y" = withr::with_seed(3456, rbinom(20, 1, 0.5)))
  
  data2 <- prep(cnv, Y, Z, rare.out = 0.03)

  lambda2 <- 2
  lambda1 <- (-2.5)
  iter.control <- .testIterControl(list("max.iter" = 10, 
                                        "tol.beta" = 1e-4, 
                                        "tol.loss" = 1e-3))  
  data2 <- .expandWTsmth(data2, weight = 'eql')
  
  X_app <- cbind(0.0, 2.0 * data2$A)
  rownames(X_app) <- rownames(data2$A)
  Y_app <- rep.int(0L, nrow(data2$A))
  names(Y_app) <- rownames(data2$A)
  
  XN = nrow(data2$XZ)
  AN = nrow(data2$A)
  XZ_app <- cbind(0.0, sqrt(2*(XN+AN)*(2^(lambda2))) * data2$A)
  rownames(XZ_app) <- rownames(data2$A)
  
  expected <- .rwlsSolution(data = data2,
                            X.app = XZ_app, Y.app = Y_app, 
                            lambda1 = lambda1,
                            iter.control = iter.control)
  
  expected <- data.frame(Vnames= c("(Intercept)", colnames(data$design), colnames(data$Z)),
                         coef=expected)
  row.names(expected) <- NULL

  expect_equal(fit_WTSMTH(data2, -2.5, 2, 
                          family = "binomial",
                          iter.control = list("max.iter" = 10, 
                                              "tol.beta" = 1e-4, 
                                              "tol.loss" = 1e-3))[, c("Vnames","coef")], 
               expected)
  expect_equal(fit_WTSMTH(prep(cnv, Y, Z, rare.out = 0.03), -2.5, 2, weight = "eql",
                          family = "binomial",
                          iter.control = list("max.iter" = 10, 
                                              "tol.beta" = 1e-4, 
                                              "tol.loss" = 1e-3))[, c("Vnames","coef")], 
               expected)
})

test_that("`fitWTSMTH()` returns expected results", {
  
  Y <- data.frame("ID" = seq(20, 1, by = -1), 
                  "Y" = withr::with_seed(3456, rbinom(20, 1, 0.5)))
  data2 <- prep(cnv, Y, Z, rare.out = 0.03)
  
  
  lambda2 <- 2
  lambda1 <- (-2.5)
  iter.control <- .testIterControl(list("max.iter" = 10, 
                                        "tol.beta" = 1e-4, 
                                        "tol.loss" = 1e-3))  
  data2 <- .expandWTsmth(data2, weight = 'eql')
  
  data2$Y <- data2$Y[c(1:16,18, 19)]
  data2$XZ <- data2$XZ[c(1:16,18, 19), , drop = FALSE]

  X_app <- cbind(0.0, 2.0 * data2$A)
  rownames(X_app) <- rownames(data2$A)
  Y_app <- rep.int(0L, nrow(data2$A))
  names(Y_app) <- rownames(data2$A)
  XN = nrow(data2$XZ)
  AN = nrow(data2$A)
  XZ_app <- cbind(0.0, sqrt(2*(XN+AN)*(2^(lambda2))) * data2$A)
  rownames(XZ_app) <- rownames(data2$A)
  
  expected <- .rwlsSolution(data = data2,
                            X.app = XZ_app, Y.app = Y_app, 
                            lambda1 = lambda1,
                            iter.control = iter.control)
  expected <- data.frame(Vnames= c("(Intercept)", colnames(data$design), colnames(data$Z)),
                         coef=expected)
  row.names(expected) <- NULL
  

  data2 <- prep(cnv, Y, Z, rare.out = 0.03)
  data2 <- .expandWTsmth(data2, weight = 'eql')
  
  expect_equal(fit_WTSMTH(data2, -2.5, 2, 
                          family = "binomial",
                          iter.control = list("max.iter" = 10, 
                                              "tol.beta" = 1e-4, 
                                              "tol.loss" = 1e-3),
                          subset = c(1:16,18, 19))[, c("Vnames","coef")], 
               expected)
  expect_equal(fit_WTSMTH(prep(cnv, Y, Z, rare.out = 0.03), -2.5, 2, weight = "eql",
                          family = "binomial",
                          iter.control = list("max.iter" = 10, 
                                              "tol.beta" = 1e-4, 
                                              "tol.loss" = 1e-3),
                          subset = c(1:16,18, 19))[, c("Vnames","coef")], 
               expected)
})