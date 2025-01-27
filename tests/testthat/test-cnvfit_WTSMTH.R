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

data <- prep(cnv, Y, Z, rare.out = 0.03)

test_that("`cvfit_WTSMTH()` returns expected errors", {
  
  expect_error(cvfit_WTSMTH(),
               "`data` must be a `WTsmth.data` object")
  expect_error(cvfit_WTSMTH(NULL),
               "`data` must be a `WTsmth.data` object")
  expect_error(cvfit_WTSMTH(unclass(data)),
               "`data` must be a `WTsmth.data` object")
  
  expect_error(cvfit_WTSMTH(data),
               "`lambda1 must be a numeric vector")
  expect_error(cvfit_WTSMTH(data, NULL),
               "`lambda1 must be a numeric vector")
  expect_error(cvfit_WTSMTH(data, numeric(0)),
               "`lambda1 must be a numeric vector")
  expect_error(cvfit_WTSMTH(data, character(3)),
               "`lambda1 must be a numeric vector")
  
  expect_error(cvfit_WTSMTH(data, lambda1 = 1),
               "`lambda2 must be a numeric vector")
  expect_error(cvfit_WTSMTH(data, lambda1 = 1, NULL),
               "`lambda2 must be a numeric vector")
  expect_error(cvfit_WTSMTH(data, lambda1 = 1, numeric(0)),
               "`lambda2 must be a numeric vector")
  expect_error(cvfit_WTSMTH(data, lambda1 = 1, character(3)),
               "`lambda2 must be a numeric vector")
  
  expect_error(cvfit_WTSMTH(data, 1, 1, family = "gaussIAN"),
               "'arg' should be one of “gaussian”, “binomial”")
  expect_error(cvfit_WTSMTH(data, 1, 1, family = c("gaussian", "Binomial")),
               "'arg' must be of length 1")
  
  expect_error(cvfit_WTSMTH(data, 1, 1, weight = NA, family = "gaussian"),
               "`weight` must be one of eql, keql, wcs, kwcs, wif, kwif")
  expect_error(cvfit_WTSMTH(data, 1, 1, family = "gaussian", weight = "EQL"),
               "`weight` must be one of eql, keql, wcs, kwcs, wif, kwif")
  
  expect_error(cvfit_WTSMTH(data, 1, 1, family = "gaussian", weight = 'eql',
                          cv.control = c("n.core" = 10, "n.fold" = 1e-4, "stratified" = FALSE)),
               "`cv.control` must be a list; allowed elements are n.fold, n.core, and stratified")
  expect_error(cvfit_WTSMTH(data, 1, 1, family = "gaussian", weight = 'eql',
                          cv.control = list("ncore" = 10, "n.fold" = 1e-4, "stratified" = FALSE)),
               "`cv.control` must be a list; allowed elements are n.fold, n.core, and stratified")
  
  expect_error(cvfit_WTSMTH(data, 1, 1, family = "gaussian", weight = 'eql',
                            iter.control = c("max.iter" = 10, "tol.beta" = 1e-4, "tol.loss" = 1e-3)),
               "`iter.control` must be a list; allowed elements are max.iter, tol.beta, and tol.loss")
  expect_error(cvfit_WTSMTH(data, 1, 1, family = "gaussian", weight = 'eql',
                            iter.control = list("maxiter" = 10, "tol.beta" = 1e-4, "tol.loss" = 1e-3)),
               "`iter.control` must be a list; allowed elements are max.iter, tol.beta, and tol.loss")
  
})


test_that("`cvfit_WTSMTH()` returns expected results", {
  
  lambda1 <- c(-2.5, -1.5)
  lambda2 <- c(2, 3)
  
  iter.control <- .testIterControl(list(max.iter = 8L, 
                                        tol.beta = 10^(-3), 
                                        tol.loss = 10^(-6)))  
  cv.control <- .testCVControl(list(n.fold = 3L, 
                                    n.core = 1L, 
                                    stratified = FALSE), "gaussian")
  
  CNV_info <- data$CNVR.info
  
  data2 <- .expandWTsmth(data, weight = "keql")
  
  withr::with_seed(25234, {

    tr <- .nfoldSplit(Y = drop(data2$Y), unique(cnv$ID), cv.control = cv.control)
  
    loss_matrix <- matrix(0.0, nrow = 2, ncol = 2)
  
    idx <- expand.grid(seq_len(cv.control$n.fold), 
                       seq_along(lambda1), 
                       seq_along(lambda2))
  
    idx_loss <- data.frame("fold" = c(1:3L, 1:3L, 1:3L, 1:3L),
                           "lambda1" = rep(rep(lambda1, each = 3), times = 2),
                           "lambda2" = rep(lambda2, each = 6))

    res <- NULL
    for (i in 1L:nrow(idx_loss)) {
    
      subset <- tr != idx_loss$fold[i]
    
      beta_lmd21 <- fit_WTSMTH(data = data2, 
                               lambda1 = idx_loss$lambda1[i], 
                               lambda2 = idx_loss$lambda2[i],
                               family = "gaussian", 
                               iter.control = iter.control,
                               subset = subset)
      head(beta_lmd21)
    
      loss <- .loss(X = data2$XZ[!subset, ], 
                    Y = data2$Y[!subset], 
                    beta = beta_lmd21$coef, 
                    family = "gaussian")
      
    
      loss_matrix[idx[i,2L], idx[i, 3L]] <- loss_matrix[idx[i,2L], idx[i, 3L]] + loss / 3.0
    }
    
    
    min_value <- Inf
    for (i in 1L:2L) {
      for (j in 1L:2L) {
        if (loss_matrix[i,j] < min_value) {
          use <- c(i,j)
          min_value <- loss_matrix[i, j]
        }
      }
    }

    beta_y_cv <- fit_WTSMTH(data = data2, 
                            lambda1 = lambda1[use[1L]], 
                            lambda2 = lambda2[use[2L]],
                            family = "gaussian", 
                            iter.control = iter.control)
    
    loss_matrix <- loss_matrix |> t() |> data.frame()
    
    colnames(loss_matrix) <- lambda1
    loss_matrix$lambda2 <- lambda2
  
    expected <-  list("Loss" = loss_matrix[, c(3,1,2)],
                      "selected.lambda" = c(lambda1[use[1L]], lambda2[use[2L]]),
                      "coef" = beta_y_cv)
  })
  
  expect_equal(withr::with_seed(
    25234,
    cvfit_WTSMTH(data, lambda1, lambda2, "keql",
                 family = "gaussian", 
                 cv.control = list(n.fold = 3L, 
                                   n.core = 1L, 
                                   stratified = FALSE),
                 iter.control = list(max.iter = 8L, 
                                     tol.beta = 10^(-3), 
                                     tol.loss = 10^(-6)))),
    expected)
})
