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
                "Y" = withr::with_seed(3456, rbinom(20, 1, 0.4)))

data <- prep(cnv, Z, Y, rare.out = 0.03)
data$weight.options[] <- withr::with_seed(3245, rnorm(length(data$weight.options)))

test_that("`rwlsSolution()` returns expected errors", {
  
  expect_error(.rwlsSolution(), "`data` must be a 'WTsth.data' object")
  expect_error(.rwlsSolution(data), "`data` must be a 'WTsth.data' object")
  data <- .expandWTsmth(data, "keql")
  
  expect_error(.rwlsSolution(data, lambda1 = 1:10),
               "`lambda1 must be a scalar numeric")
  expect_error(.rwlsSolution(data, lambda1 = NA_character_),
               "`lambda1 must be a scalar numeric")
  
  expect_error(.rwlsSolution(data, lambda1 = 0.2, 
                             iter.control = c("max.iter" = 200, "tol.beta" = 1e-4, "tol.loss" = 1e-2)),
               "`iter.control` must be a list; allowed elements are max.iter, tol.beta, and tol.loss")
  expect_error(.rwlsSolution(data, lambda1 = 0.2, 
                             iter.control = list("Max.iter" = 200, "tol.beta" = 1e-4, "tol.loss" = 1e-2)),
               "`iter.control` must be a list; allowed elements are max.iter, tol.beta, and tol.loss")
  
})

test_that("`.rwlsSolution()` returns expected results", {
  
  data <- .expandWTsmth(data, "keql")
  lambda1 <- 1e-2
  iter.control <- list("max.iter" = 0L, "tol.beta" = 1e-4, "tol.loss" = 1e-2)
  
  ZN = nrow(data$XZ)
  Xp = ncol(data$XZ)
  
  fit_yb_init <- glmnet::glmnet(x = data$XZ, y = data$Y, family = "binomial")
  
  beta_cur <- glmnet::coef.glmnet(fit_yb_init, s = 2^lambda1) |> drop()

  expect_warning(out <- .rwlsSolution(data, lambda1 = lambda1, iter.control = iter.control),
                 "maximum iterations reached")
  expect_equal(out, beta_cur)
})

test_that("`.rwlsSolution()` returns expected results", {
  
  data <- .expandWTsmth(data, "keql")
  lambda1 <- log2(0.129000)
  iter.control <- list("max.iter" = 1L, "tol.beta" = 1e-4, "tol.loss" = 1e-2)
  
  ZN = nrow(data$XZ)
  Xp = ncol(data$XZ)
  
  fit_yb_init <- glmnet::glmnet(x = data$XZ, y = data$Y, family = "binomial")

  beta_cur <- glmnet::coef.glmnet(fit_yb_init, s = 2^lambda1) |> drop()

  loss <- {-1.0 / nrow(data$XZ)} * 
    .loss(X = data$XZ, Y = drop(data$Y), beta = beta_cur, family = "binomial")

  # prepare linear prediction and probability prediction for iteration update
  linear_pred <- .linearPred(data$XZ, beta_cur)
  prob_pred <- .probPred(linear_pred)
    
  # update u*, v*
  u <- linear_pred + {data$Y - prob_pred} / {prob_pred * {1.0 - prob_pred}}
    
  sqrt_v <- sqrt(prob_pred * {1.0 - prob_pred})
    
  beta_cur <- .ctnsSolution(data = data, X.app = cbind(1.0, data$XZ) * sqrt_v, Y.app = u * sqrt_v, lambda1 = lambda1)
    
  expect_warning(out <- .rwlsSolution(data, X.app = NULL, Y.app = NULL,
                                      lambda1 = lambda1, iter.control = iter.control),
                 "maximum iterations reached")
  expect_equal(out, beta_cur)
})

test_that("`.rwlsSolution()` returns expected results", {
  
  data <- .expandWTsmth(data, "keql")
  lambda1 <- log2(0.129000)
  iter.control <- list("max.iter" = 1L, "tol.beta" = 1e-4, "tol.loss" = 1e-2)
  
  ZN = 20
  Xp = ncol(data$XZ)
  
  X.app <- withr::with_seed(1342, matrix(rnorm(Xp+1L), 1, Xp + 1L))
  Y.app <- 0.3
  
  fit_yb_init <- glmnet::glmnet(x = data$XZ, y = data$Y, family = "binomial")
  
  beta_cur <- glmnet::coef.glmnet(fit_yb_init, s = 2^lambda1) |> drop()
  
  loss <- -0.05 * .loss(X = data$XZ, Y = drop(data$Y), beta = beta_cur, family = "binomial")
  
  # prepare linear prediction and probability prediction for iteration update
  linear_pred <- .linearPred(data$XZ, beta_cur)
  prob_pred <- .probPred(linear_pred)
  
  tt <- prob_pred * {1.0 - prob_pred}
  
  X_aug <- rbind(cbind(1.0, data$XZ) * sqrt(tt), X.app)
  
  ## update Y_augmented
  y_update <- linear_pred * sqrt(tt) + {data$Y - prob_pred} / sqrt(tt)
  names(y_update) = rownames(data$XZ)
  Y_aug <- c(y_update, Y.app)
  
  beta_cur <- .ctnsSolution(data = data, X.app = X_aug, Y.app = Y_aug, lambda1 = lambda1)

  expect_warning(out <- .rwlsSolution(data, X.app = X.app, Y.app = Y.app,
                                      lambda1 = lambda1, iter.control = iter.control),
                 "maximum iterations reached")
  expect_equal(out, beta_cur)
})
