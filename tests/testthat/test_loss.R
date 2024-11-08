withr::with_seed(1234L, {
  X <- Matrix::Matrix(rnorm(100L), 100L, 1L, dimnames = list(NULL, "X1"))
  Y <- stats::rbinom(100L, 1L, 0.5)
  beta <- stats::runif(2L)
})

test_that("`.loss` returns expected errors", {
  expect_error(.loss(),
               "`X` must be provided")
  
  X <- matrix(1.0, 10L, 10L)
  expect_error(.loss(X),
               "`Y` must be a numeric vector")
  
  expect_error(.loss(X, matrix(0.0, 1L, 10L)),
               "`Y` must be a numeric vector")
  
  expect_error(.loss(X, numeric(9L)),
               "`Y` must be a numeric vector")
  
  expect_error(.loss(X, numeric(11L)),
               "`Y` must be a numeric vector")
  
  Y <- numeric(10L)
  
  expect_error(.loss(X, Y),
               "`beta` vector must be provided")
  
  beta <- numeric(11L)
  
  expect_error(.loss(X, Y, beta, binomial()),
               "'arg' must be NULL or a character vector")
  
  expect_error(.loss(X, Y, beta, c("binomial", "gaussian")),
               "'arg' must be of length 1")
  
  expect_error(.loss(X, Y, beta, c("Gaussian")),
               "'arg' should be one of “gaussian”, “binomial”")
})

test_that("`.loss()` returns expected results with family = gaussian", {
  # we've already tested Linear_Pred
  Y_hat <- .linearPred(X, beta)
  
  expected <- 0.0
  for (i in 1L:100L) {
    expected <- expected + (Y_hat[i] - Y[i]) * (Y_hat[i] - Y[i]) / 100.0
  }

  expect_equal(.loss(X, Y, beta), expected)
  expect_equal(.loss(X, Y, beta, "gaussian"), expected)
})


test_that("`.loss()` returns expected results with family = binomial", {
  # we've already tested LogLH  
  expected <- .logLH(X = X, Y = Y, beta = beta)
  
  expect_equal(.loss(X, Y, beta, "binomial"), expected)
})
