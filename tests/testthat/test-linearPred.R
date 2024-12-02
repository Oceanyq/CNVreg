test_that("`.linearPred()` returns expected errors", {
  
  expect_error(.linearPred(),
               "`X` must be a Matrix with at least 1 column")
  expect_error(.linearPred(vector(mode = "list", length = 5L),
                           vector(mode = "numeric", length = 6L)),
               "`X` must be a Matrix with at least 1 column")
  expect_error(.linearPred(data.frame("x1" = 1L:10L,
                                      "x2" = rep(1L,10L),
                                      "x3" = rep(2L,10L),
                                      "x4" = rep(1L,10L)),
                           vector(mode = "numeric", length = 6)),
               "`X` must be a Matrix with at least 1 column")
  expect_error(.linearPred(matrix("1", 10L, 5L),
                           vector(mode = "numeric", length = 6)),
               "`X` must be a Matrix with at least 1 column")
  
  X <- as(matrix(1, 10, 5L), "sparseMatrix")
  
  expect_error(.linearPred(X),
               "`beta` must be a vector with numeric values, beta has intercept")
  expect_error(.linearPred(X, matrix(1L, 5L, 1L)),
               "`beta` must be a vector with numeric values, beta has intercept")
  expect_error(.linearPred(X, data.frame("x1" = 1L:5L)),
               "`beta` must be a vector with numeric values, beta has intercept")
  expect_error(.linearPred(X, vector(mode = "character", length = 6)),
               "`beta` must be a vector with numeric values, beta has intercept")
  expect_error(.linearPred(X, vector(mode = "numeric", length = 5)),
               "`beta` must be a vector with numeric values, beta has intercept")
  
  
})


test_that("`.linearPred()` returns expected results; test 1", {
  
  X <- as(matrix(1L,10L,5L), "sparseMatrix")
  beta <- seq(1L, 6L, 1L)
  
  expected <- rep(21L, 10L)
  
  expect_equal(.linearPred(X, beta), expected)
  
})


test_that("`.linearPred()` returns expected results; test 2", {
  alternative_func <- function(X, beta) {
    
    drop(cbind(1,X)%*%beta)
    
  }
  
  X <- matrix(1L:50L, 10L, 5L)
  beta <- rep(1, 6)
  
  expected <- alternative_func(X, beta)
  
  expect_equal(.linearPred(as(matrix(1L:50L, 10L, 5L), "sparseMatrix"),beta), expected)
  
})