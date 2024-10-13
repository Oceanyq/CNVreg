
test_that("`.logLH()` returns expected errors", {
  
  #X, beta are applied/tested in function Linear_Pred()
  
  X <- as(matrix(1, 10L, 5L), "sparseMatrix")
  p <- vector(mode = "numeric", length = 6)
  
  #Checky, missingy
  expect_error(.logLH(X, p),
               "`Y` must be a vector of number (0,1)", fixed=TRUE)
  
  #Checky, not numericy
  expect_error(.logLH(X, p,
                     c("1", 1, 0, 0, 0, 0, 1, 1, 1, 0)),
               "`Y` must be a vector of number (0,1)", fixed=TRUE)
  
  #Checky, not vector
  expect_error(.logLH(X, p,
                     data.frame("x1" = 1L:10L)),
               "`Y` must be a vector of number (0,1)", fixed=TRUE)
  expect_error(.logLH(X, p,
                      matrix(1, 10L, 1L)),
               "`Y` must be a vector of number (0,1)", fixed=TRUE)
  
  #checky dimension, length(Y)=9
  expect_error(.logLH(X, p,
                      c(1, 0, 0, 0, 0, 1, 1, 1, 0)),
               "`Y` must be a vector of number (0,1)", fixed=TRUE)
  
  #warnings
  expect_warning(.logLH(X, p,
                        c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
                 "All Y are all equal")
  
})


#test_that("`.logLH()` returns expected results; test 1", {
#  X <- matrix(1L, 10L, 5L)
#  beta <- c(-2.0, seq(0.1, 0.5, 0.1))
# y <- c(rep(0,5),rep(1,5))
#
#  expected <-  prob
#
#  expect_equal(.logLH(X, beta, y), expected)
#
#})


test_that("`.logLH()` returns expected results; test 2", {
  
  alternative_func <- function(X, beta, y, loglike = TRUE) {
    
    #Get parameters
    n <- length(y)
    
    #Compute the log-likelihood
    TERMS <- rep(-Inf, n)
    for (i in 1:n) {
      XB <- beta[1]+sum(X[i, ]*beta[-1])
      TERMS[i] <- - VGAM::log1pexp(-XB) - (1-y[i])*XB }
    
    LOGLIKE <- sum(TERMS)
    
    #Return the output
    if (loglike) { LOGLIKE } else { exp(LOGLIKE) }
  }
  
  
  
  X <- as(matrix(1L, 10L, 5L), "sparseMatrix")
  beta <- c(-2.0, seq(0.1, 0.5, 0.1))
  y <- c(rep(0,5),rep(1,5))
  
  
  expected <- alternative_func(X, beta, y)
  
  expect_equal(.logLH(X, beta, y), expected)
  
})



