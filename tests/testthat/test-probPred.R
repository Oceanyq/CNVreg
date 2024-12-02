
test_that("`.probPred()` returns expected errors", {
  
  #missing X
  expect_error(.probPred(),
               "`X` must be a numeric vector")
  
  expect_error(.probPred(matrix(0, 10, 1)),
               "`X` must be a numeric vector")
  
  #warnings
  expect_warning(.probPred(rep(Inf, 10L)),
                 "probabilities are all 0 or all 1")
  
  expect_warning(.probPred(rep(-Inf, 10L)),
                 "probabilities are all 0 or all 1")
})


test_that("`.probPred()` returns expected results; test 1", {
  prob= seq(0.05, 0.85, 0.2)
  
  X <- log(prob/(1-prob))
  
  
  expected <-  prob
  
  expect_equal(.probPred(X), expected)
  
})


test_that("`.probPred()` returns expected results; test 2", {
  
  alternative_func <- function(logit){
    odds <- exp(logit)
    prob <- odds / (1 + odds)
    return(prob)
  }
  
  prob= seq(0.05, 0.85, 0.2)
  
  X <- log(prob/(1-prob))
  
  
  expected <- alternative_func(X)
  
  expect_equal(.probPred(X), expected)
  
})



