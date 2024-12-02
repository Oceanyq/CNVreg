test_that("`.nfoldSplit()` returns expected errors", {
  
  expect_error(.nfoldSplit(),
               "`Y` must be a vector")
  expect_error(.nfoldSplit(rep("A", 10)),
               "`Y` must be a vector")
  Y <- 1:10
  
  expect_error(.nfoldSplit(Y),
               "`X.id` must be a vector")
  expect_error(.nfoldSplit(Y, X.id = list("A", "B", "C")),
               "`X.id` must be a vector")
  X_id <- 1:4
  
  expect_error(.nfoldSplit(Y, X_id),
               "`cv.control` must be a list; allowed elements are n.fold, n.core, and stratified")
  expect_error(.nfoldSplit(Y, X_id, cv.control = c("n.fold" = 1, "n.core" = 2, "stratified" = FALSE)),
               "`cv.control` must be a list; allowed elements are n.fold, n.core, and stratified")
  expect_error(.nfoldSplit(Y, X_id, cv.control = list("nfold" = 1, "n.core" = 2, "stratified" = FALSE)),
               "`cv.control` must be a list; allowed elements are n.fold, n.core, and stratified")
  
})

test_that("`.nfoldSplit()` returns expected results", {
  Y <- withr::with_seed(1235, rnorm(100))
  names(Y) <- 1:100
  expect <- rep(1L, 100L)
  expect_equal(withr::with_seed(2342, 
                                .nfoldSplit(Y, 1:10, list("n.fold" = 1, "n.core" = 2, "stratified" = FALSE))),
               expect)
})


test_that("`.nfoldSplit()` returns expected results 2", {
  Y <- withr::with_seed(1235, rbinom(100, 1, 0.5))
  names(Y) <- 1:100
  expect <- withr::with_seed(2342, sample(1:2, 100, TRUE))
  expect_equal(withr::with_seed(2342, 
                                .nfoldSplit(Y, 1:100, list("n.fold" = 2, "n.core" = 2, "stratified" = FALSE))),
               expect)
})

test_that("`.nfoldSplit()` returns expected results 3", {
  Y <- withr::with_seed(1235, rbinom(100, 1, 0.5))
  names(Y) <- 1:100
  expect <- withr::with_seed(2342, sample(1:2, 100, TRUE))
  expect_equal(withr::with_seed(2342, 
                                .nfoldSplit(Y, numeric(0), list("n.fold" = 2, "n.core" = 2, "stratified" = FALSE))),
               expect)
})


test_that("`.nfoldSplit()` returns expected results with stratification", {
  Y <- withr::with_seed(1235, rbinom(100, 1, 0.5))
  names(Y) <- 1:100
  X.id <- 1:50
  expect <- integer(100)
  withr::with_seed(2342, {
    expect[51:100][Y[51:100] == 0L] <- sample(1:2, sum(Y[51:100] == 0L), TRUE)
    expect[51:100][Y[51:100] == 1L] <- sample(1:2, sum(Y[51:100] == 1L), TRUE)
    expect[1:50][Y[1:50] == 0L] <- sample(1:2, sum(Y[1:50] == 0L), TRUE)
    expect[1:50][Y[1:50] == 1L] <- sample(1:2, sum(Y[1:50] == 1L), TRUE)
  })
  expect_equal(withr::with_seed(2342, 
                                .nfoldSplit(Y, 1:50, list("n.fold" = 2, "n.core" = 2, "stratified" = TRUE))),
               expect)
})

test_that("`.nfoldSplit()` returns expected results with stratification", {
  Y <- withr::with_seed(1235, rbinom(100, 1, 0.5))
  names(Y) <- 1:100
  expect <- integer(100)
  withr::with_seed(2342, {
    expect[Y == 0L] <- sample(1:2, sum(Y == 0L), TRUE)
    expect[Y == 1L] <- sample(1:2, sum(Y == 1L), TRUE)
  })
  expect_equal(withr::with_seed(2342, 
                                .nfoldSplit(Y, 1:100, list("n.fold" = 2, "n.core" = 2, "stratified" = TRUE))),
               expect)
})

test_that("`.nfoldSplit()` returns expected results with stratification", {
  Y <- withr::with_seed(1235, rbinom(100, 1, 0.5))
  names(Y) <- 1:100
  expect <- integer(100)
  withr::with_seed(2342, {
    expect[Y == 0L] <- sample(1:2, sum(Y == 0L), TRUE)
    expect[Y == 1L] <- sample(1:2, sum(Y == 1L), TRUE)
  })
  expect_warning(withr::with_seed(2342, 
                                .nfoldSplit(Y, 1:2, list("n.fold" = 2, "n.core" = 2, "stratified" = TRUE))),
                 "cannot stratify -- too few IDs in X")
  expect_warning(withr::with_seed(2342, 
                                  .nfoldSplit(Y, integer(0), list("n.fold" = 2, "n.core" = 2, "stratified" = TRUE))),
                 "cannot stratify -- too few IDs in X")
  expect_warning(withr::with_seed(2342, 
                                  .nfoldSplit(Y, 1:3, list("n.fold" = 2, "n.core" = 2, "stratified" = TRUE))),
                 "cannot stratify -- too few IDs in X")
  
})