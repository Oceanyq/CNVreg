test_that("`.nfoldSplit()` returns expected errors", {
  
  expect_error(.nfoldSplit(),
               "`Y` must be a vector")
  expect_error(.nfoldSplit(rep("A", 10)),
               "`Y` must be a vector")
  Y <- 1:10
  
  expect_error(.nfoldSplit(Y),
               "`cv.control` must be a list; allowed elements are n.fold, n.core, and stratified")
  expect_error(.nfoldSplit(Y, cv.control = c("n.fold" = 1, "n.core" = 2, "stratified" = FALSE)),
               "`cv.control` must be a list; allowed elements are n.fold, n.core, and stratified")
  expect_error(.nfoldSplit(Y, cv.control = list("nfold" = 1, "n.core" = 2, "stratified" = FALSE)),
               "`cv.control` must be a list; allowed elements are n.fold, n.core, and stratified")
  
})

test_that("`.nfoldSplit()` returns expected results", {
  Y <- withr::with_seed(1235, rnorm(100))
  expect <- rep(1L, 100L)
  expect_equal(withr::with_seed(2342, 
                                .nfoldSplit(Y, list("n.fold" = 1, "n.core" = 2, "stratified" = FALSE))),
               expect)
})

.nfoldSplit <- function(Y, cv.control) {
  
  stopifnot(
    "`Y` must be a vector" = !missing(Y) && .isNumericVector(Y),
    "`cv.control` must be a list; allowed elements are n.fold, n.core, and stratified" = 
      .isNamedList(cv.control, c("n.fold", "n.core", "stratified")),
  )
  
  if (cv.control$stratified) {
    Y <- .confirmBinary(Y)
    tr <- rep(0L, length(Y))
    tr[Y == 0L] <- sample(seq_len(cv.control$n.fold), sum(Y == 0L), replace = TRUE)
    tr[Y != 0L] <- sample(seq_len(cv.control$n.fold), sum(Y != 0L), replace = TRUE)
  } else{
    tr <- sample(seq_len(cv.control$n.fold), length(Y), replace = TRUE)
  }
  
  sort(tr)
}
