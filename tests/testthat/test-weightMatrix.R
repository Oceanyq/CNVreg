library(Matrix)
test_that("`.weightMatrix()` returns expected errors", {
  
  #trigger missing() by not passing an empty input
  expect_error(.weightMatrix(),
               "`wide.data` must be a Matrix")
  
  expect_error(.weightMatrix(as(matrix(c(0, 2, 2, 2, 0, 0,
                                       0, 0, 0, 2, 4, 0,
                                       2, 2, 0, 0, 0, 0,
                                       0, 0, 0, 0, 4, 4), ncol = 6, nrow=4, byrow = TRUE ), "sparseMatrix")),
               "`not.rare.idx` must be a vector")
  expect_error(.weightMatrix(as(matrix(c(0, 2, 2, 2, 0, 0,
                                       0, 0, 0, 2, 4, 0,
                                       2, 2, 0, 0, 0, 0,
                                       0, 0, 0, 0, 4, 4), ncol = 6, nrow=4, byrow = TRUE ), "sparseMatrix"), 1L : 3L),
               "`freq` must be a vector of length ncol(wide.data)", fixed = TRUE)
  
  expect_error(.weightMatrix(as(matrix(c(0, 2, 2, 2, 0, 0,
                                       0, 0, 0, 2, 4, 0,
                                       2, 2, 0, 0, 0, 0,
                                       0, 0, 0, 0, 4, 4), ncol = 6, nrow=4, byrow = TRUE ), "sparseMatrix"), 1L : 3L, 1L : 3L),
               "`freq` must be a vector of length ncol(wide.data)", fixed = TRUE)
  
}
)


test_that("`.weightMatrix()` returns expected results; test 1", {
  
  cnv_wide <- as(matrix(c(0, 2, 2, 2, 0, 0,
                          0, 0, 0, 2, 4, 0,
                          2, 2, 0, 0, 0, 0,
                          0, 0, 0, 0, 4, 4), ncol = 6, nrow=4, byrow = TRUE ), "sparseMatrix")
  
  wide_idx <- c(2L, 4L, 5L)
  wide_freq <- c(1L, 2L, 1L, 2L, 2L, 1L)
  
  weight_options <- matrix(1.0, nrow = 6L, ncol = 2L, dimnames = list(c("eql", "keql", "wcs", "kwcs", "wif", "kwif")))
  weight_options[, 1L] <- 0.0
  weight_structure <- as(matrix(c(0,0,0,
                                  0,-1,1), ncol=3, nrow=2, byrow=TRUE), "sparseMatrix")
  
  expected <- list("weight.structure" = weight_structure,
                   "weight.options" = weight_options,
    "CNVR.summary"= matrix(c(1,2,2,
                          2,4,2,
                          2,5,2), ncol=3, nrow=3, byrow=TRUE, dimnames = list(NULL, c("CNV.id", "grid.id", "freq")))
    
  )
  
  expect_equal(.weightMatrix(wide.data = cnv_wide, not.rare.idx = wide_idx,  freq = wide_freq ), expected)

})

