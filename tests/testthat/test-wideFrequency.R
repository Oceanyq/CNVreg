test_that("`.wideFrequency()` returns expected errors", {
  
  #trigger missing() by not passing an empty input
  expect_error(.wideFrequency(sample.size = 10L, rare.out = 0.1),
               "`wide.raw` must be a Matrix object")
  
  tmp_wide_raw <- matrix(c(0, 2, 2, 2, 0, 0,
                           0, 0, 0, 2, 4, 0,
                           2, 2, 0, 0, 0, 0,
                           0, 0, 0, 0, 4, 4), 
                         ncol = 6, nrow = 4, byrow = TRUE)
  
  expect_error(.wideFrequency(tmp_wide_raw, sample.size = 10L, rare.out = 0.1),
               "`wide.raw` must be a Matrix object")
  
  tmp_wide_raw <- as(tmp_wide_raw, "sparseMatrix")
  
  #Matrix::Matrix(the_matrix, sparse = TRUE)
  expect_error(.wideFrequency(tmp_wide_raw),
               "`sample.size` be a scalar")
  expect_error(.wideFrequency(tmp_wide_raw, "1"),
               "`sample.size` be a scalar")
  expect_error(.wideFrequency(tmp_wide_raw, 1:2),
               "`sample.size` be a scalar")
  
  expect_error(.wideFrequency(tmp_wide_raw, 6L),
               "`rare.out` is a number in (0, 0.5)", fixed=TRUE)
  expect_error(.wideFrequency(tmp_wide_raw, 6L, "1"),
               "`rare.out` is a number in (0, 0.5)", fixed=TRUE)
  expect_error(.wideFrequency(tmp_wide_raw, 6L, 1:2),
               "`rare.out` is a number in (0, 0.5)", fixed=TRUE)
  expect_error(.wideFrequency(tmp_wide_raw, 6L, 0.55 ),
               "`rare.out` is a number in (0, 0.5)", fixed=TRUE)
  expect_error(.wideFrequency(tmp_wide_raw, 6L, 0.0 ),
               "`rare.out` is a number in (0, 0.5)", fixed=TRUE)
  
})


test_that("`.wideFrequency()` returns expected results; test 1", {
  
  cnv_wide <- as(matrix(c(0, 2, 2, 2, 0, 0,
                          0, 0, 0, 2, 4, 0,
                          2, 2, 0, 0, 0, 0,
                          0, 0, 0, 0, 4, 4), 
                        ncol = 6L, nrow = 4L, byrow = TRUE), "sparseMatrix")
  
  expected <- list("not.rare.idx" = c(2,4,5),
                   "freq" = c(1, 2, 1, 2, 2, 1))
  
  expect_equal(.wideFrequency(cnv_wide, 5, 0.3), expected)

  expected <- list("not.rare.idx" = 1:6,
                   "freq" = c(1, 2, 1, 2, 2, 1))
  
  expect_equal(.wideFrequency(cnv_wide, 5, 0.0001), expected)
})

