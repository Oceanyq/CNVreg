cnv <- data.frame("ID" = 11L:20L,
                  "CHR" = rep(22L, 10L),
                  "BP1" = seq(100L, 2000L, 200L),
                  "BP2" = seq(300L, 2200L, 200L),
                  "TYPE" = c(0L, 0L, 1L, 1L,
                             0L, 3L, 3L, 3L,
                             4L, 4L))
Z <- data.frame("ID" = seq(20, 1, by = -1),
                "X1" = withr::with_seed(1234, rnorm(20)),
                "X2" = withr::with_seed(2345, rnorm(20)),
                "X3" = factor(withr::with_seed(3456, sample(1:3, 20, replace = TRUE))))

Y <- data.frame("ID" = seq(20, 1, by = -1), 
                "Y" = withr::with_seed(3456, rnorm(20)))

test_that("`.combineXZ()` returns expected errors", {
  data <- prep(cnv, Y, Z, rare.out = 0.03)

  XZ <- Matrix::Matrix(0.0, 
                       nrow = 20L, 
                       ncol = 4L + ncol(data$design), 
                       dimnames = list(rownames(data$Z),
                                       c(colnames(data$design), "X1", "X2", "X32", "X33")))
  
  idx <- match(rownames(data$design), rownames(data$Z))
  XZ[idx, seq_len(ncol(data$design))] <- data$design
  XZ[, seq_len(4L) + ncol(data$design)] <- data$Z
  expect_equal(.combineXZ(data), XZ)

})


test_that("`.combineXZ()` returns expected errors", {
  data <- prep(cnv, Y, Z, rare.out = 0.03)
  XZ <- Matrix::Matrix(0.0, 
                       nrow = nrow(Z), 
                       ncol = 4L + ncol(data$design), 
                       dimnames = list(rownames(data$Z),
                                       c(colnames(data$design), "X1", "X2", "X32", "X33")))
  
  idx <- match(rownames(data$design), rownames(data$Z))
  XZ[idx, seq_len(ncol(data$design))] <- data$design
  XZ[, seq_len(ncol(data$Z)) + ncol(data$design)] <- data$Z
  expect_equal(.combineXZ(data), XZ)
  
})


test_that("`.createA()` returns expected errors", {
  data <- prep(cnv, Y, Z, rare.out = 0.03)
  data$XZ <- .combineXZ(data)
  
  expect_error(.createA(unclass(data)),
               "`data` must be a WTsmth.data object")
  
  expect_error(.createA(data),
               "`weight` must be a character")
  
  expect_error(.createA(data, rownames(data$weight.options)),
               "`weight` must be a character")
})

test_that("`.createA()` returns expected results", {
  
  data <- prep(cnv, Y, Z, rare.out = 1e-8)
  data$XZ <- .combineXZ(data)
  
  A <- data$weight.options[1L, ] * data$weight.structure
  rownames(A) <- paste0("A", seq_len(nrow(A)))
  A_matrix = Matrix::Matrix(0.0, 
                            nrow = nrow(data$weight.structure), ncol = ncol(data$XZ), 
                            dimnames = list(rownames(A), colnames(data$XZ)),
                            sparse = TRUE)
  A_matrix[, 1L:ncol(A)] <- A
  expect_equal(.createA(data, rownames(data$weight.options)[1L]), A_matrix)
 
  
})

test_that("`.expandWTsmth()` returns expected errors", {
  
  expect_error(.expandWTsmth(), 
               "`data` must be a WTsmth.data object")
  
  data <- prep(cnv, Y, Z, rare.out = 1e-8)
  
  expect_error(.expandWTsmth(unclass(data)), 
               "`data` must be a WTsmth.data object")

  expect_error(.expandWTsmth(data),
               "`weight` must be one of eql, keql, wcs, kwcs, wif, kwif")
  
  expect_error(.expandWTsmth(data, 1),
               "`weight` must be one of eql, keql, wcs, kwcs, wif, kwif")

  expect_error(.expandWTsmth(data, c("eql", "keql")),
               "`weight` must be one of eql, keql, wcs, kwcs, wif, kwif")
  
  expect_error(.expandWTsmth(data, c("Keql")),
               "`weight` must be one of eql, keql, wcs, kwcs, wif, kwif")

})

test_that("`.expandWTsmth()` returns expected results", {
  
  data <- prep(cnv, Y, Z, rare.out = 1e-8)
  data$XZ <- .combineXZ(data)
  data$A <- .createA(data, "keql")

  expect_equal(.expandWTsmth(data, "keql"), data)  
})

test_that("`.confirmBinary()` returns expected errors/results", {
  
  expect_error(.confirmBinary(withr::with_seed(1234, rnorm(10))),
               "Y is not integer-like")
  expect_warning(.confirmBinary(c(1,2,3,1,2,3,1,2,3)),
                 "Y does not appear to be binary")
  
  expect_equal(.confirmBinary(c(1,2,1,2,1,2,1,2)),
               c(0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L))

  expect_equal(.confirmBinary(c("A", "B", "A", "B", "A", "B", "A", "B")),
               c(0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L))
  
})