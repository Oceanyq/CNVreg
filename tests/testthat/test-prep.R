test_that("`.orderData()` returns expected errors", {
  
  expect_error(.orderData(), 
               "`CNV` must be provided")
  cnv <- data.frame("ID" = 1L:10L,
                    "CHR" = rep(22L, 10L),
                    "BP1" = seq(100L, 2000L, 200L),
                    "BP2" = seq(200L, 2000L, 200L),
                    "TYPE" = c(0L, 0L, 1L, 1L,
                               0L, 3L, 3L, 3L,
                               4L, 4L))
  expect_error(.orderData(data.matrix(cnv)),
               "`CNV` data.frame in PLINK format with columns: `ID`, `CHR`, `BP1`, `BP2`,`TYPE`")
  
  expect_error(.orderData(cnv),
               "`Z` must be a data.frame with a column named `ID`")
  expect_error(.orderData(cnv, matrix(1, 20, 2, dimnames = list(NULL, c("ID", "Z")))),
               "`Z` must be a data.frame with a column named `ID`")
  expect_error(.orderData(cnv, data.frame("id" = 1:20, "Z" = 3)),
               "`Z` must be a data.frame with a column named `ID`")
  expect_error(.orderData(cnv, data.frame("ID" = 1:5, "Z" = 3)),
               "`Z` must be a data.frame with a column named `ID`")
  
  Z <- data.frame("ID" = 1:20, "Z1" = 1, "Z2" = 2)
  
  expect_error(.orderData(cnv, Z),
               "`Y` must be a data.frame with 2 columns")
  expect_error(.orderData(cnv, Z, matrix(1, 20, 2, dimnames = list(NULL, c("ID", "Y")))),
               "`Y` must be a data.frame with 2 columns")
  expect_error(.orderData(cnv, Z, data.frame("id" = 1:20, "Y" = 3)),
               "`Y` must be a data.frame with 2 columns")
  expect_error(.orderData(cnv, Z, data.frame("ID" = 1:19, "Y" = 3)),
               "`Y` must be a data.frame with 2 columns")
  expect_error(.orderData(cnv, Z, data.frame("ID" = 1:21, "Y" = 3)),
               "`Y` must be a data.frame with 2 columns")
  expect_error(.orderData(cnv, Z, data.frame("ID" = 1:21, "Y" = 3, "Y2" = 4)),
               "`Y` must be a data.frame with 2 columns")
  Y <- data.frame("ID" = 1:20, "Y1" = 1)
  
})

test_that("`.orderData()` returns expected results", {
  
  cnv <- data.frame("ID" = as.character(11L:20L),
                    "CHR" = rep(22L, 10L),
                    "BP1" = seq(100L, 2000L, 200L),
                    "BP2" = seq(200L, 2000L, 200L),
                    "TYPE" = c(0L, 0L, 1L, 1L,
                               0L, 3L, 3L, 3L,
                               4L, 4L))
  Z <- data.frame("ID" = seq(20, 1, by = -1),
                  "X1" = withr::with_seed(1234, rnorm(20)),
                  "X2" = withr::with_seed(2345, rnorm(20)))
  
  Y <- data.frame("ID" = seq(20, 1, by = -1), 
                  "Y" = withr::with_seed(3456, rnorm(20)))
  
  # align Z and Y such that the top elements are in the same order as CNV
  IDs <- c(11:20, c(1, 10, 2:9))
  ZZ <- Z[c(10:1, 20, 11, 19:12), c("X1", "X2")] |> data.matrix()
  YY <- Y[c(10:1, 20, 11, 19:12), "Y"] |> data.matrix()
  rownames(ZZ) <- IDs
  rownames(YY) <- IDs
  colnames(YY) <- "Y"
  
  res <- list("CNV" = cnv,
              "Z" = ZZ,
              "Y" = YY)
  
  expect_equal(.orderData(cnv, Z, Y), res)
  
})

test_that(".internalStep1 returns expected results", {
  
  cnv <- data.frame("ID" = c(3L, 2L, 1L,  4L),
                    "CHR" = rep(22L, 4L),
                    "BP1" = c(100L, 400L, 200L, 500L),
                    "BP2" = c(300L, 700L, 500L, 900L),
                    "TYPE" = c(0L, 0L, 0L, 0L))
  itv_data <- .breakCNV(cnv)
  
  wide_raw_sum <- .wideDataRaw(itv_data$long.cnv)

  freq <- list("not.rare.idx" = c(2L, 4L, 5L),
               "freq" = c(1L, 2L, 1L, 2L, 2L, 1L))
  
  wide_common_data <- wide_raw_sum[, c(2L, 4L, 5L)]
  colnames(wide_common_data) <- c("del2", "del4", "del5")

  weight_deldup <- .weightMatrix(wide_raw_sum,
                                 freq$not.rare.idx,
                                 freq$freq)
  rownames(weight_deldup$CNVR.summary) <- c(2L, 4L, 5L)
  expected <- c(weight_deldup, "wide.data" = wide_common_data)
  
  expect_equal(.internalStep1(itv_data$long.cnv, 5L, 0.3, "del"),
               expected)
  
})

test_that(".mergeMatrices returns expected results", {
  old <- withr::with_seed(
    1234, 
    matrix(rnorm(100), 25, 4,
           dimnames = list(paste0("A",1:25), paste0("B", 1:4))))
  new <- withr::with_seed(
    2345, 
    matrix(rnorm(100), 20, 5,
           dimnames = list(paste0("A",1:20), paste0("b", 1:5))))
  
  merged <- matrix(0.0, 25, 9, 
                   dimnames = list(paste0("A",1:25), 
                                   c(paste0("B", 1:4), paste0("b",1:5))))
  
  merged[1:25, 1:4] <- old
  merged[1:20, 5:9] <- new
  
  expected <- as(merged, "sparseMatrix")
  expect_equal(.mergeMatrices(old, new), expected)

})

test_that("prep returns expected errors", {
  
  expect_error(prep(), 
               "`CNV` must be provided")
  cnv <- data.frame("ID" = 1L:10L,
                    "CHR" = rep(22L, 10L),
                    "BP1" = seq(100L, 2000L, 200L),
                    "BP2" = seq(200L, 2000L, 200L),
                    "TYPE" = c(0L, 0L, 1L, 1L,
                               0L, 3L, 3L, 3L,
                               4L, 4L))
  
  expect_error(prep(cnv),
               "`Y` must be provided")
  Y <- data.frame("ID" = 1:20, "Y1" = 1)
  
  expect_error(prep(cnv, Y = Y, Z = NA),
               "`Z` must be NULL or a data.frame")
  Z <- data.frame("ID" = 1:20, "Z1" = 1, "Z2" = 2)
  
  expect_error(prep(cnv, Y, Z, "1"),
               "`rare.out` is a number in (0, 0.5)", fixed=TRUE)
  expect_error(prep(cnv, Y, Z, 1:2),
               "`rare.out` is a number in (0, 0.5)", fixed=TRUE)
  expect_error(prep(cnv, Y, Z, 0.55 ),
               "`rare.out` is a number in (0, 0.5)", fixed=TRUE)
  expect_error(prep(cnv, Y, Z, 0.0 ),
               "`rare.out` is a number in (0, 0.5)", fixed=TRUE)
})

test_that("prep works as expected", {
  cnv <- data.frame("ID" = 11L:20L,
                    "CHR" = rep(22L, 10L),
                    "BP1" = seq(100L, 2000L, 200L),
                    "BP2" = seq(300L, 2200L, 200L),
                    "TYPE" = c(0L, 0L, 1L, 1L,
                               0L, 3L, 3L, 3L,
                               4L, 4L))
  Z <- data.frame("ID" = seq(20, 1, by = -1),
                  "X1" = withr::with_seed(1234, rnorm(20)),
                  "X2" = withr::with_seed(2345, rnorm(20)))
  
  Y <- data.frame("ID" = seq(20, 1, by = -1), 
                  "Y" = withr::with_seed(3456, rnorm(20)))
  
  ordered_data <- .orderData(cnv, Z, Y)
  
  out_wide_deldup <- Matrix::Matrix(seq_len(20L), 
                                    ncol  = 1L,
                                    sparse = TRUE,
                                    dimnames = list(rownames(ordered_data$Z), "mid"))
  weight_any <- matrix(0.0, nrow = 0L, ncol = 0L)
  weight_options <- matrix(0.0, nrow = 6L, ncol = 0L,
                           dimnames = list(c("eql", "keql", 
                                             "wcs", "kwcs", 
                                             "wif", "kwif"), NULL))
  
  itv_data <- .breakCNV(ordered_data$CNV)
  
  tst <- itv_data$long.cnv$deldup == "del"
  
  CNVRinfo <- NULL

  if (any(tst)) {
    res <- .internalStep1(itv_data$long.cnv[tst, ], n.samples = 20L, 
                          rare.out = 0.03, type = "del")
    
    out_wide_deldup <- .mergeMatrices(out_wide_deldup, res$wide.data)
    weight_options <- res$weight.options
    weight_any <- res$weight.structure
    
    CNVR_ITV <- merge(res$CNVR.summary, itv_data$grid.info, by = "grid.id", all.x = TRUE)
    CNVR_ITV$deldup <- "del"
    CNVRinfo <- CNVR_ITV
  }
  
  if (any(!tst)) {
    
    res <- .internalStep1(itv_data$long.cnv[!tst, ], n.samples = 20L, 
                          rare.out = 0.03, type = "dup")
    
    out_wide_deldup <- .mergeMatrices(out_wide_deldup, res$wide.data)
    weight_options <- cbind(weight_options, res$weight.options)
    weight_any <- bdiag(weight_any, res$weight.structure)
    CNVR_ITV <- merge(res$CNVR.summary, itv_data$grid.info, by = "grid.id", all.x = TRUE)
    CNVR_ITV$deldup <- "dup"
    CNVRinfo <- rbind(CNVRinfo, CNVR_ITV)
    
  }
  out_wide_deldup <- out_wide_deldup[ , colnames(out_wide_deldup) != "mid", drop = FALSE]
  
  res <- list("design" = out_wide_deldup,
              "Z" = as(ordered_data$Z, "sparseMatrix"),
              "Y" = ordered_data$Y[, "Y"] |> drop(),
              "weight.structure" = weight_any,
              "weight.options" = weight_options,
              "CNVR.info" = CNVRinfo)
  
  class(res) <- "WTsmth.data"
  expect_equal(prep(cnv, Y, Z, rare.out = 0.03), res)

})

test_that("prep works as expected without Z", {
  cnv <- data.frame("ID" = 11L:20L,
                    "CHR" = rep(22L, 10L),
                    "BP1" = seq(100L, 2000L, 200L),
                    "BP2" = seq(300L, 2200L, 200L),
                    "TYPE" = c(0L, 0L, 1L, 1L,
                               0L, 3L, 3L, 3L,
                               4L, 4L))

  Y <- data.frame("ID" = seq(20, 1, by = -1), 
                  "Y" = withr::with_seed(3456, rnorm(20)))
  
  ordered_data <- .orderData(cnv, Z = NULL, Y)
  
  out_wide_deldup <- Matrix::Matrix(seq_len(20L), 
                                    ncol  = 1L,
                                    sparse = TRUE,
                                    dimnames = list(rownames(ordered_data$Y), "mid"))
  weight_any <- matrix(0.0, nrow = 0L, ncol = 0L)
  weight_options <- matrix(0.0, nrow = 6L, ncol = 0L,
                           dimnames = list(c("eql", "keql", 
                                             "wcs", "kwcs", 
                                             "wif", "kwif"), NULL))
  
  itv_data <- .breakCNV(ordered_data$CNV)
  
  tst <- itv_data$long.cnv$deldup == "del"
  
  CNVRinfo <- NULL
  
  if (any(tst)) {
    res <- .internalStep1(itv_data$long.cnv[tst, ], n.samples = 20L, 
                          rare.out = 0.03, type = "del")
    
    out_wide_deldup <- .mergeMatrices(out_wide_deldup, res$wide.data)
    weight_options <- res$weight.options
    weight_any <- res$weight.structure
    
    CNVR_ITV <- merge(res$CNVR.summary, itv_data$grid.info, by = "grid.id", all.x = TRUE)
    CNVR_ITV$deldup <- "del"
    CNVRinfo <- CNVR_ITV
  }
  
  if (any(!tst)) {
    
    res <- .internalStep1(itv_data$long.cnv[!tst, ], n.samples = 20L, 
                          rare.out = 0.03, type = "dup")
    
    out_wide_deldup <- .mergeMatrices(out_wide_deldup, res$wide.data)
    weight_options <- cbind(weight_options, res$weight.options)
    weight_any <- bdiag(weight_any, res$weight.structure)
    CNVR_ITV <- merge(res$CNVR.summary, itv_data$grid.info, by = "grid.id", all.x = TRUE)
    CNVR_ITV$deldup <- "dup"
    CNVRinfo <- rbind(CNVRinfo, CNVR_ITV)
    
  }
  out_wide_deldup <- out_wide_deldup[ , colnames(out_wide_deldup) != "mid", drop = FALSE]
  
  res <- list("design" = out_wide_deldup,
              "Z" = as(ordered_data$Z, "sparseMatrix"),
              "Y" = ordered_data$Y[, "Y"] |> drop(),
              "weight.structure" = weight_any,
              "weight.options" = weight_options,
              "CNVR.info" = CNVRinfo)
  
  class(res) <- "WTsmth.data"
  expect_equal(prep(cnv, Y, Z = NULL, rare.out = 0.03), res)
  
})
