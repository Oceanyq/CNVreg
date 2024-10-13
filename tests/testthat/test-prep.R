test_that(".internalStep1 returns expected results", {
  
  cnv <- data.frame("ID" = c(3L, 2L, 1L,  4L),
                    "CHR" = rep(22L, 4L),
                    "BP1" = c(1L, 4L, 2L, 5L),
                    "BP2" = c(3L, 7L, 5L, 9L),
                    "TYPE" = c(0L, 0L, 0L, 0L))
  itv_data <- .breakCNV(cnv)
  
  wide_raw_sum <- .wideDataRaw(itv_data$CNV_frag_l)

  freq <- list("not.rare.idx" = c(2L, 4L, 5L),
               "freq" = c(1L, 2L, 1L, 2L, 2L, 1L))
  
  wide_common_data <- wide_raw_sum[, c(2L, 4L, 5L)]
  colnames(wide_common_data) <- c("del2", "del4", "del5")

  weight_deldup <- .weightMatrix(wide_raw_sum,
                                 freq$not.rare.idx,
                                 freq$freq)
  rownames(weight_deldup$CNVR.summary) <- c(2L, 4L, 5L)
  expected <- c(weight_deldup, "wide.data" = wide_common_data)
  
  expect_equal(.internalStep1(itv_data$CNV_frag_l, 5L, 0.3, "del"),
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
                    "BP1" = seq(1L, 20L, 2L),
                    "BP2" = seq(2L, 20L, 2L),
                    "TYPE" = c(0L, 0L, 1L, 1L,
                               0L, 3L, 3L, 3L,
                               4L, 4L))
  expect_error(prep(data.matrix(cnv)),
               "`CNV` data.frame in PLINK format with columns: `ID`, `CHR`, `BP1`, `BP2`,`TYPE`")
  
  expect_error(prep(cnv),
               "`Z` must be a data.frame with a column named `ID`")
  expect_error(prep(cnv, matrix(1, 20, 2, dimnames = list(NULL, c("ID", "Z")))),
               "`Z` must be a data.frame with a column named `ID`")
  expect_error(prep(cnv, data.frame("id" = 1:20, "Z" = 3)),
               "`Z` must be a data.frame with a column named `ID`")
  expect_error(prep(cnv, data.frame("ID" = 1:5, "Z" = 3)),
               "`Z` must be a data.frame with a column named `ID`")

  Z <- data.frame("ID" = 1:20, "Z1" = 1, "Z2" = 2)
  
  expect_error(prep(cnv, Z),
               "`Y` must be a data.frame with 2 columns")
  expect_error(prep(cnv, Z, matrix(1, 20, 2, dimnames = list(NULL, c("ID", "Y")))),
               "`Y` must be a data.frame with 2 columns")
  expect_error(prep(cnv, Z, data.frame("id" = 1:20, "Y" = 3)),
               "`Y` must be a data.frame with 2 columns")
  expect_error(prep(cnv, Z, data.frame("ID" = 1:19, "Y" = 3)),
               "`Y` must be a data.frame with 2 columns")
  expect_error(prep(cnv, Z, data.frame("ID" = 1:21, "Y" = 3)),
               "`Y` must be a data.frame with 2 columns")
  expect_error(prep(cnv, Z, data.frame("ID" = 1:21, "Y" = 3, "Y2" = 4)),
               "`Y` must be a data.frame with 2 columns")
  Y <- data.frame("ID" = 1:20, "Y1" = 1)
  
  expect_error(prep(cnv, Z, Y),
               "`rare.out` is a number in (0, 0.5)", fixed=TRUE)
  expect_error(prep(cnv, Z, Y, "1"),
               "`rare.out` is a number in (0, 0.5)", fixed=TRUE)
  expect_error(prep(cnv, Z, Y, 1:2),
               "`rare.out` is a number in (0, 0.5)", fixed=TRUE)
  expect_error(prep(cnv, Z, Y, 0.55 ),
               "`rare.out` is a number in (0, 0.5)", fixed=TRUE)
  expect_error(prep(cnv, Z, Y, 0.0 ),
               "`rare.out` is a number in (0, 0.5)", fixed=TRUE)
})

test_that("prep works as expected", {
  cnv <- data.frame("ID" = 11L:20L,
                    "CHR" = rep(22L, 10L),
                    "BP1" = seq(1L, 20L, 2L),
                    "BP2" = seq(2L, 20L, 2L),
                    "TYPE" = c(0L, 0L, 1L, 1L,
                               0L, 3L, 3L, 3L,
                               4L, 4L))
  Z <- data.frame("ID" = seq(20, 1, by = -1),
                  "X1" = withr::with_seed(1234, rnorm(20)),
                  "X2" = withr::with_seed(2345, rnorm(20)))
  
  Y <- data.frame("ID" = seq(20, 1, by = -1), 
                  "Y" = withr::with_seed(3456, rnorm(20)))
  
  # align Z and Y such that the top elements are in the same order as CNV
  IDs <- c(11:20, 1:10)
  ZZ <- Z[c(10:1, 20:11), c("X1", "X2")] |> data.matrix()
  YY <- Y[c(10:1, 20:11), "Y"] |> data.matrix()
  rownames(ZZ) <- IDs
  rownames(YY) <- IDs
  colnames(YY) <- "Y"
  
  out_wide_deldup <- Matrix::Matrix(seq_len(20L), 
                                    ncol  = 1L,
                                    sparse = TRUE,
                                    dimnames = list(IDs, "mid"))
  weight_any <- matrix(0.0, nrow = 0L, ncol = 0L)
  weight_options <- matrix(0.0, nrow = 6L, ncol = 0L,
                           dimnames = list(c("eql", "keql", 
                                             "wcs", "kwcs", 
                                             "wif", "kwif"), NULL))
  
  itv_data <- .breakCNV(cnv)
  
  tst <- itv_data$CNV_frag_l$deldup == "del"
  print(sum(tst))
  if (any(tst)) {
    res <- .internalStep1(itv_data$CNV_frag_l[tst, ], n.samples = 20L, 
                          rare.out = 0.03, type = "del")
    
    out_wide_deldup <- .mergeMatrices(out_wide_deldup, res$wide.data)
    weight_options <- res$weight.options
    weight_any <- res$weight.structure
    
    CNVR_ITV <- merge(res$CNVR.summary, itv_data$ITV_info, by = "idx", all.x = TRUE)
    CNVR_ITV$deldup <- "del"
    CNVRinfo <- CNVR_ITV
  }
  
  if (any(!tst)) {
    
    res <- .internalStep1(itv_data$CNV_frag_l[!tst, ], n.samples = 20L, 
                          rare.out = 0.03, type = "dup")
    
    out_wide_deldup <- .mergeMatrices(out_wide_deldup, res$wide.data)
    weight_options <- cbind(weight_options, res$weight.options)
    weight_any <- bdiag(weight_any, res$weight.structure)
    CNVR_ITV <- merge(res$CNVR.summary, itv_data$ITV_info, by = "idx", all.x = TRUE)
    CNVR_ITV$deldup <- "dup"
    CNVRinfo <- rbind(CNVRinfo, CNVR_ITV)
    
  }
  out_wide_deldup <- out_wide_deldup[ , colnames(out_wide_deldup) != "mid", drop = FALSE]
  
  res <- list("design" = out_wide_deldup,
              "Z" = as(ZZ, "sparseMatrix"),
              "Y" = YY,
              "weight.structure" = weight_any,
              "weight.options" = weight_options,
              "CNVR.info" = CNVRinfo)
  
  class(res) <- "WTsmth.data"
  
  expect_equal(prep(cnv, Z, Y, rare.out = 0.03), res)
  
  
  
})