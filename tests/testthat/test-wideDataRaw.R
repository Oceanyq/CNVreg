test_that("`.wideDataRaw()` returns expected errors", {
  
  expect_error(.wideDataRaw(),
               "`itv.long` data.frame is missing or not appropriately defined")
  
  expect_error(.wideDataRaw(list("ID" = c(1L : 5L),
                                 "CHR" = rep(22L, 5L),
                                 "idx" = 1:5,
                                 "TYPE" = c(0L, 1L, 0L, 0L, 0L),
                                 "deldup" = 1:5,
                                 "AUC" = 1:5)),
               "`itv.long` data.frame is missing or not appropriately defined")
  
  expect_error(.wideDataRaw(data.frame("ID" = c(1L : 5L),
                                       "CHR" = rep(22L, 5L),
                                       "TYPE" = c(0L, 1L, 0L, 0L, 0L),
                                       "deldup" = 1:5,
                                       "AUC" = 1:5)),
               "`itv.long` data.frame is missing or not appropriately defined")
  
  expect_warning(.wideDataRaw(data.frame("ID" = c(1L, 1L:4L),
                                         "CHR" = rep(22L, 5L),
                                         "idx" = c(1, 1:4),
                                         "TYPE" = c(0L, 1L, 0L, 0L, 0L),
                                         "deldup" = 1:5,
                                         "AUC" = 1:5)),
                 "One or more samples have more than 1 record for a CNV del/dup event at the same location, kept the highest dosage")
  
})



test_that("`.wideDataRaw()` returns expected results; test 1", {
  
  cnv <- data.frame("ID" = c(3L, 2L, 1L,  4L),
                    "CHR" = rep(22L, 4L),
                    "BP1" = c(1L, 4L, 2L, 5L),
                    "BP2" = c(3L, 7L, 5L, 9L),
                    "TYPE" = c(0L, 0L, 0L, 0L))
  itv_data <- .breakCNV(cnv)

  expected <- as(matrix(c(0, 2, 2, 2, 0, 0,
                          0, 0, 0, 2, 4, 0,
                          2, 2, 0, 0, 0, 0,
                          0, 0, 0, 0, 4, 4), 
                        ncol = 6L, nrow = 4L, byrow = TRUE, 
                        dimnames = list(1L:4L, 1L:6L) ), "sparseMatrix")
  
  expect_equal(.wideDataRaw(itv_data$CNV_frag_l), expected)
})

