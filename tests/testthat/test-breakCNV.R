test_that("`.breakCNV()` returns expected errors", {
  
  # trigger missing(CNV) by not passing an input
  expect_error(.breakCNV(),
               "`CNV` data.frame in PLINK format with columns: `ID`, `CHR`, `BP1`, `BP2`,`TYPE`")
  
})

test_that(".createITVData returns expected results", {
  
  cnv <- data.frame("ID" = 1L:10L,
                    "CHR" = rep(22L, 10L),
                    "BP1" = seq(1L, 20L, 2L),
                    "BP2" = seq(2L, 20L, 2L),
                    "TYPE" = c(0L, 0L, 1L, 1L,
                               0L, 3L, 3L, 3L,
                               4L, 4L))
  expected <- data.frame("CHR" = rep(22L, 19L),
                         "grid.id" = 1L:19L,
                         "CNV.start" = 1L:19L,
                         "CNV.end" = 2L:20L)
  
  expect_equal(.createGrid(cnv), expected)
})

test_that(".createLongData returns expected results", {
  
  cnv <- data.frame("ID" = 1L:10L,
                    "CHR" = rep(22L, 10L),
                    "BP1" = seq(1L, 20L, 2L),
                    "BP2" = seq(2L, 20L, 2L),
                    "TYPE" = c(0L, 0L, 1L, 1L,
                               0L, 3L, 3L, 3L,
                               4L, 4L))

  itv_data <- .createGrid(cnv)
  
  expected <- data.frame("ID" = 1L:10L,
                         "CHR" = rep(22L, 10L),
                         "grid.id" = seq(1L, 20L, 2L),
                         "TYPE" = c(0L, 0L, 1L, 1L,
                                    0L, 3L, 3L, 3L,
                                    4L, 4L),
                         "deldup" = c(rep("del", 5), rep("dup", 5)),
                         "AUC" = c(2.0, 2.0, 1.0, 1.0,
                                   2.0, 1.0, 1.0, 1.0,
                                   2.0, 2.0))

  expect_equal(.createLongData(cnv, itv_data), expected)

})

test_that("`Break_CNV_interval()` returns expected results; test 1", {
  
  cnv <- data.frame("ID" = 1L:10L,
                    "CHR" = rep(22L, 10L),
                    "BP1" = seq(100L, 2000L, 200L),
                    "BP2" = seq(200L, 2000L, 200L),
                    "TYPE" = c(0L, 0L, 1L, 1L,
                               0L, 3L, 3L, 3L,
                               4L, 4L))
  expected <- list("long.cnv" = data.frame("ID" = 1L:10L,
                                             "CHR" = rep(22L, 10L),
                                             "grid.id" = seq(1L, 20L, 2L),
                                             "TYPE" = c(0L, 0L, 1L, 1L,
                                                        0L, 3L, 3L, 3L,
                                                        4L, 4L),
                                             "deldup" = c(rep("del", 5), rep("dup", 5)),
                                             "AUC" = c(200.0, 200.0, 100.0, 100.0,
                                                       200.0, 100.0, 100.0, 100.0,
                                                       200.0, 200.0)),
                   "grid.info" = data.frame("CHR" = rep(22L, 19L),
                                           "grid.id" = 1L:19L,
                                           "CNV.start" = seq(100L, 1900L, 100L),
                                           "CNV.end" = seq(200L, 2000L, 100L)))

  expect_equal(.breakCNV(cnv), expected)
})


test_that("`Break_CNV_interval()` returns expected results; test 2", {
  
  cnv <- data.frame("ID" = 1L:4L,
                    "CHR" = rep(22L, 4L),
                    "BP1" = c(100L, 1000L, 1100L, 2000L),
                    "BP2" = c(500L, 1200L, 1500L, 3500L),
                    "TYPE" = c(0L, 0L, 1L, 1L))
  
  expected <- list("long.cnv" = data.frame("ID" = c(1L, 2L, 2L, 3L, 3L, 4L),
                                             "CHR" = rep(22L, 6L),
                                             "grid.id" = c(1L, 3L, 4L, 4L, 5L, 7L),
                                             "TYPE" = c(0L, 0L, 0L, 1L, 1L, 1L),
                                             "deldup" = c(rep("del", 6)),
                                             "AUC" = c(2*400, 2*100, 2*100, 1*100, 1*300, 1*1500)),
                   "grid.info" = data.frame("CHR" = rep(22L, 7L),
                                           "grid.id" = 1L:7L,
                                           "CNV.start" = c(100L, 500L, 1000L, 1100L, 1200L, 1500L, 2000L),
                                           "CNV.end" = c(500L, 1000L, 1100L, 1200L, 1500L, 2000L, 3500L)))
  expect_equal(.breakCNV(cnv), expected)
})


test_that("`Break_CNV_interval()` returns expected results; test 3", {
  
  alternative_func <- function(cnv) {
    BPs <- sort(unique(c(cnv$BP1, cnv$BP2)))
    long_form <- NULL
    for (i in 1L:nrow(cnv)) {
      for (j in 1L:{length(BPs) - 1}) {
        if (BPs[j] >= cnv$BP1[i] && BPs[j+1] <= cnv$BP2[i]) {
          if (is.null(long_form)) {
            long_form <- data.frame("ID" = cnv$ID[i],
                                    "CHR" = cnv$CHR[i],
                                    "grid.id" = j,
                                    "TYPE" = cnv$TYPE[i],
                                    "deldup" = if(cnv$TYPE[i] < 2) {"del"} else {"dup"},
                                    "AUC" = abs(2.0 - cnv$TYPE[i]) * (BPs[j+1] - BPs[j]))
          } else {
            long_form <- rbind(long_form,
                               data.frame("ID" = cnv$ID[i],
                                          "CHR" = cnv$CHR[i],
                                          "grid.id" = j,
                                          "TYPE" = cnv$TYPE[i],
                                          "deldup" = if(cnv$TYPE[i] < 2) {"del"} else {"dup"},
                                          "AUC" = abs(2.0 - cnv$TYPE[i]) * (BPs[j+1] - BPs[j])))
          }
        }
      }
    }
    
    ITV_info <- data.frame("CHR" = rep(cnv$CHR[1L], length(BPs) - 1L),
                           "grid.id" = 1L:{length(BPs) - 1L},
                           "CNV.start" = utils::head(BPs, -1L),
                           "CNV.end" = utils::tail(BPs, -1L))
    
    
    list("long.cnv" = long_form[order(long_form$grid.id), ],
         "grid.info" = ITV_info)
  }
  
  cnv <- data.frame("ID" = 1L:4L,
                    "CHR" = rep(22L, 4L),
                    "BP1" = c(100L, 1000L, 1100L, 2000L),
                    "BP2" = c(500L, 1200L, 1500L, 3500L),
                    "TYPE" = c(0L, 0L, 1L, 1L))
  
  expected <- alternative_func(cnv)
  expect_equal(.breakCNV(cnv), expected)
})
