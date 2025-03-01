test_that("`.isCNV returns expected errors", {
  
  # trigger missing(CNV) by not passing an input
  expect_error(.isCNV(),
               "`CNV` must be provided")
  
  # trigger is.data.frame by passing objects that resemble data.frames; such as
  # a list with 5 elements and a matrix with 5 columns
  expect_error(.isCNV(vector(mode = "list", length = 5L)),
               "`CNV` data.frame in PLINK format with columns: `ID`, `CHR`, `BP1`, `BP2`,`TYPE`")
  
  expect_error(.isCNV(matrix(1, 10, 5L)),
               "`CNV` data.frame in PLINK format with columns: `ID`, `CHR`, `BP1`, `BP2`,`TYPE`")
  
  # trigger the column header test -- test for "close"headers and missing headers
  expect_error(.isCNV(data.frame("Id"= 1L:10L,
                                 "CHR"= rep(22L, 10L),
                                 "BP1"= seq(1L, 20L, 2L),
                                 "BP2"= seq(2L, 20L, 2L),
                                 "TYPE"= c(0L, 0L, 1L, 1L,
                                           0L, 3L, 3L, 3L,
                                           4L, 4L))),
               "`CNV` data.frame in PLINK format with columns: `ID`, `CHR`, `BP1`, `BP2`,`TYPE`")
  
  expect_error(.isCNV(data.frame("ID"= 1L:10L,
                                 "ChR"= rep(22L, 10L),
                                 "BP1"= seq(1L, 20L, 2L),
                                 "BP2"= seq(2L, 20L, 2L),
                                 "TYPE"= c(0L, 0L, 1L, 1L,
                                           0L, 3L, 3L, 3L,
                                           4L, 4L))),
               "`CNV` data.frame in PLINK format with columns: `ID`, `CHR`, `BP1`, `BP2`,`TYPE`")
  
  expect_error(.isCNV(data.frame("ID"= 1L:10L,
                                 "CHR"= rep(22L, 10L),
                                 "BP1"= seq(1L, 20L, 2L),
                                 "BP2"= seq(2L, 20L, 2L))),
               "`CNV` data.frame in PLINK format with columns: `ID`, `CHR`, `BP1`, `BP2`,`TYPE`")
  
  expect_error(.isCNV(data.frame("ID"= 1L:10L,
                                 "CHR"= rep(22L, 10L),
                                 "BP1"= seq(1L, 20L, 2L),
                                 "TYPE"= c(0L, 0L, 1L, 1L,
                                           0L, 3L, 3L, 3L,
                                           4L, 4L))),
               "`CNV` data.frame in PLINK format with columns: `ID`, `CHR`, `BP1`, `BP2`,`TYPE`")
  
  
  
  # trigger integer vector test
  expect_error(.isCNV(data.frame("ID"= 1L:10L,
                                 "CHR"= rep(22.0, 10L),
                                 "BP1"= seq(1L, 20L, 2L),
                                 "BP2"= seq(2L, 20L, 2L),
                                 "TYPE"= c(0L, 0L, 1L, 1L,
                                           0L, 3L, 3L, 3L,
                                           4L, 4L))),
               "`CNV$CHR` must take a single value in [1, 22]",
               fixed = TRUE)
  
  expect_error(.isCNV(data.frame("ID"= 1L:10L,
                                 "CHR"= rep("22L", 10L),
                                 "BP1"= seq(1L, 20L, 2L),
                                 "BP2"= seq(2L, 20L, 2L),
                                 "TYPE"= c(0L, 0L, 1L, 1L,
                                           0L, 3L, 3L, 3L,
                                           4L, 4L))),
               "`CNV$CHR` must take a single value in [1, 22]",
               fixed = TRUE)
  
  # trigger only 1 value
  expect_error(.isCNV(data.frame("ID"= 1L:10L,
                                 "CHR"= c(rep(22L, 9L), 21L),
                                 "BP1"= seq(1L, 20L, 2L),
                                 "BP2"= seq(2L, 20L, 2L),
                                 "TYPE"= c(0L, 0L, 1L, 1L,
                                           0L, 3L, 3L, 3L,
                                           4L, 4L))),
               "`CNV$CHR` must take a single value in [1, 22]",
               fixed = TRUE)
  
  # trigger out of expected range [0, 22]
  expect_error(.isCNV(data.frame("ID"= 1L:10L,
                                 "CHR"= rep(23L, 10L),
                                 "BP1"= seq(1L, 20L, 2L),
                                 "BP2"= seq(2L, 20L, 2L),
                                 "TYPE"= c(0L, 0L, 1L, 1L,
                                           0L, 3L, 3L, 3L,
                                           4L, 4L))),
               "`CNV$CHR` must take a single value in [1, 22]",
               fixed = TRUE)
  
  
  # trigger integer vector test
  expect_error(.isCNV(data.frame("ID"= 1L:10L,
                                 "CHR"= rep(22L, 10L),
                                 "BP1"= as.numeric(seq(1L, 20L, 2L)) -3.1,
                                 "BP2"= seq(2L, 20L, 2L),
                                 "TYPE"= c(0L, 0L, 1L, 1L,
                                           0L, 3L, 3L, 3L,
                                           4L, 4L))),
               "`CNV$BP1 must contain non-negative numbers",
               fixed = TRUE)
  # trigger all > 0 test
  expect_error(.isCNV(data.frame("ID"= 1L:10L,
                                 "CHR"= rep(22L, 10L),
                                 "BP1"= seq(1L, 20L, 2L) - 2L,
                                 "BP2"= seq(2L, 20L, 2L),
                                 "TYPE"= c(0L, 0L, 1L, 1L,
                                           0L, 3L, 3L, 3L,
                                           4L, 4L))),
               "`CNV$BP1 must contain non-negative numbers",
               fixed = TRUE)
  
  # "`CNV$BP2 must contain non-negative numbers"=
  #      is.vector(CNV$BP2, mode = "integer") &&
  #      all(CNV$BP2 >= 0L),
  # trigger integer vector test
  expect_error(.isCNV(data.frame("ID"= 1L:10L,
                                 "CHR"= rep(22L, 10L),
                                 "BP1"= seq(1L, 20L, 2L),
                                 "BP2"= as.numeric(seq(2L, 20L, 2L)) -2.1,
                                 "TYPE"= c(0L, 0L, 1L, 1L,
                                           0L, 3L, 3L, 3L,
                                           4L, 4L))),
               "`CNV$BP2 must contain non-negative numbers",
               fixed = TRUE)
  
  # trigger all > 0 test
  expect_error(.isCNV(data.frame("ID"= 1L:10L,
                                 "CHR"= rep(22L, 10L),
                                 "BP1"= seq(1L, 20L, 2L),
                                 "BP2"= seq(2L, 20L, 2L) - 3L,
                                 "TYPE"= c(0L, 0L, 1L, 1L,
                                           0L, 3L, 3L, 3L,
                                           4L, 4L))),
               "`CNV$BP2 must contain non-negative numbers",
               fixed = TRUE)
  
  #  "`CNV$BP2 must be great than CNV$BP1 for all records"=
  #        !any(CNV$BP1 > CNV$BP2),
  
  # trigger for all greater
  expect_error(.isCNV(data.frame("ID"= 1L:10L,
                                 "CHR"= rep(22L, 10L),
                                 "BP1"= seq(2L, 20L, 2L),
                                 "BP2"= seq(1L, 20L, 2L),
                                 "TYPE"= c(0L, 0L, 1L, 1L,
                                           0L, 3L, 3L, 3L,
                                           4L, 4L))),
               "`CNV$BP2 must be great than CNV$BP1 for all records",
               fixed = TRUE)
  # trigger for 1 greater
  expect_error(.isCNV(data.frame("ID"= 1L:10L,
                                 "CHR"= rep(22L, 10L),
                                 "BP1"= seq(1L, 20L, 2L),
                                 "BP2"= c( 2L,  2L,  6L,  8L, 10L,
                                           12L, 14L, 16L, 18L, 20L),
                                 "TYPE"= c(0L, 0L, 1L, 1L,
                                           0L, 3L, 3L, 3L,
                                           4L, 4L))),
               "`CNV$BP2 must be great than CNV$BP1 for all records",
               fixed = TRUE)
  
  # "`CNV$Type` must be integer in [0,1, 3, ..., 9999]"=
  #        is.vector(CNV$TYPE, mode = "integer") &&
  #        all(CNV$TYPE %in% c(0L, 1L, 2L, 3L, 4L))
  
  # trigger integer vector test
  expect_error(.isCNV(data.frame("ID"= 1L:10L,
                                 "CHR"= rep(22L, 10L),
                                 "BP1"= seq(1L, 20L, 2L),
                                 "BP2"= seq(2L, 20L, 2L),
                                 "TYPE"= as.numeric(c(0L, 0L, 1L, 1L,
                                                      0L, 3L, 3L, 3L,
                                                      4L, 4L)))),
               "`CNV$Type` must be integer in [0, 1, 3, ..., 9999]",
               fixed = TRUE)
  #out of range
  
  expect_error(.isCNV(data.frame("ID"= 1L:10L,
                                 "CHR"= rep(22L, 10L),
                                 "BP1"= seq(1L, 20L, 2L),
                                 "BP2"= seq(2L, 20L, 2L),
                                 "TYPE"= c(-1L, 0L, 1L, 1L,
                                           0L, 3L, 3L, 3L,
                                           4L, 4L))),
               "`CNV$Type` must be integer in [0, 1, 3, ..., 9999]",
               fixed = TRUE)
  
  expect_error(.isCNV(data.frame("ID"= 1L:10L,
                                 "CHR"= rep(22L, 10L),
                                 "BP1"= seq(1L, 20L, 2L),
                                 "BP2"= seq(2L, 20L, 2L),
                                 "TYPE"= c(2L, 0L, 1L, 1L,
                                           0L, 3L, 3L, 3L,
                                           4L, 4L))),
               "`CNV$Type` must be integer in [0, 1, 3, ..., 9999]",
               fixed = TRUE)

  expect_equal(.isCNV(data.frame("ID"= 1L:10L,
                                 "CHR"= rep(22L, 10L),
                                 "BP1"= seq(100L, 2000L, 200L),
                                 "BP2"= seq(200L, 2000L, 200L),
                                 "TYPE"= c(0L, 0L, 1L, 1L,
                                           0L, 3L, 3L, 3L,
                                           4L, 4L))), TRUE)
})

test_that("`.testIterControl` returns expected errors", {
  
  expect_error(.testIterControl(),
               "`iter.control` must be provided")
  
  expect_error(.testIterControl(c("max.iter" = 3)),
               "`iter.control` must be a list")
  
  expect_warning(.testIterControl(list("maxiter" = 5)),
                 "ignored iter.control elements maxiter")
  
  expect_warning(.testIterControl(list("maxiter" = 5, "tolbeta" = 0.01)),
                 "ignored iter.control elements maxiter, tolbeta")
  
  expect_error(.testIterControl(list("max.iter" = c(1, 2))),
               "iter.control$max.iter is not appropriately defined",
               fixed = TRUE)
  expect_error(.testIterControl(list("max.iter" = 1)),
               "iter.control$max.iter is not appropriately defined",
               fixed = TRUE)  
  
  expect_error(.testIterControl(list("max.iter" = 100,
                                     "tol.beta" = 0)),
               "iter.control$tol.beta is not appropriately defined (0 < tol.beta <= 1e-2)",
               fixed = TRUE)  
  expect_error(.testIterControl(list("max.iter" = 100,
                                     "tol.beta" = 0.1)),
               "iter.control$tol.beta is not appropriately defined (0 < tol.beta <= 1e-2)",
               fixed = TRUE)  
  expect_error(.testIterControl(list("max.iter" = 100,
                                     "tol.beta" = c(0.0001, 0.0001))),
               "iter.control$tol.beta is not appropriately defined (0 < tol.beta <= 1e-2)",
               fixed = TRUE)  
  
  expect_error(.testIterControl(list("max.iter" = 100,
                                     "tol.beta" = 0.01,
                                     "tol.loss" = 0)),
               "iter.control$tol.loss is not appropriately defined (0 < tol.loss <= 1e-3)",
               fixed = TRUE)  
  expect_error(.testIterControl(list("max.iter" = 100,
                                     "tol.beta" = 0.01,
                                     "tol.loss" = 0.1)),
               "iter.control$tol.loss is not appropriately defined (0 < tol.loss <= 1e-3)",
               fixed = TRUE)  
  expect_error(.testIterControl(list("max.iter" = 100,
                                     "tol.beta" = 0.01,
                                     "tol.loss" = c(0.001, 0.001))),
               "iter.control$tol.loss is not appropriately defined (0 < tol.loss <= 1e-3)",
               fixed = TRUE)  
})

test_that("`.testIterControl` returns expected results", {
  
  default_iter <- list(max.iter = 8L, tol.beta = 1e-3, tol.loss = 1e-6)
  
  expect_equal(.testIterControl(NULL), default_iter)
  expect_equal(.testIterControl(list("max.iter" = 100)),
               list(max.iter = 100L, tol.beta = 1e-3, tol.loss = 1e-6))
  expect_equal(.testIterControl(list("tol.beta" = 0.002)),
               list(max.iter = 8L, tol.beta = 0.002, tol.loss = 1e-6))
  expect_equal(.testIterControl(list("tol.loss" = 0.0004)),
               list(max.iter = 8L, tol.beta = 1e-3, tol.loss = 0.0004))
  expect_equal(.testIterControl(list("tol.loss" = 0.0004, "tol.beta" = 0.002)),
               list(max.iter = 8L, tol.beta = 0.002, tol.loss = 0.0004))
  expect_equal(.testIterControl(list("tol.loss" = 0.0004, "tol.beta" = 0.002, "max.iter" = 5)),
               list(max.iter = 5L, tol.beta = 0.002, tol.loss = 0.0004))
  expect_warning(out <- .testIterControl(list("tolloss" = 0.0004, "tolbeta" = 0.002, "maxiter" = 5)))
  expect_equal(out, default_iter)
})

test_that("`.testCVControl returns expected errors", {
  
  expect_error(.testCVControl(),
               "`cv.control` must be provided")
  
  expect_error(.testCVControl(c("n.fold" = 3)),
               "`cv.control` must be a list")
  
  expect_warning(.testCVControl(list("nfold" = 5), family = "gaussian"),
                 "ignored cv.control elements nfold")
  
  expect_warning(.testCVControl(list("nfold" = 5, "ncore" = 1), family = "gaussian"),
                 "ignored cv.control elements nfold, ncore")
  
  expect_error(.testCVControl(list("n.fold" = c(1, 2))),
               "cv.control$n.fold is not appropriately defined",
               fixed = TRUE)
  expect_error(.testCVControl(list("n.fold" = 1)),
               "cv.control$n.fold is not appropriately defined",
               fixed = TRUE)  
  
  expect_error(.testCVControl(list("n.fold" = 4,
                                   "n.core" = 1:2)),
               "cv.control$n.core is not appropriately defined",
               fixed = TRUE)  
  expect_error(.testCVControl(list("n.fold" = 4,
                                   "n.core" = "1")),
               "cv.control$n.core is not appropriately defined",
               fixed = TRUE)  
  
  expect_error(.testCVControl(list("n.fold" = 4,
                                   "n.core" = 2,
                                   "stratified" = 0)),
               "cv.control$stratified is not appropriately defined",
               fixed = TRUE)  
  expect_error(.testCVControl(list("n.fold" = 4,
                                   "n.core" = 2,
                                   "stratified" = c(FALSE, TRUE))),
               "cv.control$stratified is not appropriately defined",
               fixed = TRUE)  

})

test_that("`.testCVControl` returns expected results", {
  
  default_cv <- list("n.fold" = 5L, "n.core" = 1L, "stratified" = FALSE)
  
  expect_equal(.testCVControl(NULL, family = "gaussian"), default_cv)
  expect_equal(.testCVControl(list("n.fold" = 100), family = "gaussian"),
               list("n.fold" = 100L, "n.core" = 1L, "stratified" = FALSE))
  expect_equal(.testCVControl(list("n.core" = 16), family = "gaussian"),
               list("n.fold" = 5L, "n.core" = 16L, "stratified" = FALSE))
  expect_equal(.testCVControl(list("stratified" = TRUE), family = "gaussian"),
               list("n.fold" = 5L, "n.core" = 1L, "stratified" = FALSE))
  expect_equal(.testCVControl(list("stratified" = TRUE), family = "binomial"),
               list("n.fold" = 5L, "n.core" = 1L, "stratified" = TRUE))
  expect_equal(.testCVControl(list("stratified" = TRUE, "n.core" = 8), family = "binomial"),
               list("n.fold" = 5L, "n.core" = 8L, "stratified" = TRUE))
  expect_equal(.testCVControl(list("stratified" = TRUE, "n.core" = 8, "n.fold" = 15), family = "binomial"),
               list("n.fold" = 15L, "n.core" = 8L, "stratified" = TRUE))
  expect_warning(out <- .testCVControl(list("stratifieD" = TRUE, "ncore" = 8, "nfold" = 15), family = "gaussian"))
  expect_equal(out, default_cv)
})