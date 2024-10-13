XZ_test =as(matrix(c( 2, 2, 0, 0, 6,	10, 1,	3,	0,
                      0, 2, 4, 1,	3,	0, 0, 2, 4,
                      0, 0, 0, 0, 0, 0,1,	0,	0,
                      0, 0, 0, 0, 0, 0,1,	0,	0,
                      2, 0, 0, 1,	0,	0,4, 0,	0,
                      0, 0, 0, 0, 0, 0,1,	3,	0,
                      0, 0, 0, 0, 0, 0,0, 0,	0,
                      0, 0, 4, 0,	0, 10, 0,1,	3), ncol = 9, nrow=8, byrow = TRUE, dimnames = list(paste0("id", c(1L:8L)), c(paste0("del",c(2L,4L,5L)),paste0("dup",c(2L,4L,5L)),paste0("pc", c(1L:3L))))), "sparseMatrix")
y_test0 = vector(mode = "numeric", length = 8L)
y_test1= vector(mode = "numeric", length = 8L)
names(y_test1)<-paste0("id", c(1L:8L))
A_test=matrix(1, 7L, 9L)

test_that("`.ctnsSolution()` returns expected errors", {
  
  
  # trigger missing by not passing an input
  expect_error(.ctnsSolution(),
               "`data` must be a 'WTsth.data' object")
  
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
  data <- prep(cnv, Z, Y, rare.out = 0.03)  
  
})

#STH YOU ARE HERE!!!!!!