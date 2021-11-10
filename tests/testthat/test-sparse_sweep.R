test_that("sweep for dgCMatrices works", {
  mat <- matrix(0, nrow =100, ncol = 17)
  mat[sample.int(length(mat), 20)] <- rnorm(n = 20)
  sp_mat <- as(mat, "dgCMatrix")

  sparse_divide_out_size_factor(sp_mat, 1:17)
  expect_equal(sparse_divide_out_size_factor(sp_mat, 8), DelayedArray::sweep(sp_mat, 2, 8, FUN = "/"))
  expect_equal(sparse_divide_out_size_factor(sp_mat, 1:17), DelayedArray::sweep(sp_mat, 2, 1:17, FUN = "/"))
})






