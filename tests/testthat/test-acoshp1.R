
test_that("acoshp1 works", {

  expect_equal(acoshp1(17), acosh(17 + 1))

})

test_that("sparse acosh plus 1 works", {

  mat <- matrix(rpois(n = 10 * 4, lambda = 0.3), nrow = 10, ncol = 4)
  expect_equal(acoshp1(mat), acosh(mat + 1))

  sp_mat <- as(mat, "dgCMatrix")
  expect_s4_class(acoshp1(sp_mat), "CsparseMatrix")
  expect_equal(acoshp1(sp_mat), as(acosh(mat + 1), "dgCMatrix"))
})
