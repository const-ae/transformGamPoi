

test_that("acosh_transformation works", {


  mat <- matrix(rpois(n = 10 * 4, lambda = 0.3), nrow = 10, ncol = 4)
  expect_equal(acosh_transform(mat), .acosh_trans_impl(mat, 0.05))

  # Setting the overdispersion works
  alpha <- rnorm(10)^2
  expect_equal(acosh_transform(mat, overdispersion = alpha), .acosh_trans_impl(mat, alpha))

  # Exact zeros in overdispersion work
  sel <- which.max(rowSums(mat))
  alpha[sel] <- 0
  res <- acosh_transform(mat, overdispersion = alpha)
  expect_equal(res[-sel, ], .acosh_trans_impl(mat[-sel, ], alpha = alpha[-sel]))
  expect_equal(res[sel, ], .sqrt_trans_impl(mat[sel, ]))


  alpha <- matrix(rnorm(10 * 4)^2, nrow = 10, ncol = 4)
  expect_equal(acosh_transform(mat, overdispersion = alpha), .acosh_trans_impl(mat, alpha))

  # Check vector
  expect_equal(acosh_transform(1:3), .acosh_trans_impl(1:3, alpha = 0.05))
})

test_that("acosh_transformation works for sparse input", {
  mat <- matrix(rpois(n = 10 * 4, lambda = 0.3), nrow = 10, ncol = 4)
  sp_mat <- as(mat, "dgCMatrix")
  expect_equal(acosh_transform(sp_mat), as(acosh_transform(mat), "dgCMatrix"))
  expect_s4_class(acosh_transform(sp_mat), "CsparseMatrix")

  # Setting the overdispersion works
  alpha <- rnorm(10)^2
  expect_equal(acosh_transform(sp_mat, overdispersion = alpha),
               as(acosh_transform(mat, overdispersion = alpha), "dgCMatrix"))

  # Exact zeros in overdispersion work
  sel <- which.max(rowSums(mat))
  alpha[sel] <- 0
  expect_equal(acosh_transform(sp_mat, overdispersion = alpha),
               as(acosh_transform(mat, overdispersion = alpha), "dgCMatrix"))

})


test_that("the transition from acosh to sqrt is smooth", {

  expect_lt(acosh_transform(3, overdispersion = 1e-5), acosh_transform(3, overdispersion = 0))
  expect_equal(acosh_transform(3, overdispersion = 1e-5), acosh_transform(3, overdispersion = 0), tolerance = 1e-4)


})


test_that("shifted log is correct", {
  # Values between 1e-10 and 1e6
  xg <- rchisq(100, df = 0.3) * 1e6
  alpha <- 0.3
  expect_equal(.log_plus_alpha_impl(xg, alpha = alpha), 1/sqrt(alpha) * (log(xg + 1/(4 * alpha)) + log(4 * alpha)))

  # Zero stays zero
  expect_equal(.log_plus_alpha_impl(0, alpha = alpha), 1/sqrt(alpha) * (log(0 + 1/(4 * alpha)) + log(4 * alpha)))
})


test_that("shifted_log_transform errors if overdispersion and pseudo_count are specified", {
  expect_silent(shifted_log_transform(3, overdispersion = 0.01))
  expect_silent(shifted_log_transform(3, pseudo_count = 1))
  expect_silent(shifted_log_transform(3, overdispersion = 0.01, pseudo_count = 1/(4 * 0.01)))
  expect_error(shifted_log_transform(3, overdispersion = 0.01, pseudo_count = 1))
})


test_that("acosh, sqrt, and shifted log converge to each other", {

  expect_equal(shifted_log_transform(1e5),
               acosh_transform(1e5), tolerance = 1e-4)

  expect_equal(acosh_transform(5, overdispersion = 1e-6),
               2 * sqrt(5), tolerance = 1e-4)

})

