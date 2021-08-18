

test_that("acosh_transformation works", {


  mat <- matrix(rpois(n = 10 * 4, lambda = 0.3), nrow = 10, ncol = 4)
  expect_equal(acosh_transform(mat, size_factors = 1), .acosh_trans_impl(mat, 0.05))

  # Setting the overdispersion works
  alpha <- rnorm(10)^2
  expect_equal(acosh_transform(mat, overdispersion = alpha, size_factors = 1), .acosh_trans_impl(mat, alpha))

  # Exact zeros in overdispersion work
  sel <- which.max(rowSums(mat))
  alpha[sel] <- 0
  res <- acosh_transform(mat, overdispersion = alpha, size_factors = 1)
  expect_equal(res[-sel, ], .acosh_trans_impl(mat[-sel, ], alpha = alpha[-sel]))
  expect_equal(res[sel, ], .sqrt_trans_impl(mat[sel, ]))


  alpha <- matrix(rnorm(10 * 4)^2, nrow = 10, ncol = 4)
  expect_equal(acosh_transform(mat, overdispersion = alpha, size_factors = 1), .acosh_trans_impl(mat, alpha))

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

test_that("acosh_transform calling to glm_gp works as expected", {
  mat <- matrix(rpois(n = 10 * 4, lambda = 0.3), nrow = 10, ncol = 4)
  sp_mat <- as(mat, "dgCMatrix")

  expect_equal(acosh_transform(sp_mat, overdispersion = TRUE, on_disk = FALSE),
               acosh_transform(mat, overdispersion = TRUE), ignore_attr = TRUE)

  fit <- glmGamPoi::glm_gp(mat, overdispersion_shrinkage = TRUE)
  expect_equal(acosh_transform(mat, overdispersion = TRUE),
               acosh_transform(mat, overdispersion = fit$overdispersion_shrinkage_list$dispersion_trend), ignore_attr = TRUE)

  fit <- glmGamPoi::glm_gp(mat, overdispersion_shrinkage = FALSE)
  expect_equal(acosh_transform(mat, overdispersion = TRUE, overdispersion_shrinkage = FALSE),
               acosh_transform(mat, overdispersion = fit$overdispersions), ignore_attr = TRUE)
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
})


test_that("acosh, sqrt, and shifted log converge to each other", {

  expect_equal(shifted_log_transform(1e5),
               acosh_transform(1e5), tolerance = 1e-4)

  expect_equal(acosh_transform(5, overdispersion = 1e-6),
               2 * sqrt(5), tolerance = 1e-4)

})


test_that("different input types work", {

  n_genes <- 100
  n_cells <- 500

  beta0 <- rnorm(n = n_genes, mean = 2, sd = 0.3)
  sf <- rchisq(n = n_cells, df = 100)
  sf <- sf / mean(sf)

  Mu <- exp( beta0 %*% t(log(sf)) )

  Y <- matrix(rnbinom(n = n_genes * n_cells, mu = Mu, size = 0.1), nrow = n_genes, ncol = n_cells)

  fit <- glmGamPoi::glm_gp(Y, design = ~ 1)

  # matrix
  res <- acosh_transform(Y, verbose = TRUE)
  # glmGamPoi
  res2 <- acosh_transform(fit)
  # SummarizedExperiment
  res3 <- acosh_transform(fit$data)

  expect_equal(res, res2)
  expect_equal(res, res3)

})


test_that("overdispersion handling works", {

  n_genes <- 100
  n_cells <- 500
  beta0 <- rnorm(n = n_genes, mean = 2, sd = 0.3)
  sf <- rchisq(n = n_cells, df = 100)
  sf <- sf / mean(sf)
  Mu <- exp( beta0 %*% t(log(sf)) )
  Y <- matrix(rnbinom(n = n_genes * n_cells, mu = Mu, size = 0.1), nrow = n_genes, ncol = n_cells)

  fit <- glmGamPoi::glm_gp(Y, design = ~ 1)
  res1 <- acosh_transform(Y, overdispersion = TRUE)
  res2 <- acosh_transform(Y, overdispersion = fit$overdispersion_shrinkage_list$dispersion_trend)

  expect_equal(res1, res2)

  res3 <- shifted_log_transform(Y, overdispersion = TRUE)
  res4 <- shifted_log_transform(Y, overdispersion = fit$overdispersion_shrinkage_list$dispersion_trend)
  expect_equal(res3, res4)

  fit <- glmGamPoi::glm_gp(Y, design = ~ 1, overdispersion_shrinkage = FALSE)
  res5 <- acosh_transform(Y, overdispersion = TRUE)
  res6 <- acosh_transform(Y, overdispersion = fit$overdispersion_shrinkage_list$dispersion_trend)

  expect_equal(res5, res6)

  res7 <- shifted_log_transform(Y, overdispersion = TRUE)
  res8 <- shifted_log_transform(Y, overdispersion = fit$overdispersion_shrinkage_list$dispersion_trend)
  expect_equal(res7, res8)

})



test_that("on_disk works", {
  n_genes <- 100
  n_cells <- 30
  beta0 <- rnorm(n = n_genes, mean = 2, sd = 0.3)
  sf <- rchisq(n = n_cells, df = 100)
  sf <- sf / mean(sf)
  Mu <- exp( beta0 %*% t(log(sf)) )
  Y <- matrix(rnbinom(n = n_genes * n_cells, mu = Mu, size = 0.1), nrow = n_genes, ncol = n_cells)

  Y_hdf5 <- HDF5Array::writeHDF5Array(Y)

  res <- acosh_transform(Y)
  res1 <- acosh_transform(Y_hdf5)
  res2 <- acosh_transform(Y, on_disk = TRUE)

  expect_s4_class(res1, "DelayedMatrix")
  expect_s4_class(res2, "DelayedMatrix")

  expect_equal(res, as.matrix(res1))
  expect_equal(res, as.matrix(res2))

  Y_sp <- as(Y, "dgCMatrix")
  Y_sp_hdf5 <- HDF5Array::writeHDF5Array(Y_sp)
  res3 <- acosh_transform(Y_sp)
  res4 <- acosh_transform(Y_sp_hdf5)
  res5 <- acosh_transform(Y_sp, on_disk = TRUE)

  expect_true(HDF5Array::is_sparse(res3))
  expect_true(HDF5Array::is_sparse(res4))
  expect_true(HDF5Array::is_sparse(res5))
  expect_equal(res3, as(res4, "dgCMatrix"), ignore_attr = TRUE)
  expect_equal(res3, as(res5, "dgCMatrix"), ignore_attr = TRUE)
})


test_that("Clipping works for Pearson residuals", {

  Y <- matrix(rnbinom(n = 10 * 30, mu = 3, size = 1/0.15), nrow = 10, ncol = 5)
  resid1 <- residual_transform(Y, "pearson", clipping = FALSE)
  resid2 <- residual_transform(Y, "pearson", clipping = TRUE)

  clip <- sqrt(5)
  large_clip <- resid1 > clip
  small_clip <- resid1 < -clip

  expect_equal(resid2[large_clip], rep(clip, sum(large_clip)))
  expect_equal(resid2[small_clip], rep(-clip, sum(small_clip)))
  expect_equal(resid2[! (large_clip | small_clip)],
               resid1[! (large_clip | small_clip)])


  resid3 <- residual_transform(Y, "pearson", clipping = 0.234)
  clip <- 0.234
  large_clip <- resid1 > clip
  small_clip <- resid1 < -clip

  expect_equal(resid3[large_clip], rep(clip, sum(large_clip)))
  expect_equal(resid3[small_clip], rep(-clip, sum(small_clip)))
  expect_equal(resid3[! (large_clip | small_clip)],
               resid1[! (large_clip | small_clip)])

})

