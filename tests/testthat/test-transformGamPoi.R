test_that("residual_transform works", {
    n_genes <- 100
    n_cells <- 500

    beta0 <- rnorm(n = n_genes, mean = 2, sd = 0.3)
    sf <- rchisq(n = n_cells, df = 100)
    sf <- sf / mean(sf)

    Mu <- exp( beta0 %*% t(log(sf)) )

    Y <- matrix(rnbinom(n = n_genes * n_cells, mu = Mu, size = 0.1), nrow = n_genes, ncol = n_cells)

    summary(MatrixGenerics::colMeans2(Y))
    summary(MatrixGenerics::rowMeans2(Y))


    resids <- residual_transform(Y, verbose = FALSE)
    res2 <- residual_transform(Y, offset_model = FALSE, verbose = FALSE, return_fit = TRUE)
    expect_true(all(abs(res2$fit$Beta[,2] - 1) < 0.1))

})


test_that("different input types work", {
    n_genes <- 100
    n_cells <- 500

    beta0 <- rnorm(n = n_genes, mean = 2, sd = 0.3)
    sf <- rchisq(n = n_cells, df = 100)
    sf <- sf / mean(sf)

    Mu <- exp( beta0 %*% t(log(sf)) )

    Y <- matrix(rnbinom(n = n_genes * n_cells, mu = Mu, size = 0.1), nrow = n_genes, ncol = n_cells)

    # matrix
    res <- residual_transform(Y, verbose = FALSE, return_fit = TRUE, residual_type = "pearson")
    # glmGamPoi
    res2 <- residual_transform(res$fit, residual_type = "pearson")
    # SummarizedExperiment
    res3 <- residual_transform(res$fit$data, residual_type = "pearson")

    expect_equal(res$Residuals, res2)
    expect_equal(res$Residuals, res3)
})


test_that("overdisperion = 'global' works", {
    set.seed(1)
    n_genes <- 100
    n_cells <- 500

    beta0 <- rnorm(n = n_genes, mean = 2, sd = 0.3)
    sf <- rchisq(n = n_cells, df = 100)
    sf <- sf / mean(sf)

    Mu <- exp( beta0 %*% t(log(sf)) )

    Y <- matrix(rnbinom(n = n_genes * n_cells, mu = Mu, size = 1/0.1), nrow = n_genes, ncol = n_cells)

    tmp <- transformGamPoi(Y, "rand", overdispersion = "global", verbose = FALSE, on_disk = FALSE, return_fit = TRUE)
    expect_equal(tmp$fit$overdispersions, rep(0.1, n_genes), tolerance = 0.1)
})


test_that("transformGamPoi errors if 'residual_type' is specified", {

    set.seed(1)
    Y <- matrix(rnbinom(n = 24, mu = 3, size = 1/0.1), nrow = 3, ncol = 8)
    expect_error(transformGamPoi(Y, residual_type = "pearson"))

})

