test_that("transformGamPoi works", {
    n_genes <- 100
    n_cells <- 500

    beta0 <- rnorm(n = n_genes, mean = 2, sd = 0.3)
    sf <- rchisq(n = n_cells, df = 100)
    sf <- sf / mean(sf)

    Mu <- exp( beta0 %*% t(log(sf)) )

    Y <- matrix(rnbinom(n = n_genes * n_cells, mu = Mu, size = 0.1), nrow = n_genes, ncol = n_cells)

    summary(matrixStats::colMeans2(Y))
    summary(matrixStats::rowMeans2(Y))


    resids <- transformGamPoi(Y, verbose = TRUE)
    res2 <- transformGamPoi(Y, offset_model = FALSE, verbose = TRUE, return_fit = TRUE)
    expect_true(all(abs(res2$fit$Beta[,2] - 1) < 0.1))

})





