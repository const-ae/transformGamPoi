

#'  Residual-based Variance Stabilizing Transformation
#'
#'  Fit an intercept Gamma-Poisson model that corrects for sequencing depth and return the residuals
#'  as variance stabilized results for further downstream application, for which no proper count-based
#'  method exist or is performant enough (e.g., clustering, dimensionality reduction).
#'
#' @param data any matrix-like object (e.g. matrix, DelayedArray, HDF5Matrix)
#'   with one column per sample and row per gene. It can also be an object of type `glmGamPoi`,
#'   in which case it is directly used to calculate the variance-stabilized values.
#' @param offset_model boolean to decide if \eqn{\beta_1} in \eqn{y = \beta_0 + \beta_1 log(sf)},
#'   is set to 1 (i.e., treating the log of the size factors as an offset ) or is estimated per gene.
#'   From a theoretical point, it should be fine to treat \eqn{\beta_1} as an offset, because a cell that is
#'   twice as big, should have twice as many counts per gene (without any gene-specific effects).
#'   However, `sctransform` suggested that it would be advantageous to nonetheless estimate
#'   \eqn{\beta_0} as it may counter data artifacts. On the other side, Lause et al. (2020)
#'   demonstrated that the estimating \eqn{\beta_0} and \eqn{\beta_1} together can be difficult. If
#'   you still want to fit `sctransform`'s model, you can set the `ridge_penalty` argument to a
#'   non-zero value, which shrinks \eqn{\beta_1} towards 1 and resolved the degeneracy. \cr
#'   Default: `TRUE`.
#' @param residual_type a string that specifies what kind of residual is returned as variance stabilized-value.
#'   \describe{
#'     \item{`"randomized_quantile"`}{The discrete nature of count distribution stops simple transformations from
#'     obtaining a truly standard normal residuals. The trick of of quantile randomized residuals is to match the
#'     cumulative density function of the Gamma-Poisson and the Normal distribution. Due to the discrete nature of
#'     Gamma-Poisson distribution, a count does not correspond to a single quantile of the Normal distribution, but
#'     to a range of possible value. This is resolved by randomly choosing one of the mapping values from the
#'     Normal distribution as the residual. This ensures perfectly normal distributed
#'     residuals, for the cost of introducing randomness. More details are available in the documentation
#'     of [`statmod::qresiduals()`] and the corresponding publication by Dunn and Smyth (1996).}
#'     \item{`"pearson"`}{The Pearson residuals are defined as \eqn{res = (y - m) / sqrt(m + m^2 * theta)}.}
#'   }
#'   Default: `"randomized_quantile"`
#' @param overdispersion,overdispersion_shrinkage,size_factors arguments that are passed to the underlying
#'   call to [`glmGamPoi::glm_gp()`]. Default for each: `TRUE`.
#' @param ridge_penalty another argument that is passed to [`glmGamPoi::glm_gp()`]. It is ignored if
#'   `offset_model = TRUE`. Default: `2`.
#' @param return_fit boolean to decide if the matrix of residuals is returned directly (`return_fit = FALSE`)
#'   or if in addition the `glmGamPoi`-fit is returned (`return_fit = TRUE`) . Default: `FALSE`.
#' @param verbose boolean that decides if information about the individual steps are printed.
#'   Default: `FALSE`
#' @param ... additional parameters passed to [`glmGamPoi::glm_gp()`].
#'
#'
#'
#' @seealso [`glmGamPoi::glm_gp()`], [`glmGamPoi::residuals.glmGamPoi()`], [`sctransform::vst()`],
#'   [`statmod::qresiduals()`]
#'
#' @references
#'   Ahlmann-Eltze, Constantin, and Wolfgang Huber. "glmGamPoi: Fitting Gamma-Poisson Generalized Linear
#'   Models on Single Cell Count Data." Bioinformatics (2020)
#'
#'   Dunn, Peter K., and Gordon K. Smyth. "Randomized quantile residuals." Journal of Computational and
#'   Graphical Statistics 5.3 (1996): 236-244.
#'
#'   Hafemeister, Christoph, and Rahul Satija. "Normalization and variance stabilization of single-cell
#'   RNA-seq data using regularized negative binomial regression." Genome biology 20.1 (2019): 1-15.
#'
#'   Hafemeister, Christoph, and Rahul Satija. "Analyzing scRNA-seq data with the sctransform and offset
#'   models." (2020)
#'
#'   Lause, Jan, Philipp Berens, and Dmitry Kobak. "Analytic Pearson residuals for normalization of
#'   single-cell RNA-seq UMI data." bioRxiv (2020).
#'
#' @examples
#'  # Make example data
#'  n_genes <- 100
#'  n_cells <- 200
#'  beta0 <- rnorm(n = n_genes, mean = 2, sd = 0.3)
#'  beta1 <- rnorm(n = n_genes, mean = 0, sd = 2.5)
#'  X <- cbind(1, rep_len(c(-1,1), n_cells))
#'  sf <- rchisq(n = n_cells, df = 100)
#'  sf <- sf / mean(sf)
#'  Mu <- exp( cbind(beta0, beta1) %*% t(X) +
#'        matrix(log(sf), nrow = n_genes, ncol = n_cells, byrow = TRUE) )
#'
#'  # Generate count data
#'  Y <- matrix(rnbinom(n = n_genes * n_cells, mu = Mu, size = 0.1),
#'              nrow = n_genes, ncol = n_cells)
#'
#'  # Apply VST
#'  vst <- gamPoiTransform(Y, verbose = TRUE)
#'
#'  # Plot the PCA of the result
#'  vst_pca <- prcomp(t(vst), rank. = 2)
#'  plot(vst_pca$x, col = as.factor((X[,2] == -1)))
#'
#' @export
gamPoiTransform <- function(data,
                            offset_model = TRUE,
                            residual_type = c("randomized_quantile", "pearson"),
                            size_factors = TRUE,
                            overdispersion = TRUE,
                            overdispersion_shrinkage = TRUE,
                            ridge_penalty = 2,
                            return_fit = FALSE,
                            verbose = FALSE, ...){
  residual_type <- residual_type[1]

  if(inherits(data, "glmGamPoi")){
    fit <- data
  }else if(offset_model){
    fit <- glmGamPoi::glm_gp(data, design = ~ 1, size_factors = size_factors,
                             overdispersion = overdispersion,
                             overdispersion_shrinkage = overdispersion_shrinkage,
                             verbose = verbose, ...)
  }else{
    # Calculate size factors
    if(isTRUE(size_factors)){
      size_factors <- matrixStats::colSums2(data)
      size_factors <- size_factors / mean(size_factors)
    }
    log_sf <- log(size_factors)
    attr(ridge_penalty, "target") <- c(0, 1)

    fit <- glmGamPoi::glm_gp(data, design = ~ log_sf + 1, size_factors = 1,
                             overdispersion = overdispersion,
                             overdispersion_shrinkage = overdispersion_shrinkage,
                             ridge_penalty = ridge_penalty,
                             verbose = verbose, ...)
  }

  if(overdispersion_shrinkage){
    # Use the dispersion trend when calculating the residuals
    fit$overdispersion_shrinkage_list$original_overdispersions <- fit$overdispersions
    fit$overdispersions <- fit$overdispersion_shrinkage_list$dispersion_trend
  }


  if(verbose){message("Calculate ", residual_type, " residuals")}

  if(! return_fit){
    stats::residuals(fit, type = residual_type)
  }else{
    list(Residuals = stats::residuals(fit, type = residual_type),
         fit = fit)
  }
}


