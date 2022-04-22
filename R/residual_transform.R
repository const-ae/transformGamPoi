

#'  Residual-based Variance Stabilizing Transformation
#'
#'  Fit an intercept Gamma-Poisson model that corrects for sequencing depth and return the residuals
#'  as variance stabilized results for further downstream application, for which no proper count-based
#'  method exist or is performant enough (e.g., clustering, dimensionality reduction).
#'
#' @param offset_model boolean to decide if \eqn{\beta_1} in \eqn{y = \beta_0 + \beta_1 log(sf)},
#'   is set to 1 (i.e., treating the log of the size factors as an offset ) or is estimated per gene.
#'   From a theoretical point, it should be fine to treat \eqn{\beta_1} as an offset, because a cell that is
#'   twice as big, should have twice as many counts per gene (without any gene-specific effects).
#'   However, `sctransform` suggested that it would be advantageous to nonetheless estimate
#'   \eqn{\beta_0} as it may counter data artifacts. On the other side, Lause et al. (2020)
#'   demonstrated that the estimating \eqn{\beta_0} and \eqn{\beta_1} together can be difficult. If
#'   you still want to fit `sctransform`'s model, you can set the `ridge_penalty` argument to a
#'   non-zero value, which shrinks \eqn{\beta_1} towards 1 and resolves the degeneracy. \cr
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
#'     \item{`"analytic_pearson"`}{Similar to the method above, however, instead of estimating \eqn{m} using a
#'     GLM model fit, \eqn{m} is approximated by \eqn{m_ij = (\sum_j y_{ij}) (\sum_i y_{ij}) / (\sum_{i,j} y_{ij})}.
#'     For all details, see Lause et al. (2021).
#'     Note that `overdispersion_shrinkage` and `ridge_penalty` are ignored when fitting analytic Pearson residuals.}
#'   }
#'   The two above options are the most common choices, however you can use any `residual_type` supported by
#'   [`glmGamPoi::residuals.glmGamPoi()`]. Default: `"randomized_quantile"`
#' @param clipping a single boolean or numeric value specifying that all residuals are in the range
#'   `[-clipping, +clipping]`. If `clipping = TRUE`, we use the default of `clipping = sqrt(ncol(data))`
#'   which is the default behavior for `sctransform`. Default: `FALSE`, which means no clipping is applied.
#' @param overdispersion_shrinkage,size_factors arguments that are passed to the underlying
#'   call to [`glmGamPoi::glm_gp()`]. Default for each: `TRUE`.
#' @param ridge_penalty another argument that is passed to [`glmGamPoi::glm_gp()`]. It is ignored if
#'   `offset_model = TRUE`. Default: `2`.
#' @param return_fit boolean to decide if the matrix of residuals is returned directly (`return_fit = FALSE`)
#'   or if in addition the `glmGamPoi`-fit is returned (`return_fit = TRUE`) . Default: `FALSE`.
#' @param ... additional parameters passed to [`glmGamPoi::glm_gp()`].
#' @inherit transformGamPoi
#'
#' @details
#'  Internally, this method uses the `glmGamPoi` package. The function goes through the following steps
#'  \enumerate{
#'    \item fit model using [`glmGamPoi::glm_gp()`]
#'    \item plug in the trended overdispersion estimates
#'    \item call [`glmGamPoi::residuals.glmGamPoi()`] to calculate the residuals.
#'  }
#'
#' @return a matrix (or a vector if the input is a vector) with the transformed values. If `return_fit = TRUE`,
#'   a list is returned with two elements: `fit` and `Residuals`.
#'
#' @seealso [`glmGamPoi::glm_gp()`], [`glmGamPoi::residuals.glmGamPoi()`], `sctransform::vst()`,
#'   `statmod::qresiduals()`
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
#'   single-cell RNA-seq UMI data." Genome Biology (2021).
#'
#' @examples
#'   # Load a single cell dataset
#'   sce <- TENxPBMCData::TENxPBMCData("pbmc4k")
#'   # Reduce size for this example
#'   set.seed(1)
#'   sce_red <- sce[sample(which(rowSums2(counts(sce)) > 0), 1000),
#'                  sample(ncol(sce), 200)]
#'   counts(sce_red) <- as.matrix(counts(sce_red))
#'
#'   # Residual Based Variance Stabilizing Transformation
#'   rq <- residual_transform(sce_red, residual_type = "randomized_quantile",
#'                            verbose = TRUE)
#'   pearson <- residual_transform(sce_red, residual_type = "pearson", verbose = TRUE)
#'
#'   # Plot first two principal components
#'   pearson_pca <- prcomp(t(pearson), rank. = 2)
#'   rq_pca <- prcomp(t(rq), rank. = 2)
#'   plot(rq_pca$x, asp = 1)
#'   points(pearson_pca$x, col = "red")
#'
#' @export
residual_transform <- function(data,
                               residual_type = c("randomized_quantile", "pearson", "analytic_pearson"),
                               clipping = FALSE,
                               overdispersion = 0.05,
                               size_factors = TRUE,
                               offset_model = TRUE,
                               overdispersion_shrinkage = TRUE,
                               ridge_penalty = 2,
                               on_disk = NULL,
                               return_fit = FALSE,
                               verbose = FALSE, ...){

  # Allow any valid argument from glmGamPoi::residual.glmGamPoi()
  residual_type <- match.arg(residual_type[1], c("deviance", "pearson", "randomized_quantile",
                                                 "working", "response", "quantile", "analytic_pearson"))

  if(residual_type == "analytic_pearson"){
    return(analytic_pearson_residual_transform(data = data, clipping = clipping, overdispersion = overdispersion, size_factors = size_factors, on_disk = on_disk, return_fit = return_fit, verbose = verbose))
  }

  if(inherits(data, "glmGamPoi")){
    fit <- data
  }else if(offset_model){
    counts <- .handle_data_parameter(data, on_disk, allow_sparse = FALSE  )
    size_factors <- .handle_size_factors(size_factors, counts)

    fit <- glmGamPoi::glm_gp(counts, design = ~ 1, size_factors = size_factors,
                             overdispersion = overdispersion,
                             overdispersion_shrinkage = overdispersion_shrinkage,
                             verbose = verbose, ...)
  }else{
    counts <- .handle_data_parameter(data, on_disk, allow_sparse = FALSE  )
    size_factors <- .handle_size_factors(size_factors, counts)

    log_sf <- log(size_factors)
    attr(ridge_penalty, "target") <- c(0, 1)

    fit <- glmGamPoi::glm_gp(counts, design = ~ log_sf + 1, size_factors = 1,
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

  resid <- stats::residuals(fit, type = residual_type)
  resid <- clip_residuals(resid, clipping)
  resid <- .convert_to_output(resid, data)

  if(! return_fit){
    resid
  }else{
    list(Residuals = resid, fit = fit)
  }
}

# Original implementation in scanpy by Jan Lause
# https://github.com/scverse/scanpy/blob/bd06cc3d1e0bd990f6994e54414512fa0b25fea0/scanpy/experimental/pp/_normalization.py
# Translated to R by Constantin Ahlmann-Eltze
analytic_pearson_residual_transform <- function(data,
                                                clipping = FALSE,
                                                overdispersion = 0.05,
                                                size_factors = TRUE,
                                                on_disk = NULL,
                                                return_fit = FALSE,
                                                verbose = FALSE){
  if(isFALSE(overdispersion)){
    overdispersion <- 0
  }
  if(! is.numeric(overdispersion)){
    stop("For 'analytic_pearson' residuals, the overdispersion must be a non-negative double.")
  }

  counts <- .handle_data_parameter(data, on_disk, allow_sparse = TRUE)
  size_factors <- .handle_size_factors(size_factors, counts)

  make_offset_hdf5_mat <- is(counts, "DelayedMatrix") && is(DelayedArray::seed(counts), "HDF5ArraySeed")

  sum_genes <- MatrixGenerics::rowSums2(counts)
  size_factors <- size_factors / sum(size_factors)

  Mu <- if(make_offset_hdf5_mat){
    delayed_matrix_multiply(DelayedArray::DelayedArray(matrix(sum_genes, ncol = 1)), DelayedArray::DelayedArray(matrix(size_factors, nrow = 1)))
  }else{
    tcrossprod(sum_genes, size_factors)
  }
  resid <- (counts - Mu) / sqrt(Mu + Mu^2 * overdispersion)
  resid <- clip_residuals(resid, clipping)
  resid <- .convert_to_output(resid, data)

  if(! return_fit){
    resid
  }else{
    list(Residuals = resid, fit = NULL)
  }
}



clip_residuals <- function(resid, clipping){
  if(isFALSE(clipping)){
    # Do nothing
  }else{
    if(isTRUE(clipping)){
      clipping <- sqrt(ncol(resid))
    }
    if(! is.numeric(clipping) || length(clipping) != 1){
      stop("Clipping has to be 'TRUE'/'FALSE' or a single numeric value.")
    }
    if(clipping < 0){
      stop("'clipping = ", clipping, "' is negative. Only positive values are allowed.")
    }
    resid[resid > clipping] <- clipping
    resid[resid < -clipping] <- -clipping
  }
  resid
}


