
#'  Variance Stabilizing Transformation for Gamma Poisson Data
#'
#' @param data any matrix-like object (e.g. matrix, dgCMatrix, DelayedArray, HDF5Matrix)
#'   with one column per sample and row per gene. It can also be an object of type `glmGamPoi`,
#'   in which case it is directly used to calculate the variance-stabilized values.
#' @param transformation one of `c("acosh", "shifted_log", "randomized_quantile_residuals", "pearson_residuals")`.
#'   See [`acosh_transform`], [`shifted_log_transform`], or [`residual_transform`] for more information.
#' @param overdispersion the simplest count model is the Poisson model. However, the Poisson model
#'   assumes that \eqn{variance = mean}. For many applications this is too rigid and the Gamma-Poisson
#'   allows a more flexible mean-variance relation (\eqn{variance = mean + mean^2 * overdispersion}). \cr
#'   `overdispersion` can either be
#'   \itemize{
#'      \item a single boolean that indicates if an overdispersion is estimated for each gene.
#'      \item a numeric vector of length `nrow(data)` fixing the overdispersion to those values.
#'      \item the string `"global"` to indicate that one dispersion is fit across all genes.
#'   }
#'   Note that `overdispersion = 0` and `overdispersion = FALSE` are equivalent and both reduce
#'   the Gamma-Poisson to the classical Poisson model. Default: `0.05` which is roughly the
#'   overdispersion observed on ostensibly homogeneous cell lines.
#' @param verbose boolean that decides if information about the individual steps are printed.
#'   Default: `FALSE`
#' @param ... additional parameters passed to [`acosh_transform`], [`shifted_log_transform`], or [`residual_transform`]
#' @inheritParams glmGamPoi::glm_gp
#'
#' @return a matrix (or a vector if the input is a vector) with the transformed values.
#'
#' @seealso [`acosh_transform`], [`shifted_log_transform`], and [`residual_transform`]
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
#'   # Load a single cell dataset
#'   sce <- TENxPBMCData::TENxPBMCData("pbmc4k")
#'   # Reduce size for this example
#'   set.seed(1)
#'   sce_red <- sce[sample(which(rowSums2(counts(sce)) > 0), 1000),
#'                  sample(ncol(sce), 200)]
#'
#'   assay(sce_red, "acosh") <- transformGamPoi(sce_red, "acosh")
#'   assay(sce_red, "shifted_log") <- transformGamPoi(sce_red, "shifted_log")
#'
#'   # Residual Based Variance Stabilizing Transformation
#'   rq <- transformGamPoi(sce_red, transformation = "randomized_quantile", on_disk = FALSE,
#'                         verbose = TRUE)
#'   pearson <- transformGamPoi(sce_red, transformation = "pearson", on_disk = FALSE, verbose = TRUE)
#'
#'   plot(rowMeans2(counts(sce_red)), rowVars(assay(sce_red, "acosh")), log = "x")
#'   points(rowMeans2(counts(sce_red)), rowVars(assay(sce_red, "shifted_log")), col = "red")
#'   points(rowMeans2(counts(sce_red)), rowVars(rq), col = "blue")
#'
#'
#'   # Plot first two principal components
#'   acosh_pca <- prcomp(t(assay(sce_red, "acosh")), rank. = 2)
#'   rq_pca <- prcomp(t(rq), rank. = 2)
#'   pearson_pca <- prcomp(t(pearson), rank. = 2)
#'
#'   plot(acosh_pca$x, asp = 1)
#'   points(rq_pca$x, col = "blue")
#'   points(pearson_pca$x, col = "green")
#'
#' @export
transformGamPoi <- function(data,
                            transformation = c("acosh", "shifted_log", "randomized_quantile_residuals", "pearson_residuals"),
                            overdispersion = 0.05, size_factors = TRUE, ..., on_disk = NULL, verbose = FALSE){

  transformation <- match.arg(transformation)

  dots <- list(...)
  if(! is.null(dots$residual_type)){
    stop("You specified 'residual_type = \"", dots$residual_type, "\"'. However, ",
         "for the 'transformGamPoi()' function rather specify 'transformation = \"",
         "randomized_quantile_residuals\"' or 'transformation = \"pearson_residuals\"'.")
  }

  if(transformation == "acosh"){
    acosh_transform(data, overdispersion = overdispersion, size_factors = size_factors,
                    on_disk = on_disk, verbose = verbose)
  }else if(transformation == "shifted_log"){
    shifted_log_transform(data, overdispersion = overdispersion, size_factors = size_factors,
                    on_disk = on_disk, verbose = verbose, ...)
  }else if(transformation == "randomized_quantile_residuals"){
    residual_transform(data, residual_type = "randomized_quantile",
                       overdispersion = overdispersion, size_factors = size_factors,
                       on_disk = on_disk, verbose = verbose, ...)
  }else if(transformation == "pearson_residuals"){
    residual_transform(data, residual_type = "pearson", overdispersion = overdispersion, size_factors = size_factors,
                       on_disk = on_disk, verbose = verbose, ...)
  }else{
    stop("Unsupported transformation of type: ", transformation, ". The available options are: \n",
         paste0(c("acosh", "shifted_log", "randomized_quantile_residuals", "pearson_residuals"), collapse = ", "))
  }

}

