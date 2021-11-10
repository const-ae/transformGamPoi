




.acosh_trans_impl <- function(x, alpha){
  1/sqrt(alpha) * acoshp1(2 * alpha * x)
}

.sqrt_trans_impl <- function(x){
  2 * sqrt(x)
}


#' Delta method-based variance stabilizing transformation
#'
#'
#' @inherit transformGamPoi
#' @param pseudo_count instead of specifying the overdispersion, the
#'   `shifted_log_transform` is commonly parameterized with a pseudo-count
#'   (\eqn{pseudo-count = 1/(4 * overdispersion)}). If both the `pseudo-count`
#'   and `overdispersion` is specified, the `overdispersion` is ignored.
#'   Default: `1/(4 * overdispersion)`
#' @param minimum_overdispersion the `acosh_transform` converges against
#'   \eqn{2 * sqrt(x)} for `overdispersion == 0`. However, the `shifted_log_transform`
#'   would just become `0`, thus here we apply the `minimum_overdispersion` to avoid
#'   this behavior.
#' @param ... additional parameters for `glmGamPoi::glm_gp()` which is called in
#'   case `overdispersion = TRUE`.
#'
#' @describeIn  acosh_transform \eqn{1/sqrt(alpha)} acosh(2 * alpha * x + 1)
#'
#' @return a matrix (or a vector if the input is a vector) with the transformed values.
#'
#' @examples
#'   # Load a single cell dataset
#'   sce <- TENxPBMCData::TENxPBMCData("pbmc4k")
#'   # Reduce size for this example
#'   set.seed(1)
#'   sce_red <- sce[sample(which(rowSums2(counts(sce)) > 0), 1000),
#'                  sample(ncol(sce), 200)]
#'
#'   assay(sce_red, "acosh") <- acosh_transform(sce_red)
#'   assay(sce_red, "shifted_log") <- shifted_log_transform(sce_red)
#'   plot(rowMeans2(assay(sce_red, "acosh")), rowVars(assay(sce_red, "acosh")), log = "x")
#'   points(rowMeans2(assay(sce_red, "shifted_log")), rowVars(assay(sce_red, "shifted_log")),
#'          col = "red")
#'
#'   # Sqrt transformation
#'   sqrt_dat <- acosh_transform(sce_red, overdispersion = 0, size_factor = 1)
#'   plot(2 * sqrt(assay(sce_red))[,1], sqrt_dat[,1]); abline(0,1)
#'
#' @export
acosh_transform <- function(data, overdispersion = 0.05,
                            size_factors = TRUE,
                            ...,
                            on_disk = NULL,
                            verbose = FALSE){

  counts <- .handle_data_parameter(data, on_disk, allow_sparse = TRUE)

  if(inherits(data, "glmGamPoi")){
    size_factors <- data$size_factors
  }else{
    size_factors <- .handle_size_factors(size_factors, counts)
  }

  if(all(isTRUE(overdispersion)) || all(overdispersion == "global")){
    if(HDF5Array::is_sparse(counts)){
      counts <- .handle_data_parameter(data, on_disk, allow_sparse = FALSE)
    }
    dots <- list(...)
    overdispersion_shrinkage <- if("overdispersion_shrinkage" %in% names(dots)){
      dots[["overdispersion_shrinkage"]]
    }else{
      TRUE
    }
    fit <- glmGamPoi::glm_gp(counts, design = ~ 1, size_factors = size_factors,
                             overdispersion = overdispersion,
                             overdispersion_shrinkage = overdispersion_shrinkage,
                             verbose = verbose)
    if(overdispersion_shrinkage){
      # Use the dispersion trend when calculating the residuals
      fit$overdispersion_shrinkage_list$original_overdispersions <- fit$overdispersions
      fit$overdispersions <- fit$overdispersion_shrinkage_list$dispersion_trend
    }
    overdispersion <- fit$overdispersions
  }else{
    overdispersion <- .handle_overdispersion(overdispersion, counts)
  }

  if(is(counts, "dgCMatrix")){
    norm_counts <- sparse_divide_out_size_factor(counts, size_factors)
  }else{
    norm_counts <- DelayedArray::sweep(counts, 2, size_factors, FUN = "/")
  }

  overdispersion_near_zero <- .near(overdispersion, 0)

  result <- if(! any(overdispersion_near_zero)){
    # no overdispersion is zero
    .acosh_trans_impl(norm_counts, overdispersion)
  }else if(all(overdispersion_near_zero)){
    # all overdispersion is zero
    .sqrt_trans_impl(norm_counts)
  }else{
    # overdispersion is a mix of zeros and other values.
    if(is.matrix(overdispersion)){
      norm_counts[overdispersion_near_zero] <- .sqrt_trans_impl(norm_counts[overdispersion_near_zero])
      norm_counts[!overdispersion_near_zero] <- .acosh_trans_impl(norm_counts[!overdispersion_near_zero],
                                                                  overdispersion[! overdispersion_near_zero])
      norm_counts
    }else{
      norm_counts[overdispersion_near_zero, ] <- .sqrt_trans_impl(norm_counts[overdispersion_near_zero, ])
      norm_counts[!overdispersion_near_zero, ] <- .acosh_trans_impl(norm_counts[!overdispersion_near_zero, ],
                                                                    overdispersion[! overdispersion_near_zero])
      norm_counts
    }
  }

  .convert_to_output(result, data)
}



.log_plus_alpha_impl <- function(x, alpha){
  1/sqrt(alpha) * log1p(4 * alpha * x)
}


#' @describeIn acosh_transform \eqn{1/sqrt(alpha) log(4 * alpha * x + 1)}
#' @export
shifted_log_transform <- function(data,
                                  overdispersion = 0.05,
                                  pseudo_count = 1/(4 * overdispersion),
                                  size_factors = TRUE,
                                  minimum_overdispersion = 0.001,
                                  ...,
                                  on_disk = NULL,
                                  verbose = FALSE){

  counts <- .handle_data_parameter(data, on_disk, allow_sparse = TRUE)

  if(inherits(data, "glmGamPoi")){
    size_factors <- data$size_factors
  }else{
    size_factors <- .handle_size_factors(size_factors, counts)
  }

  if(all(isTRUE(overdispersion)) || all(overdispersion == "global")){
    if(HDF5Array::is_sparse(counts)){
      counts <- .handle_data_parameter(data, on_disk, allow_sparse = FALSE)
    }
    dots <- list(...)
    overdispersion_shrinkage <- if("overdispersion_shrinkage" %in% names(dots)){
      dots[["overdispersion_shrinkage"]]
    }else{
      TRUE
    }
    fit <- glmGamPoi::glm_gp(counts, design = ~ 1, size_factors = size_factors,
                             overdispersion = overdispersion,
                             overdispersion_shrinkage = overdispersion_shrinkage,
                             verbose = verbose)
    if(overdispersion_shrinkage){
      fit$overdispersion_shrinkage_list$original_overdispersions <- fit$overdispersions
      fit$overdispersions <- fit$overdispersion_shrinkage_list$dispersion_trend
    }
    overdispersion <- fit$overdispersions
  }else{
    overdispersion <- 1/(4 * pseudo_count)
    overdispersion <- .handle_overdispersion(overdispersion, counts)
  }

  if(is(counts, "dgCMatrix")){
    norm_counts <- sparse_divide_out_size_factor(counts, size_factors)
  }else{
    norm_counts <- DelayedArray::sweep(counts, 2, size_factors, FUN = "/")
  }


  overdispersion[overdispersion < minimum_overdispersion] <- minimum_overdispersion

  result <- .log_plus_alpha_impl(norm_counts, overdispersion)
  .convert_to_output(result, data)
}




.convert_to_output <- function(result, data){
  if(is.vector(data)){
    as.vector(result)
  }else{
    result
  }
}


.near <- function (x, y, tol = .Machine$double.eps^0.5){
  abs(x - y) < tol
}
