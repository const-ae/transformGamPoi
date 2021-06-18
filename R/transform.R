



variance_stabilization <- function(data, overdispersion = 0.05, method = c("acosh", "logp")){



}


.acosh_trans_impl <- function(x, alpha){
  1/sqrt(alpha) * acoshp1(2 * alpha * x)
}

.sqrt_trans_impl <- function(x){
  2 * sqrt(x)
}

acosh_transform <- function(data, overdispersion = 0.05, size_factors = TRUE, on_disk = NULL, verbose = FALSE){

  counts <- .handle_data_parameter(data, on_disk, allow_sparse = TRUE)

  if(inherits(data, "glmGamPoi")){
    size_factors <- data$size_factors
  }else{
    size_factors <- .handle_size_factors(size_factors, counts)
  }

  if(length(overdispersion) == 1 && overdispersion == "fitted"){
    fit <- glmGamPoi::glm_gp(counts, design = ~ 1, size_factors = size_factors,
                             overdispersion = TRUE,
                             overdispersion_shrinkage = FALSE,
                             verbose = verbose)
    overdispersion <- fit$overdispersions
  }else{
    overdispersion <- .handle_overdispersion(overdispersion, counts)
  }

  norm_counts <- DelayedArray::sweep(counts, 2, size_factors, FUN = "/")

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

shifted_log_transform <- function(data, overdispersion = 0.05, pseudo_count = 1/(4 * overdispersion),
                                  size_factors = TRUE, minimum_overdispersion = 0.001, on_disk = NULL, verbose = FALSE){

  counts <- .handle_data_parameter(data, on_disk, allow_sparse = TRUE)

  if(inherits(data, "glmGamPoi")){
    size_factors <- data$size_factors
  }else{
    size_factors <- .handle_size_factors(size_factors, counts)
  }

  if(length(overdispersion) == 1 && overdispersion == "fitted"){
    fit <- glmGamPoi::glm_gp(counts, design = ~ 1, size_factors = size_factors,
                             overdispersion = TRUE,
                             overdispersion_shrinkage = FALSE,
                             verbose = verbose)
    data <- fit
    overdispersion <- fit$overdispersions
  }else{
    overdispersion <- 1/(4 * pseudo_count)
    overdispersion <- .handle_overdispersion(overdispersion, counts)
  }




  norm_counts <- DelayedArray::sweep(counts, 2, size_factors, FUN = "/")

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
