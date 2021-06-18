



variance_stabilization <- function(data, overdispersion = 0.05, method = c("acosh", "logp")){



}


.acosh_trans_impl <- function(x, alpha){
  1/sqrt(alpha) * acoshp1(2 * alpha * x)
}

.sqrt_trans_impl <- function(x){
  2 * sqrt(x)
}

acosh_transform <- function(data, overdispersion = 0.05, size_factors = 1, verbose = FALSE){
  counts <- .get_count_matrix(data)
  norm_counts <- sweep(counts, 2, size_factors, FUN = "/")
  stopifnot(length(overdispersion) == 1 || length(overdispersion) == nrow(norm_counts) || all(dim(overdispersion) == dim(norm_counts)))


  if(verbose && any(overdispersion > 1)){
    warning("The overdispersion at position", paste0(head(which(overdispersion > 1)), collapse = ", "), " seems unusually large.")
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

shifted_log_transform <- function(data, overdispersion = 0.05, pseudo_count = 1/(4 * overdispersion), size_factors = 1, minimum_overdispersion = 0.001, verbose = FALSE){
  if((length(overdispersion) != 1 || overdispersion != 0.05) && !all(pseudo_count == 1/(4 * overdispersion))){
    stop("You can only change the default value of 'overdispersion' or 'pseudo_count'")
  }
  overdispersion <- 1/(4 * pseudo_count)

  counts <- .get_count_matrix(data)
  norm_counts <- sweep(counts, 2, size_factors, FUN = "/")
  stopifnot(length(overdispersion) == 1 || length(overdispersion) == nrow(norm_counts) || all(dim(overdispersion) == dim(norm_counts)))


  if(verbose && any(overdispersion > 1)){
    warning("The overdispersion at position", paste0(head(which(overdispersion > 1)), collapse = ", "), " seems unusually large.")
  }

  overdispersion[overdispersion < minimum_overdispersion] <- minimum_overdispersion

  result <- .log_plus_alpha_impl(norm_counts, overdispersion)
  .convert_to_output(result, data)
}




.get_count_matrix <- function(data){
  if(is(data, "SummarizedExperiment")){
    BiocGenerics::counts(data)
  }else if(is.vector(data)){
    matrix(data, nrow = 1)
  }else {
    data
  }
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
