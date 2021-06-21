
#' Take any kind of data and extract the matrix
#'
#'
#' Adapted from glmGamPoi:::handle_data_parameter
#'
#'
#' @keywords internal
.handle_data_parameter <- function(data, on_disk, allow_sparse = TRUE){
  if(is.vector(data)){
    data <- matrix(data, nrow = 1)
  }
  if(is.matrix(data)){
    if(! is.numeric(data)){
      stop("The data argument must consist of numeric values and not of ", mode(data), " values")
    }
    if(is.null(on_disk) || isFALSE(on_disk)){
      data_mat <- data
    }else if(isTRUE(on_disk)){
      data_mat <- HDF5Array::writeHDF5Array(data)
    }else{
      stop("Illegal argument type for on_disk. Can only handle 'NULL', 'TRUE', or 'FALSE'")
    }
  }else if(is(data, "DelayedArray")){
    if(is.null(on_disk) || isTRUE(on_disk)){
      data_mat <- data
    }else if(isFALSE(on_disk)){
      data_mat <- as.matrix(data)
    }else{
      stop("Illegal argument type for on_disk. Can only handle 'NULL', 'TRUE', or 'FALSE'")
    }
  }else if(is(data, "SummarizedExperiment")){
    data_mat <- .handle_data_parameter(SummarizedExperiment::assay(data, "counts"), on_disk, allow_sparse)
  }else if(canCoerce(data, "SummarizedExperiment")){
    se <- as(data, "SummarizedExperiment")
    data_mat <- .handle_data_parameter(SummarizedExperiment::assay(se, "counts"), on_disk, allow_sparse)
  }else if(is(data, "dgCMatrix") || is(data, "dgTMatrix")) {
    if(isTRUE(on_disk)){
      data_mat <- HDF5Array::writeHDF5Array(data)
    }else if(isFALSE(on_disk)){
      if(allow_sparse){
        data_mat <- data
      }else{
        data_mat <- as.matrix(data)
      }
    }else{
      if(allow_sparse){
        data_mat <- data
      }else{
        stop("transformGamPoi does not yet support sparse input data of type '", class(data),"'. ",
             "Please explicitly set the 'on_disk' parameter to force a conversion to a dense format either in-memory ('on_disk = FALSE') ",
             "or on-disk ('on_disk = TRUE')")
      }
    }
  }else if(inherits(data, "glmGamPoi")){
    data_mat <- .handle_data_parameter(data$data, on_disk, allow_sparse)
  }else{
    stop("Cannot handle data of class '", class(data), "'.",
         "It must be of a matrix object (i.e., a base matrix or DelayedArray),",
         " or a container for such a matrix (for example: SummarizedExperiment).")
  }
  data_mat
}




.handle_size_factors <- function(size_factors, Y, verbose = FALSE){
  n_samples <- ncol(Y)

  if(isTRUE(size_factors) || is.character(size_factors)){
    method <- if(isTRUE(size_factors)){
      "normed_sum"
    }else if(all(size_factors == c("normed_sum", "deconvolution", "poscounts"))){
      "normed_sum"
    }else if(length(size_factors) == 1 && size_factors %in% c("normed_sum", "deconvolution", "poscounts")){
      size_factors
    }else{
      stop("Cannot handle size factor ", paste0(size_factors, collapse = ", "))
    }
    if(verbose){ message("Calculate Size Factors (", method, ")") }
    sf <- estimate_size_factors(Y, method = method, verbose = verbose)
  }else if(isFALSE(size_factors)){
    sf <- rep(1, n_samples)
  }else{
    stopifnot(is.numeric(size_factors) && (length(size_factors) == 1 || length(size_factors) == n_samples))
    if(any(size_factors < 0)){
      stop("size factor 'size_factors' must be larger than 0")
    }
    if(length(size_factors) == 1){
      sf <- rep(size_factors, n_samples)
    }else{
      sf <- size_factors
    }
  }
  sf
}







#' Estimate the Size Factors
#'
#' @param Y any matrix-like object (\code{base::matrix()}, \code{DelayedArray}, \code{HDF5Matrix},
#'   \code{Matrix::Matrix()}, etc.)
#' @param method one of \code{c("normed_sum", "deconvolution", "poscounts")}
#'
#'
#' @return a vector with one size factor per column of `Y`
#'
#' @keywords internal
estimate_size_factors <- function(Y, method, verbose = FALSE){
  if(nrow(Y) <= 1){
    if(verbose) {
      message("Skipping size factor estimation, because there is only a single gene.",
              call. = FALSE)
    }
    return(rep(1, ncol(Y)))
  }
  stopifnot(length(method) == 1 && is.character(method))

  if(method == "poscounts"){
    # Accept any matrix-like object
    log_geometric_means <- MatrixGenerics::rowMeans2(log(Y + 0.5))
    Y2 <- Y
    Y2[Y2 == 0] <- NA
    sf <- exp(MatrixGenerics::colMedians(subtract_vector_from_each_column(log(Y2), log_geometric_means), na.rm = TRUE))
  }else if(method == "deconvolution"){
    if(requireNamespace("scran", quietly = TRUE)){
      tryCatch({
        sf <- scran::calculateSumFactors(Y)
      }, error = function(err){
        stop("Error in size factor estimation with 'size_factors = \"deconvolution\"'.\n",
             "'scran::calculateSumFactors(Y)' threw the following error: \n\t",
             err, call. = FALSE)
      })
    }else{
      stop("To use the \"deconvolution\" method for size factor calculation, you need to install the ",
           "'scran' package from Bioconductor. Alternatively, you can use \"normed_sum\" or \"poscounts\"",
           "to calculate the size factors.", call. = FALSE)
    }
  }else if(method == "normed_sum"){
    sf <- colSums2(Y)
    # sf <- matrixStats::colSums2(as.matrix(Y))
  }else{
    stop("Unknown size factor estimation method: ", method)
  }



  # stabilize size factors to have geometric mean of 1
  all_zero_column <- is.nan(sf) | sf <= 0
  sf[all_zero_column] <- NA
  if(any(all_zero_column)){
    sf <- sf/exp(mean(log(sf), na.rm=TRUE))
    sf[all_zero_column] <- 0.001
  }else{
    sf <- sf/exp(mean(log(sf)))
  }
  sf
}


.handle_overdispersion <- function(overdispersion, counts, verbose = FALSE){
  stopifnot(length(overdispersion) == 1 || length(overdispersion) == nrow(counts) || all(dim(overdispersion) == dim(counts)))
  if(length(overdispersion) == 1){
    rep(overdispersion, nrow(counts))
  }else{
    overdispersion
  }

  if(verbose && any(overdispersion > 1)){
    warning("The overdispersion at position", paste0(head(which(overdispersion > 1)), collapse = ", "), " seems unusually large.")
  }
  overdispersion
}


subtract_vector_from_each_column <- function(matrix, vector){
  stopifnot(length(vector) == 1 || length(vector) == nrow(matrix))
  matrix - vector
}
