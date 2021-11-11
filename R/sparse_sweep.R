

sparse_divide_out_size_factor <- function(sp_mat, sf){
  if(length(sf) == 1){
    sp_mat@x <- sp_mat@x / sf
  }else if(length(sf) == ncol(sp_mat)){
    sp_mat@x <- sparse_divide_out_size_factor_impl(sp_mat@x, sp_mat@p, sf)
  }else{
    stop("Length of sf does not match the number of columns in sp_mat")
  }
  sp_mat
}


