

sparse_divide_out_size_factor <- function(sp_mat, sf){
  stopifnot(is(sp_mat, "dgCMatrix"))

  if(length(sf) == 1){
    sp_mat@x <- sp_mat@x / sf
  }else if(length(sf) == ncol(sp_mat)){
    for(col_idx in seq_along(sp_mat@p[-1])){
      if(sp_mat@p[col_idx] < sp_mat@p[col_idx + 1]){
        sel <- (sp_mat@p[col_idx]+1):(sp_mat@p[col_idx+1])
        sp_mat@x[sel] <- sp_mat@x[sel] / sf[col_idx]
      }
    }
  }else{
    stop("Length of sf does not match the number of columns in sp_mat")
  }
  sp_mat
}
