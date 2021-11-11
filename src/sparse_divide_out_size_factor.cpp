#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector sparse_divide_out_size_factor_impl(const NumericVector& x, const IntegerVector& p, const NumericVector& s){
  const int n_elem = x.size();
  NumericVector res(n_elem);
  int col_idx = 0;
  int pointer_loc = p[1];
  for(int idx = 0; idx < n_elem; ++idx){
    while(idx >= pointer_loc){
      pointer_loc = p[++col_idx + 1];
    }
    res[idx] = x[idx] / s[col_idx];
  }
  return res;
}

