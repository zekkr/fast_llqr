#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// Even more efficient version
// [[Rcpp::export]]
arma::rowvec fast_weighted_colsums_v2(arma::mat A, arma::vec w) {
  arma::rowvec result(A.n_cols, arma::fill::zeros);
  
  // Calculate weighted sums directly
  for (int j = 0; j < A.n_cols; j++) {
    for (int i = 0; i < A.n_rows; i++) {
      result(j) += A(i, j) * w(i);
    }
  }
  
  return result;
}

// [[Rcpp::export]]
arma::vec fast_mat_mult(arma::mat A, arma::vec x) {
  if (A.n_cols != x.n_elem) {
    stop("Dimension mismatch: A.n_cols != x.length()");
  }
  arma::vec result(A.n_rows, arma::fill::zeros);
  for (int i = 0; i < A.n_rows; i++) {
    for (int j = 0; j < A.n_cols; j++) {
      result(i) += A(i, j) * x(j);
    }
  }
  return result;
}