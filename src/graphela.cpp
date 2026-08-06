#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
double energy(
    const arma::rowvec& state,
    const arma::rowvec& alpha,
    const arma::mat& beta
) {
  mat res = -state * alpha.t() - (state * (state * beta).t() ) / 2;
  return as_scalar(res);
}
