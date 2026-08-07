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
  mat e = - state * alpha.t() - (state * (state * beta).t() ) / 2;
  return as_scalar(e);
}

// [[Rcpp::export]]
arma::rowvec stpd(
    const arma::rowvec& state,
    const arma::rowvec& alpha,
    const arma::mat& beta,
    bool print = false
) {

  // {Declare variables}
  // eflip: vector of energy after a flip of species "i"
  // sign: sign vector for a flip of species "i"
  // emin: candidate energy minimum after flipping
  // idx: species index that reduces the energy most
  rowvec eflip;
  rowvec sign;
  double emin;
  uword idx;

  // {Initialize}
  // s: temporary state vector - subject to updates
  // e: temporary energy scalar - subject to updates
  rowvec s = state;
  double e = energy(state, alpha, beta);

  // {Steepest descent}
  while(true) {
    sign = (1 - 2 * s);
    eflip = ((-alpha - s * beta) % sign) + e;
    emin = eflip.min();

    if (print)
      std::cout << e << " -> " << emin << std::endl;

    // break if no improvement in energy
    if (emin >= e)
      break;

    // identify index to flip
    idx = eflip.index_min();

    // flip 0/1
    s(idx) = 1 - s(idx);
    e = emin;
  }

  return join_rows(s, rowvec({e}));
}

// [[Rcpp::export]]
arma::mat rss(
    arma::rowvec alpha,
    arma::mat beta,
    int n = 20000
) {
  // {Declare}
  // m: matrix of stable states (each row represents random initial state)
  // s: temporary state
  // ss: temporary stable state
  mat m = zeros(n, beta.n_cols + 1);
  rowvec s;
  rowvec ss;

  // {Steepest descent}
  for (int k = 0; k < n; ++k) {
    s = randi<rowvec>(1, beta.n_cols, distr_param(0, 1));
    ss = stpd(s, alpha, beta);
    m.row(k) = ss;
  }

  return m;
}
