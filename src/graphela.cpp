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

// [[Rcpp::export]]
arma::mat ridge(
    const arma::rowvec s0,
    const arma::rowvec s1,
    const arma::rowvec alpha,
    const arma::mat beta,
    double temp = 1,
    const double r = 0.01,
    const int n = 10000
) {

  // {set path}
  // flip: indices of species flipped (index)
  // nf: number of flips, or steps from s0 to s1
  // ns: number of species
  // path: initial path sequence from s0 to s1, shuffled
  uvec flip = find(abs(s0 - s1) == 1);
  const int nf = flip.n_elem;
  const int ns = alpha.n_elem;
  uvec path = shuffle(flip);
  uvec u = path;

  // {declare}
  // idx: index for swapping
  // ms0: initialized matrix
  // ms: matrix for intermediate states
  // e: energy
  // omega:
  // etop:
  // stip: state vector of a tipping point

  // scalar
  double etop = datum::inf;
  double pr;

  // vectors
  uvec idx(2);
  vec e0(nf + 1);
  vec omega(n + 1);
  rowvec stip(ns + 1);

  // matrices
  mat ms0(nf + 1, ns);

  // initialize
  ms0.row(0) = s0;
  e0(0) = energy(s0, alpha, beta);

  // ---- initial path ----
  mat ms = ms0;
  vec e = e0;

  // state sequence from s0 to s1
  for (int i = 0; i < nf; ++i) {
    // flip one species, update state
    ms.row(i + 1) = ms.row(i);
    ms(i + 1, u(i)) = 1 - ms(i + 1, u(i));

    // energy of current state
    e(i + 1) = energy(ms.row(i + 1), alpha, beta);
  }

  omega(0) = e.max();

  // ---- simulated annealing ----
  for (int t = 1; t < n; ++t) {
    // idx: indices for swap
    u = path;
    idx = randperm(u.n_elem, 2);

    // update path sequence by swapping indices
    u.swap_rows(idx(0), idx(1));

    // reset the initial state
    ms = ms0;

    // reset energy
    e = e0;

    // state sequence from s0 to s1
    for (int i = 0; i < nf; ++i) {
      // flip one species, update state
      ms.row(i + 1) = ms.row(i);
      ms(i + 1, u(i)) = 1 - ms(i + 1, u(i));

      // energy of current state
      e(i + 1) = energy(ms.row(i + 1), alpha, beta);
    }

    etop = e.max();

    // define temperature and acceptance formula
    pr = std::min(1.0, std::exp((omega(t) - etop) / temp));
    temp *= (1 - r);
    std::cout << temp << std::endl;

    if (randu<double>() < pr) {
      path = u;
      omega(t + 1) = etop;
      stip.cols(0, ns - 1) = ms.row(e.index_max());
      stip.col(ns) = etop;
    } else {
      omega(t + 1) = omega(t);
    }
  }

  return omega;
}

