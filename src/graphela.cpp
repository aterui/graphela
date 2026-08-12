#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

// [[Rcpp::export]]
double energy(
    const arma::rowvec& state,
    const arma::rowvec& alpha,
    const arma::mat& beta
) {
  return -(arma::dot(state, alpha) + arma::dot(state, state * beta) * 0.5);
}

// [[Rcpp::export]]
arma::rowvec stpd(
    const arma::rowvec& state,
    const arma::rowvec& alpha,
    const arma::mat& beta,
    bool print = false
) {
  // ---- declare ----
  // {scalar}
  // e: current energy
  // emin: minimum energy after one possible flip
  // idx: index of the species that gives the largest energy decrease
  double e = energy(state, alpha, beta);
  double emin;
  arma::uword idx;

  // {vector}
  // eflip: energy after flipping each species
  // sign: sign vector associated with each possible flip
  // s: current state
  arma::rowvec eflip, sign;
  arma::rowvec s = state;

  // ---- steepest descent ----
  while (true) {
    sign = (1 - 2 * s);
    eflip = ((-alpha - s * beta) % sign) + e;
    emin = eflip.min();

    if (print)
      Rcpp::Rcout << e << " -> " << emin << "\n";

    // stop if no flip decreases energy
    if (emin >= e)
      break;

    // identify the species that gives the largest energy decrease
    idx = eflip.index_min();

    // flip 0/1
    s(idx) = 1 - s(idx);
    e = emin;
  }

  return arma::join_rows(s, arma::rowvec({e}));
}

// [[Rcpp::export]]
arma::mat rss(
    arma::rowvec alpha,
    arma::mat beta,
    const arma::uword n = 20000
) {
  // ---- declare ----
  // m: matrix of stable states (each row represents random initial state)
  // s: temporary state
  // ss: temporary stable state
  arma::mat m = arma::zeros(n, beta.n_cols + 1);
  arma::rowvec s;
  arma::rowvec ss;

  // ---- steepest descent ----
  for (arma::uword k = 0; k < n; ++k) {

    s = arma::randi<arma::rowvec>(
      1, beta.n_cols,
      arma::distr_param(0, 1)
    );

    ss = stpd(s, alpha, beta);
    m.row(k) = ss;
  }

  return m;
}

arma::rowvec findpath_cpp(
    const arma::rowvec& s0,
    const arma::rowvec& s1,
    const arma::rowvec& alpha,
    const arma::mat& beta,
    double temp = 1,
    const double r = 0.01,
    const arma::uword n = 10000,
    arma::mat* mse_out = nullptr,
    arma::vec* omega_out = nullptr
) {
  // ---- set path ----
  // flip: indices of species flipped (index)
  // nf: number of flips, or steps from s0 to s1
  // ns: number of species
  // path: initial path sequence from s0 to s1, shuffled
  arma::uvec flip = arma::find(abs(s0 - s1) == 1);
  arma::uvec path = arma::shuffle(flip);
  const arma::uword nf = flip.n_elem;
  const arma::uword ns = alpha.n_elem;

  // ---- declare ----
  // {scalar}
  // de: energy increment by single flipping
  // etop: highest energy for the path
  // pr: acceptance probability
  // k: index for flipping
  double de, etop, pr;
  arma::uword k;

  // {vector}
  // idx: index for swapping
  // e0: energy vector
  // omega: lowest ridge energy, dynamic updates
  // stip: state vector of a tipping point
  arma::uvec idx(2);
  arma::vec e0(nf + 1);
  arma::vec omega(n);
  arma::rowvec stip(ns + 1);

  // {matrix}
  // ms0: initialized matrix for state vectors
  arma::mat ms0(nf + 1, ns);

  // initialize
  ms0.row(0) = s0;
  e0(0) = energy(s0, alpha, beta);

  // ---- initial path ----
  // temporary intermediate objects
  // s: state vector, dynamic updates
  // u: temporary path vector, dynamic updates
  // ms: matrix for state sequence, dynamic updates
  // e: current energy, dynamic updates
  arma::rowvec s = s0;
  arma::uvec u = path;
  arma::mat ms = ms0;
  arma::vec e = e0;

  // state sequence from s0 to s1
  for (arma::uword i = 0; i < nf; ++i) {
    k = u(i);

    // flip one species, update state
    ms.row(i + 1) = ms.row(i);
    ms(i + 1, k) = 1 - ms(i + 1, k);

    // update energy
    de = -(1 - 2 * s(k)) * (alpha(k) + arma::dot(beta.row(k), s));
    e(i + 1) = de + e(i);
    s(k) = 1 - s(k);
  }

  arma::mat mse = arma::join_rows(ms, e);
  stip = mse.row(e.index_max());
  omega(0) = e.max();

  // ---- simulated annealing ----
  for (arma::uword t = 1; t < n; ++t) {

    // idx: indices for swap
    u = path;
    idx = arma::randperm(u.n_elem, 2);

    // update path sequence by swapping indices
    u.swap_rows(idx(0), idx(1));

    // reset state and energy
    s = s0;
    e = e0;

    // state sequence from s0 to s1
    for (arma::uword i = 0; i < nf; ++i) {
      k = u(i);
      de = -(1 - 2 * s(k)) * (alpha(k) + arma::dot(beta.row(k), s));
      e(i + 1) = de + e(i);
      s(k) = 1 - s(k);
    }

    // record barrier energy for the current path
    etop = e.max();

    // define temperature and acceptance formula
    pr = std::min(
      1.0,
      std::exp((omega(t - 1) - etop) / temp)
    );

    // update temperature
    temp *= (1 - r);

    // update if the new value is accepted
    if (arma::randu<double>() < pr) {
      // update path and barrier
      path = u;
      omega(t) = etop;

      // reset matrix
      ms = ms0;

      // re-construct matrix
      for (arma::uword i = 0; i < nf; ++i) {
        // flip one species, update state
        k = path(i);
        ms.row(i + 1) = ms.row(i);
        ms(i + 1, k) = 1 - ms(i + 1, k);
      }

      mse = arma::join_rows(ms, e);
      stip = mse.row(e.index_max());
    } else {
      omega(t) = omega(t - 1);
    }
  }

  // output
  if (mse_out != nullptr)
    *mse_out = mse;

  if (omega_out != nullptr)
    *omega_out = omega;

  return stip;
}

// [[Rcpp::export]]
Rcpp::List findpath(
    const arma::rowvec& s0,
    const arma::rowvec& s1,
    const arma::rowvec& alpha,
    const arma::mat& beta,
    double temp = 1,
    const double r = 0.01,
    const arma::uword n = 10000
) {
  arma::mat mse;
  arma::vec omega;

  arma::rowvec stip = findpath_cpp(
    s0, s1, alpha, beta,
    temp, r, n,
    &mse, &omega
  );

  return Rcpp::List::create(
    Rcpp::Named("state") = stip,
    Rcpp::Named("path") = mse,
    Rcpp::Named("omega") = omega
  );
}

// [[Rcpp::export]]
arma::mat ridge(
    const arma::mat& sse,
    const arma::rowvec& alpha,
    const arma::mat& beta,
    double temp = 1,
    const double r = 0.01,
    const arma::uword n = 10000
) {
  // ---- declare ----
  // {scalar}
  // ess0, ess1: stable state energy
  // etip: energy at a tipping point
  // cost: cumulative energy costs
  // ed: exp(energy[i+1] - energy[i])
  // ns: number of species
  // nss: number of stable states
  // nr: number of combinations
  // k: index
  double ess0, ess1, etip, cost, ed;
  const arma::uword ns = alpha.n_elem;
  const arma::uword nss = sse.n_rows;
  const arma::uword nr = nss * (nss - 1) / 2;
  arma::uword k = 0;

  if (sse.n_cols != ns + 1)
    Rcpp::stop("sse must contain ns species columns plus one energy column.");

  // {vector}
  // stip: state vector of a tipping point
  // e: vector of stable state energy values
  // epath: vector of path energy values
  arma::rowvec stip;
  arma::vec e = sse.col(sse.n_cols - 1);
  arma::vec epath;

  // {matrix}
  // ss: matrix of stable states
  // combn: output matrix
  // mse: matrix for energy path
  arma::mat ss = sse.cols(0, ns - 1);
  arma::mat combn(nr, 7);
  arma::mat mse;

  for (arma::uword i = 0; i < nss - 1; ++i) {
    for (arma::uword j = i + 1; j < nss; ++j) {

      // calculate state i -> state j
      stip = findpath_cpp(
        ss.row(i), ss.row(j),
        alpha, beta,
        temp, r, n,
        &mse
      );

      ess0 = e(i);
      ess1 = e(j);

      etip = stip.back();
      epath = mse.col(mse.n_cols - 1);

      cost = 0;

      for (arma::uword m = 0; m < epath.n_elem - 1; ++m) {
        ed = epath(m + 1) - epath(m);
        cost += std::exp(ed) - 1.0;
      }

      if (ess0 > ess1) {
        combn(k, 0) = i + 1;  // higher-energy state
        combn(k, 1) = j + 1;  // lower-energy state
      } else {
        combn(k, 0) = j + 1;  // higher-energy state
        combn(k, 1) = i + 1;  // lower-energy state
      }

      combn(k, 2) = std::max(ess0, ess1); // ss energy higher
      combn(k, 3) = std::min(ess0, ess1); // ss energy lower
      combn(k, 4) = etip; // tipping point
      combn(k, 5) = cost; // cumulative energy cost
      combn(k, 6) = etip - combn(k, 2);  // energy barrier

      ++k;
    }
  }

  return combn;
}

// [[Rcpp::export]]
arma::uword find_shallow(
    const arma::mat& pem
) {

  // ---- key column indices ----
  // idx_depth: column index for basin depth
  // idx_cost: column index for path energy cost
  arma::uword idx_depth = pem.n_cols - 1;
  arma::uword idx_cost = pem.n_cols - 2;

  // indices for candidate rows
  arma::uvec idx_sh =
    arma::find(pem.col(idx_depth) == pem.col(idx_depth).min());

  // select the minimum path cost amongst the candidates
  arma::mat candid = pem.rows(idx_sh);
  arma::uword r = idx_sh(candid.col(idx_cost).index_min());

  return r;
}

// [[Rcpp::export]]
Rcpp::List prune(
    arma::mat& pem,
    const double th = 0.2
) {
  // ---- declare ----
  // idx_depth: column index for basin depth
  // idx_cost: column index for path energy cost
  // idx_rm: index for basin to be removed
  // dmax, dmin: basin depth max, min
  arma::mat map(pem.n_rows, 2);
  arma::uword idx_depth = pem.n_cols - 1;
  arma::uword idx_rm;
  arma::uword k = 0;
  double dmax, dmin;

  while (true) {
    // find deepest & shallowest
    arma::uword r = find_shallow(pem);
    dmax = pem.col(idx_depth).max();
    dmin = pem(r, idx_depth);

    if (dmin >= th * dmax)
      break;

    // basin to be removed
    idx_rm = static_cast<arma::uword>(pem(r, 0));

    // retain pairs not involving shallow basin
    arma::uvec keep =
      arma::find(
        (pem.col(0) != idx_rm) && (pem.col(1) != idx_rm)
      );

    // basin mapping for merge
    map(k, 0) = pem(r, 0); // shallower
    map(k, 1) = pem(r, 1); // deeper
    pem = pem.rows(keep);

    k++;
  }

  if (k == 0) {
    return Rcpp::List::create(
      Rcpp::Named("pem") = pem,
      Rcpp::Named("map") = R_NilValue
    );
  }

  return Rcpp::List::create(
    Rcpp::Named("pem") = pem,
    Rcpp::Named("map") = map.resize(k, 2)
  );
}

