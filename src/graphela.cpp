#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

// [[Rcpp::export]]
double energy_cpp(
    const arma::rowvec& state,
    const arma::rowvec& alpha,
    const arma::mat& beta
) {
  return -(arma::dot(state, alpha) + arma::dot(state, state * beta) * 0.5);
}

// [[Rcpp::export]]
arma::rowvec stpd_cpp(
    const arma::rowvec& state,
    const arma::rowvec& alpha,
    const arma::mat& beta
) {
  // ---- declare ----
  // {scalar}
  // e: current energy
  // emin: minimum energy after one possible flip
  // idx: index of the species that gives the largest energy decrease
  double e = energy_cpp(state, alpha, beta);
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
    sign = (1.0 - 2.0 * s);
    eflip = (-(alpha + s * beta) % sign) + e;
    emin = eflip.min();

    // if (print)
    //   Rcpp::Rcout << e << " -> " << emin << "\n";

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
arma::mat rss_cpp(
    arma::rowvec alpha,
    arma::mat beta,
    const arma::uword n = 10000,
    const bool replace = true
) {
  // ---- declare ----
  // m: matrix of stable states
  // s: temporary initial state
  // ss: stable state reached from the initial state
  // ns: number of species
  // nstate: number of possible binary states when ns < 32
  // nr: number of initial states to evaluate

  const arma::uword ns = alpha.n_elem;
  const arma::uword nstate = (ns < 32) ? (1ULL << ns) : 0;
  const arma::uword nr = (ns < 32 && !replace) ? std::min(nstate, n) : n;

  arma::mat m(nr, ns + 1);
  arma::rowvec s(ns);
  arma::rowvec ss(ns + 1);

  // initial states represented as integers when ns < 32
  arma::uvec init;

  if (ns < 32 && !replace) {

    if (nstate <= n) {
      // enumerate all possible binary states
      init = arma::regspace<arma::uvec>(0, nstate - 1);
    } else {
      // randomly sample n unique binary states
      init = arma::randperm(nstate, n);
    }

  }

  // ---- steepest descent ----
  for (arma::uword k = 0; k < nr; ++k) {

    if (ns < 32 && !replace) {

      // convert integer state to binary state vector
      const arma::uword z = init(k);

      for (arma::uword i = 0; i < ns; ++i)
        s(i) = (z >> i) & 1ULL;

    } else {

      // randomly sample an initial binary state
      s = arma::randi<arma::rowvec>(
        1, ns,
        arma::distr_param(0, 1)
      );
    }

    ss = stpd_cpp(s, alpha, beta);
    m.row(k) = ss;
  }

  return m;
}

arma::rowvec findpath_inline(
    const arma::rowvec& s0,
    const arma::rowvec& s1,
    const arma::rowvec& alpha,
    const arma::mat& beta,
    double temp = 1,
    const double r = 0.01,
    const arma::uword iter = 10000,
    arma::mat* mse_out = nullptr,
    arma::vec* omega_out = nullptr
) {
  // ---- set path ----
  // flip: indices of species that differ between s0 and s1
  // path: shuffled sequence of species to flip
  // nf: number of flips, or steps from s0 to s1
  // ns: number of species
  arma::uvec flip = arma::find(abs(s0 - s1) == 1);
  arma::uvec path = arma::shuffle(flip);
  const arma::uword nf = flip.n_elem;
  const arma::uword ns = alpha.n_elem;

  // ---- declare ----
  // {scalar}
  // de: energy increment from a single flip
  // etop: maximum energy along the candidate path
  // pr: acceptance probability
  // k: index of the species to flip
  double de, etop, pr;
  arma::uword k;

  // {vector}
  // idx: indices of the two positions to swap
  // e0: energy vector initialized at s0
  arma::uvec idx(2);
  arma::vec e0(nf + 1);

  // {matrix}
  // ms0: state sequence initialized at s0
  arma::mat ms0(nf + 1, ns);

  // initialize
  ms0.row(0) = s0;
  e0(0) = energy_cpp(s0, alpha, beta);

  // ---- initial path ----
  // s: current state
  // u: temporary path sequence
  // ms: state sequence for the current path
  // e: energy sequence for the current path
  arma::rowvec s = s0;
  arma::uvec u = path;
  arma::mat ms = ms0;
  arma::vec e = e0;

  // state sequence from s0 to s1
  for (arma::uword i = 0; i < nf; ++i) {
    k = u(i);

    // flip one species and update the state
    ms.row(i + 1) = ms.row(i);
    ms(i + 1, k) = 1 - ms(i + 1, k);

    // update energy
    de = -(1.0 - 2.0 * s(k)) * (alpha(k) + arma::dot(beta.row(k), s));
    e(i + 1) = de + e(i);
    s(k) = 1 - s(k);
  }

  // mse: state sequence for the current path with energy
  // stip: state vector and energy at the highest-energy point
  // omega: barrier energy of the current accepted path
  arma::mat mse = arma::join_rows(ms, e);
  arma::rowvec stip = mse.row(e.index_max());
  arma::vec omega(iter, arma::fill::value(e.max()));

  // ---- simulated annealing ----
  if (nf > 1) {

    for (arma::uword t = 1; t < iter; ++t) {

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

        // update energy
        de = -(1.0 - 2.0 * s(k)) * (alpha(k) + arma::dot(beta.row(k), s));
        e(i + 1) = de + e(i);
        s(k) = 1 - s(k);
      }

      // record the maximum energy along the candidate path
      etop = e.max();

      // calculate acceptance probability
      pr = std::min(
        1.0,
        std::exp((omega(t - 1) - etop) / temp)
      );

      // update temperature
      temp *= (1 - r);

      // accept the candidate path with probability pr
      if (arma::randu<double>() < pr) {
        // update path and barrier
        path = u;
        omega(t) = etop;

        // reset matrix
        ms = ms0;

        // reconstruct state sequence
        for (arma::uword i = 0; i < nf; ++i) {
          k = path(i);

          // flip one species and update the state
          ms.row(i + 1) = ms.row(i);
          ms(i + 1, k) = 1 - ms(i + 1, k);
        }

        mse = arma::join_rows(ms, e);
        stip = mse.row(e.index_max());

      } else {
        omega(t) = omega(t - 1);
      }
    }// for loop t

  }// if nf

  // output
  if (mse_out != nullptr)
    *mse_out = mse;

  if (omega_out != nullptr)
    *omega_out = omega;

  return stip;
}

// [[Rcpp::export]]
Rcpp::List findpath_cpp(
    const arma::rowvec& s0,
    const arma::rowvec& s1,
    const arma::rowvec& alpha,
    const arma::mat& beta,
    double temp = 1,
    const double r = 0.01,
    const arma::uword iter = 10000
) {
  arma::mat mse;
  arma::vec omega;

  arma::rowvec stip = findpath_inline(
    s0, s1, alpha, beta,
    temp, r, iter,
    &mse, &omega
  );

  return Rcpp::List::create(
    Rcpp::Named("state") = stip,
    Rcpp::Named("path") = mse,
    Rcpp::Named("omega") = omega
  );
}

// [[Rcpp::export]]
Rcpp::List ridge_cpp(
    const arma::mat& sse,
    const arma::rowvec& alpha,
    const arma::mat& beta,
    double temp = 1,
    const double r = 0.01,
    const arma::uword iter = 10000,
    Rcpp::Nullable<arma::uvec> index = R_NilValue
) {
  // ---- declare ----
  // {scalar}
  // ess0, ess1: energy of the two stable states
  // etip: energy at a tipping point
  // cost: cumulative energy cost along the path
  // ed: energy difference between consecutive states
  // ns: number of species
  // nss: number of stable states
  // nr: number of state pairs
  // k: index of the state pair
  double ess0, ess1, etip, cost, ed;
  const arma::uword ns = alpha.n_elem;
  const arma::uword nss = sse.n_rows;
  const arma::uword nr = nss * (nss - 1) / 2;
  arma::uword k = 0;

  if (sse.n_cols != ns + 1)
    Rcpp::stop("sse must contain ns species columns plus one energy column.");

  // {vector}
  // idx: index for stable states (if any)
  // stip: state vector and energy at a tipping point
  // e: energy of each stable state
  // epath: energy along the selected path
  arma::uvec idx;
  arma::rowvec stip;
  arma::vec e = sse.col(sse.n_cols - 1);
  arma::vec epath;

  if (index.isNull()) {
    idx = arma::regspace<arma::uvec>(1, nss);
  } else {
    idx = index;
  }

  // {matrix}
  // ss: matrix of stable states
  // barrier: output matrix for energy/barrier summary
  // rs: output matrix for "tipping point" states of all stable state pairs
  // mse: state sequence and energy along the selected path
  arma::mat ss = sse.cols(0, ns - 1);
  arma::mat barrier(nr, 8);
  arma::mat rs(nr, ns + 3);
  arma::mat mse;

  for (arma::uword i = 0; i < nss - 1; ++i) {
    for (arma::uword j = i + 1; j < nss; ++j) {

      // calculate path from state i to state j
      stip = findpath_inline(
        ss.row(i), ss.row(j),
        alpha, beta,
        temp, r, iter,
        &mse
      );

      // stack state vectors
      rs.row(k).cols(0, ns) = stip;

      // barrier summary
      ess0 = e(i);
      ess1 = e(j);

      etip = stip.back();
      epath = mse.col(mse.n_cols - 1);

      cost = 0;

      for (arma::uword m = 0; m < epath.n_elem - 1; ++m) {
        ed = epath(m + 1) - epath(m);
        cost += std::max(0.0, ed);
      }

      if (ess0 > ess1) {
        // higher-energy state
        barrier(k, 0) = idx(i);
        rs(k, ns + 1) = idx(i);

        // lower-energy state
        barrier(k, 1) = idx(j);
        rs(k, ns + 2) = idx(j);
      } else {
        // higher-energy state
        barrier(k, 0) = idx(j);
        rs(k, ns + 1) = idx(j);

        // lower-energy state
        barrier(k, 1) = idx(i);
        rs(k, ns + 2) = idx(i);
      }

      barrier(k, 2) = std::max(ess0, ess1); // higher stable-state energy
      barrier(k, 3) = std::min(ess0, ess1); // lower stable-state energy
      barrier(k, 4) = etip; // tipping-point energy
      barrier(k, 5) = arma::accu(ss.row(i) != ss.row(j)); // L1 distance between two states
      barrier(k, 6) = cost; // cumulative energy cost
      barrier(k, 7) = etip - barrier(k, 2); // energy barrier

      ++k;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("barrier") = barrier,
    Rcpp::Named("state") = rs
  );
}

// [[Rcpp::export]]
arma::uword find_shallow_cpp(
    const arma::mat& barrier
) {
  // ---- key column indices ----
  // idx_depth: column index for basin depth
  // idx_cost: column index for path energy cost
  const arma::uword idx_depth = barrier.n_cols - 1;
  const arma::uword idx_cost = barrier.n_cols - 2;

  // identify rows with the minimum basin depth
  arma::uvec idx_sh =
    arma::find(barrier.col(idx_depth) == barrier.col(idx_depth).min());

  // select the minimum path cost among the shallowest candidates
  arma::mat candid = barrier.rows(idx_sh);
  arma::uword r = idx_sh(candid.col(idx_cost).index_min());

  return r;
}

// [[Rcpp::export]]
Rcpp::List prune_cpp(
    arma::mat& barrier,
    const double th = 0.2
) {
  // ---- declare ----
  // idx_depth: column index for basin depth
  // idx_rm: identifier of the basin to be removed
  // dmax, dmin: maximum and minimum basin depth
  arma::mat map(barrier.n_rows, 2);
  arma::uword idx_depth = barrier.n_cols - 1;
  arma::uword idx_dist = barrier.n_cols - 3;
  arma::uword idx_rm;
  arma::uword k = 0;
  double dmax, dmin;

  while (true) {
    // stop if no barriers remain
    if (barrier.n_rows == 0)
      break;

    // find pairs with dist == 1
    arma::uvec nei = arma::find(barrier.col(idx_dist) == 1);

    // find deepest and shallowest basin
    arma::uword r = find_shallow_cpp(barrier);
    dmax = barrier.col(idx_depth).max();
    dmin = barrier(r, idx_depth);

    // stop if the shallowest basin is sufficiently deep
    if (nei.is_empty() && dmin > th * dmax)
      break;

    // identify basin to be removed
    idx_rm = static_cast<arma::uword>(barrier(r, 0));

    // retain pairs not involving the shallow basin
    arma::uvec keep =
      arma::find(
        (barrier.col(0) != idx_rm) && (barrier.col(1) != idx_rm)
      );

    // record basin mapping for merge
    map(k, 0) = barrier(r, 0); // shallower basin
    map(k, 1) = barrier(r, 1); // deeper basin
    barrier = barrier.rows(keep);

    k++;
  }

  // remove excess rows
  map.resize(k, 2);

  if (k == 0) {
    return Rcpp::List::create(
      Rcpp::Named("barrier") = barrier,
      Rcpp::Named("map") = R_NilValue
    );
  }

  if (barrier.n_rows == 0) {
    return Rcpp::List::create(
      Rcpp::Named("barrier") = R_NilValue,
      Rcpp::Named("map") = map
    );
  }

  return Rcpp::List::create(
    Rcpp::Named("barrier") = barrier,
    Rcpp::Named("map") = map
  );
}
