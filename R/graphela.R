#' Calculate energy
#'
#' Calculate community energy of a given state.
#'
#' @param state A binary row vector specifying the community state.
#' @param alpha A numeric vector of model parameters controlling the intrinsic
#'   contribution of each species to system energy.
#' @param beta A numeric matrix of pairwise interaction parameters among
#'   species.
#'
#' @useDynLib graphela, .registration = TRUE
#' @importFrom Rcpp evalCpp
#'
#' @return Numeric value of community energy.
#'
#' @export

energy <- function(
    state,
    alpha,
    beta
) {
  ## validate input
  check_dim(
    state = state,
    alpha = alpha,
    beta = beta
  )

  ## run cpp function
  energy_cpp(
    state = state,
    alpha = alpha,
    beta = beta
  )
}


#' Steepest descent method
#'
#' Identifies a local minimum in system energy using the steepest descent
#' algorithm, starting from a specified initial state.
#'
#' @param state A binary row vector specifying the initial state.
#' @param alpha A numeric vector of model parameters controlling the intrinsic
#'   contribution of each species to system energy.
#' @param beta A numeric matrix of pairwise interaction parameters among
#'   species.
#'
#' @useDynLib graphela, .registration = TRUE
#' @importFrom Rcpp evalCpp
#'
#' @return A binary state vector corresponding to the stable state reached by
#'   the steepest descent algorithm, along with its energy in the last element.
#'
#' @export

stpd <- function(
    state,
    alpha,
    beta
) {
  ## validate input
  check_dim(
    state = state,
    alpha = alpha,
    beta = beta
  )

  ## run cpp function
  stpd_cpp(
    state = state,
    alpha = alpha,
    beta = beta
  )
}


#' Identify stable states
#'
#' Identifies stable states by applying the steepest descent algorithm to
#' randomly sampled initial states.
#'
#' @param alpha A numeric vector of model parameters controlling the intrinsic
#'   contribution of each species to system energy.
#' @param beta A numeric matrix of pairwise interaction parameters among
#'   species.
#' @param n An integer specifying the number of initial states to sample.
#'   Defaults to `10000`.
#' @param replace A logical value indicating whether identical initial states can be
#'   sampled more than once. Defaults to `TRUE`.
#'
#' @useDynLib graphela, .registration = TRUE
#' @importFrom Rcpp evalCpp
#'
#' @return A matrix of stable states and their corresponding energy values.
#'
#' @export

rss <- function(
    alpha,
    beta,
    n = 10000,
    replace = TRUE
) {
  ## validate input
  check_dim(
    state = NULL,
    alpha = alpha,
    beta = beta
  )

  ## run cpp function
  rss_cpp(
    alpha = alpha,
    beta = beta,
    n = n,
    replace = replace
  )

}


#' Identify a transition path between stable states
#'
#' Identifies a transition path between two stable states using simulated
#' annealing.
#'
#' @param s0 A binary numeric vector specifying the initial stable state.
#' @param s1 A binary numeric vector specifying the destination stable state.
#' @param alpha A numeric vector of model parameters controlling the intrinsic
#'   contribution of each species to system energy.
#' @param beta A numeric matrix of pairwise interaction parameters among
#'   species.
#' @param temp A numeric value specifying the initial temperature of simulated
#'   annealing.
#' @param r A numeric value specifying the cooling rate of simulated annealing.
#'   Must be between 0 and 1.
#' @param iter An integer specifying the number of simulated annealing
#'   iterations.
#'
#' @return A matrix representing the transition path from `s0` to `s1`, with
#'   the energy of each state in the last column.
#'
#' @export

findpath <- function(
    s0,
    s1,
    alpha,
    beta,
    temp,
    r,
    iter
) {
  ## validate input
  lapply(
    list(s0, s1),
    FUN = function(s) {
      check_dim(
        state = s,
        alpha = alpha,
        beta = beta
      )
    })

  if (temp < 0)
    stop("`temp` must be positive numeric.")

  if (r < 0 || r > 1)
    stop("`r` must be between 0 and 1.")

  ## run cpp function
  path <- findpath_cpp(
    s0 = s0,
    s1 = s1,
    alpha = alpha,
    beta = beta,
    temp = temp,
    r = r,
    n = iter
  )

  path$omega <- drop(path$omega)

  return(path)
}


#' Identify energy ridges between stable states
#'
#' Identifies the energy ridge separating pairs of stable states using
#' simulated annealing. For each pair of stable states, the function identifies
#' a tipping point along the transition path and calculates the associated
#' path cost and energy barrier.
#'
#' @param m A matrix of stable states. The last column must contain the
#'   energy of each stable state.
#' @param alpha A numeric vector of model parameters controlling the intrinsic
#'   contribution of each species to system energy.
#' @param beta A numeric matrix of pairwise interaction parameters among
#'   species.
#' @param temp A numeric value specifying the initial temperature for simulated
#'   annealing. Defaults to `10`.
#' @param r A numeric value specifying the cooling rate of simulated annealing.
#'   Defaults to `0.001`.
#' @param iter An integer specifying the number of simulated annealing
#'   iterations. Defaults to `10000`.
#'
#' @useDynLib graphela, .registration = TRUE
#' @importFrom Rcpp evalCpp
#'
#' @return A matrix with one row for each pair of stable states and seven
#'   columns:
#'   \describe{
#'     \item{ss1}{Index of the shallower stable state.}
#'     \item{ss2}{Index of the deeper stable state.}
#'     \item{e1}{Energy of the shallower stable state.}
#'     \item{e2}{Energy of the deeper stable state.}
#'     \item{tp}{Energy at the tipping point along the transition path.}
#'     \item{cost}{Energy cost of the transition path.}
#'     \item{barrier}{Energy barrier separating the two stable states.}
#'   }
#'
#' @export

ridge <- function(
    m,
    alpha,
    beta,
    temp = 10,
    r = 0.001,
    iter = 10000
) {
  ## validate input
  check_dim(
    state = NULL,
    alpha = alpha,
    beta = beta
  )

  ## validate matrix
  if (nrow(m) == 1)
    stop("Only one stable state: no ridge can be defined.")

  s <- ncol(m) - 1
  if (s != length(alpha) || !all(dim(beta) == c(s, s)))
    stop("Invalid matrix dimension: `m' must contain `length(alpha) + 1` columns")

  if (!all(m[, seq_len(s)] %in% c(0, 1)))
    stop("`state` must be a binary numeric vector.")

  ## run analysis
  res <- ridge_cpp(
    sse = m,
    alpha = alpha,
    beta = beta,
    temp = temp,
    r = r,
    n = iter
  )

  colnames(res) <- c("ss1",
                     "ss2",
                     "e1",
                     "e2",
                     "tp",
                     "cost",
                     "barrier")

  return(res)
}


#' Prune shallow energy basins
#'
#' Removes stable states associated with shallow energy basins based on an
#' energy-barrier threshold.
#'
#' @param m A matrix of pairwise stable-state relationships returned by
#'   [ridge()]. The matrix must conform to the output format of [ridge()].
#' @param th A numeric value between 0 and 1 specifying the threshold used
#'   to prune shallow basins. Defaults to `0.2`.
#'
#' @useDynLib graphela, .registration = TRUE
#' @importFrom Rcpp evalCpp
#'
#' @return A matrix containing the stable-state relationships remaining after
#'   shallow basins have been pruned.
#'
#' @export

prune <- function(m, th = 0.2) {

  ## validate input
  cnm <- c("ss1",
           "ss2",
           "e1",
           "e2",
           "tp",
           "cost",
           "barrier")

  if (any(colnames(m) != cnm))
    stop("The matrix `m` must conform to the output format of `ridge()`")

  ## run cpp function
  res <- prune_cpp(
    pem = m,
    th = th
  )

  colnames(res$pem) <- cnm
  return(res)
}


#' Identify ecological basins from stable states and transition dynamics
#'
#' Identifies stable states from random or exhaustive sampling, estimates
#' transitions among stable states, and prunes shallow basins. Summarizes
#' the resulting basins by their stable-state configuration, energy, depth,
#' and width.
#'
#' @param alpha A numeric vector of model parameters controlling the intrinsic
#'   contribution of each species to system energy.
#' @param beta A numeric matrix of pairwise interaction parameters among
#'   species.
#' @param n An integer specifying the number of initial random states to estimate stable states.
#'   Defaults to `10000`.
#' @param replace A logical value indicating whether initial random states are sampled
#'   with replacement. Defaults to `TRUE`.
#' @param temp Temperature of simulated annealing in the ridge
#'   search. Defaults to `10`.
#' @param r Cooling rate of simulated annealing the ridge search. Defaults to `0.001`.
#' @param iter An integer specifying the number of iterations used in the
#'   ridge search. Defaults to `10000`.
#' @param th A numeric threshold used to prune shallow basins. Defaults to `0.2`.
#'
#' @useDynLib graphela, .registration = TRUE
#' @importFrom Rcpp evalCpp
#'
#' @return A list containing:
#'   \describe{
#'     \item{state}{A matrix containing the species-state configuration of
#'       each basin.}
#'     \item{energy}{A data frame summarizing each basin, including its
#'       stable-state ID (`ss`), energy, basin depth, and basin width.}
#'   }
#'
#' @details
#' Stable states are first identified using `rss()` and sorted by energy.
#' Unique stable states are then used to search tipping points between stable states with `ridge()`.
#' Shallow basins are pruned using `prune()`.
#'
#' Basin depth is calculated from the minimum energy barrier among transitions
#' originating from each stable state. Both directions of each transition are
#' considered so that the barrier is evaluated relative to the energy of the
#' starting state.
#'
#' Basin width is calculated as the proportion of initial random states assigned to
#' each final basin.
#'
#' @export

basin <- function(
    alpha,
    beta,
    n = 10000,
    replace = TRUE,
    temp = 10,
    r = 0.001,
    iter = 10000,
    th = 0.2
) {

  ## stable states
  m_ss <- rss(
    alpha = alpha,
    beta = beta,
    n = n,
    replace = replace
  )

  ## sort stable states by energy
  m_ss <- m_ss[order(m_ss[, ncol(m_ss)]), ]

  ## assign unique integer IDs to stable states based on energy
  v_ss <- as.numeric(factor(m_ss[, ncol(m_ss)]))
  rownames(m_ss) <- v_ss

  ## retain unique stable states
  m_uss <- unique(m_ss)

  if (nrow(m_uss) == 1) {
  ## if only one stable state
    return(
      list(
        ## stable-state configurations for the final basins
        state = m_uss,

        ## summary of energy, depth, and width for each basin
        energy = data.frame(
          ss = 1,
          energy = m_uss[, ncol(m_uss)],
          depth = NA,
          width = 1.0,
          row.names = NULL
        )
      )
    )

  }

  ## ridge and pruning
  list_p <- ridge(
    m = m_uss,
    alpha = alpha,
    beta = beta,
    temp = temp,
    r = r,
    iter = iter
  ) |>
    prune(th = th)

  ## basin depth
  ## each row represents a transition between two stable states:
  ## ss1 -> ss2, with energies e1 and e2 and tipping-point energy tp.
  pem <- list_p$pem[, 1:5, drop = FALSE]

  ## include both directions of each transition so that each stable
  ## state can be evaluated as the starting (shallower) state.
  m_depth <- rbind(
    pem,
    pem[, c(2, 1, 4, 3, 5)]
  ) |>
    transform(b = tp - e1)

  ## minimum basin depth among all transitions originating from each state
  v_depth <- tapply(
    m_depth[, "b"],
    m_depth[, "ss1"],
    min
  )

  ## basin width
  ## merge the stable-state IDs through the merging map.
  v_merge <- v_ss

  if (!is.null(list_p$map)) {

    for (i in seq_len(nrow(list_p$map))) {
      v_merge[v_merge == list_p$map[i, 1]] <- list_p$map[i, 2]
    }

  }

  ## IDs of the final merged basins
  idx_mss <- unique(v_merge)
  m_mss <- m_uss[idx_mss, , drop = FALSE]

  list(
    ## stable-state configurations for the final basins
    state = m_mss,

    ## summary of energy, depth, and width for each basin
    energy = data.frame(
      ss = idx_mss,
      energy = m_mss[, ncol(m_mss)],
      depth = v_depth[as.character(idx_mss)],
      width = tabulate(v_merge)[idx_mss] / n,
      row.names = NULL
    )
  )
}

