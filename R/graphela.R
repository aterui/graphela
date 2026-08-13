#' Steepest descent method
#'
#' @param state A binary row vector of initial state
#' @param alpha A numeric vector of model parameters controlling the intrinsic
#'   contribution of each species to system energy.
#' @param beta A numeric matrix of pairwise interaction parameters among
#'   species.
#'
#' @export

stpd <- function(
    state,
    alpha,
    beta
) {

  stpd_cpp(
    state = state,
    alpha = alpha,
    beta = beta
  )

}

#' Identify stable states
#'
#' @param alpha A numeric vector of model parameters controlling the intrinsic
#'   contribution of each species to system energy.
#' @param beta A numeric matrix of pairwise interaction parameters among
#'   species.
#' @param n An integer specifying the number of initial random states.
#'   Defaults to `10000`.
#'
#' @export

rss <- function(
    alpha,
    beta,
    n = 10000,
    replace = TRUE
) {

  rss_cpp(
    alpha = alpha,
    beta = beta,
    n = n,
    replace = replace
  )

}

#' Identify energy ridge between stable states
#'
#' @param m A matrix of stable state vectors. The last column must contain energy of each stable state.
#' @param alpha A numeric vector of model parameters controlling the intrinsic
#'   contribution of each species to system energy.
#' @param beta A numeric matrix of pairwise interaction parameters among
#'   species.
#' @param temp
#' @param r A scalar specifing the cooling rate of the simulated annealing.
#'   Defaults to `0.001`.
#' @param iter An integer specifying the number of simulated annealing.
#'   Defaults to `10000`.
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

  ridge_cpp(
    sse = m,
    alpha = alpha,
    beta = beta,
    temp = temp,
    r = r,
    n = iter
  )

}

#' Identify ecological basins from stable states and transition dynamics
#'
#' Identifies stable states from random or exhaustive sampling, estimates
#' transitions among stable states, prunes weak transitions, and summarizes
#' the resulting basins by their stable-state configuration, energy, depth,
#' and width.
#'
#' @param alpha A numeric vector of model parameters controlling the intrinsic
#'   contribution of each species to system energy.
#' @param beta A numeric matrix of pairwise interaction parameters among
#'   species.
#' @param n An integer specifying the number of states sampled by `rss`.
#'   Defaults to `10000`.
#' @param replace A logical value indicating whether initial states are sampled
#'   with replacement. Defaults to `TRUE`.
#' @param temp A numeric value controlling the temperature used in the ridge
#'   search. Defaults to `10`.
#' @param r A numeric value controlling the ridge search. Defaults to `0.01`.
#' @param iter An integer specifying the number of iterations used in the
#'   ridge search. Defaults to `5000`.
#' @param th A numeric threshold used to prune transitions. Defaults to `0.2`.
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
#' Stable states are first identified using `rss` and sorted by energy.
#' Unique stable states are then used as starting points for `ridge`, and
#' transitions are pruned using `prune`.
#'
#' Basin depth is calculated from the minimum energy barrier among transitions
#' originating from each stable state. Both directions of each transition are
#' considered so that the barrier is evaluated relative to the energy of the
#' starting state.
#'
#' Basin width is calculated as the proportion of sampled states assigned to
#' each final basin after applying the transition map.
#'
#' @export

basin <- function(
    alpha,
    beta,
    n = 10000,
    replace = TRUE,
    temp = 10,
    r = 0.01,
    iter = 5000,
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

  ## ridge and pruning
  list_p <- ridge(
    sse = m_uss,
    alpha = alpha,
    beta = beta,
    temp = temp,
    r = r,
    n = iter
  ) |>
    prune(th = th)

  ## basin depth
  ## each row represents a transition between two stable states:
  ## ss1 -> ss2, with energies e1 and e2 and tipping-point energy tp.
  pem <- list_p$pem[, 1:5]
  colnames(pem) <- c("ss1", "ss2", "e1", "e2", "tp")

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

  for (i in seq_len(nrow(list_p$map))) {
    v_merge[v_merge == list_p$map[i, 1]] <- list_p$map[i, 2]
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
