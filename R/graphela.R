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
  check_dim(state, alpha, beta)

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
  check_dim(state, alpha, beta)

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
#' @param seed An optional integer used to control random-number
#'   generation. If `NULL`, the current random-number state is used.
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
    replace = TRUE,
    seed = NULL
) {
  ## validate input
  check_dim(
    state = NULL,
    alpha = alpha,
    beta = beta
  )

  if (length(n) != 1 ||
      !is.numeric(n) ||
      !is.finite(n) ||
      n < 1 ||
      n != as.integer(n))
    stop("`n` must be a positive integer.")

  if (length(replace) != 1 ||
      !is.logical(replace))
    stop("`replace` must be a single logical value.")

  ## run cpp function
  run <- function() {
    rss_cpp(
      alpha = alpha,
      beta = beta,
      n = n,
      replace = replace
    )
  }

  if (is.null(seed)) {
    run()
  } else {
    withr::with_seed(seed, run())
  }

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
  check_dim(rbind(s0, s1), alpha, beta)
  check_sa(temp, r, iter)

  ## run cpp function
  path <- findpath_cpp(
    s0 = s0,
    s1 = s1,
    alpha = alpha,
    beta = beta,
    temp = temp,
    r = r,
    iter = iter
  )

  path$omega <- drop(path$omega)

  ## return
  path
}


#' Identify energy ridges between stable states
#'
#' Identifies energy ridges separating pairs of stable states using
#' simulated annealing. For each pair of stable states, the function
#' identifies a tipping point along the transition path and calculates
#' the associated path cost and energy barrier.
#'
#' @param m A matrix of stable states. Each row represents a stable state,
#'   and the last column must contain its energy.
#' @param alpha A numeric vector of model parameters controlling the
#'   intrinsic contribution of each species to system energy.
#' @param beta A numeric matrix of pairwise interaction parameters among
#'   species.
#' @param focus A character string specifying which results to return.
#'   `"barrier"` returns energy-ridge information for each pair of stable
#'   states, `"state"` returns the tipping-point states and associated
#'   information, and `"all"` returns both as a list. Defaults to
#'   `"barrier"`.
#' @param temp A numeric value specifying the initial temperature for
#'   simulated annealing. Defaults to `10`.
#' @param r A numeric value specifying the cooling rate of simulated
#'   annealing. Defaults to `0.001`.
#' @param iter An integer specifying the number of simulated annealing
#'   iterations. Defaults to `10000`.
#' @param seed An optional integer used to control random-number
#'   generation. If `NULL`, the current random-number state is used.
#'
#' @useDynLib graphela, .registration = TRUE
#' @importFrom Rcpp evalCpp
#'
#' @return If `focus = "barrier"`, a matrix with one row for each pair
#'   of stable states and seven columns:
#'   \describe{
#'     \item{ss1}{Index of the shallower stable state.}
#'     \item{ss2}{Index of the deeper stable state.}
#'     \item{e1}{Energy of the shallower stable state.}
#'     \item{e2}{Energy of the deeper stable state.}
#'     \item{tp}{Energy at the tipping point along the transition path.}
#'     \item{cost}{Energy cost of the transition path.}
#'     \item{barrier}{Energy barrier separating the two stable states.}
#'   }
#'   If `focus = "state"`, a matrix containing the tipping-point state
#'   vectors, their energies, and the corresponding stable-state indices.
#'   If `focus = "all"`, a list containing both `"barrier"` and `"state"`
#'   matrices.
#'
#' @export

ridge <- function(
    m,
    alpha,
    beta,
    focus = c("barrier", "state", "all"),
    temp = 10,
    r = 0.001,
    iter = 10000,
    seed = NULL
) {
  ## validate input
  focus <- match.arg(focus)
  s <- check_sse(m, alpha, beta)
  check_sa(temp, r, iter)

  ## run analysis
  run <- function() {
    ridge_cpp(
      sse = m,
      alpha = alpha,
      beta = beta,
      temp = temp,
      r = r,
      iter = iter
    )
  }

  if (is.null(seed)) {
    res <- run()
  } else {
    res <- withr::with_seed(seed, run())
  }

  ## format output
  ## - barrier matrix
  colnames(res$barrier) <- c(
    "ss1", "ss2", "e1", "e2", "tp", "cost", "barrier"
  )

  ## - state matrix
  state_names <- colnames(m)[seq_len(s)]

  if (is.null(state_names))
    state_names <- as.character(seq_len(s))

  colnames(res$state) <- c(state_names, "energy", "ss1", "ss2")

  ## return
  if (focus == "all")
    res
  else
    res[[focus]]

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

  ## expected output format from ridge()
  cnm <- c(
    "ss1", "ss2", "e1", "e2", "tp", "cost", "barrier"
  )

  ## validate input
  if (!is.matrix(m) || !is.numeric(m))
    stop("`m` must be a numeric matrix.")

  if (ncol(m) != length(cnm) ||
      !identical(colnames(m), cnm))
    stop("`m` must conform to the 'barrier' format of `ridge()`.")

  if (!is.numeric(th) || length(th) != 1L || is.na(th) || th < 0)
    stop("`th` must be a single non-negative numeric value.")

  ## run cpp function
  res <- prune_cpp(
    barrier = m,
    th = th
  )

  colnames(res$barrier) <- cnm

  res
}


#' Identify basins of attraction
#'
#' Identifies stable states from random or exhaustive sampling, estimates
#' transition barriers and tipping points among stable states, and prunes
#' shallow basins.
#' The resulting basins are summarized by their energy, depth, and width.
#' Both the pruned results and the underlying unpruned results are retained.
#'
#' @param alpha A numeric vector of model parameters controlling the intrinsic
#'   contribution of each species to system energy.
#' @param beta A numeric matrix of pairwise interaction parameters among
#'   species.
#' @param cnm An optional character vector of names for the species-state
#'   variables. If `NULL`, variables are named sequentially from `1` to
#'   the number of species.
#' @param n An integer specifying the number of initial states sampled by
#'   `rss()` to identify stable states. Defaults to `10000`.
#' @param replace A logical value indicating whether initial states are
#'   sampled with replacement by `rss()`. Defaults to `TRUE`.
#' @param temp A numeric value specifying the initial temperature used for
#'   simulated annealing in the ridge search. Defaults to `10`.
#' @param r A numeric value specifying the cooling rate used for simulated
#'   annealing in the ridge search. Defaults to `0.001`.
#' @param iter An integer specifying the maximum number of iterations used
#'   for each ridge search. Defaults to `10000`.
#' @param th A numeric threshold used by `prune()` to remove shallow
#'   transitions and merge the corresponding stable states. Defaults to `0.2`.
#' @param seed An optional integer used to control random-number generation
#'   in `rss()` and `ridge()`. If `NULL`, the current random-number state is
#'   used.
#'
#' @useDynLib graphela, .registration = TRUE
#' @importFrom Rcpp evalCpp
#'
#' @return A list containing two components:
#'   \describe{
#'     \item{pruned}{Results after pruning shallow basins. Contains:
#'       \describe{
#'         \item{state}{A matrix containing the species-state configuration
#'           and energy of each final basin.}
#'         \item{summary}{A data frame summarizing each final basin, including
#'           its stable-state ID (`ss`), energy, basin depth, and basin width.}
#'         \item{tps}{A matrix containing the tipping-point states associated
#'           with transitions among the final basins.}
#'       }
#'     }
#'     \item{raw}{Results before pruning. Contains:
#'       \describe{
#'         \item{state}{A matrix containing the unique stable-state
#'           configurations and their energies identified by `rss()`.}
#'         \item{barrier}{A data frame containing the transition barriers
#'           identified by `ridge()`.}
#'         \item{tps}{A matrix containing the tipping-point states identified
#'           by `ridge()`.}
#'         \item{map}{A matrix describing the mapping of stable-state IDs
#'           before pruning to IDs after merging. `NULL` if no merging occurs.}
#'       }
#'     }
#'   }
#'
#' @details
#' Stable states are first identified using `rss()` and sorted by energy.
#' Duplicate stable-state configurations are then removed. If only one unique
#' stable state is identified, that state is returned directly without ridge
#' searching or pruning.
#'
#' When multiple stable states are present, `ridge()` is used to identify
#' transition barriers and tipping-point states among all unique stable states.
#' Shallow basins are subsequently pruned using `prune()`. Tipping points are
#' retained only for transitions involving stable states that remain after
#' pruning.
#'
#' Basin depth is calculated as the minimum energy barrier among transitions
#' originating from each stable state.
#'
#' Basin width is calculated as the proportion of the initial stable-state
#' assignments from `rss()` that belong to each final basin. Stable states
#' merged during pruning are therefore combined when calculating basin width.
#'
#' The returned object also stores `alpha`, `beta`, and `seed` as attributes.
#'
#' @export

basin <- function(
    alpha,
    beta,
    cnm = NULL,
    n = 10000,
    replace = TRUE,
    temp = 10,
    r = 0.001,
    iter = 10000,
    th = 0.2,
    seed = NULL
) {

  ## stable states
  m_ss <- rss(
    alpha = alpha,
    beta = beta,
    n = n,
    replace = replace,
    seed = seed
  )

  s <- length(alpha)

  if (is.null(cnm))
    state_names <- as.character(seq_len(s))
  else
    state_names <- cnm

  colnames(m_ss) <- c(state_names, "energy")

  ## sort stable states by energy
  idx <- order(m_ss[, ncol(m_ss), drop = TRUE])
  m_ss <- m_ss[idx, , drop = FALSE]

  ## assign unique integer IDs to stable states based on energy
  label <- apply(
    m_ss[, seq_len(s), drop = FALSE],
    MARGIN = 1,
    \(x) paste0(x, collapse = "")
  )

  v_ss <- factor(label, levels = unique(label)) |>
    as.numeric()

  rownames(m_ss) <- v_ss

  ## retain unique stable states
  m_uss <- unique(m_ss)

  if (nrow(m_uss) == 1) {
    ## if only one stable state
    return(
      structure(
        ## main output
        list(
          pruned = list(
            ## pruned stable states
            state = m_uss,

            ## summary of energy, depth, and width for each basin
            summary = data.frame(
              ss = 1,
              energy = m_uss[, ncol(m_uss)],
              depth = NA,
              width = 1.0,
              row.names = NULL
            ),

            ## tipping point state matrix
            tps = NULL
          ),

          raw = list(
            state = m_uss,
            barrier = NULL,
            tps = NULL,
            map = NULL
          )
        ),

        ## attributes
        alpha = alpha,
        beta = beta,
        temp = temp,
        r = r,
        iter = iter,
        th = th,
        seed = seed
      )
    )
  }

  ## ridge and pruning
  list_r <- ridge(
    m = m_uss,
    alpha = alpha,
    beta = beta,
    focus = "all",
    temp = temp,
    r = r,
    iter = iter,
    seed = seed
  )

  list_ss <- prune(
    m = list_r$barrier,
    th = th
  )

  ## keep tipping points for basins not pruned
  ss_keep <- unique(c(list_ss$barrier[, 1:2]))
  tp_keep <- apply(
    X = list_r$state[, c("ss1", "ss2")],
    MARGIN = 1,
    \(x) all(x %in% ss_keep)
  )

  m_tps <- list_r$state[tp_keep, , drop = FALSE]

  ## basin depth
  ## each row represents a transition between two stable states:
  ## ss1 -> ss2, with energies e1 and e2 and tipping-point energy tp.
  m_tpe <- list_ss$barrier[, 1:5, drop = FALSE]

  ## include both directions of each transition so that each stable
  ## state can be evaluated as the starting (shallower) state.
  m_depth <- rbind(
    m_tpe,
    m_tpe[, c("ss2", "ss1", "e2", "e1", "tp")]
  ) |>
    transform(depth = tp - e1)

  ## minimum basin depth among all transitions originating from each state
  v_depth <- tapply(
    m_depth[, "depth"],
    m_depth[, "ss1"],
    min
  )

  ## basin width
  ## merge the stable-state IDs through the merging map.
  ## note: this code is sensitive to the row order of the `map` object
  ## validity affirmed by `prune_cpp()` implementation
  v_merge <- v_ss

  if (!is.null(list_ss$map)) {

    for (i in seq_len(nrow(list_ss$map))) {
      v_merge[v_merge == list_ss$map[i, 1]] <- list_ss$map[i, 2]
    }

  }

  ## IDs of the final merged basins
  idx_mss <- unique(v_merge)
  m_mss <- m_uss[idx_mss, , drop = FALSE]

  structure(
    ## main output
    list(
      pruned = list(
        ## pruned stable states
        state = m_mss,

        ## summary of energy, depth, and width for each basin
        summary = data.frame(
          ss = idx_mss,
          energy = m_mss[, ncol(m_mss)],
          depth = v_depth[as.character(idx_mss)],
          width = tabulate(v_merge)[idx_mss] / nrow(m_ss),
          row.names = NULL
        ),

        ## tipping point state matrix
        tps = m_tps
      ),

      raw = list(
        ## raw stable states
        state = m_uss,

        ## ridge information
        barrier = list_r$barrier,

        ## tipping point state matrix
        tps = list_r$state,

        ## mapping from raw ss to merged ss
        map = list_ss$map
      )
    ),

    ## attributes
    alpha = alpha,
    beta = beta,
    temp = temp,
    r = r,
    iter = iter,
    th = th,
    seed = seed
  )
}

#' @export

egap <- function(
    b,
    obs,
    temp = NULL,
    r = NULL,
    iter = NULL,
    th = NULL,
    seed = NULL
) {

  ## validate input
  if (is.null(alpha))
    alpha = attr(b, "alpha")

  if (is.null(beta))
    beta = attr(b, "beta")

  s <- check_dim(obs, alpha, beta)

  ## energy of observed states
  v_e <- apply(
    matrix(obs, ncol = s),
    MARGIN = 1,
    FUN = \(x) energy(x, alpha, beta)
  )

  ## stable states to which observed states belong
  m_oss <- t(
    apply(
      matrix(obs, ncol = s),
      MARGIN = 1,
      FUN = \(x) stpd(x, alpha, beta)
    )
  )

  v_match <- with(b$raw, {
    v_match <- match(
      apply(m_oss[, seq_len(s), drop = FALSE], 1, paste0, collapse = ""),
      apply(state[, seq_len(s), drop = FALSE], 1, paste0, collapse = "")
    )

    for (i in 1:nrow(map))
      v_match[v_match == map[i, 1]] <- map[i, 2]

    v_match
  })

  idx <- sapply(v_match, \(x) which(x == b$pruned$summary$ss))

  ## output
  with(b$pruned$summary, {
    data.frame(
      gap = v_e - energy[idx],
      energy = v_e,
      ss = v_match,
      bottom = energy[idx]
    )
  })

}
