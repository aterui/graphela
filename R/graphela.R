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
    "ss1", "ss2", "e1", "e2", "tp", "dist", "cost", "barrier"
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
    "ss1", "ss2", "e1", "e2", "tp", "dist", "cost", "barrier"
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

  if (!is.null(res$barrier))
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

  ## assign unique integer IDs to stable states
  label <- apply(
    m_ss[, seq_len(s), drop = FALSE],
    MARGIN = 1,
    \(x) paste0(x, collapse = "")
  )

  v_ss <- factor(label, levels = unique(label)) |>
    as.numeric()

  rownames(m_ss) <- v_ss

  ## retain unique stable states
  m_uss <- m_ss[!duplicated(v_ss), , drop = FALSE]

  ## if only one stable state
  if (nrow(m_uss) == 1) {

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

  ## if all states but one are pruned
  if (is.null(list_ss$barrier)) {

    ## merge ss indices
    v_merge <- v_ss

    if (!is.null(list_ss$map)) {

      for (i in seq_len(nrow(list_ss$map))) {
        v_merge[v_merge == list_ss$map[i, 1]] <- list_ss$map[i, 2]
      }

    }

    idx_ss <- unique(v_merge)

    return(
      structure(
        ## main output
        list(
          pruned = list(
            ## pruned stable states
            state = m_uss[idx_ss, , drop = FALSE],

            ## summary of energy, depth, and width for each basin
            summary = data.frame(
              ss = idx_ss,
              energy = m_uss[idx_ss, ncol(m_uss)],
              depth = NA,
              width = 1.0,
              row.names = NULL
            ),

            ## tipping point state matrix
            tps = NULL
          ),

          raw = list(
            state = m_uss,
            barrier = list_r$barrier,
            tps = list_r$state,
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
    )

  }

  ## keep tipping points for basins not pruned
  ss_keep <- unique(c(list_ss$barrier[, c("ss1", "ss2")]))
  tp_keep <- apply(
    X = list_r$state[, c("ss1", "ss2"), drop = FALSE],
    MARGIN = 1,
    \(x) all(x %in% ss_keep)
  )

  m_tps <- list_r$state[tp_keep, , drop = FALSE]

  ## basin depth
  ## each row represents a transition between two stable states:
  ## ss1 -> ss2, with energies e1 and e2 and tipping-point energy tp.
  m_tpe <- list_ss$barrier[, c("ss1", "ss2", "e1", "e2", "tp"), drop = FALSE]

  ## include both directions of each transition so that each stable
  ## state can be evaluated as the starting (shallower) state.
  m_depth <- rbind(
    m_tpe,
    m_tpe[, c("ss2", "ss1", "e2", "e1", "tp"), drop = FALSE]
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


#' Calculate energy gaps between observations and stable states
#'
#' Calculates the energy gap between each observed community state and the
#' stable state reached by steepest descent. Observed states are first
#' assigned to stable states identified during the original basin analysis.
#' If an observation leads to a previously unidentified stable state, the
#' stable-state set is expanded, transition barriers are recalculated, and
#' shallow basins are pruned.
#'
#' The function returns the energy of each observed state, the energy of its
#' associated stable state, and the resulting energy gap. It also provides
#' an updated set of retained stable states and their basin summaries.
#'
#' @param b A basin object returned by [basin()]. The model parameters
#'   `alpha` and `beta`, as well as the basin-analysis settings, are retrieved
#'   from attributes of this object.
#' @param obs A numeric matrix or vector of observed community states, with
#'   rows representing observations and columns representing species.
#'   States must be compatible with the dimensions of `alpha` and `beta`.
#' @param temp Optional temperature parameter used when recalculating
#'   transition barriers for newly discovered stable states. If `NULL`, the
#'   value stored in `b` is used.
#' @param r Optional parameter controlling the transition-barrier search.
#'   If `NULL`, the value stored in `b` is used.
#' @param iter Optional number of iterations used for transition-barrier
#'   estimation. If `NULL`, the value stored in `b` is used.
#' @param th Optional threshold for pruning shallow basins. The value stored
#'   in `b` is always used when newly discovered stable states require
#'   re-evaluation.
#' @param seed Optional random seed used when recalculating transition
#'   barriers. If `NULL`, the seed stored in `b` is used.
#'
#' @return A list with three components:
#' \describe{
#'   \item{gap}{A data frame containing the energy gap (`gap`), observed-state
#'   energy (`energy`), associated stable-state index (`ss`), and stable-state
#'   energy (`bottom`) for each observation.}
#'   \item{state}{The stable states retained after incorporating the observed
#'   states and re-evaluating transition barriers and basin pruning.}
#'   \item{summary}{A data frame summarizing the retained stable states,
#'   including their energy, basin depth, and basin width. Newly discovered
#'   stable states have `NA` for depth and width because these quantities are
#'   not estimated by `egap()`.}
#' }
#'
#' @details
#' For each observation, [stpd()] is used to identify the stable state reached
#' by steepest descent. The energy gap is calculated as the difference between
#' the observed-state energy and the energy of this associated stable state.
#'
#' If all observed states lead to stable states already identified by
#' [basin()], the original basin results are retained. If new stable states
#' are discovered, they are added to the original stable-state set, duplicate
#' states are removed, and transition barriers are recalculated using
#' [ridge()]. The expanded set is then pruned using [prune()] before the
#' observed states are reassigned to the resulting stable states.
#'
#' @seealso [basin()], [stpd()], [ridge()], [prune()]
#'
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

  ## validate input and retrieve model parameters
  ## alpha and beta are stored as attributes of the basin object
  alpha <- attr(b, "alpha")
  beta  <- attr(b, "beta")

  ## check dimensions and validity of the observed states
  ## returns the number of species/states (s)
  s <- check_dim(obs, alpha, beta)

  ## calculate energy of each observed state
  v_e <- apply(
    X = matrix(obs, ncol = s),
    MARGIN = 1,
    FUN = \(x) energy(x, alpha, beta)
  )

  ## identify stable states associated with observations
  ## each observed state is assigned to the stable state reached
  ## by steepest descent
  m_oss <- t(
    apply(
      X = matrix(obs, ncol = s),
      MARGIN = 1,
      FUN = \(x) stpd(x, alpha, beta)
    )
  )

  ## match observed stable states to original stable states
  idx0 <- with(b$raw, {

    ## match observed stable states to the stable states identified
    ## in the original basin analysis
    v_match <- match(
      apply(
        m_oss[, seq_len(s), drop = FALSE],
        1,
        paste0,
        collapse = ""
      ),
      apply(
        state[, seq_len(s), drop = FALSE],
        1,
        paste0,
        collapse = ""
      )
    )

    ## update indices if stable states were merged during pruning
    if (!is.null(map)) {

      ## skip if no merging occurred in prune()
      for (i in 1:nrow(map))
        v_match[v_match == map[i, 1]] <- map[i, 2]
    }

    v_match
  })

  ## return if all observed stable states are already known
  if (!any(is.na(idx0))) {

    ## no new stable states were found
    return(
      list(
        gap = with(b$raw, {
          data.frame(
            gap = v_e - state[idx0, "energy"],
            energy = v_e,
            ss = idx0,
            bottom = state[idx0, "energy"]
          )
        }),
        state = b$pruned$state,
        summary = b$pruned$summary
      )
    )
  }

  message("New states were found; re-evaluate pruning")

  ## combine original and observed stable states
  ## append observed stable states to the original set and remove
  ## duplicate states
  m_uss <- rbind(
    b$raw$state,
    m_oss
  ) |>
    unique()

  ## assign sequential row names for stable-state indexing
  rownames(m_uss) <- seq_len(nrow(m_uss))

  ## identify barriers and prune the expanded stable-state set
  ## calculate transition barriers among original and newly observed
  ## stable states, then remove shallow basins
  list_ss <- ridge(
    m = m_uss,
    alpha = alpha,
    beta = beta,
    focus = "barrier",
    temp = if(is.null(temp)) attr(b, "temp") else temp,
    r = if(is.null(r)) attr(b, "r") else r,
    iter = if(is.null(iter)) attr(b, "iter") else iter,
    seed = if(is.null(seed)) attr(b, "seed") else seed
  ) |>
    prune(th = attr(b, "th"))

  ## identify stable states retained after pruning
  ## collect stable states that participate in at least one retained
  ## barrier
  ss_keep <- c(list_ss$barrier[, c("ss1", "ss2")]) |>
    unique() |>
    sort()

  ## match observed stable states to expanded stable-state set
  idx1 <- with(list_ss, {

    ## match observed stable states to the combined set of original
    ## and newly observed stable states
    v_match <- match(
      apply(
        m_oss[, seq_len(s), drop = FALSE],
        1,
        paste0,
        collapse = ""
      ),
      apply(
        m_uss[, seq_len(s), drop = FALSE],
        1,
        paste0,
        collapse = ""
      )
    )

    ## update indices if stable states were merged during pruning
    if (!is.null(map)) {

      ## skip if no merging occurred in prune()
      for (i in 1:nrow(map))
        v_match[v_match == map[i, 1]] <- map[i, 2]
    }

    v_match
  })

  ## identify newly discovered stable states
  ## stable states retained after pruning that were not in the
  ## original stable-state set
  ss_new <- setdiff(
    ss_keep,
    unique(idx0)
  )

  ## return energy gaps and updated stable-state summary
  list(
    ## energy gap between each observation and its associated
    ## stable state
    gap = data.frame(
      gap = v_e - m_uss[idx1, "energy"],
      energy = v_e,
      ss = idx1,
      bottom = m_uss[idx1, "energy"]
    ),

    ## stable states retained after incorporating observations
    state = m_uss[ss_keep, ],

    ## append newly discovered stable states to the original summary
    ## depth and width are not estimated for these new states
    summary = rbind(
      b$pruned$summary,
      data.frame(
        ss = ss_new,
        energy = m_uss[ss_new, "energy"],
        depth = NA,
        width = NA
      )
    )
  )
}


#' Fit regularized regressions for a multivariate response (conditional Markov random fields)
#'
#' Fits a cross-validated regularized regression model for each response
#' variable in \code{Y}, using the predictors in \code{X} and the remaining
#' response variables as predictors. Coefficients are separated into effects
#' of the predictors in \code{X} (\code{alpha}) and pairwise effects among
#' response variables (\code{beta}).
#'
#' @param Y A matrix or data frame of response variables. Must have the same
#'   number of rows as \code{X}.
#' @param X A matrix or data frame of predictor variables.
#' @param family A character string specifying the response distribution used
#'   by \code{glmnet::cv.glmnet()}, such as \code{"gaussian"},
#'   \code{"binomial"}, or \code{"poisson"}.
#' @param sym.method A character string specifying the method used to
#'   symmetrize pairwise coefficients. Passed to \code{symmetrize()}.
#' @param lambda.method A character string specifying the criterion used to
#'   select the regularization parameter when extracting coefficients.
#'   Must be either \code{"lambda.1se"} or \code{"lambda.min"}.
#' @param type.measure A character string specifying the loss used to evaluate
#'   models during cross-validation. Passed to
#'   \code{glmnet::cv.glmnet()}.
#' @param nfolds Number of folds used for cross-validation. Defaults to the
#'   number of rows in \code{Y}.
#' @param grouped Logical; whether to use grouped cross-validation statistics.
#'   Passed to \code{glmnet::cv.glmnet()}.
#' @param control A named list of algorithm control parameters,
#'   providing per-call overrides of session defaults set by
#'   [glmnet::glmnet.control()]. See [glmnet::glmnet()] for details.
#' @param future.seed Logical; whether to generate reproducible random-number
#'   streams for parallel computation.
#' @param progress Logical; whether to show a progress bar.
#'
#' @return A list with two components: \code{alpha}, a matrix of coefficients
#'   for predictors in \code{X}, and \code{beta}, a symmetric matrix of
#'   pairwise coefficients among response variables. The diagonal of
#'   \code{beta} is set to zero. Model-fitting warnings are stored as a
#'   \code{"warning"} attribute containing a data frame with the response
#'   index, response name, and warning message.
#'
#' @export

cmrf <- function(
    Y,
    X = NULL,
    family,
    sym.method = "min",
    lambda.method = c("lambda.1se", "lambda.min"),
    type.measure = "deviance",
    nfolds = nrow(Y),
    grouped = FALSE,
    control = list(),
    future.seed = TRUE,
    progress = TRUE
) {

  ## validate input
  if (is.null(X)) {
    X <- matrix(numeric(0), nrow = nrow(Y), ncol = 0)
  }

  if (nrow(X) != nrow(Y))
    stop("X and Y must have the same number of rows.")

  ## select the lambda criterion used to extract coefficients
  lambda.method <- match.arg(lambda.method)

  ## assign default names to predictors if column names are absent
  if (ncol(X) > 0 && is.null(colnames(X)))
    colnames(X) <- paste0("x", seq_len(ncol(X)))

  ## assign default names to responses if column names are absent
  if (is.null(colnames(Y)))
    colnames(Y) <- paste0("y", seq_len(ncol(Y)))

  ## fit regularized regressions
  fit <- function(i, p = NULL) {

    ## response variable and remaining biotic factors
    y <- Y[, i, drop = TRUE]
    Y_minus_i <- Y[, -i, drop = FALSE]

    ## full predictor matrix, combine abiotic and biotic factors
    ## then remove intercept column
    if (ncol(X) > 0) {

      # w/ predictors
      Z <- stats::model.matrix(
        ~.,
        data = data.frame(X, Y_minus_i)
      )[, -1, drop = FALSE]

    } else {

      # w/o predictors
      Z <- as.matrix(Y_minus_i)

    }

    ## collect warning messages generated during model fitting
    warn <- character()

    ## fit cross-validated regularized regression while capturing warnings
    m <-
      withCallingHandlers(
        glmnet::cv.glmnet(
          x = Z,
          y = y,
          family = family,
          type.measure = type.measure,
          nfolds = nfolds,
          grouped = grouped,
          control = control
        ),
        warning = function(w) {
          warn <<- c(warn, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )

    if (!is.null(p))
      p()

    ## return both the fitted model and any warnings generated
    list(
      model = m,
      warning = warn
    )
  }

  ## fit one regularized regression for each response variable in parallel
  if (progress) {

    res <- progressr::with_progress({

      p <- progressr::progressor(steps = ncol(Y))

      future.apply::future_lapply(
        seq_len(ncol(Y)),
        function(i) {
          fit(i, p = p)
        },
        future.seed = future.seed
      )
    })

  } else {

    res <- future.apply::future_lapply(
      seq_len(ncol(Y)),
      fit,
      future.seed = future.seed
    )

  }

  ## extract fitted models from the results
  list_m <- lapply(res, FUN = `[[`, "model")

  ## warnings
  list_warn <- lapply(
    seq_along(res),
    function(i) {

      ## remove duplicate warning messages for each response
      warn <- unique(res[[i]]$warning)

      ## return NULL when the model produced no warnings
      if (!length(warn)) {
        return(NULL)
      }

      ## associate each warning with its response variable
      data.frame(
        i = i,
        response = colnames(Y)[i],
        warning = warn,
        row.names = NULL
      )
    }
  )

  ## remove responses that produced no warnings
  list_warn <- Filter(Negate(is.null), list_warn)

  ## combine warning records into a single data frame
  ## and return an empty data frame when no warnings were generated
  if (length(list_warn)) {
    df_warn <- do.call(rbind, list_warn)
    rownames(df_warn) <- NULL
  } else {
    df_warn <- NULL
  }

  ## abiotic factors
  l_a <- lapply(seq_len(ncol(Y)),
                function(j) {

                  ## extract coefficients selected by the specified lambda
                  beta <- stats::coef(list_m[[j]], lambda.method)
                  nm <- rownames(beta)

                  ## retain coefficients corresponding to abiotic predictors
                  alpha <- beta[!(nm %in% colnames(Y))]
                  names(alpha) <- nm[!(nm %in% colnames(Y))]

                  alpha
                })

  m_a <- do.call(cbind, l_a)

  ## biotic factors (co-occurrence components)
  l_b <- lapply(seq_len(ncol(Y)),
                function(j) {

                  ## initialize coefficients for all biotic factors
                  b <- numeric(ncol(Y))

                  ## extract coefficients selected by the specified lambda
                  beta <- stats::coef(list_m[[j]], lambda.method)

                  ## exclude the response itself and retain coefficients
                  ## corresponding to other biotic factors
                  b[-j] <- beta[rownames(beta) %in% colnames(Y)]

                  b
                })

  ## enforce symmetry of pairwise biotic effects
  m_b <- do.call(cbind, l_b) |>
    symmetrize(
      method = sym.method,
      diagonal = FALSE
    )

  ## name dimensions
  colnames(m_a) <- colnames(Y)
  dimnames(m_b) <- list(
    colnames(Y),
    colnames(Y)
  )

  ## output list
  cout <- list(
    alpha = m_a,
    beta = m_b
  )

  ## attach model-fitting warnings as an attribute
  attr(cout, "warning") <- df_warn

  ## export
  cout
}


#' Calculate constant terms from coefficients and new predictor values
#'
#' Calculates constant terms for a given set of predictor values
#' using a single set of estimated model coefficients. Predictor columns in
#' \code{newdata} are matched to the coefficient names in \code{parm}.
#' An intercept, if present in \code{parm}, is included automatically.
#'
#' @param parm A matrix or vector of estimated model coefficients. Coefficient
#'   names must be provided as row names when \code{parm} is a matrix. An
#'   intercept should be named \code{"(Intercept)"}.
#' @param newdata A data frame containing new values of the predictors.
#'   Predictor names must match the row names of \code{parm}, excluding
#'   \code{"(Intercept)"}.
#'
#' @return A numeric matrix containing the calculated constant terms for
#'   each predictor set in \code{newdata}.
#'
#' @export

const <- function(
    parm,
    newdata
) {

  ## predictor names in coefficient vector
  pname <- setdiff(rownames(parm), "(Intercept)")

  ## check that all predictors are available
  missing <- setdiff(pname, colnames(newdata))

  if (length(missing) > 0) {
    stop(
      "Missing predictors in `newdata`: ",
      paste(missing, collapse = ", ")
    )
  }

  ## match newdata columns to coefficient order
  X <- stats::model.matrix(
    ~.,
    data = newdata[, pname, drop = FALSE]
  )

  ## calculate linear predictor
  as.matrix(X %*% parm)
}
