
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
    state = m_mss[, seq_len(length(alpha)), drop = FALSE],

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
