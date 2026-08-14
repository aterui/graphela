#' Check dimensions and validity of model inputs
#'
#' Internal utility for validating the dimensions and values of `state`,
#' `alpha`, and `beta`. The function accepts either a single state vector
#' or a matrix containing multiple states.
#'
#' @param state An optional binary numeric vector or matrix. If supplied as
#'   a vector, its length must equal the number of elements in `alpha`.
#'   If supplied as a matrix, the number of columns must equal the number
#'   of elements in `alpha`.
#' @param alpha A numeric vector of model parameters controlling the
#'   intrinsic contribution of each species to system energy.
#' @param beta A square numeric matrix of pairwise interaction parameters
#'   among species. Its dimensions must equal `length(alpha)` by
#'   `length(alpha)`.
#'
#' @return Invisibly returns the number of state variables, i.e.,
#'   `length(alpha)`, if all checks pass.
#'
#' @noRd

check_dim <- function(
    state = NULL,
    alpha,
    beta
) {

  if (!is.numeric(alpha) || is.null(alpha))
    stop("`alpha` must be a numeric vector.")

  if (!is.matrix(beta) || !is.numeric(beta))
    stop("`beta` must be a numeric matrix.")

  if (!all(diag(beta) == 0))
    stop("The diagonal elements of `beta` must be zero")

  if (!isTRUE(all.equal(beta, t(beta), tolerance = 1e-10)))
    stop("`beta` must be symmetric.")

  s <- length(alpha)

  if (!all(dim(beta) == c(s, s)))
    stop("`beta` must have dimensions `length(alpha)` x `length(alpha)`.")

  if (!is.null(state)) {

    if (!is.numeric(state))
      stop("`state` must be numeric.")

    if (is.null(dim(state))) {
      ## vector
      if (length(state) != s)
        stop("`state` must have the same length as `alpha`.")

    } else {
      ## matrix
      if (ncol(state) != s)
        stop("The state matrix must have `length(alpha)` columns.")
    }

    if (anyNA(state) || !all(state %in% c(0, 1)))
      stop("`state` must contain only 0 and 1.")
  }

  invisible(s)
}


#' Check stable-state and energy matrix
#'
#' Internal utility for validating a matrix containing binary stable-state
#' configurations and their corresponding energy values. The state
#' variables are expected in all columns except the last, which must
#' contain energy values.
#'
#' @param sse A matrix containing binary state configurations in the first
#'   columns and the corresponding energy values in the last column.
#' @param alpha A numeric vector of model parameters controlling the
#'   intrinsic contribution of each species to system energy.
#' @param beta A square numeric matrix of pairwise interaction parameters
#'   among species.
#' @param min_rows An integer specifying the minimum number of rows
#'   required in `sse`. Defaults to `1`.
#'
#' @return Invisibly returns the number of state variables, i.e.,
#'   `length(alpha)`, if all checks pass.
#'
#' @noRd

check_sse <- function(
    sse,
    alpha,
    beta,
    min_rows = 2
) {

  if (!is.matrix(sse))
    stop("The stable state input must be a matrix.")

  if (nrow(sse) < min_rows)
    stop("The stable state matrix must contain at least ", min_rows, " row(s).")

  s <- check_dim(
    state = sse[, -ncol(sse), drop = FALSE],
    alpha = alpha,
    beta = beta
  )

  if (ncol(sse) != s + 1)
    stop("The stable state matrix must contain `length(alpha) + 1` columns.")

  if (!is.numeric(sse[, ncol(sse)]))
    stop("The last column of the stable state matrix must contain numeric energy values.")

  if (anyNA(sse) || any(!is.finite(sse[, ncol(sse)])))
    stop("The stable state matrix cannot contain missing or non-finite values.")

  invisible(s)
}


#' Check simulated annealing parameters
#'
#' Internal utility for validating the parameters controlling simulated
#' annealing. Checks that the initial temperature and cooling rate are
#' positive finite numeric values and that the number of iterations is a
#' positive integer.
#'
#' @param temp A positive numeric value specifying the initial temperature
#'   of simulated annealing.
#' @param r A positive numeric value specifying the cooling rate of
#'   simulated annealing.
#' @param iter A positive integer specifying the number of simulated
#'   annealing iterations.
#'
#' @return Invisibly returns `NULL` if all checks pass.
#'
#' @noRd

check_sa <- function(temp, r, iter) {

  if (!is.numeric(temp) || length(temp) != 1 ||
      !is.finite(temp) || temp <= 0)
    stop("`temp` must be a positive numeric value.")

  if (!is.numeric(r) || length(r) != 1 ||
      !is.finite(r) || r <= 0)
    stop("`r` must be a positive numeric value.")

  if (length(iter) != 1 ||
      !is.numeric(iter) ||
      !is.finite(iter) ||
      iter < 1 ||
      iter != as.integer(iter))
    stop("`iter` must be a positive integer.")

  invisible(NULL)
}
