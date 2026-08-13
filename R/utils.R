#' Check dimensions and validity of model inputs
#'
#' Internal utility for validating `state`, `alpha`, and `beta`.
#'
#' @param state An optional binary numeric vector.
#' @param alpha A numeric vector of model parameters.
#' @param beta A square numeric matrix of pairwise interaction parameters.
#'
#' @return Invisibly returns `NULL` if all checks pass.
#'
#' @noRd

check_dim <- function(state = NULL, alpha, beta) {

  ## check `state`
  if (!is.null(state)) {
    if (!is.numeric(state) || !all(state %in% c(0, 1))) {
      stop("`state` must be a binary numeric vector.")
    }

    s <- length(state)

    if (length(alpha) != s) {
      stop("`state` and `alpha` must have the same length.")
    }

  } else {
    s <- length(alpha)
  }

  ## check `beta`
  if (!is.matrix(beta) || !all(dim(beta) == c(s, s))) {
    stop(
      "`beta` must be a square matrix with dimensions equal to the length of",
      ifelse(is.null(state), "`alpha`", "`state` and `alpha`"),
      "."
    )
  }

  if (anyNA(beta)) {
    stop("`beta` must not contain missing values.")
  }

  if (any(diag(beta) != 0)) {
    stop("The diagonal elements of `beta` must be zero.")
  }

  invisible(NULL)
}
