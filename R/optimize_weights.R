#' Combine gene rankings with a weighted linear combination
#'
#' Merges any number of named numeric ranking vectors using weights that
#' sum to 1. Genes present in any input but absent from another receive a
#' score of 0 for the missing metric.
#'
#' @param rankings A list of two or more named numeric vectors (gene
#'   rankings), typically one per omics layer / metric.
#' @param weights Numeric vector of the same length as `rankings`. Values
#'   should be non-negative; they are used as-is (caller controls the
#'   constraint that they sum to 1).
#'
#' @return A named numeric vector of combined scores sorted in decreasing
#'   order.
#' @export
#'
#' @examples
#' r1 <- c(A = 0.5, B = 0.3, C = 0.2)
#' r2 <- c(B = 0.6, C = 0.3, D = 0.1)
#' combine_rankings(list(r1, r2), weights = c(0.7, 0.3))
combine_rankings <- function(rankings, weights) {
  if (!is.list(rankings) || length(rankings) < 2L) {
    stop("rankings must be a list of two or more named numeric vectors")
  }
  if (length(weights) != length(rankings)) {
    stop("weights must have the same length as rankings")
  }

  all_genes <- unique(unlist(lapply(rankings, names)))
  scores <- numeric(length(all_genes))
  names(scores) <- all_genes

  for (i in seq_along(rankings)) {
    matched <- match(names(rankings[[i]]), all_genes)
    scores[matched] <- scores[matched] + weights[i] * rankings[[i]]
  }

  sort(scores, decreasing = TRUE)
}


#' Grid-search weight optimisation for gene rankings
#'
#' Evaluates weight combinations on a simplex grid (weights non-negative and
#' summing to 1, in steps of `step_size`) for any number of ranking metrics.
#' For two metrics the weights are `(w, 1 - w)`; for three metrics all valid
#' triplets `(w1, w2, w3)` are tested; and the same generalises to any number
#' of metrics.
#'
#' @param rankings A list of two or more named numeric vectors.
#' @param step_size Numeric step between 0 and 1 for the weight grid
#'   (default 0.1).
#' @param max_combos Integer safety cap on the number of weight combinations
#'   (default 20000). Exceeding it is an error asking for a larger
#'   `step_size`, since a fine grid over many metrics explodes combinatorially.
#'
#' @return A data frame with one weight column per metric (`w1`, `w2`, ...)
#'   and a `ranking` list-column of the combined named numeric vectors.
#' @export
#'
#' @examples
#' r1 <- c(A = 0.5, B = 0.3, C = 0.2)
#' r2 <- c(B = 0.6, C = 0.3, D = 0.1)
#' grid <- optimize_weights(list(r1, r2), step_size = 0.5)
optimize_weights <- function(rankings, step_size = 0.1, max_combos = 20000L) {
  n <- length(rankings)
  if (!is.list(rankings) || n < 2L) {
    stop("rankings must contain two or more vectors")
  }
  if (step_size <= 0 || step_size > 1) {
    stop("step_size must be in (0, 1]")
  }

  grid <- simplex_grid(n, step_size)
  if (nrow(grid) > max_combos) {
    stop(
      "Weight grid has ", nrow(grid), " combinations for ", n, " metrics, ",
      "exceeding max_combos (", max_combos, "). Use a larger step_size."
    )
  }

  results <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    results[[i]] <- combine_rankings(rankings, grid[i, ])
  }

  out <- as.data.frame(grid)
  colnames(out) <- paste0("w", seq_len(n))
  out$ranking <- I(results)
  out
}


#' Enumerate a simplex grid of weights
#'
#' Returns all length-`n` weight vectors whose entries are non-negative
#' multiples of `step_size` and sum to 1. Used by [optimize_weights()].
#'
#' @param n Number of metrics (>= 1).
#' @param step_size Grid step in (0, 1]; `1 / step_size` should be (close to)
#'   an integer.
#'
#' @return A numeric matrix with `n` columns; each row sums to 1.
#' @keywords internal
simplex_grid <- function(n, step_size) {
  k <- round(1 / step_size)
  comps <- integer_compositions(k, n)
  comps / k
}

# All ways to write `total` as an ordered sum of `parts` non-negative
# integers, returned as a matrix (one composition per row).
integer_compositions <- function(total, parts) {
  if (parts == 1L) {
    return(matrix(total, nrow = 1L))
  }
  rows <- list()
  for (i in 0:total) {
    rest <- integer_compositions(total - i, parts - 1L)
    rows[[length(rows) + 1L]] <- cbind(i, rest, deparse.level = 0L)
  }
  do.call(rbind, rows)
}
