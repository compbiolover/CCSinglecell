#' Combine gene rankings with weighted linear combination
#'
#' Merges two or three named numeric ranking vectors using weights that
#' sum to 1. Genes present in any input but absent from another receive a
#' score of 0 for the missing metric.
#'
#' @param rankings A list of 2 or 3 named numeric vectors (gene rankings).
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
  if (!is.list(rankings) || length(rankings) < 2L || length(rankings) > 3L) {
    stop("rankings must be a list of 2 or 3 named numeric vectors")
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
#' Evaluates all weight combinations (in steps of `step_size`) for 2 or 3
#' ranking metrics. For 2 metrics the weights are `(w, 1-w)`. For 3 metrics
#' all triplets `(w1, w2, w3)` with `w1 + w2 + w3 ≈ 1` are tested.
#'
#' @param rankings A list of 2 or 3 named numeric vectors.
#' @param step_size Numeric step between 0 and 1 for weight grid
#'   (default 0.1).
#'
#' @return A data frame with columns `w1`, `w2`, (optionally `w3`), and
#'   `ranking` (a list-column of named numeric vectors).
#' @export
#'
#' @examples
#' r1 <- c(A = 0.5, B = 0.3, C = 0.2)
#' r2 <- c(B = 0.6, C = 0.3, D = 0.1)
#' grid <- optimize_weights(list(r1, r2), step_size = 0.5)
optimize_weights <- function(rankings, step_size = 0.1) {
  n <- length(rankings)
  if (n < 2L || n > 3L) stop("rankings must contain 2 or 3 vectors")

  steps <- seq(0, 1, by = step_size)

  if (n == 2L) {
    results <- vector("list", length(steps))
    w1_out <- numeric(length(steps))
    w2_out <- numeric(length(steps))

    for (i in seq_along(steps)) {
      w <- steps[i]
      w1_out[i] <- w
      w2_out[i] <- 1 - w
      results[[i]] <- combine_rankings(rankings, c(w, 1 - w))
    }

    data.frame(w1 = w1_out, w2 = w2_out, ranking = I(results))
  } else {
    combos <- expand.grid(w1 = steps, w2 = steps)
    combos$w3 <- 1 - combos$w1 - combos$w2
    combos <- combos[combos$w3 >= -1e-10 & combos$w3 <= 1 + 1e-10, ]
    combos$w3 <- pmax(combos$w3, 0)
    rownames(combos) <- NULL

    results <- vector("list", nrow(combos))
    for (i in seq_len(nrow(combos))) {
      w <- as.numeric(combos[i, 1:3])
      results[[i]] <- combine_rankings(rankings, w)
    }

    combos$ranking <- I(results)
    combos
  }
}
