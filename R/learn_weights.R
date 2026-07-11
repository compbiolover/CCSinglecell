#' Learn integration weights from the data
#'
#' Learns how much weight to give each metric / omics layer when blending their
#' gene rankings, instead of grid-searching with [optimize_weights()]. The
#' default `"pca"` method is unsupervised and dependency-free: it runs a
#' principal component analysis on the gene-by-metric score matrix and takes the
#' first principal component — the dominant axis of shared variation across
#' metrics — as the integration axis. Each metric's weight is its (absolute)
#' loading on that axis, normalised to sum to 1. Metrics that co-vary with the
#' shared signal get more weight; a flat or idiosyncratic metric gets less.
#'
#' This is a lightweight, always-available stand-in for full latent-variable
#' integration (MOFA) and supervised multi-block selection (DIABLO); see
#' `ROADMAP.md` for those planned backends.
#'
#' @param rankings A **named** list of named numeric vectors (one per metric /
#'   omics layer), as consumed by [combine_rankings()].
#' @param method Integration method. Currently `"pca"` (unsupervised first
#'   principal component).
#'
#' @return A named numeric vector of non-negative weights summing to 1, one per
#'   metric (in the order of `rankings`). The proportion of variance explained
#'   by the integration axis is attached as attribute `"variance_explained"`.
#' @export
#'
#' @examples
#' r1 <- c(g1 = 0.4, g2 = 0.3, g3 = 0.2, g4 = 0.1)
#' r2 <- c(g1 = 0.35, g2 = 0.30, g3 = 0.25, g4 = 0.10) # correlated with r1
#' r3 <- c(g1 = 0.1, g2 = 0.4, g3 = 0.1, g4 = 0.4)     # different structure
#' learn_weights(list(mad = r1, switchde = r2, mirna = r3))
learn_weights <- function(rankings, method = c("pca")) {
  method <- match.arg(method)

  if (!is.list(rankings) || length(rankings) == 0L) {
    stop("`rankings` must be a non-empty named list of gene rankings")
  }
  if (is.null(names(rankings)) || any(!nzchar(names(rankings)))) {
    stop("`rankings` must be named (one name per metric / omics layer)")
  }
  metrics <- names(rankings)
  if (anyDuplicated(metrics)) {
    stop("metric names must be unique; got: ", paste(metrics, collapse = ", "))
  }
  for (i in seq_along(rankings)) {
    r <- rankings[[i]]
    if (!is.numeric(r) || is.null(names(r))) {
      stop("ranking '", metrics[i], "' must be a named numeric vector")
    }
    nm <- names(r)
    if (any(!nzchar(nm)) || anyNA(nm)) {
      stop("ranking '", metrics[i], "' has empty or missing gene names")
    }
    if (anyDuplicated(nm)) {
      stop("ranking '", metrics[i], "' has duplicated gene names")
    }
  }
  n <- length(rankings)

  # One metric: it gets all the weight.
  if (n == 1L) {
    w <- stats::setNames(1, metrics)
    attr(w, "variance_explained") <- NA_real_
    return(w)
  }

  # Build the gene-by-metric score matrix (union of genes, 0 for missing).
  genes <- unique(unlist(lapply(rankings, names)))
  mat <- vapply(rankings, function(r) {
    v <- r[genes]
    v[is.na(v)] <- 0
    v
  }, numeric(length(genes)))
  colnames(mat) <- metrics

  # Metrics with no variance across genes carry no information: weight 0.
  col_sd <- apply(mat, 2L, stats::sd)
  informative <- col_sd > 0

  if (sum(informative) < 2L) {
    # Not enough varying metrics to learn structure: fall back to equal weights
    # across informative metrics (or all, if none vary).
    keep <- if (any(informative)) informative else rep(TRUE, n)
    w <- ifelse(keep, 1, 0)
    w <- w / sum(w)
    w <- stats::setNames(w, metrics)
    attr(w, "variance_explained") <- NA_real_
    return(w)
  }

  pca <- stats::prcomp(mat[, informative, drop = FALSE],
                       center = TRUE, scale. = TRUE)
  # Weight each informative metric by its absolute loading on the first
  # principal component; the sign of the axis does not affect the magnitude.
  loading <- pca$rotation[, 1L]

  weights <- stats::setNames(numeric(n), metrics)
  weights[informative] <- abs(loading)
  total <- sum(weights)
  if (total > 0) weights <- weights / total

  var_explained <- (pca$sdev[1L]^2) / sum(pca$sdev^2)
  attr(weights, "variance_explained") <- var_explained
  weights
}
