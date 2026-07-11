#' Score an arbitrary set of genes
#'
#' High-level wrapper that scores a user-supplied set of genes using one or
#' more of the package's ranking metrics and combines them into a single
#' score. This is the main entry point for the command-line tool and for
#' interactive "score my genes" workflows.
#'
#' Each metric is computed across the full `expression_matrix` (so scores
#' reflect a genome-wide ranking), then restricted to the requested `genes`.
#' The available metrics are:
#'
#' \describe{
#'   \item{`mad`}{Expression variability via [calculate_mad()]. Always
#'     computed.}
#'   \item{`switchde`}{Switch-like dynamics along a trajectory via
#'     [calculate_switchde()]. Computed only when `pseudotime` is supplied
#'     and the \pkg{switchde} package is available.}
#'   \item{`mirna`}{Cancer-miRNA targeting via [calculate_mirna()]. Computed
#'     only when `mirna_matrix` is supplied.}
#' }
#'
#' Present metrics are blended with [combine_rankings()] using `weights`
#' (equal weights by default). When only one metric is available the combined
#' score equals that metric.
#'
#' @param genes A character vector of gene symbols, or a path to a gene file
#'   (see [read_gene_list()]).
#' @param expression_matrix Numeric matrix with genes as rows and cells as
#'   columns. Must have row names (gene identifiers).
#' @param pseudotime Optional numeric vector of pseudotime values, one per
#'   cell (length `ncol(expression_matrix)`). Enables the `switchde` metric.
#' @param mirna_matrix Optional numeric matrix/data frame (genes x miRNAs) of
#'   interaction counts. Enables the `mirna` metric.
#' @param weights Optional numeric vector of metric weights, one per available
#'   metric, in the order (mad, switchde, mirna) for the metrics that are
#'   present. Defaults to equal weights. Rescaled to sum to 1.
#' @param renormalize Logical; if `TRUE` (default) each metric's scores are
#'   renormalised to sum to 1 *across the requested gene set* before
#'   combining, so the combined score is a proportion within the set.
#'
#' @return A data frame with one row per requested gene, sorted by `combined`
#'   score (decreasing), containing:
#'   \describe{
#'     \item{gene}{Gene symbol.}
#'     \item{one column per available metric}{`mad`, and optionally
#'       `switchde` / `mirna` — the (renormalised) per-metric score.}
#'     \item{combined}{Weighted combination of the available metrics.}
#'     \item{rank}{Integer rank (1 = highest combined score).}
#'   }
#'   The `weights` used and the metric names are attached as attributes
#'   `"weights"` and `"metrics"`.
#' @export
#'
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(200), nrow = 20,
#'   dimnames = list(paste0("g", 1:20), paste0("c", 1:10)))
#' score_gene_set(c("g1", "g5", "g10"), mat)
score_gene_set <- function(
    genes,
    expression_matrix,
    pseudotime = NULL,
    mirna_matrix = NULL,
    weights = NULL,
    renormalize = TRUE
) {
  genes <- read_gene_list(genes)
  if (length(genes) == 0L) stop("`genes` is empty")

  if (!is.matrix(expression_matrix) || !is.numeric(expression_matrix)) {
    stop("expression_matrix must be a numeric matrix")
  }
  if (is.null(rownames(expression_matrix))) {
    stop("expression_matrix must have row names (gene identifiers)")
  }

  # --- Compute each available metric genome-wide -----------------------------
  rankings <- list()

  rankings$mad <- calculate_mad(expression_matrix, normalize = TRUE)

  if (!is.null(pseudotime)) {
    if (!requireNamespace("switchde", quietly = TRUE)) {
      warning(
        "`pseudotime` supplied but the 'switchde' package is not installed; ",
        "skipping the switchde metric."
      )
    } else {
      rankings$switchde <- calculate_switchde(
        expression_matrix, pseudotime, normalize = TRUE
      )
    }
  }

  if (!is.null(mirna_matrix)) {
    rankings$mirna <- calculate_mirna(mirna_matrix, normalize = TRUE)
  }

  metrics <- names(rankings)

  # --- Restrict every metric to the requested gene set -----------------------
  missing_all <- setdiff(genes, unique(unlist(lapply(rankings, names))))
  if (length(missing_all) > 0L) {
    warning(
      length(missing_all), " of ", length(genes),
      " requested gene(s) were not found in any metric and scored 0 ",
      "(e.g. ", paste(utils::head(missing_all, 3L), collapse = ", "), ")."
    )
  }

  subset_metric <- function(r) {
    v <- r[genes]
    v[is.na(v)] <- 0
    names(v) <- genes
    if (renormalize) {
      total <- sum(v)
      if (total > 0) v <- v / total
    }
    v
  }
  rankings <- lapply(rankings, subset_metric)

  # --- Resolve weights -------------------------------------------------------
  n_metrics <- length(rankings)
  if (is.null(weights)) {
    weights <- rep(1 / n_metrics, n_metrics)
  } else {
    if (length(weights) != n_metrics) {
      stop(
        "`weights` must have one value per available metric (",
        n_metrics, ": ", paste(metrics, collapse = ", "), ")"
      )
    }
    if (any(weights < 0)) stop("`weights` must be non-negative")
    total <- sum(weights)
    if (total <= 0) stop("`weights` must not all be zero")
    weights <- weights / total
  }

  # --- Combine ---------------------------------------------------------------
  if (n_metrics == 1L) {
    combined <- rankings[[1L]]
  } else {
    combined <- combine_rankings(rankings, weights)
  }

  # --- Assemble tidy output --------------------------------------------------
  out <- data.frame(gene = genes, stringsAsFactors = FALSE)
  for (m in metrics) {
    out[[m]] <- unname(rankings[[m]][genes])
  }
  out$combined <- unname(combined[genes])
  out <- out[order(out$combined, decreasing = TRUE), , drop = FALSE]
  out$rank <- seq_len(nrow(out))
  rownames(out) <- NULL

  attr(out, "weights") <- stats::setNames(weights, metrics)
  attr(out, "metrics") <- metrics
  out
}
