#' Calculate SwitchDE gene rankings
#'
#' Ranks genes by switch-like differential expression along pseudotime
#' trajectories using the switchde package (Campbell & Yau, 2017).
#'
#' @param expression_matrix Numeric matrix with genes as rows and cells as
#'   columns.
#' @param pseudotime Numeric vector of pseudotime values (one per cell), or a
#'   data.frame with a `"Pseudotime"` column.
#' @param zero_inflated Logical; use zero-inflated model (default `FALSE`).
#' @param q_threshold Numeric; q-value cutoff for significance (default 0.05).
#' @param normalize Logical; if `TRUE` (default), scores are normalized to
#'   sum to 1.
#'
#' @return A named numeric vector of SwitchDE k values sorted by absolute
#'   magnitude. If `normalize = TRUE`, values are absolute and sum to 1.
#'   Returns `numeric(0)` if no genes pass the threshold.
#' @references Campbell KR, Yau C (2017). "switchde: inference of switch-like
#'   differential expression along single-cell trajectories."
#'   *Bioinformatics*, 33(8), 1241-1242.
#' @export
#'
#' @examples
#' \dontrun{
#' rankings <- calculate_switchde(expr_mat, pseudotime_vec)
#' }
calculate_switchde <- function(
    expression_matrix,
    pseudotime,
    zero_inflated = FALSE,
    q_threshold = 0.05,
    normalize = TRUE
) {
  if (!is.matrix(expression_matrix) || !is.numeric(expression_matrix)) {
    stop("expression_matrix must be a numeric matrix")
  }
  if (length(expression_matrix) == 0L) {
    stop("expression_matrix must not be empty")
  }

  # Accept data.frame with Pseudotime column
  if (is.data.frame(pseudotime)) {
    if ("Pseudotime" %in% colnames(pseudotime)) {
      pseudotime <- as.numeric(pseudotime$Pseudotime)
    } else {
      stop("pseudotime data.frame must contain a 'Pseudotime' column")
    }
  }
  if (!is.numeric(pseudotime)) {
    stop("pseudotime must be numeric")
  }
  if (length(pseudotime) != ncol(expression_matrix)) {
    stop("length(pseudotime) must equal ncol(expression_matrix)")
  }

  if (!requireNamespace("switchde", quietly = TRUE)) {
    stop("Package 'switchde' is required. Install from Bioconductor.")
  }

  sde <- switchde::switchde(
    expression_matrix,
    pseudotime,
    verbose = FALSE,
    zero_inflated = zero_inflated
  )

  # Filter by q-value
  sde <- sde[sde$qval < q_threshold, , drop = FALSE]

  if (nrow(sde) == 0L) {
    warning("No genes passed the q-value threshold of ", q_threshold)
    return(numeric(0))
  }

  # Sort by absolute k value descending
  sde <- sde[order(abs(sde$k), decreasing = TRUE), ]

  ranking <- sde$k
  names(ranking) <- sde$gene

  if (normalize) {
    total <- sum(abs(ranking))
    if (total > 0) ranking <- abs(ranking) / total
  }

  ranking
}
