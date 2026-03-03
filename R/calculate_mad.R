#' Calculate Median Absolute Deviation (MAD) gene rankings
#'
#' Ranks genes by their expression variability across cells using the median
#' absolute deviation. Higher MAD indicates more variable expression.
#'
#' @param expression_matrix Numeric matrix with genes as rows and cells as
#'   columns. Must have row names (gene identifiers).
#' @param normalize Logical; if `TRUE` (default), scores are normalized to
#'   sum to 1.
#'
#' @return A named numeric vector of MAD values sorted in decreasing order.
#'   If `normalize = TRUE`, values sum to 1.
#' @export
#'
#' @examples
#' mat <- matrix(rnorm(30), nrow = 3, dimnames = list(paste0("g", 1:3), NULL))
#' calculate_mad(mat)
calculate_mad <- function(expression_matrix, normalize = TRUE) {
  if (!is.matrix(expression_matrix) || !is.numeric(expression_matrix)) {
    stop("expression_matrix must be a numeric matrix")
  }
  if (length(expression_matrix) == 0L) {
    stop("expression_matrix must not be empty")
  }

  gene_mads <- apply(expression_matrix, 1L, stats::mad)
  gene_mads <- sort(gene_mads, decreasing = TRUE)

  if (normalize) {
    total <- sum(abs(gene_mads))
    if (total > 0) {
      gene_mads <- gene_mads / total
    }
  }

  gene_mads
}
