#' Calculate miRNA interaction gene rankings
#'
#' Builds a gene ranking based on the number of cancer-associated miRNA
#' interactions per gene. Uses a pre-computed interaction matrix (genes x
#' miRNAs) where each cell counts the number of targeting interactions.
#'
#' @param mirna_matrix A numeric matrix or data frame with genes as rows
#'   and miRNAs as columns, where values represent interaction counts.
#'   Typically produced by querying TargetScan for cancer-associated miRNAs.
#' @param normalize Logical; if `TRUE` (default), scores are normalized to
#'   sum to 1.
#'
#' @return A named numeric vector of miRNA interaction scores sorted in
#'   decreasing order. If `normalize = TRUE`, values sum to 1.
#' @export
#'
#' @examples
#' mat <- matrix(c(3, 0, 1, 2, 1, 0), nrow = 3,
#'   dimnames = list(c("TP53", "KRAS", "BRAF"), c("miR-1", "miR-2")))
#' calculate_mirna(mat)
calculate_mirna <- function(mirna_matrix, normalize = TRUE) {
  if (!is.matrix(mirna_matrix) && !is.data.frame(mirna_matrix)) {
    stop("mirna_matrix must be a matrix or data.frame")
  }
  if (nrow(mirna_matrix) == 0L || ncol(mirna_matrix) == 0L) {
    stop("mirna_matrix must not be empty")
  }

  scores <- rowSums(as.matrix(mirna_matrix))
  scores <- sort(scores, decreasing = TRUE)

  if (normalize) {
    total <- sum(abs(scores))
    if (total > 0) scores <- abs(scores) / total
  }

  scores
}
