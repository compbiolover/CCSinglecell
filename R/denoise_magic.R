#' Denoise single-cell RNA-seq data with MAGIC
#'
#' Applies the MAGIC (Markov Affinity-based Graph Imputation of Cells)
#' algorithm to denoise a single-cell expression matrix.
#'
#' @param expression_matrix Numeric matrix with genes as rows and cells as
#'   columns.
#' @param seed Integer random seed for reproducibility (default 123).
#' @param knn Integer number of nearest neighbours (default 10).
#' @param solver Character; `"approximate"` (default) or `"exact"`.
#' @param n_jobs Integer number of parallel jobs (default 1).
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{denoised}{Numeric matrix (genes x cells) of denoised values.}
#'     \item{gene_metadata}{Data frame with `gene_short_name` column, suitable
#'       for Monocle3 cell data set construction.}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' result <- denoise_magic(raw_counts)
#' denoised_mat <- result$denoised
#' }
denoise_magic <- function(
    expression_matrix,
    seed = 123L,
    knn = 10L,
    solver = "approximate",
    n_jobs = 1L
) {
  if (!is.matrix(expression_matrix) || !is.numeric(expression_matrix)) {
    stop("expression_matrix must be a numeric matrix")
  }
  if (length(expression_matrix) == 0L) {
    stop("expression_matrix must not be empty")
  }
  if (!requireNamespace("Rmagic", quietly = TRUE)) {
    stop("Package 'Rmagic' is required. Install with: install.packages('Rmagic')")
  }

  # MAGIC expects cells x genes
  denoised <- Rmagic::magic(
    t(expression_matrix),
    seed = seed,
    solver = solver,
    knn = knn,
    n.jobs = n_jobs,
    verbose = FALSE
  )

  denoised_mat <- t(as.matrix(denoised[["result"]]))

  gene_metadata <- data.frame(
    gene_short_name = rownames(denoised_mat),
    row.names = rownames(denoised_mat),
    stringsAsFactors = FALSE
  )

  list(denoised = denoised_mat, gene_metadata = gene_metadata)
}
