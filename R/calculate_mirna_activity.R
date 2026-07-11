#' Calculate expression-anchored miRNA activity gene rankings
#'
#' An expression-anchored alternative to [calculate_mirna()]. Instead of
#' counting how many cancer miRNAs are *predicted* to target a gene, this ranks
#' genes by how strongly they are *observed* to be repressed by their targeting
#' miRNAs in the user's own samples. For each predicted gene-miRNA interaction,
#' the score accumulates the anti-correlation between the gene's and the
#' miRNA's expression across matched samples (the classic miRNA-repression /
#' competing-endogenous-RNA signal). Genes whose targeting miRNAs actually move
#' opposite to them score highest.
#'
#' @param mirna_expr Numeric matrix of miRNA expression (miRNAs as rows,
#'   samples as columns). Must have row and column names.
#' @param gene_expr Numeric matrix of gene expression (genes as rows, samples
#'   as columns). Must have row and column names. Samples are matched to
#'   `mirna_expr` by column name.
#' @param target_matrix Predicted gene-miRNA interactions (genes as rows,
#'   miRNAs as columns), e.g. from TargetScan for cancer-associated miRNAs.
#'   Non-zero entries mark candidate interactions; their magnitude is used as a
#'   weight. This is the same kind of matrix consumed by [calculate_mirna()].
#' @param method Either `"anticorrelation"` (default; only repressive, i.e.
#'   negative, correlations contribute evidence) or `"signed"` (positive
#'   correlations subtract, capturing net regulatory direction).
#' @param normalize Logical; if `TRUE` (default), scores are divided by the sum
#'   of their absolute values. For `method = "anticorrelation"` scores are
#'   non-negative and therefore sum to 1; for `method = "signed"` scores may be
#'   negative, so they are put on a comparable scale but are not guaranteed to
#'   sum to 1.
#'
#' @return A named numeric vector of miRNA-activity scores sorted in decreasing
#'   order. When `normalize = TRUE`, non-negative (`"anticorrelation"`) scores
#'   sum to 1; signed scores are scaled by the sum of absolute values and may be
#'   negative.
#' @export
#'
#' @examples
#' set.seed(1)
#' samples <- paste0("s", 1:20)
#' genes <- paste0("g", 1:5)
#' mirs <- paste0("mir", 1:3)
#' gene_expr <- matrix(rnorm(5 * 20), nrow = 5, dimnames = list(genes, samples))
#' mirna_expr <- matrix(rnorm(3 * 20), nrow = 3, dimnames = list(mirs, samples))
#' # Make g1 strongly anti-correlated with mir1
#' mirna_expr["mir1", ] <- -gene_expr["g1", ] + rnorm(20, sd = 0.1)
#' target_matrix <- matrix(1, nrow = 5, ncol = 3, dimnames = list(genes, mirs))
#' calculate_mirna_activity(mirna_expr, gene_expr, target_matrix)
calculate_mirna_activity <- function(
    mirna_expr,
    gene_expr,
    target_matrix,
    method = c("anticorrelation", "signed"),
    normalize = TRUE
) {
  method <- match.arg(method)

  for (nm in c("mirna_expr", "gene_expr")) {
    m <- get(nm)
    if (!is.matrix(m) || !is.numeric(m)) {
      stop(nm, " must be a numeric matrix")
    }
    if (is.null(rownames(m)) || is.null(colnames(m))) {
      stop(nm, " must have row names (features) and column names (samples)")
    }
  }
  if (!is.matrix(target_matrix) && !is.data.frame(target_matrix)) {
    stop("target_matrix must be a matrix or data.frame (genes x miRNAs)")
  }
  target_matrix <- as.matrix(target_matrix)
  if (!is.numeric(target_matrix)) {
    stop("target_matrix must be numeric (genes x miRNAs interaction weights); ",
         "did a non-numeric column slip in, e.g. gene names not set as row ",
         "names when reading a CSV?")
  }
  if (is.null(rownames(target_matrix)) || is.null(colnames(target_matrix))) {
    stop("target_matrix must have gene row names and miRNA column names")
  }

  # Match samples across the two expression matrices.
  common_samples <- intersect(colnames(gene_expr), colnames(mirna_expr))
  if (length(common_samples) < 3L) {
    stop("mirna_expr and gene_expr must share at least 3 samples (by column ",
         "name); found ", length(common_samples))
  }
  ge <- gene_expr[, common_samples, drop = FALSE]
  me <- mirna_expr[, common_samples, drop = FALSE]

  # Restrict to the genes/miRNAs that appear in both the interaction map and
  # the expression matrices.
  genes <- intersect(rownames(target_matrix), rownames(ge))
  mirs <- intersect(colnames(target_matrix), rownames(me))
  if (length(genes) == 0L || length(mirs) == 0L) {
    stop("target_matrix shares no genes/miRNAs with the expression matrices")
  }
  ge <- ge[genes, , drop = FALSE]
  me <- me[mirs, , drop = FALSE]
  weights <- abs(target_matrix[genes, mirs, drop = FALSE])

  # Correlation between every gene and every miRNA across samples.
  cor_mat <- suppressWarnings(stats::cor(t(ge), t(me)))
  cor_mat[is.na(cor_mat)] <- 0 # constant rows -> undefined -> no evidence

  repression <- if (method == "anticorrelation") {
    pmax(0, -cor_mat)
  } else {
    -cor_mat
  }

  # Only predicted interactions contribute, weighted by interaction strength.
  scores <- rowSums(repression * weights)
  scores <- sort(scores, decreasing = TRUE)

  if (normalize) {
    total <- sum(abs(scores))
    if (total > 0) scores <- scores / total
  }

  scores
}
