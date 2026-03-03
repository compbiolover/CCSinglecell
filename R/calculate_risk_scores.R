#' Calculate patient risk scores from Cox model coefficients
#'
#' Binarises gene expression relative to per-gene medians, weights by
#' coefficient sign, and sums across genes to produce a risk score per
#' patient. Patients are then stratified into "high" and "low" risk groups.
#'
#' @param cox_df Data frame with gene expression columns plus
#'   `days_to_last_follow_up` and `vital_status`.
#' @param active_genes Character vector of gene names selected by the Cox
#'   model.
#' @param coefficients Named numeric vector of Cox model coefficients
#'   (same order / names as `active_genes`).
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{risk_score}{Numeric sum of binarised expression x coefficient
#'       sign.}
#'     \item{risk_group}{`"high"` or `"low"` relative to the median absolute
#'       risk.}
#'     \item{vital_status}{From `cox_df`.}
#'     \item{time}{Survival time from `cox_df`.}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' scores <- calculate_risk_scores(patient_data, model$active_genes, model$coefficients)
#' table(scores$risk_group)
#' }
calculate_risk_scores <- function(cox_df, active_genes, coefficients) {
  if (!is.data.frame(cox_df)) stop("cox_df must be a data.frame")
  if (length(active_genes) == 0L) stop("active_genes must not be empty")

  # Sanitise column names to match fit_cox_model convention
  clean <- function(x) gsub("[-_/]", ".", x)
  active_genes <- clean(active_genes)
  colnames(cox_df) <- clean(colnames(cox_df))
  names(coefficients) <- clean(names(coefficients))

  missing_genes <- setdiff(active_genes, colnames(cox_df))
  if (length(missing_genes) > 0L) {
    stop("Genes not found in cox_df: ", paste(head(missing_genes, 5), collapse = ", "))
  }

  expr_mat <- as.matrix(cox_df[, active_genes, drop = FALSE])

  # Propagate coefficient sign via matrix multiplication
  gene_sign <- ifelse(coefficients[active_genes] > 0, 1, -1)
  expr_signed <- expr_mat %*% diag(gene_sign)
  colnames(expr_signed) <- active_genes

  # Binarise: compare each gene to its column median
  col_medians <- apply(expr_signed, 2L, stats::median)
  binarised <- vapply(seq_along(active_genes), function(j) {
    med <- col_medians[j]
    if (med > 0) {
      ifelse(expr_signed[, j] > med, 1, 0)
    } else {
      ifelse(expr_signed[, j] < med, -1, 0)
    }
  }, numeric(nrow(expr_signed)))

  risk_score <- rowSums(binarised)
  median_risk <- stats::median(abs(binarised))

  data.frame(
    risk_score = risk_score,
    risk_group = ifelse(risk_score > median_risk, "high", "low"),
    vital_status = cox_df$vital_status,
    time = cox_df$days_to_last_follow_up,
    row.names = rownames(cox_df),
    stringsAsFactors = FALSE
  )
}
