#' Fit a penalized Cox proportional hazards model
#'
#' Selects survival-predictive genes using elastic-net regularized Cox
#' regression via [glmnet::cv.glmnet()].
#'
#' @param cox_df Data frame containing gene expression columns plus
#'   `days_to_last_follow_up` (numeric) and `vital_status` (0/1).
#' @param gene_names Character vector of gene column names to use as
#'   predictors. Genes absent from `cox_df` are silently dropped.
#' @param alpha Elastic-net mixing parameter: 1 = lasso (default),
#'   0 = ridge.
#' @param nfolds Integer number of cross-validation folds (default 10).
#' @param seed Integer random seed for reproducibility (default 1).
#' @param max_genes Maximum number of genes to include (default `NULL`,
#'   uses all).
#'
#' @return A list with:
#'   \describe{
#'     \item{cv_fit}{The [glmnet::cv.glmnet()] object.}
#'     \item{active_genes}{Character vector of genes with non-zero
#'       coefficients.}
#'     \item{coefficients}{Named numeric vector of non-zero coefficients.}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' result <- fit_cox_model(patient_data, gene_names = top_genes)
#' result$active_genes
#' }
fit_cox_model <- function(
    cox_df,
    gene_names,
    alpha = 1,
    nfolds = 10L,
    seed = 1L,
    max_genes = NULL
) {
  if (!is.data.frame(cox_df)) stop("cox_df must be a data.frame")
  if (!is.character(gene_names) || length(gene_names) == 0L) {
    stop("gene_names must be a non-empty character vector")
  }

  required <- c("days_to_last_follow_up", "vital_status")
  missing_cols <- setdiff(required, colnames(cox_df))
  if (length(missing_cols) > 0L) {
    stop("cox_df is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Sanitise gene names (match legacy behaviour: replace - _ / with .)
  clean <- function(x) gsub("[-_/]", ".", x)
  gene_names <- clean(gene_names)
  colnames(cox_df) <- clean(colnames(cox_df))

  # Subset to available genes
  gene_names <- intersect(gene_names, colnames(cox_df))
  if (length(gene_names) == 0L) stop("No gene_names found in cox_df columns")

  if (!is.null(max_genes) && length(gene_names) > max_genes) {
    gene_names <- gene_names[seq_len(max_genes)]
  }

  x <- as.matrix(cox_df[, gene_names, drop = FALSE])
  y <- survival::Surv(
    time = cox_df$days_to_last_follow_up,
    event = cox_df$vital_status
  )

  set.seed(seed)
  fold_ids <- sample(seq_len(nfolds), size = nrow(x), replace = TRUE)

  cv_fit <- glmnet::cv.glmnet(
    x = x, y = y,
    family = "cox",
    alpha = alpha,
    nfolds = nfolds,
    foldid = fold_ids,
    type.measure = "C",
    maxit = 100000L
  )

  coefs <- as.numeric(glmnet::coef.glmnet(cv_fit, s = "lambda.min"))
  names(coefs) <- gene_names
  active <- coefs[coefs != 0]

  list(
    cv_fit = cv_fit,
    active_genes = names(active),
    coefficients = active
  )
}
