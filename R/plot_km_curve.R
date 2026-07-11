#' Plot a Kaplan-Meier survival curve for risk groups
#'
#' Wraps [survminer::ggsurvplot()] to draw the paper's primary result: how
#' well a risk score stratifies patients into high- and low-risk groups. This
#' consumes the output of [calculate_risk_scores()].
#'
#' @param risk_df A data frame with `time`, `vital_status`, and `risk_group`
#'   columns (as returned by [calculate_risk_scores()]).
#' @param title Optional plot title.
#' @param pval Logical; annotate the log-rank p-value (default `TRUE`).
#' @param conf_int Logical; draw confidence bands (default `TRUE`).
#' @param risk_table Logical; add a numbers-at-risk table (default `TRUE`).
#' @param palette Colours for the high/low risk groups.
#'
#' @return A `ggsurvplot` object (a list containing a \pkg{ggplot2} plot).
#' @export
#'
#' @examples
#' \dontrun{
#' risk <- calculate_risk_scores(patient_df, model$active_genes, model$coefficients)
#' plot_km_curve(risk, title = "TCGA-COAD")
#' }
plot_km_curve <- function(
    risk_df,
    title = NULL,
    pval = TRUE,
    conf_int = TRUE,
    risk_table = TRUE,
    palette = c("#E7B800", "#2E9FDF")
) {
  stopifnot(is.data.frame(risk_df))
  required <- c("time", "vital_status", "risk_group")
  missing_cols <- setdiff(required, colnames(risk_df))
  if (length(missing_cols) > 0L) {
    stop("risk_df is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  if (!requireNamespace("survminer", quietly = TRUE)) {
    stop("Package 'survminer' is required for plot_km_curve(). ",
         "Install with: install.packages('survminer')")
  }

  fit <- survival::survfit(
    survival::Surv(time, vital_status) ~ risk_group,
    data = risk_df
  )

  survminer::ggsurvplot(
    fit,
    data = risk_df,
    pval = pval,
    conf.int = conf_int,
    risk.table = risk_table,
    palette = palette,
    xlab = "Time (days)",
    ylab = "Survival probability",
    legend.title = "Risk group",
    title = title
  )
}
