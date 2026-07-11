#' Plot metric-weight optimisation results
#'
#' Visualises how model performance (e.g. concordance / C-index) varies across
#' the grid of metric weights explored during optimisation. For two metrics a
#' line chart of performance versus the first weight is drawn; for three
#' metrics a `w1` x `w2` heatmap is drawn (the third weight is
#' `1 - w1 - w2`).
#'
#' @param results A data frame with weight columns (`w1`, `w2`, and optionally
#'   `w3`) and a performance column (default name `cindex`). This can be built
#'   by evaluating each row of an [optimize_weights()] grid.
#' @param metric_col Name of the performance column (default `"cindex"`).
#' @param title Optional plot title.
#'
#' @return A \pkg{ggplot2} object.
#' @importFrom ggplot2 .data
#' @export
#'
#' @examples
#' grid <- data.frame(w1 = seq(0, 1, 0.1), w2 = seq(1, 0, -0.1))
#' grid$cindex <- 0.5 + 0.2 * grid$w1 - 0.15 * grid$w1^2
#' plot_weight_optimization(grid)
plot_weight_optimization <- function(results, metric_col = "cindex",
                                     title = "Weight optimisation") {
  stopifnot(is.data.frame(results))
  if (!metric_col %in% colnames(results)) {
    stop("`results` must contain a '", metric_col, "' column")
  }
  if (!"w1" %in% colnames(results)) {
    stop("`results` must contain weight columns (w1, w2[, w3])")
  }

  has_w3 <- "w3" %in% colnames(results)

  if (!has_w3) {
    ggplot2::ggplot(
      results, ggplot2::aes(x = .data$w1, y = .data[[metric_col]])
    ) +
      ggplot2::geom_line(colour = "grey60") +
      ggplot2::geom_point(colour = "#2E9FDF", size = 2.5) +
      ggplot2::labs(
        x = "Weight on metric 1 (w1)", y = toupper(metric_col), title = title
      ) +
      ggplot2::theme_minimal(base_size = 12)
  } else {
    ggplot2::ggplot(
      results,
      ggplot2::aes(x = .data$w1, y = .data$w2, fill = .data[[metric_col]])
    ) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_viridis_c(name = toupper(metric_col)) +
      ggplot2::labs(
        x = "Weight on metric 1 (w1)",
        y = "Weight on metric 2 (w2)",
        title = title
      ) +
      ggplot2::theme_minimal(base_size = 12)
  }
}
