#' Plot the top-scoring genes
#'
#' Draws a horizontal lollipop chart of the highest combined-score genes from
#' a [score_gene_set()] result.
#'
#' @param scores A data frame returned by [score_gene_set()] (must contain
#'   `gene` and `combined` columns).
#' @param top_n Integer; number of top genes to display (default 20).
#' @param title Optional plot title.
#'
#' @return A \pkg{ggplot2} object.
#' @importFrom ggplot2 .data
#' @export
#'
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(400), nrow = 40,
#'   dimnames = list(paste0("g", 1:40), paste0("c", 1:10)))
#' sc <- score_gene_set(paste0("g", 1:40), mat)
#' plot_gene_scores(sc, top_n = 10)
plot_gene_scores <- function(scores, top_n = 20, title = "Top scoring genes") {
  stopifnot(is.data.frame(scores), all(c("gene", "combined") %in% names(scores)))

  df <- utils::head(scores[order(scores$combined, decreasing = TRUE), ], top_n)
  df$gene <- factor(df$gene, levels = rev(df$gene))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$combined, y = .data$gene)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = .data$combined, yend = .data$gene),
      colour = "grey70"
    ) +
    ggplot2::geom_point(colour = "#2E9FDF", size = 3) +
    ggplot2::labs(
      x = "Combined score", y = NULL, title = title
    ) +
    ggplot2::theme_minimal(base_size = 12)
}

#' Plot per-metric score contributions
#'
#' Draws a stacked bar chart showing how much each metric (MAD, switchDE,
#' miRNA) contributes to the combined score of the top genes. Contributions
#' are the per-metric score multiplied by that metric's weight, so the bars
#' sum to the combined score.
#'
#' @param scores A data frame returned by [score_gene_set()]. Uses the
#'   `"weights"` and `"metrics"` attributes attached by that function.
#' @param top_n Integer; number of top genes to display (default 20).
#' @param title Optional plot title.
#'
#' @return A \pkg{ggplot2} object.
#' @importFrom ggplot2 .data
#' @export
#'
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(400), nrow = 40,
#'   dimnames = list(paste0("g", 1:40), paste0("c", 1:10)))
#' mir <- matrix(sample(0:3, 120, replace = TRUE), nrow = 40,
#'   dimnames = list(paste0("g", 1:40), paste0("mir", 1:3)))
#' sc <- score_gene_set(paste0("g", 1:40), mat, mirna_matrix = mir)
#' plot_score_contributions(sc, top_n = 10)
plot_score_contributions <- function(scores, top_n = 20,
                                      title = "Score contributions by metric") {
  stopifnot(is.data.frame(scores), all(c("gene", "combined") %in% names(scores)))

  metrics <- attr(scores, "metrics")
  weights <- attr(scores, "weights")
  if (is.null(metrics) || is.null(weights)) {
    stop("`scores` must come from score_gene_set() (missing weight attributes)")
  }

  df <- utils::head(scores[order(scores$combined, decreasing = TRUE), ], top_n)
  gene_levels <- rev(df$gene)

  long <- do.call(rbind, lapply(metrics, function(m) {
    data.frame(
      gene = df$gene,
      metric = m,
      contribution = df[[m]] * weights[[m]],
      stringsAsFactors = FALSE
    )
  }))
  long$gene <- factor(long$gene, levels = gene_levels)
  long$metric <- factor(long$metric, levels = metrics)

  ggplot2::ggplot(
    long,
    ggplot2::aes(x = .data$contribution, y = .data$gene, fill = .data$metric)
  ) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_brewer(palette = "Set2", name = "Metric") +
    ggplot2::labs(x = "Weighted contribution", y = NULL, title = title) +
    ggplot2::theme_minimal(base_size = 12)
}
