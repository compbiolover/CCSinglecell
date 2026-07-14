#' Plot a permutation-null assessment
#'
#' Histogram of the permuted-block `combined - clinical` gains with the observed
#' gain marked, so the reader can see at a glance whether the real genomic block
#' lands inside or beyond the noise band.
#'
#' @param assessment A `null_assessment` from [assess_null()].
#' @return A \pkg{ggplot2} object.
#' @importFrom ggplot2 .data
#' @export
plot_null_assessment <- function(assessment) {
  if (!inherits(assessment, "null_assessment")) {
    stop("`assessment` must come from assess_null()")
  }
  null <- data.frame(delta = assessment$null[is.finite(assessment$null)])
  ggplot2::ggplot(null, ggplot2::aes(x = .data$delta)) +
    ggplot2::geom_histogram(bins = 30, fill = "grey75", colour = "white") +
    ggplot2::geom_vline(xintercept = 0, linetype = 3, colour = "grey50") +
    ggplot2::geom_vline(xintercept = assessment$observed, colour = "#d7301f",
      linewidth = 1) +
    ggplot2::annotate("text", x = assessment$observed, y = Inf,
      label = sprintf("observed  p=%.3f", assessment$p_value),
      colour = "#d7301f", hjust = -0.05, vjust = 1.5, size = 3.5) +
    ggplot2::labs(
      x = "Permuted-block gain (combined - clinical C-index)",
      y = "Permutations",
      title = "Genomic increment vs permutation null"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}

#' Plot a detectable-effect (power) curve
#'
#' Detection rate versus injected effect size, with the target-power line and
#' the minimum detectable effect marked. Shows what incremental signal the
#' cohort could actually have found.
#'
#' @param pc A `power_curve` from [power_curve()].
#' @return A \pkg{ggplot2} object.
#' @importFrom ggplot2 .data
#' @export
plot_power_curve <- function(pc) {
  if (!inherits(pc, "power_curve")) stop("`pc` must come from power_curve()")
  g <- ggplot2::ggplot(pc$curve,
      ggplot2::aes(x = .data$effect, y = .data$detect_rate)) +
    ggplot2::geom_hline(yintercept = pc$power_target, linetype = 2,
      colour = "grey50") +
    ggplot2::geom_line(colour = "grey60") +
    ggplot2::geom_point(colour = "#2E9FDF", size = 2.5) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      x = "Injected incremental effect (hazard-log scale)",
      y = "Detection rate",
      title = "Detectable-effect (power) curve"
    ) +
    ggplot2::theme_minimal(base_size = 12)
  if (!is.na(pc$mde)) {
    g <- g + ggplot2::geom_vline(xintercept = pc$mde, colour = "#d7301f",
      linetype = 3)
  }
  g
}

#' Plot a learning curve
#'
#' Out-of-fold C-index against training-set size for the signature, clinical,
#' and combined models. A combined curve still rising at full n signals a
#' power-limited null; a plateau signals a signal-limited one.
#'
#' @param lc A `learning_curve` from [learning_curve()].
#' @return A \pkg{ggplot2} object.
#' @importFrom ggplot2 .data
#' @export
plot_learning_curve <- function(lc) {
  if (!inherits(lc, "learning_curve")) stop("`lc` must come from learning_curve()")
  long <- do.call(rbind, lapply(c("signature", "clinical", "combined"),
    function(m) data.frame(n = lc$curve$n, model = m, cindex = lc$curve[[m]])))
  ggplot2::ggplot(long,
      ggplot2::aes(x = .data$n, y = .data$cindex, colour = .data$model)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_colour_manual(values = c(
      signature = "#2E9FDF", clinical = "grey50", combined = "#d7301f")) +
    ggplot2::labs(
      x = "Training-set size (subjects)",
      y = "Out-of-fold C-index", colour = NULL,
      title = "Learning curve"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}
