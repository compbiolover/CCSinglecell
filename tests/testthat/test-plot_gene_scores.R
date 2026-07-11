make_scores <- function(with_mirna = TRUE) {
  set.seed(1)
  expr <- matrix(
    stats::rpois(40 * 12, lambda = 3), nrow = 40,
    dimnames = list(paste0("g", 1:40), paste0("c", 1:12))
  )
  mir <- if (with_mirna) {
    matrix(sample(0:3, 40 * 3, replace = TRUE), nrow = 40,
      dimnames = list(paste0("g", 1:40), paste0("mir", 1:3)))
  } else {
    NULL
  }
  score_gene_set(paste0("g", 1:40), expr, mirna_matrix = mir)
}

test_that("plot_gene_scores returns a ggplot", {
  p <- plot_gene_scores(make_scores(), top_n = 10)
  expect_s3_class(p, "ggplot")
})

test_that("plot_score_contributions returns a ggplot", {
  p <- plot_score_contributions(make_scores(), top_n = 10)
  expect_s3_class(p, "ggplot")
})

test_that("plot_score_contributions needs score_gene_set attributes", {
  df <- data.frame(gene = c("a", "b"), combined = c(0.6, 0.4))
  expect_error(plot_score_contributions(df), "score_gene_set")
})

test_that("plot_weight_optimization handles 2- and 3-metric grids", {
  g2 <- data.frame(w1 = seq(0, 1, 0.25), w2 = seq(1, 0, -0.25))
  g2$cindex <- 0.5 + 0.1 * g2$w1
  expect_s3_class(plot_weight_optimization(g2), "ggplot")

  g3 <- expand.grid(w1 = c(0, 0.5, 1), w2 = c(0, 0.5))
  g3$w3 <- 1 - g3$w1 - g3$w2
  g3$cindex <- 0.5 + 0.05 * g3$w1
  expect_s3_class(plot_weight_optimization(g3), "ggplot")
})

test_that("plot_weight_optimization validates the performance column", {
  expect_error(
    plot_weight_optimization(data.frame(w1 = 0:1), metric_col = "cindex"),
    "cindex"
  )
})
