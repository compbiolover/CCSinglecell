test_that("learn_weights returns a valid weight vector", {
  r1 <- c(g1 = 0.4, g2 = 0.3, g3 = 0.2, g4 = 0.1)
  r2 <- c(g1 = 0.35, g2 = 0.30, g3 = 0.25, g4 = 0.10)
  r3 <- c(g1 = 0.1, g2 = 0.4, g3 = 0.1, g4 = 0.4)
  w <- learn_weights(list(mad = r1, switchde = r2, mirna = r3))

  expect_named(w, c("mad", "switchde", "mirna"))
  expect_true(all(w >= 0))
  expect_equal(sum(w), 1, tolerance = 1e-8)
  expect_true(is.numeric(attr(w, "variance_explained")))
})

test_that("learn_weights gives a non-informative (constant) metric zero weight", {
  r1 <- c(g1 = 0.4, g2 = 0.3, g3 = 0.2, g4 = 0.1)
  r2 <- c(g1 = 0.35, g2 = 0.30, g3 = 0.25, g4 = 0.10)
  flat <- c(g1 = 0.25, g2 = 0.25, g3 = 0.25, g4 = 0.25) # zero variance
  w <- learn_weights(list(a = r1, b = r2, flat = flat))

  expect_equal(unname(w["flat"]), 0)
  expect_gt(w["a"], 0)
  expect_gt(w["b"], 0)
  expect_equal(sum(w), 1, tolerance = 1e-8)
})

test_that("learn_weights handles single and all-constant metrics", {
  expect_equal(as.numeric(learn_weights(list(only = c(g1 = 1, g2 = 2)))), 1)
  # Both metrics constant across genes -> fall back to equal weights
  w <- learn_weights(list(a = c(g1 = 0.5, g2 = 0.5), b = c(g1 = 0.3, g2 = 0.3)))
  expect_equal(as.numeric(w), c(0.5, 0.5))
})

test_that("score_rankings supports weights = 'learn'", {
  r1 <- c(g1 = 0.4, g2 = 0.3, g3 = 0.2, g4 = 0.1)
  r2 <- c(g1 = 0.35, g2 = 0.30, g3 = 0.25, g4 = 0.10)
  sc <- score_rankings(c("g1", "g2", "g3", "g4"),
    list(a = r1, b = r2), weights = "learn")
  w <- attr(sc, "weights")
  expect_equal(sum(w), 1, tolerance = 1e-8)
  expect_true(all(w >= 0))
  expect_error(
    score_rankings(c("g1", "g2"),
      list(a = c(g1 = 0.5, g2 = 0.5), b = c(g1 = 0.3, g2 = 0.7)),
      weights = "bogus"),
    "must be"
  )
})

test_that("weights = 'learn' preserves the variance-explained metadata", {
  r1 <- c(g1 = 0.4, g2 = 0.3, g3 = 0.2, g4 = 0.1)
  r2 <- c(g1 = 0.35, g2 = 0.30, g3 = 0.25, g4 = 0.10)
  sc <- score_rankings(c("g1", "g2", "g3", "g4"),
    list(a = r1, b = r2), weights = "learn")
  ve <- attr(attr(sc, "weights"), "variance_explained")
  expect_true(is.numeric(ve) && !is.na(ve))
})

test_that("learn_weights validates its rankings", {
  expect_error(
    learn_weights(list(a = c(g1 = 1, g2 = 2), b = c(1, 2))),
    "named numeric vector"
  )
  expect_error(
    learn_weights(list(a = c(g1 = 1, g1 = 2), b = c(g1 = 1, g2 = 2))),
    "duplicated gene"
  )
  dup <- list(x = c(g1 = 1, g2 = 2), x = c(g1 = 2, g2 = 1))
  expect_error(learn_weights(dup), "must be unique")
})

test_that("score_multiomics and score_gene_set accept weights = 'learn'", {
  set.seed(1)
  genes <- paste0("g", 1:30)
  expr <- matrix(stats::rpois(30 * 12, 3), nrow = 30,
    dimnames = list(genes, paste0("c", 1:12)))
  mir <- matrix(sample(0:3, 30 * 3, replace = TRUE), nrow = 30,
    dimnames = list(genes, paste0("mir", 1:3)))

  sc1 <- score_gene_set(genes, expr, mirna_matrix = mir, weights = "learn")
  expect_equal(sum(attr(sc1, "weights")), 1, tolerance = 1e-8)

  sc2 <- score_multiomics(genes, list(
    omics_block("expr", expr, "mad"),
    omics_block("mir", mir, "mirna")
  ), weights = "learn")
  expect_equal(sum(attr(sc2, "weights")), 1, tolerance = 1e-8)
})
