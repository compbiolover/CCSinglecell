make_expr <- function() {
  set.seed(1)
  matrix(
    stats::rpois(40 * 12, lambda = 3),
    nrow = 40, ncol = 12,
    dimnames = list(paste0("g", 1:40), paste0("c", 1:12))
  )
}

test_that("score_gene_set returns tidy scores for MAD-only input", {
  expr <- make_expr()
  sc <- score_gene_set(paste0("g", 1:10), expr)

  expect_s3_class(sc, "data.frame")
  expect_equal(names(sc), c("gene", "mad", "combined", "rank"))
  expect_equal(nrow(sc), 10L)
  expect_equal(sc$rank, seq_len(10L))
  # Sorted by combined, decreasing
  expect_true(all(diff(sc$combined) <= 0))
  expect_equal(attr(sc, "metrics"), "mad")
})

test_that("score_gene_set adds the miRNA metric when supplied", {
  expr <- make_expr()
  mir <- matrix(
    sample(0:3, 40 * 3, replace = TRUE), nrow = 40,
    dimnames = list(paste0("g", 1:40), paste0("mir", 1:3))
  )
  sc <- score_gene_set(paste0("g", 1:10), expr, mirna_matrix = mir)

  expect_true(all(c("mad", "mirna", "combined") %in% names(sc)))
  expect_equal(sort(attr(sc, "metrics")), c("mad", "mirna"))
  # Equal weights by default
  expect_equal(unname(attr(sc, "weights")), c(0.5, 0.5))
})

test_that("score_gene_set honours custom weights", {
  expr <- make_expr()
  mir <- matrix(1, nrow = 40, ncol = 2,
    dimnames = list(paste0("g", 1:40), c("mir1", "mir2")))
  sc <- score_gene_set(paste0("g", 1:10), expr, mirna_matrix = mir,
                       weights = c(3, 1))
  expect_equal(unname(attr(sc, "weights")), c(0.75, 0.25))
})

test_that("score_gene_set warns and zero-scores genes absent from the matrix", {
  expr <- make_expr()
  expect_warning(
    sc <- score_gene_set(c("g1", "not_a_gene"), expr),
    "not found"
  )
  expect_equal(sc$combined[sc$gene == "not_a_gene"], 0)
})

test_that("score_gene_set validates inputs", {
  expr <- make_expr()
  expect_error(score_gene_set(character(0), expr), "empty")
  expect_error(score_gene_set("g1", as.data.frame(expr)), "numeric matrix")
  expect_error(
    score_gene_set("g1", matrix(1:4, 2)),
    "row names"
  )
  mir <- matrix(1, nrow = 40, ncol = 2,
    dimnames = list(paste0("g", 1:40), c("mir1", "mir2")))
  expect_error(
    score_gene_set("g1", expr, mirna_matrix = mir, weights = c(1)),
    "one value per metric"
  )
})
