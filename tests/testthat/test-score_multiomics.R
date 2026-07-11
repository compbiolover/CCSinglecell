make_layers <- function() {
  set.seed(1)
  genes <- paste0("g", 1:40)
  expr <- matrix(
    stats::rpois(40 * 12, lambda = 3), nrow = 40,
    dimnames = list(genes, paste0("c", 1:12))
  )
  mir <- matrix(
    sample(0:3, 40 * 3, replace = TRUE), nrow = 40,
    dimnames = list(genes, paste0("mir", 1:3))
  )
  meth <- matrix(
    stats::runif(40 * 12), nrow = 40,
    dimnames = list(genes, paste0("c", 1:12))
  )
  list(genes = genes, expr = expr, mir = mir, meth = meth)
}

test_that("omics_block accepts built-in and custom metrics", {
  d <- make_layers()
  b1 <- omics_block("expression", d$expr, "mad")
  b2 <- omics_block("mean_expr", d$expr, function(m) {
    sort(rowMeans(m), decreasing = TRUE)
  })
  expect_s3_class(b1, "omics_block")
  expect_s3_class(b2, "omics_block")
  expect_true(is.function(b1$metric))
  expect_error(omics_block("x", d$expr, "not_a_metric"), "Unknown built-in metric")
  expect_error(omics_block("", d$expr, "mad"), "non-empty string")
})

test_that("score_multiomics blends an arbitrary number of layers", {
  d <- make_layers()
  blocks <- list(
    omics_block("expression", d$expr, "mad"),
    omics_block("mirna", d$mir, "mirna"),
    omics_block("methylation", d$meth, function(m) {
      sort(apply(m, 1L, stats::sd), decreasing = TRUE)
    })
  )
  sc <- score_multiomics(d$genes, blocks)

  expect_s3_class(sc, "data.frame")
  expect_equal(
    names(sc),
    c("gene", "expression", "mirna", "methylation", "combined", "rank")
  )
  expect_equal(nrow(sc), 40L)
  expect_true(all(diff(sc$combined) <= 0))
  expect_equal(sort(attr(sc, "metrics")),
    c("expression", "methylation", "mirna"))
  expect_equal(unname(attr(sc, "weights")), rep(1 / 3, 3))
})

test_that("score_multiomics honours custom weights and validates blocks", {
  d <- make_layers()
  blocks <- list(
    omics_block("a", d$expr, "mad"),
    omics_block("b", d$mir, "mirna")
  )
  sc <- score_multiomics(d$genes, blocks, weights = c(3, 1))
  expect_equal(unname(attr(sc, "weights")), c(0.75, 0.25))

  expect_error(score_multiomics(d$genes, list()), "non-empty list")
  expect_error(score_multiomics(d$genes, list(1, 2)), "must be an omics_block")
  expect_error(
    score_multiomics(d$genes, list(
      omics_block("dup", d$expr, "mad"),
      omics_block("dup", d$mir, "mirna")
    )),
    "names must be unique"
  )
})

test_that("score_rankings works directly on precomputed rankings", {
  r1 <- c(g1 = 0.5, g2 = 0.3, g3 = 0.2)
  r2 <- c(g2 = 0.6, g3 = 0.3, g4 = 0.1)
  sc <- score_rankings(c("g1", "g2", "g3"),
    list(expression = r1, methylation = r2))
  expect_equal(names(sc),
    c("gene", "expression", "methylation", "combined", "rank"))
  expect_error(
    score_rankings("g1", list(r1)),
    "must be named"
  )
})

test_that("score_gene_set remains a special case of the multi-omics engine", {
  d <- make_layers()
  via_wrapper <- score_gene_set(d$genes, d$expr, mirna_matrix = d$mir)
  via_engine <- score_multiomics(d$genes, list(
    omics_block("mad", d$expr, "mad"),
    omics_block("mirna", d$mir, "mirna")
  ))
  # Same genes, same combined scores (column names differ by design).
  expect_equal(via_wrapper$gene, via_engine$gene)
  expect_equal(via_wrapper$combined, via_engine$combined)
})
