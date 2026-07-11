make_activity_inputs <- function() {
  set.seed(1)
  samples <- paste0("s", 1:30)
  genes <- paste0("g", 1:6)
  mirs <- paste0("mir", 1:3)
  gene_expr <- matrix(stats::rnorm(6 * 30), nrow = 6,
    dimnames = list(genes, samples))
  mirna_expr <- matrix(stats::rnorm(3 * 30), nrow = 3,
    dimnames = list(mirs, samples))
  # g1 is strongly repressed by mir1 (anti-correlated); g2 is co-expressed with
  # mir2 (positive correlation, not repression).
  mirna_expr["mir1", ] <- -gene_expr["g1", ] + stats::rnorm(30, sd = 0.1)
  mirna_expr["mir2", ] <- gene_expr["g2", ] + stats::rnorm(30, sd = 0.1)
  target_matrix <- matrix(1, nrow = 6, ncol = 3,
    dimnames = list(genes, mirs))
  list(mirna_expr = mirna_expr, gene_expr = gene_expr,
       target_matrix = target_matrix)
}

test_that("calculate_mirna_activity ranks anti-correlated genes highest", {
  d <- make_activity_inputs()
  scores <- calculate_mirna_activity(d$mirna_expr, d$gene_expr, d$target_matrix)

  expect_type(scores, "double")
  expect_named(scores)
  expect_true(all(diff(scores) <= 0))
  expect_equal(sum(scores), 1, tolerance = 1e-8)
  # g1 (repressed) should score above g2 (co-expressed, no repression evidence)
  expect_gt(scores["g1"], scores["g2"])
  expect_equal(unname(which.max(scores)), unname(which(names(scores) == "g1")))
})

test_that("signed method penalises positive correlation", {
  d <- make_activity_inputs()
  signed <- calculate_mirna_activity(
    d$mirna_expr, d$gene_expr, d$target_matrix,
    method = "signed", normalize = FALSE
  )
  # g2 is positively correlated with mir2 -> negative signed evidence
  expect_lt(signed["g2"], 0)
})

test_that("calculate_mirna_activity validates inputs", {
  d <- make_activity_inputs()
  expect_error(
    calculate_mirna_activity(as.data.frame(d$mirna_expr), d$gene_expr, d$target_matrix),
    "must be a numeric matrix"
  )
  # No shared samples
  bad <- d$gene_expr
  colnames(bad) <- paste0("x", seq_len(ncol(bad)))
  expect_error(
    calculate_mirna_activity(d$mirna_expr, bad, d$target_matrix),
    "at least 3 samples"
  )
  # No overlap between target map and expression
  tm <- d$target_matrix
  rownames(tm) <- paste0("z", seq_len(nrow(tm)))
  expect_error(
    calculate_mirna_activity(d$mirna_expr, d$gene_expr, tm),
    "shares no genes"
  )
})

test_that("mirna_activity works as a built-in omics_block metric", {
  d <- make_activity_inputs()
  block <- omics_block(
    "mirna_activity",
    list(mirna_expr = d$mirna_expr, gene_expr = d$gene_expr,
         target_matrix = d$target_matrix),
    "mirna_activity"
  )
  sc <- score_multiomics(rownames(d$gene_expr), list(block))
  expect_true("mirna_activity" %in% names(sc))
  expect_equal(sc$gene[1], "g1")
})
