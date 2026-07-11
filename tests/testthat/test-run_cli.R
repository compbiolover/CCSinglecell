test_that("run_cli help prints usage and returns NULL invisibly", {
  expect_output(res <- run_cli("--help"), "ccscore")
  expect_null(res)
})

test_that("run_cli rejects unknown commands", {
  expect_error(run_cli("frobnicate"), "Unknown command")
})

test_that("run_cli score mode writes a scores CSV and plots", {
  skip_if_not_installed("ggplot2")

  dir <- withr::local_tempdir()
  expr_path <- file.path(dir, "expr.csv")
  genes_path <- file.path(dir, "genes.csv")
  out_dir <- file.path(dir, "out")

  set.seed(1)
  expr <- matrix(
    stats::rpois(40 * 12, lambda = 3), nrow = 40,
    dimnames = list(paste0("g", 1:40), paste0("c", 1:12))
  )
  utils::write.csv(expr, expr_path)
  utils::write.csv(data.frame(x = paste0("g", 1:20)), genes_path)

  res <- suppressMessages(run_cli(c(
    "score",
    "--genes", genes_path,
    "--expression", expr_path,
    "--out", out_dir
  )))

  expect_true(file.exists(file.path(out_dir, "scores.csv")))
  expect_true(file.exists(file.path(out_dir, "gene_scores.png")))
  scores <- utils::read.csv(file.path(out_dir, "scores.csv"))
  expect_equal(nrow(scores), 20L)
  expect_true(all(c("gene", "mad", "combined", "rank") %in% names(scores)))
})

test_that("run_cli score mode errors without required options", {
  expect_error(run_cli(c("score", "--genes", "g.csv")), "requires --expression")
})

test_that("read_cohort_csv preserves patient IDs from a leading column", {
  dir <- withr::local_tempdir()
  path <- file.path(dir, "cohort.csv")
  df <- data.frame(
    TP53 = stats::rnorm(3),
    vital_status = c(0L, 1L, 0L),
    days_to_last_follow_up = c(10, 20, 30),
    row.names = c("P1", "P2", "P3")
  )
  utils::write.csv(df, path) # row.names = TRUE by default

  got <- read_cohort_csv(path)
  expect_equal(rownames(got), c("P1", "P2", "P3"))
  expect_false("X" %in% names(got))
  expect_true(all(
    c("TP53", "vital_status", "days_to_last_follow_up") %in% names(got)
  ))
})

test_that("read_pseudotime aligns by cell ID regardless of row order", {
  dir <- withr::local_tempdir()
  path <- file.path(dir, "pt.csv")
  # Rows are shuffled relative to the expression columns.
  utils::write.csv(
    data.frame(cell = c("c2", "c1", "c3"), pseudotime = c(0.2, 0.1, 0.3)),
    path, row.names = FALSE
  )
  got <- read_pseudotime(path, c("c1", "c2", "c3"))
  expect_equal(got, c(0.1, 0.2, 0.3))
})

test_that("read_pseudotime errors on length mismatch without a cell-ID column", {
  dir <- withr::local_tempdir()
  path <- file.path(dir, "pt.csv")
  utils::write.csv(data.frame(pseudotime = c(0.1, 0.2)), path, row.names = FALSE)
  expect_error(read_pseudotime(path, c("c1", "c2", "c3")), "no cell-ID column")
})
