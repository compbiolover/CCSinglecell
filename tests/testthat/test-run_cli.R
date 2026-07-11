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
