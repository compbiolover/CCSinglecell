test_that("fit_cox_model validates inputs", {
  expect_error(fit_cox_model("not_df", "gene"), "data.frame")
  expect_error(fit_cox_model(data.frame(x = 1), character(0)), "non-empty")
  expect_error(
    fit_cox_model(data.frame(x = 1), "x"),
    "missing required columns"
  )
})

test_that("fit_cox_model drops missing genes gracefully", {
  df <- data.frame(
    gene1 = rnorm(50), gene2 = rnorm(50),
    days_to_last_follow_up = rexp(50, 0.01),
    vital_status = sample(0:1, 50, replace = TRUE)
  )
  # gene3 doesn't exist; should still work with gene1 and gene2
  result <- fit_cox_model(df, c("gene1", "gene2", "gene3"), seed = 42)

  expect_type(result, "list")
  expect_named(result, c("cv_fit", "active_genes", "coefficients"))
  expect_s3_class(result$cv_fit, "cv.glmnet")
  expect_true(is.character(result$active_genes))
  expect_true(is.numeric(result$coefficients))
})

test_that("fit_cox_model respects max_genes", {
  set.seed(99)
  df <- data.frame(
    g1 = rnorm(50), g2 = rnorm(50), g3 = rnorm(50),
    days_to_last_follow_up = rexp(50, 0.01),
    vital_status = sample(0:1, 50, replace = TRUE)
  )
  result <- fit_cox_model(df, c("g1", "g2", "g3"), max_genes = 2, seed = 42)
  # Active genes should be a subset of the first 2
  expect_true(all(result$active_genes %in% c("g1", "g2")))
})

test_that("fit_cox_model sanitises gene names with hyphens", {
  df <- data.frame(
    `gene-1` = rnorm(50), gene2 = rnorm(50),
    days_to_last_follow_up = rexp(50, 0.01),
    vital_status = sample(0:1, 50, replace = TRUE),
    check.names = FALSE
  )
  # gene-1 should become gene.1 internally
  result <- fit_cox_model(df, c("gene-1", "gene2"), seed = 42)
  expect_type(result, "list")
})
