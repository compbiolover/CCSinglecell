test_that("calculate_risk_scores validates inputs", {
  expect_error(calculate_risk_scores("not_df", "g1", c(g1 = 0.5)), "data.frame")
  expect_error(
    calculate_risk_scores(data.frame(g1 = 1), character(0), numeric(0)),
    "must not be empty"
  )
  expect_error(
    calculate_risk_scores(
      data.frame(x = 1, days_to_last_follow_up = 1, vital_status = 0),
      "missing_gene", c(missing_gene = 0.5)
    ),
    "not found"
  )
})

test_that("calculate_risk_scores produces correct output structure", {
  set.seed(42)
  n <- 30
  df <- data.frame(
    geneA = rnorm(n, 5, 2),
    geneB = rnorm(n, 3, 1),
    days_to_last_follow_up = rexp(n, 0.005),
    vital_status = sample(0:1, n, replace = TRUE)
  )
  coefs <- c(geneA = 0.3, geneB = -0.2)
  result <- calculate_risk_scores(df, c("geneA", "geneB"), coefs)

  expect_s3_class(result, "data.frame")
  expect_named(result, c("risk_score", "risk_group", "vital_status", "time"))
  expect_equal(nrow(result), n)
  expect_true(all(result$risk_group %in% c("high", "low")))
  expect_true(is.numeric(result$risk_score))
})

test_that("calculate_risk_scores handles single gene", {
  set.seed(1)
  df <- data.frame(
    g1 = rnorm(20, 10, 3),
    days_to_last_follow_up = rexp(20, 0.01),
    vital_status = sample(0:1, 20, replace = TRUE)
  )
  result <- calculate_risk_scores(df, "g1", c(g1 = 1.5))
  expect_equal(nrow(result), 20L)
})
