sim_cv <- function(n = 250, beta_g = 0.8, seed = 1) {
  set.seed(seed)
  g <- matrix(stats::rnorm(n * 20), n, dimnames = list(NULL, paste0("g", 1:20)))
  age <- stats::rnorm(n, 60, 10)
  lp <- beta_g * g[, 1] - beta_g * g[, 2] + 0.03 * (age - 60)
  data.frame(time = stats::rexp(n, exp(lp - mean(lp))),
    event = stats::rbinom(n, 1, 0.7), age = age, g)
}

test_that("cv_validate_survival detects real genomic signal over clinical", {
  skip_if_not_installed("glmnet")
  d <- sim_cv(beta_g = 1.0)
  v <- cv_validate_survival(d, predictors = paste0("g", 1:20),
    clinical = "age", n_repeats = 5, seed = 1)
  expect_s3_class(v, "cv_survival_validation")
  expect_equal(v$cindex$model, c("signature", "clinical", "combined"))
  expect_true(all(v$cindex$cindex >= 0 & v$cindex$cindex <= 1))
  cidx <- stats::setNames(v$cindex$cindex, v$cindex$model)
  expect_gt(cidx["signature"], cidx["clinical"])   # genes carry the signal
  expect_gt(v$delta$mean, 0)                        # combined beats clinical
  expect_gt(v$delta$win_rate, 0.5)
})

test_that("pure-noise genomics does not beat clinical", {
  skip_if_not_installed("glmnet")
  d <- sim_cv(beta_g = 0)          # genes are noise; only age matters
  d$time <- stats::rexp(nrow(d), exp(0.04 * (d$age - 60)))
  v <- cv_validate_survival(d, predictors = paste0("g", 1:20),
    clinical = "age", n_repeats = 5, seed = 2)
  # combined should not reliably beat clinical: interval includes 0
  expect_lte(v$delta$lower, 0)
})

test_that("cv_validate_survival validates inputs", {
  skip_if_not_installed("glmnet")
  d <- sim_cv()
  expect_error(cv_validate_survival(d, predictors = paste0("g", 1:20)),
    "clinical.*required")
  expect_error(
    cv_validate_survival(d, predictors = "nope", clinical = "age", n_repeats = 2),
    "missing columns")
})

test_that("group holds out each cohort once and scores it within that cohort", {
  skip_if_not_installed("glmnet")
  d <- sim_cv(n = 300, beta_g = 1.0)
  d$cohort <- rep(c("a", "b", "c"), length.out = nrow(d))
  v <- cv_validate_survival(d, predictors = paste0("g", 1:20),
    clinical = "age", group = "cohort", seed = 1)
  expect_equal(rownames(v$per_repeat), c("a", "b", "c"))
  expect_true(all(v$per_repeat >= 0 & v$per_repeat <= 1))
  expect_gt(v$delta$mean, 0)
  expect_output(print(v), "Leave-one-cohort-out")

  d$cohort <- "a"
  expect_error(cv_validate_survival(d, paste0("g", 1:20), "age", group = "cohort"),
    "at least 2")
})
