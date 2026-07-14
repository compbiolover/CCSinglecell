sim_null <- function(n = 200, seed = 1) {
  set.seed(seed)
  age <- stats::rnorm(n, 60, 10)
  g <- matrix(stats::rnorm(n * 12), n, dimnames = list(NULL, paste0("g", 1:12)))
  # genes are pure noise; only age drives survival
  data.frame(time = stats::rexp(n, exp(0.05 * (age - 60))),
    event = stats::rbinom(n, 1, 0.7), age = age, g)
}

sim_signal <- function(n = 220, seed = 1) {
  set.seed(seed)
  age <- stats::rnorm(n, 60, 10)
  g <- matrix(stats::rnorm(n * 12), n, dimnames = list(NULL, paste0("g", 1:12)))
  lp <- 0.9 * g[, 1] - 0.7 * g[, 2] + 0.05 * (age - 60)
  data.frame(time = stats::rexp(n, exp(lp - mean(lp))),
    event = stats::rbinom(n, 1, 0.7), age = age, g)
}

test_that("assess_null does not flag a genuine null", {
  skip_if_not_installed("glmnet")
  d <- sim_null()
  a <- assess_null(d, predictors = paste0("g", 1:12), clinical = "age",
    n_perm = 40, n_repeats = 4, seed = 1)
  expect_s3_class(a, "null_assessment")
  expect_length(a$null, 40)
  expect_true(a$p_value > 0.1)                 # noise gain sits inside the null
  expect_true(a$p_value >= 0 && a$p_value <= 1)
})

test_that("assess_null detects a real incremental signal", {
  skip_if_not_installed("glmnet")
  d <- sim_signal()
  a <- assess_null(d, predictors = paste0("g", 1:12), clinical = "age",
    n_perm = 40, n_repeats = 4, seed = 1)
  expect_gt(a$observed, 0.02)                  # real gain
  expect_lt(a$p_value, 0.1)                     # beyond the permutation null
})

test_that("assess_null requires a clinical baseline and valid columns", {
  skip_if_not_installed("glmnet")
  d <- sim_null()
  expect_error(assess_null(d, predictors = paste0("g", 1:12)), "clinical.*required")
  expect_error(
    assess_null(d, predictors = "nope", clinical = "age", n_perm = 2, n_repeats = 2),
    "missing columns")
})

test_that("power_curve traces a monotone-ish detectable floor", {
  skip_if_not_installed("glmnet")
  d <- sim_null()
  pc <- power_curve(d, predictors = paste0("g", 1:12), clinical = "age",
    effect_sizes = c(0, 0.4, 0.8), n_sim = 6, n_repeats = 4, seed = 1)
  expect_s3_class(pc, "power_curve")
  expect_equal(nrow(pc$curve), 3L)
  expect_true(all(pc$curve$detect_rate >= 0 & pc$curve$detect_rate <= 1))
  # a strong injected effect is detected more often than none
  expect_gte(pc$curve$detect_rate[pc$curve$effect == 0.8],
             pc$curve$detect_rate[pc$curve$effect == 0])
})

test_that("power_curve rejects a colliding column name", {
  skip_if_not_installed("glmnet")
  d <- sim_null()
  d$planted_signal <- 0
  expect_error(
    power_curve(d, predictors = paste0("g", 1:12), clinical = "age"),
    "planted_signal")
})

test_that("learning_curve reports C-index by training size", {
  skip_if_not_installed("glmnet")
  d <- sim_signal(n = 260)
  lc <- learning_curve(d, predictors = paste0("g", 1:12), clinical = "age",
    fractions = c(0.5, 0.75, 1.0), n_subsample = 3, seed = 1)
  expect_s3_class(lc, "learning_curve")
  expect_equal(nrow(lc$curve), 3L)
  expect_true(all(c("signature", "clinical", "combined") %in% colnames(lc$curve)))
  expect_true(all(lc$curve$combined >= 0 & lc$curve$combined <= 1))
  expect_true(is.finite(lc$combined_slope))
})

test_that("null-diagnostic plots return ggplot objects", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("ggplot2")
  d <- sim_null()
  a <- assess_null(d, predictors = paste0("g", 1:12), clinical = "age",
    n_perm = 20, n_repeats = 3, seed = 1)
  pc <- power_curve(d, predictors = paste0("g", 1:12), clinical = "age",
    effect_sizes = c(0, 0.5), n_sim = 4, n_repeats = 3, seed = 1)
  lc <- learning_curve(d, predictors = paste0("g", 1:12), clinical = "age",
    fractions = c(0.6, 1.0), n_subsample = 2, seed = 1)
  expect_s3_class(plot_null_assessment(a), "ggplot")
  expect_s3_class(plot_power_curve(pc), "ggplot")
  expect_s3_class(plot_learning_curve(lc), "ggplot")
})
