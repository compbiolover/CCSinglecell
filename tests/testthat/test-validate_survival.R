sim_surv <- function(n = 300, seed = 1) {
  set.seed(seed)
  g1 <- stats::rnorm(n); g2 <- stats::rnorm(n)
  age <- stats::rnorm(n, 60, 10)
  noise <- stats::rnorm(n)
  lp <- 0.9 * g1 - 0.6 * g2 + 0.02 * (age - 60)
  time <- stats::rexp(n, rate = exp(lp - mean(lp)))
  data.frame(time = time, event = stats::rbinom(n, 1, 0.7),
    g1 = g1, g2 = g2, age = age, noise = noise)
}

test_that("cindex_ci returns a valid, ordered interval and detects signal", {
  set.seed(1)
  risk <- stats::rnorm(300)
  time <- stats::rexp(300, rate = exp(0.8 * risk))
  event <- stats::rbinom(300, 1, 0.7)
  ci <- cindex_ci(time, event, risk, n_boot = 200)
  expect_named(ci, c("cindex", "lower", "upper"))
  expect_true(ci["lower"] <= ci["cindex"] && ci["cindex"] <= ci["upper"])
  expect_true(all(ci >= 0 & ci <= 1))
  expect_gt(ci["cindex"], 0.6) # real signal -> C well above 0.5
})

test_that("validate_survival reports signature, clinical, and combined", {
  d <- sim_surv()
  v <- validate_survival(d, predictors = c("g1", "g2"), clinical = "age",
    n_boot = 200)
  expect_s3_class(v, "survival_validation")
  expect_equal(v$cindex$model, c("signature", "clinical", "combined"))
  expect_true(all(v$cindex$lower <= v$cindex$cindex))
  expect_true(all(v$cindex$cindex <= v$cindex$upper))
  # informative genes beat the uninformative age baseline
  cidx <- stats::setNames(v$cindex$cindex, v$cindex$model)
  expect_gt(cidx["signature"], cidx["clinical"])
})

test_that("pure-noise predictor lands near 0.5", {
  d <- sim_surv()
  v <- validate_survival(d, predictors = "noise", n_boot = 200)
  expect_lt(abs(v$cindex$cindex[1] - 0.5), 0.12)
})

test_that("external evaluation and calibration work", {
  tr <- sim_surv(n = 300, seed = 1)
  te <- sim_surv(n = 300, seed = 2)
  v <- validate_survival(tr, predictors = c("g1", "g2"), clinical = "age",
    test = te, horizon = 1, n_groups = 4, n_boot = 200)
  expect_true(v$meta$external)
  expect_s3_class(v$calibration, "data.frame")
  expect_equal(nrow(v$calibration), 4L)
  expect_true(all(c("predicted", "observed", "obs_lower", "obs_upper") %in%
    names(v$calibration)))
  expect_s3_class(plot_calibration(v), "ggplot")
})

test_that("validate_survival warns without a clinical baseline and validates cols", {
  d <- sim_surv()
  expect_warning(validate_survival(d, predictors = "g1", n_boot = 50),
    "clinical")
  expect_error(validate_survival(d, predictors = "nope", n_boot = 50),
    "missing columns")
})
