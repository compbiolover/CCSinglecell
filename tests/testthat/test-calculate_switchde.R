test_that("calculate_switchde validates inputs", {
  expect_error(calculate_switchde(c(1, 2), c(1)), "numeric matrix")
  expect_error(
    calculate_switchde(matrix(0, 0, 0), c()), "must not be empty"
  )
  expect_error(
    calculate_switchde(matrix(1, 2, 3), c(1, 2)),
    "ncol"
  )
  expect_error(
    calculate_switchde(matrix(1, 2, 3), "a"),
    "numeric"
  )
})

test_that("calculate_switchde accepts data.frame pseudotime", {
  mat <- matrix(rnorm(20), nrow = 2, ncol = 10)
  pt_df <- data.frame(Pseudotime = runif(10))
  # Should not error on input validation (will fail at switchde call)
  expect_error(
    calculate_switchde(mat, data.frame(other = 1:10)),
    "Pseudotime"
  )
})

test_that("calculate_switchde requires switchde package", {
  skip_if(requireNamespace("switchde", quietly = TRUE))
  mat <- matrix(rnorm(20), nrow = 2, ncol = 10)
  expect_error(calculate_switchde(mat, runif(10)), "switchde")
})
