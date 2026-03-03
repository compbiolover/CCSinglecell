test_that("combine_rankings merges two rankings", {
  r1 <- c(A = 0.5, B = 0.3, C = 0.2)
  r2 <- c(B = 0.6, C = 0.3, D = 0.1)
  result <- combine_rankings(list(r1, r2), c(0.5, 0.5))

  expect_type(result, "double")
  expect_named(result)
  # All genes from both inputs should be present

  expect_true(all(c("A", "B", "C", "D") %in% names(result)))
  # Sorted decreasing
  expect_true(all(diff(result) <= 0))
  # A only in r1: 0.5*0.5 + 0.5*0 = 0.25
  expect_equal(unname(result["A"]), 0.25)
  # B in both: 0.5*0.3 + 0.5*0.6 = 0.45
  expect_equal(unname(result["B"]), 0.45)
})

test_that("combine_rankings merges three rankings", {
  r1 <- c(A = 1, B = 0)
  r2 <- c(B = 1, C = 0)
  r3 <- c(C = 1, A = 0)
  result <- combine_rankings(list(r1, r2, r3), c(1/3, 1/3, 1/3))

  expect_length(result, 3L)
  expect_true(all(c("A", "B", "C") %in% names(result)))
  # Each gene appears once with weight 1/3
  expect_equal(unname(sort(result)), rep(1/3, 3), tolerance = 1e-10)
})

test_that("combine_rankings validates inputs", {
  expect_error(combine_rankings(list(c(A = 1)), c(1)), "2 or 3")
  expect_error(
    combine_rankings(list(c(A = 1), c(B = 1)), c(1)),
    "same length"
  )
})

test_that("optimize_weights returns correct grid for 2 metrics", {
  r1 <- c(A = 0.5, B = 0.3, C = 0.2)
  r2 <- c(B = 0.6, C = 0.3, D = 0.1)
  grid <- optimize_weights(list(r1, r2), step_size = 0.5)

  expect_s3_class(grid, "data.frame")
  expect_true("w1" %in% names(grid))
  expect_true("w2" %in% names(grid))
  expect_true("ranking" %in% names(grid))
  # step_size 0.5 gives weights 0, 0.5, 1 = 3 combos
  expect_equal(nrow(grid), 3L)
  expect_equal(grid$w1 + grid$w2, rep(1, 3))
})

test_that("optimize_weights returns correct grid for 3 metrics", {
  r1 <- c(A = 0.5, B = 0.5)
  r2 <- c(B = 0.5, C = 0.5)
  r3 <- c(A = 0.3, C = 0.7)
  grid <- optimize_weights(list(r1, r2, r3), step_size = 0.5)

  expect_s3_class(grid, "data.frame")
  expect_true(all(c("w1", "w2", "w3") %in% names(grid)))
  # All weight triplets should sum to 1
  sums <- grid$w1 + grid$w2 + grid$w3
  expect_true(all(abs(sums - 1) < 1e-10))
  # Each ranking should be a named numeric vector
  expect_true(is.numeric(grid$ranking[[1]]))
})

test_that("optimize_weights validates input", {
  expect_error(optimize_weights(list(c(A = 1)), step_size = 0.5), "2 or 3")
})
