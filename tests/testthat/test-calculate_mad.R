test_that("calculate_mad returns correct format", {
  mat <- matrix(
    c(1, 2, 3, 4, 5, 6, 7, 8, 9),
    nrow = 3,
    dimnames = list(c("g1", "g2", "g3"), c("c1", "c2", "c3"))
  )
  result <- calculate_mad(mat)

  expect_type(result, "double")
  expect_true(is.vector(result))
  expect_named(result)
  expect_length(result, 3L)
  expect_true(all(diff(result) <= 0))
  expect_equal(sum(result), 1)
})

test_that("calculate_mad matches manual calculation", {
  mat <- matrix(
    c(1, 2, 3, 4, 5, 6, 7, 8, 9),
    nrow = 3,
    dimnames = list(c("g1", "g2", "g3"), c("c1", "c2", "c3"))
  )
  manual <- c(
    g1 = stats::mad(c(1, 4, 7)),
    g2 = stats::mad(c(2, 5, 8)),
    g3 = stats::mad(c(3, 6, 9))
  )
  manual <- sort(manual, decreasing = TRUE)
  manual <- manual / sum(manual)

  expect_equal(calculate_mad(mat), manual)
})

test_that("calculate_mad normalize = FALSE returns raw values", {
  mat <- matrix(
    c(1, 10, 3, 4, 50, 6),
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("c1", "c2", "c3"))
  )
  result <- calculate_mad(mat, normalize = FALSE)
  expect_true(all(result >= 0))
  # Not normalized, so shouldn't sum to 1 in general
  manual <- c(g1 = stats::mad(c(1, 3, 50)), g2 = stats::mad(c(10, 4, 6)))
  manual <- sort(manual, decreasing = TRUE)
  expect_equal(result, manual)
})

test_that("calculate_mad handles single-row matrix", {
  mat <- matrix(1:5, nrow = 1, dimnames = list("g1", paste0("c", 1:5)))
  result <- calculate_mad(mat)
  expect_length(result, 1L)
  expect_equal(unname(result), 1)
})

test_that("calculate_mad handles constant rows (MAD = 0)", {
  mat <- matrix(5, nrow = 3, ncol = 4,
    dimnames = list(paste0("g", 1:3), paste0("c", 1:4)))
  result <- calculate_mad(mat)
  # All MADs are 0, total is 0, so no division occurs
  expect_true(all(result == 0))
})

test_that("calculate_mad rejects invalid input", {
  expect_error(calculate_mad(c(1, 2, 3)), "numeric matrix")
  expect_error(calculate_mad(list(a = 1)), "numeric matrix")
  expect_error(calculate_mad(NULL), "numeric matrix")
  expect_error(
    calculate_mad(matrix(character(0), nrow = 0, ncol = 0)), "numeric matrix"
  )
  expect_error(
    calculate_mad(matrix(0, nrow = 0, ncol = 0)), "must not be empty"
  )
})
