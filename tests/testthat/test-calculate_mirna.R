test_that("calculate_mirna returns correct ranking", {
  mat <- matrix(
    c(3, 0, 1, 2, 1, 0),
    nrow = 3,
    dimnames = list(c("TP53", "KRAS", "BRAF"), c("miR-1", "miR-2"))
  )
  result <- calculate_mirna(mat)

  expect_type(result, "double")
  expect_named(result)
  expect_length(result, 3L)
  # TP53 has 3+2=5, KRAS has 0+1=1, BRAF has 1+0=1
  expect_equal(names(result)[1], "TP53")
  expect_equal(sum(result), 1)
  expect_true(all(diff(result) <= 0))
})

test_that("calculate_mirna normalize = FALSE returns raw sums", {
  mat <- matrix(c(3, 1, 2, 0), nrow = 2,
    dimnames = list(c("A", "B"), c("m1", "m2")))
  result <- calculate_mirna(mat, normalize = FALSE)
  expect_equal(unname(result), c(5, 1))
  expect_equal(names(result), c("A", "B"))
})

test_that("calculate_mirna accepts data.frame input", {
  df <- data.frame(m1 = c(2, 1), m2 = c(0, 3), row.names = c("X", "Y"))
  result <- calculate_mirna(df)
  expect_length(result, 2L)
  expect_equal(sum(result), 1)
})

test_that("calculate_mirna validates input", {
  expect_error(calculate_mirna(c(1, 2)), "matrix or data.frame")
  expect_error(calculate_mirna(matrix(0, 0, 0)), "must not be empty")
  expect_error(
    calculate_mirna(matrix(0, nrow = 1, ncol = 0)), "must not be empty"
  )
})
