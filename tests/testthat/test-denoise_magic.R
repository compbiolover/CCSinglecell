test_that("denoise_magic validates inputs", {
  expect_error(denoise_magic(c(1, 2, 3)), "numeric matrix")
  expect_error(denoise_magic(matrix(0, 0, 0)), "must not be empty")
  expect_error(
    denoise_magic(matrix(character(0), 0, 0)), "numeric matrix"
  )
})

test_that("denoise_magic requires Rmagic package", {
  skip_if(requireNamespace("Rmagic", quietly = TRUE))
  mat <- matrix(rnorm(20), nrow = 4, ncol = 5)
  expect_error(denoise_magic(mat), "Rmagic")
})
