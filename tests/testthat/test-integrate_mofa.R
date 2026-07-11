make_mofa_blocks <- function(n = 40) {
  z <- stats::rnorm(n)
  rna <- matrix(stats::rnorm(n * 8), n,
    dimnames = list(paste0("s", seq_len(n)), paste0("g", 1:8)))
  mir <- matrix(stats::rnorm(n * 5), n,
    dimnames = list(paste0("s", seq_len(n)), paste0("mir", 1:5)))
  # g1 and mir1 share a latent axis
  rna[, "g1"] <- 3 * z + stats::rnorm(n, 0, 0.3)
  mir[, "mir1"] <- -3 * z + stats::rnorm(n, 0, 0.3)
  list(rna = rna, mir = mir)
}

test_that("integrate_mofa returns a normalized per-feature ranking", {
  skip_if_not_installed("MOFA2")
  skip_if_not_installed("reticulate")
  skip_if_not(reticulate::py_module_available("mofapy2"), "mofapy2 not available")
  set.seed(1)
  b <- make_mofa_blocks()
  sc <- integrate_mofa(b, nfactors = 3)

  expect_true(is.numeric(sc) && !is.null(names(sc)))
  expect_length(sc, 8 + 5)
  expect_true(all(sc >= 0))
  expect_equal(sum(sc), 1, tolerance = 1e-8)
  expect_false(is.unsorted(rev(as.numeric(sc))))
  # the shared-latent features should dominate their views
  expect_equal(names(which.max(sc[paste0("g", 1:8)])), "g1")
  expect_equal(names(which.max(sc[paste0("mir", 1:5)])), "mir1")
  expect_s4_class(attr(sc, "model"), "MOFA")
})

test_that("integrate_mofa output plugs into score_rankings", {
  skip_if_not_installed("MOFA2")
  skip_if_not_installed("reticulate")
  skip_if_not(reticulate::py_module_available("mofapy2"), "mofapy2 not available")
  set.seed(2)
  b <- make_mofa_blocks()
  sc <- integrate_mofa(b, nfactors = 3)
  genes <- paste0("g", 1:8)
  scored <- score_rankings(genes, list(mofa = sc[genes]))
  expect_equal(nrow(scored), length(genes))
  expect_true("combined" %in% names(scored))
})

test_that("integrate_mofa validates its inputs before training", {
  skip_if_not_installed("MOFA2")
  b <- make_mofa_blocks()
  expect_error(integrate_mofa(list()), "non-empty named list")
  expect_error(integrate_mofa(unname(b)), "must be named")
  bad <- b
  bad$mir <- bad$mir[1:10, , drop = FALSE]
  expect_error(integrate_mofa(bad), "same number of samples")
  expect_error(integrate_mofa(b, nfactors = 0), "positive integer")
})
