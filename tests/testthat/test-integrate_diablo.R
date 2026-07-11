make_blocks <- function(n = 30) {
  grp <- factor(rep(c("low", "high"), each = n / 2))
  rna <- matrix(stats::rnorm(n * 8), n,
    dimnames = list(paste0("s", seq_len(n)), paste0("g", 1:8)))
  mir <- matrix(stats::rnorm(n * 5), n,
    dimnames = list(paste0("s", seq_len(n)), paste0("mir", 1:5)))
  # g1 and mir1 carry the class signal
  rna[, "g1"] <- rna[, "g1"] + as.integer(grp) * 3
  mir[, "mir1"] <- mir[, "mir1"] - as.integer(grp) * 3
  list(blocks = list(rna = rna, mir = mir), grp = grp)
}

test_that("integrate_diablo returns a normalized per-feature ranking", {
  skip_if_not_installed("mixOmics")
  set.seed(1)
  d <- make_blocks()
  sc <- integrate_diablo(d$blocks, d$grp)

  expect_true(is.numeric(sc) && !is.null(names(sc)))
  expect_length(sc, 8 + 5) # all genes + miRNAs
  expect_true(all(sc >= 0))
  expect_equal(sum(sc), 1, tolerance = 1e-8)
  expect_false(is.unsorted(rev(sc))) # sorted decreasing
  # the discriminating features should rank above pure-noise ones
  expect_gt(sc["g1"], sc["g8"])
  expect_gt(sc["mir1"], sc["mir5"])
  expect_s3_class(attr(sc, "model"), "block.splsda")
})

test_that("integrate_diablo output plugs into score_rankings", {
  skip_if_not_installed("mixOmics")
  set.seed(2)
  d <- make_blocks()
  sc <- integrate_diablo(d$blocks, d$grp)
  genes <- paste0("g", 1:8)
  scored <- score_rankings(genes, list(diablo = sc[genes]))
  expect_equal(nrow(scored), length(genes))
  expect_true("combined" %in% names(scored))
})

test_that("integrate_diablo validates its inputs", {
  skip_if_not_installed("mixOmics")
  d <- make_blocks()
  expect_error(integrate_diablo(list(), d$grp), "non-empty named list")
  expect_error(integrate_diablo(unname(d$blocks), d$grp), "must be named")
  bad <- d$blocks
  bad$mir <- bad$mir[1:10, , drop = FALSE] # mismatched samples
  expect_error(integrate_diablo(bad, d$grp), "same number of samples")
  expect_error(integrate_diablo(d$blocks, d$grp[1:10]), "one label per sample")
  expect_error(
    integrate_diablo(d$blocks, factor(rep("only", nrow(d$blocks$rna)))),
    "at least 2 classes"
  )
})
