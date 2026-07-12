#' Learn per-gene attributions with MOFA (unsupervised latent factors)
#'
#' The unsupervised, latent-factor counterpart to [integrate_diablo()]. Where
#' DIABLO learns against a response, `integrate_mofa()` needs no outcome: it fits
#' MOFA (`MOFA2::run_mofa`), an unsupervised multi-omics factor model, and turns
#' the learned per-view feature weights into a single per-feature attribution
#' ranking. Genes (and miRNAs) that load strongly on the shared latent factors —
#' the axes of coordinated variation across omics layers — score highest.
#'
#' This is the MOFA backend from `ROADMAP.md` (Phase 3b). It is gated behind the
#' optional `MOFA2` package, which in turn drives a Python `mofapy2` backend via
#' `reticulate`. If no Python binding is configured, the function looks for a
#' virtualenv named `r-mofapy2` (created with `mofapy2` installed); otherwise it
#' respects `RETICULATE_PYTHON` / any Python you have already initialised, or set
#' `use_basilisk = TRUE` to let MOFA2 manage its own environment.
#'
#' @param blocks A **named** list of numeric matrices, one per omics layer, with
#'   **samples as rows** and features (genes / miRNAs) as columns. Every block
#'   must have the same samples in the same order (matched by row name). (This is
#'   the same orientation [integrate_diablo()] takes; MOFA's own
#'   features-by-samples layout is handled internally.)
#' @param nfactors Number of latent factors to fit (default 5). MOFA prunes
#'   factors that explain no variance, so the final attribution may pool over
#'   fewer.
#' @param normalize Logical; if `TRUE` (default) the returned attributions are
#'   divided by their sum so they sum to 1.
#' @param use_basilisk Passed to [MOFA2::run_mofa()]; if `TRUE`, MOFA2 manages
#'   its own Python environment via `basilisk` instead of the active
#'   `reticulate` binding. Default `FALSE`.
#' @param ... Further arguments merged into the MOFA **model** options (e.g.
#'   `likelihoods`, `spikeslab_weights`); see [MOFA2::get_default_model_options()].
#'
#' @return A named numeric vector of per-feature attributions (mean absolute
#'   factor weight across factors, pooled across views), sorted decreasing and
#'   summing to 1 when `normalize = TRUE`. The trained model is attached as
#'   attribute `"model"` and the raw per-view weights as `"weights"`. The result
#'   is a ready-to-use ranking: feed it to [score_rankings()] as one metric, or
#'   use it directly as a learned, unsupervised integration.
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("MOFA2", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 40
#'   rna <- matrix(rnorm(n * 8), n, dimnames = list(paste0("s", 1:n), paste0("g", 1:8)))
#'   mir <- matrix(rnorm(n * 5), n, dimnames = list(paste0("s", 1:n), paste0("mir", 1:5)))
#'   # a shared latent axis that couples g1 and mir1
#'   z <- rnorm(n)
#'   rna[, "g1"] <- 3 * z + rnorm(n, 0, 0.3)
#'   mir[, "mir1"] <- -3 * z + rnorm(n, 0, 0.3)
#'   attr_scores <- integrate_mofa(list(rna = rna, mir = mir))
#'   head(attr_scores)
#' }
#' }
integrate_mofa <- function(blocks, nfactors = 5, normalize = TRUE,
                           use_basilisk = FALSE, ...) {
  if (!requireNamespace("MOFA2", quietly = TRUE)) {
    stop("integrate_mofa() requires the 'MOFA2' package. Install it with ",
         "BiocManager::install('MOFA2'), plus a Python 'mofapy2' backend.")
  }
  blocks <- validate_omics_blocks(blocks) # shared with integrate_diablo()

  nfactors <- as.integer(nfactors)
  if (is.na(nfactors) || nfactors < 1L) stop("`nfactors` must be a positive integer")

  if (!use_basilisk) mofa_bind_python()

  # MOFA wants features-by-samples views; our blocks are samples-by-features.
  views <- lapply(blocks, t)

  object <- MOFA2::create_mofa(views)
  model_opts <- MOFA2::get_default_model_options(object)
  model_opts$num_factors <- nfactors
  extra <- list(...)
  for (nm in names(extra)) model_opts[[nm]] <- extra[[nm]]
  train_opts <- MOFA2::get_default_training_options(object)
  train_opts$verbose <- FALSE
  object <- MOFA2::prepare_mofa(object,
    model_options = model_opts, training_options = train_opts)

  outfile <- tempfile(fileext = ".hdf5")
  on.exit(unlink(outfile), add = TRUE)
  model <- suppressWarnings(
    MOFA2::run_mofa(object, outfile = outfile, use_basilisk = use_basilisk)
  )

  # Aggregate absolute factor weights into one score per feature, pooled across
  # views (get_weights() returns a features-by-factors matrix per view).
  weights <- MOFA2::get_weights(model)
  per_view <- lapply(weights, function(W) rowMeans(abs(as.matrix(W))))
  scores <- unlist(unname(per_view))
  if (anyDuplicated(names(scores))) {
    scores <- vapply(split(scores, names(scores)), sum, numeric(1L))
  }
  scores <- sort(scores, decreasing = TRUE)

  if (normalize) {
    total <- sum(scores)
    if (total > 0) scores <- scores / total
  }

  attr(scores, "model") <- model
  attr(scores, "weights") <- weights
  scores
}

# Opportunistically bind reticulate to a local mofapy2 virtualenv, yielding to
# any Python the user has already configured. Only acts when nothing is set and
# the conventional `r-mofapy2` env exists — a no-op otherwise.
mofa_bind_python <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(invisible())
  if (reticulate::py_available(initialize = FALSE)) return(invisible())
  if (nzchar(Sys.getenv("RETICULATE_PYTHON"))) return(invisible())
  if (reticulate::virtualenv_exists("r-mofapy2")) {
    reticulate::use_virtualenv("r-mofapy2", required = FALSE)
  }
  invisible()
}
