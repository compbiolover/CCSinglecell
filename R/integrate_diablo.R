#' Learn per-gene attributions with DIABLO (supervised multi-block)
#'
#' The heavyweight, supervised counterpart to [learn_weights()]. Where
#' `learn_weights()` blends already-computed per-metric rankings without an
#' outcome, `integrate_diablo()` learns directly from the raw omics blocks
#' *against a response*: it fits DIABLO (`mixOmics::block.splsda`), a sparse
#' multi-block PLS discriminant model, and turns the learned block loadings into
#' a single per-feature attribution ranking. Genes (and miRNAs) that the model
#' relies on to separate the response classes score highest.
#'
#' This is the DIABLO backend planned in `ROADMAP.md` (Phase 3b). It is gated
#' behind the optional `mixOmics` dependency. The unsupervised latent-factor
#' backend (MOFA2) is not included: its Python/reticulate backend is not
#' installable in a plain R environment, so it cannot be exercised end-to-end;
#' add it the same way once `MOFA2` is available.
#'
#' @param blocks A **named** list of numeric matrices, one per omics layer, with
#'   **samples as rows** and features (genes / miRNAs) as columns. Every block
#'   must have the same samples in the same order (matched by row name).
#' @param response A per-sample class label (factor or coercible vector), length
#'   equal to the number of rows in each block. For a survival use case, pass a
#'   binarised outcome (e.g. high- vs low-risk).
#' @param ncomp Number of DIABLO components to fit (default 2). Attributions
#'   aggregate the absolute loadings across all components.
#' @param keepX Optional per-block feature-selection budget passed to
#'   [mixOmics::block.splsda()] — a named list (matching `blocks`) of integer
#'   vectors of length `ncomp`. `NULL` (default) keeps all features.
#' @param normalize Logical; if `TRUE` (default) the returned attributions are
#'   divided by their sum so they sum to 1.
#' @param ... Further arguments forwarded to [mixOmics::block.splsda()] (e.g. a
#'   custom `design`).
#'
#' @return A named numeric vector of per-feature attributions (mean absolute
#'   loading across components, pooled across blocks), sorted decreasing and
#'   summing to 1 when `normalize = TRUE`. The fitted model is attached as
#'   attribute `"model"` and the raw per-block loadings as `"loadings"`. The
#'   result is a ready-to-use ranking: feed it to [score_rankings()] as one
#'   metric, or use it directly as a learned, outcome-aware integration.
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("mixOmics", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 30
#'   rna <- matrix(rnorm(n * 8), n, dimnames = list(paste0("s", 1:n), paste0("g", 1:8)))
#'   mir <- matrix(rnorm(n * 5), n, dimnames = list(paste0("s", 1:n), paste0("mir", 1:5)))
#'   grp <- factor(rep(c("low", "high"), each = n / 2))
#'   # make g1 and mir1 discriminate the groups
#'   rna[, "g1"] <- rna[, "g1"] + as.integer(grp) * 2
#'   mir[, "mir1"] <- mir[, "mir1"] - as.integer(grp) * 2
#'   attr_scores <- integrate_diablo(list(rna = rna, mir = mir), grp)
#'   head(attr_scores)
#' }
#' }
integrate_diablo <- function(blocks, response, ncomp = 2, keepX = NULL,
                             normalize = TRUE, ...) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("integrate_diablo() requires the 'mixOmics' package. Install it with ",
         "BiocManager::install('mixOmics').")
  }
  blocks <- validate_omics_blocks(blocks)
  block_names <- names(blocks)
  n_samples <- nrow(blocks[[1L]])

  response <- as.factor(response)
  if (length(response) != n_samples) {
    stop("`response` must have one label per sample (", n_samples, "); got ",
         length(response))
  }
  if (nlevels(response) < 2L) {
    stop("`response` must have at least 2 classes; got ", nlevels(response))
  }

  ncomp <- as.integer(ncomp)
  if (is.na(ncomp) || ncomp < 1L) stop("`ncomp` must be a positive integer")

  fit <- suppressMessages(
    mixOmics::block.splsda(X = blocks, Y = response, ncomp = ncomp,
                           keepX = keepX, ...)
  )

  # Aggregate absolute loadings across components into one score per feature,
  # then pool across blocks (dropping the response's own "Y" loadings).
  load_list <- fit$loadings[block_names]
  per_block <- lapply(load_list, function(L) rowMeans(abs(as.matrix(L))))
  scores <- unlist(unname(per_block))
  if (anyDuplicated(names(scores))) {
    scores <- vapply(split(scores, names(scores)), sum, numeric(1L))
  }
  scores <- sort(scores, decreasing = TRUE)

  if (normalize) {
    total <- sum(scores)
    if (total > 0) scores <- scores / total
  }

  attr(scores, "model") <- fit
  attr(scores, "loadings") <- load_list
  scores
}

# Validate a named list of sample-by-feature omics blocks (shared by
# integrate_diablo() and integrate_mofa()). Returns the coerced-to-matrix
# blocks; errors unless every block has the same samples in the same row order.
validate_omics_blocks <- function(blocks) {
  if (!is.list(blocks) || length(blocks) == 0L) {
    stop("`blocks` must be a non-empty named list of sample x feature matrices")
  }
  block_names <- names(blocks)
  if (is.null(block_names) || any(!nzchar(block_names))) {
    stop("`blocks` must be named (one name per omics layer)")
  }
  if (anyDuplicated(block_names)) {
    stop("block names must be unique; got: ", paste(block_names, collapse = ", "))
  }

  blocks <- lapply(blocks, function(m) {
    m <- as.matrix(m)
    if (!is.numeric(m)) stop("every block must be a numeric matrix")
    if (is.null(rownames(m)) || is.null(colnames(m))) {
      stop("every block must have sample row names and feature column names")
    }
    m
  })
  names(blocks) <- block_names

  ref_rows <- rownames(blocks[[1L]])
  for (nm in block_names) {
    if (nrow(blocks[[nm]]) != length(ref_rows)) {
      stop("all blocks must have the same number of samples (rows); ",
           "block '", nm, "' has ", nrow(blocks[[nm]]), ", expected ",
           length(ref_rows))
    }
    if (!identical(rownames(blocks[[nm]]), ref_rows)) {
      stop("all blocks must have the same samples in the same order (row ",
           "names); block '", nm, "' does not match '", block_names[1L], "'")
    }
  }
  blocks
}
