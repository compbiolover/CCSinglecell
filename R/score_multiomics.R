#' Score a gene set from an arbitrary set of gene rankings
#'
#' The generic engine behind [score_gene_set()] and [score_multiomics()].
#' Given a named list of per-metric gene rankings (one per omics layer /
#' metric), it restricts each to the requested genes, optionally renormalizes
#' within the gene set, and blends them into a single combined score with
#' [combine_rankings()].
#'
#' @param genes A character vector of gene symbols, or a path to a gene file
#'   (see [read_gene_list()]).
#' @param rankings A **named** list of named numeric vectors. Each element is
#'   one metric's genome-wide ranking; the list name becomes the metric's
#'   column name in the output.
#' @param weights Optional metric weights. Either a numeric vector (one per
#'   element of `rankings`, in list order; rescaled to sum to 1), the string
#'   `"learn"` to learn data-driven weights with [learn_weights()], or `NULL`
#'   (default) for equal weights.
#' @param renormalize Logical; if `TRUE` (default) each metric is renormalised
#'   to sum to 1 across the requested gene set before combining.
#'
#' @return A data frame with one row per gene, sorted by `combined`
#'   (decreasing): `gene`, one column per metric, `combined`, and `rank`. The
#'   weights used and metric names are attached as attributes `"weights"` and
#'   `"metrics"`.
#' @export
#'
#' @examples
#' r1 <- c(g1 = 0.5, g2 = 0.3, g3 = 0.2)
#' r2 <- c(g2 = 0.6, g3 = 0.3, g4 = 0.1)
#' score_rankings(c("g1", "g2", "g3"), list(expression = r1, methylation = r2))
score_rankings <- function(genes, rankings, weights = NULL, renormalize = TRUE) {
  genes <- read_gene_list(genes)
  if (length(genes) == 0L) stop("`genes` is empty")

  if (!is.list(rankings) || length(rankings) == 0L) {
    stop("`rankings` must be a non-empty named list of gene rankings")
  }
  if (is.null(names(rankings)) || any(!nzchar(names(rankings)))) {
    stop("`rankings` must be named (one name per metric / omics layer)")
  }
  metrics <- names(rankings)
  if (anyDuplicated(metrics)) {
    stop("metric names must be unique; got: ", paste(metrics, collapse = ", "))
  }
  reserved <- c("gene", "combined", "rank")
  if (any(metrics %in% reserved)) {
    stop("metric names must not collide with reserved output columns (",
         paste(reserved, collapse = ", "), ")")
  }
  for (m in metrics) {
    r <- rankings[[m]]
    if (!is.numeric(r) || is.null(names(r))) {
      stop("ranking '", m, "' must be a named numeric vector")
    }
    if (any(!nzchar(names(r))) || anyNA(names(r))) {
      stop("ranking '", m, "' has empty or missing gene names")
    }
  }

  missing_all <- setdiff(genes, unique(unlist(lapply(rankings, names))))
  if (length(missing_all) > 0L) {
    warning(
      length(missing_all), " of ", length(genes),
      " requested gene(s) were not found in any metric and scored 0 ",
      "(e.g. ", paste(utils::head(missing_all, 3L), collapse = ", "), ")."
    )
  }

  subset_metric <- function(r, metric) {
    v <- r[genes]
    v[is.na(v)] <- 0
    names(v) <- genes
    if (renormalize) {
      if (any(v < 0)) {
        stop("renormalize = TRUE requires non-negative scores, but metric '",
             metric, "' has negative values. Set renormalize = FALSE for ",
             "signed metrics.")
      }
      total <- sum(v)
      if (total > 0) v <- v / total
    }
    v
  }
  rankings <- Map(subset_metric, rankings, metrics)

  n_metrics <- length(rankings)
  learned_var_explained <- NULL
  if (is.null(weights)) {
    weights <- rep(1 / n_metrics, n_metrics)
  } else if (is.character(weights)) {
    if (length(weights) != 1L || weights != "learn") {
      stop("character `weights` must be \"learn\"")
    }
    learned <- learn_weights(rankings, method = "pca")
    learned_var_explained <- attr(learned, "variance_explained")
    weights <- as.numeric(learned[metrics])
  } else {
    if (length(weights) != n_metrics) {
      stop(
        "`weights` must have one value per metric (",
        n_metrics, ": ", paste(metrics, collapse = ", "), ")"
      )
    }
    if (any(weights < 0)) stop("`weights` must be non-negative")
    total <- sum(weights)
    if (total <= 0) stop("`weights` must not all be zero")
    weights <- weights / total
  }

  combined <- if (n_metrics == 1L) {
    rankings[[1L]]
  } else {
    combine_rankings(rankings, weights)
  }

  out <- data.frame(gene = genes, stringsAsFactors = FALSE)
  for (m in metrics) {
    out[[m]] <- unname(rankings[[m]][genes])
  }
  out$combined <- unname(combined[genes])
  out <- out[order(out$combined, decreasing = TRUE), , drop = FALSE]
  out$rank <- seq_len(nrow(out))
  rownames(out) <- NULL

  weights_attr <- stats::setNames(weights, metrics)
  if (!is.null(learned_var_explained)) {
    attr(weights_attr, "variance_explained") <- learned_var_explained
  }
  attr(out, "weights") <- weights_attr
  attr(out, "metrics") <- metrics
  out
}


#' Describe one omics layer for multi-omics scoring
#'
#' Bundles an omics data object with a *metric* that turns it into a
#' genome-wide, per-gene ranking. Pass a list of these to [score_multiomics()].
#' This is the extension point that lets any omics layer participate in scoring
#' without changing the integration core.
#'
#' @param name Character label for the layer (becomes the metric's column name
#'   in the scored output).
#' @param data The omics data the metric consumes (e.g. an expression matrix,
#'   a miRNA interaction matrix, or a list of inputs for metrics that need
#'   several — such as `list(expression = ..., pseudotime = ...)`).
#' @param metric Either a built-in metric name (`"mad"`, `"mirna"`, or
#'   `"switchde"`) or a function taking `data` and returning a **named numeric
#'   vector** of gene scores (higher = more important).
#'
#' @return An `omics_block` object.
#' @export
#'
#' @examples
#' set.seed(1)
#' expr <- matrix(rpois(200, 3), nrow = 20,
#'   dimnames = list(paste0("g", 1:20), paste0("c", 1:10)))
#' blk <- omics_block("expression", expr, metric = "mad")
#' # A custom metric: rank genes by mean expression
#' blk2 <- omics_block("mean_expr", expr, metric = function(m) sort(rowMeans(m), decreasing = TRUE))
omics_block <- function(name, data, metric) {
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    stop("`name` must be a non-empty string")
  }
  if (is.character(metric)) {
    metric <- builtin_metric(metric)
  }
  if (!is.function(metric)) {
    stop("`metric` must be a function or a built-in metric name")
  }
  structure(
    list(name = name, data = data, metric = metric),
    class = "omics_block"
  )
}

# Resolve a built-in metric name to a scoring function `data -> named ranking`.
builtin_metric <- function(name) {
  switch(
    name,
    mad = function(data) calculate_mad(data, normalize = TRUE),
    mirna = function(data) calculate_mirna(data, normalize = TRUE),
    mirna_activity = function(data) {
      calculate_mirna_activity(
        data$mirna_expr, data$gene_expr, data$target_matrix,
        normalize = TRUE
      )
    },
    switchde = function(data) {
      calculate_switchde(data$expression, data$pseudotime, normalize = TRUE)
    },
    stop("Unknown built-in metric: '", name,
         "'. Known metrics: mad, mirna, mirna_activity, switchde")
  )
}

#' Score a gene set across an arbitrary number of omics layers
#'
#' Truly multi-omics scoring: evaluate each [omics_block()] to a genome-wide
#' gene ranking, then blend them into one combined score. Unlike
#' [score_gene_set()] (which is limited to the built-in expression / pseudotime
#' / miRNA metrics), this accepts any number of layers with any metrics.
#'
#' @param genes A character vector of gene symbols, or a path to a gene file
#'   (see [read_gene_list()]).
#' @param blocks A list of [omics_block()] objects (one per omics layer).
#' @param weights Optional per-block weights: a numeric vector (block order),
#'   the string `"learn"` for data-driven weights (see [learn_weights()]), or
#'   `NULL` (default) for equal weights.
#' @param renormalize Logical; passed to [score_rankings()].
#'
#' @return A scored data frame, identical in shape to [score_gene_set()] but
#'   with one column per block.
#' @export
#'
#' @examples
#' set.seed(1)
#' expr <- matrix(rpois(400, 3), nrow = 40,
#'   dimnames = list(paste0("g", 1:40), paste0("c", 1:10)))
#' mir <- matrix(sample(0:3, 120, replace = TRUE), nrow = 40,
#'   dimnames = list(paste0("g", 1:40), paste0("mir", 1:3)))
#' meth <- matrix(runif(400), nrow = 40,
#'   dimnames = list(paste0("g", 1:40), paste0("c", 1:10)))
#' blocks <- list(
#'   omics_block("expression", expr, "mad"),
#'   omics_block("mirna", mir, "mirna"),
#'   omics_block("methylation", meth, function(m) sort(apply(m, 1, stats::sd), decreasing = TRUE))
#' )
#' score_multiomics(paste0("g", 1:40), blocks)
score_multiomics <- function(genes, blocks, weights = NULL, renormalize = TRUE) {
  if (!is.list(blocks) || length(blocks) == 0L) {
    stop("`blocks` must be a non-empty list of omics_block objects")
  }
  is_block <- vapply(blocks, inherits, logical(1L), what = "omics_block")
  if (!all(is_block)) {
    stop("every element of `blocks` must be an omics_block (see omics_block())")
  }

  block_names <- vapply(blocks, function(b) b$name, character(1L))
  if (anyDuplicated(block_names)) {
    stop("omics_block names must be unique; got: ",
         paste(block_names, collapse = ", "))
  }

  rankings <- lapply(blocks, function(b) {
    r <- b$metric(b$data)
    if (!is.numeric(r) || is.null(names(r))) {
      stop("metric for block '", b$name,
           "' must return a named numeric vector")
    }
    r
  })
  names(rankings) <- block_names

  score_rankings(genes, rankings, weights = weights, renormalize = renormalize)
}
