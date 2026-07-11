#' Command-line interface for gene scoring
#'
#' Programmatic entry point behind the `ccscore` command-line tool (see
#' `inst/bin/ccscore`). Parses arguments, runs the requested analysis, and
#' writes tables and plots to an output directory. Exposed as a function so
#' the CLI logic is unit-testable.
#'
#' Supported commands:
#' \describe{
#'   \item{`score`}{Score a gene set against an expression matrix and write
#'     `scores.csv` plus ranking and contribution plots.}
#'   \item{`survival`}{Fit a penalised Cox model on a gene set, compute risk
#'     scores, and write `risk_scores.csv` plus a Kaplan-Meier plot.}
#' }
#'
#' @param args Character vector of command-line arguments (defaults to
#'   [commandArgs()] with `trailingOnly = TRUE`).
#'
#' @return Invisibly, a list of the paths written (or `NULL` for help).
#' @export
#'
#' @examples
#' \dontrun{
#' run_cli(c("score",
#'   "--genes", "genes.csv",
#'   "--expression", "expr.csv",
#'   "--out", "results"))
#' }
run_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) == 0L || args[1L] %in% c("-h", "--help", "help")) {
    cat(cli_usage())
    return(invisible(NULL))
  }

  command <- args[1L]
  opts <- parse_cli_opts(args[-1L])

  result <- switch(
    command,
    score = cli_score(opts),
    survival = cli_survival(opts),
    stop("Unknown command: '", command, "'\n\n", cli_usage())
  )
  invisible(result)
}

cli_usage <- function() {
  paste0(
    "ccscore - score gene sets from single-cell data (CcSinglecell)\n\n",
    "Usage:\n",
    "  ccscore score    --genes FILE --expression FILE --out DIR \\\n",
    "                   [--pseudotime FILE] [--mirna FILE] \\\n",
    "                   [--weights w1,w2,...] [--top N]\n",
    "  ccscore survival --genes FILE --survival FILE --out DIR [--alpha A]\n",
    "  ccscore --help\n\n",
    "score options:\n",
    "  --genes FILE       Gene set (CSV with an 'x' column, or one gene/line)\n",
    "  --expression FILE  Expression matrix CSV (genes in rows, cells in cols,\n",
    "                     first column = gene names)\n",
    "  --pseudotime FILE  Optional pseudotime per cell (enables switchDE)\n",
    "  --mirna FILE       Optional miRNA interaction matrix (enables miRNA)\n",
    "  --weights LIST     Comma-separated metric weights (default: equal)\n",
    "  --top N            Number of genes to show in plots (default 20)\n",
    "  --out DIR          Output directory (created if missing)\n\n",
    "survival options:\n",
    "  --genes FILE       Gene set to use as predictors\n",
    "  --survival FILE    Cox data frame CSV: gene columns plus\n",
    "                     days_to_last_follow_up and vital_status\n",
    "  --alpha A          Elastic-net mixing (1 = lasso, default; 0 = ridge)\n",
    "  --out DIR          Output directory (created if missing)\n"
  )
}

# Parse "--key value" pairs and bare "--flag" switches into a named list.
parse_cli_opts <- function(args) {
  opts <- list()
  i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    if (!grepl("^--", a)) {
      stop("Expected an option starting with '--' but got: '", a, "'")
    }
    key <- sub("^--", "", a)
    if (i + 1L <= length(args) && !grepl("^--", args[i + 1L])) {
      opts[[key]] <- args[i + 1L]
      i <- i + 2L
    } else {
      opts[[key]] <- TRUE
      i <- i + 1L
    }
  }
  opts
}

# Require an option to be present.
require_opt <- function(opts, name, command) {
  if (is.null(opts[[name]])) {
    stop("'", command, "' requires --", name, "\n\n", cli_usage())
  }
  opts[[name]]
}

# Read an expression matrix CSV (first column = gene names).
read_expression_csv <- function(path) {
  df <- utils::read.csv(path, row.names = 1L, check.names = FALSE)
  mat <- as.matrix(df)
  if (!is.numeric(mat)) {
    storage.mode(mat) <- "double"
  }
  mat
}

# Read a pseudotime vector aligned to the expression matrix columns. Accepts a
# file with a numeric pseudotime column and, optionally, a cell-ID column used
# to align to `cell_ids`. Falls back to positional order only when no ID column
# matches, and then requires the lengths to agree. Errors on NA so switchDE is
# never computed on silently misaligned pseudotime.
read_pseudotime <- function(path, cell_ids) {
  pt <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)

  # Choose the value column: prefer one named "pseudotime", else the last
  # numeric column.
  lower <- tolower(names(pt))
  if ("pseudotime" %in% lower) {
    value_col <- which(lower == "pseudotime")[1L]
  } else {
    numeric_cols <- which(vapply(pt, is.numeric, logical(1L)))
    if (length(numeric_cols) == 0L) {
      stop("No numeric pseudotime column found in ", path)
    }
    value_col <- numeric_cols[length(numeric_cols)]
  }
  values <- as.numeric(pt[[value_col]])

  # Align by any non-value column whose entries cover all expression cell IDs.
  id_col <- NULL
  for (j in seq_along(pt)) {
    if (j == value_col) next
    if (all(cell_ids %in% as.character(pt[[j]]))) {
      id_col <- j
      break
    }
  }
  if (!is.null(id_col)) {
    values <- values[match(cell_ids, as.character(pt[[id_col]]))]
  } else if (length(values) != length(cell_ids)) {
    stop("Pseudotime has ", length(values), " rows but the expression matrix ",
         "has ", length(cell_ids), " cells, and no cell-ID column matched to ",
         "align them. Provide a cell-ID column or matching row order.")
  }
  if (anyNA(values)) {
    stop("Pseudotime contains missing (NA) values after alignment to cells.")
  }
  values
}

# Read a survival cohort CSV, using a leading ID column as row names when
# present (e.g. patient IDs written by `write.csv(row.names = TRUE)`), so risk
# outputs keep patient identifiers instead of row numbers.
read_cohort_csv <- function(path) {
  df <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  first <- names(df)[1L]
  looks_like_ids <- first %in% c("", "X") ||
    (is.character(df[[1L]]) && anyDuplicated(df[[1L]]) == 0L)
  if (looks_like_ids) {
    rownames(df) <- as.character(df[[1L]])
    df[[1L]] <- NULL
  }
  df
}

cli_score <- function(opts) {
  genes_path <- require_opt(opts, "genes", "score")
  expr_path <- require_opt(opts, "expression", "score")
  out_dir <- require_opt(opts, "out", "score")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  genes <- read_gene_list(genes_path)
  expr <- read_expression_csv(expr_path)

  pseudotime <- NULL
  if (!is.null(opts[["pseudotime"]])) {
    pseudotime <- read_pseudotime(opts[["pseudotime"]], colnames(expr))
  }

  mirna <- NULL
  if (!is.null(opts[["mirna"]])) {
    mirna <- as.matrix(utils::read.csv(
      opts[["mirna"]], row.names = 1L, check.names = FALSE
    ))
  }

  weights <- NULL
  if (!is.null(opts[["weights"]])) {
    weights <- as.numeric(strsplit(opts[["weights"]], ",")[[1L]])
  }

  top_n <- if (!is.null(opts[["top"]])) as.integer(opts[["top"]]) else 20L

  scores <- score_gene_set(
    genes, expr,
    pseudotime = pseudotime,
    mirna_matrix = mirna,
    weights = weights
  )

  scores_path <- file.path(out_dir, "scores.csv")
  utils::write.csv(scores, scores_path, row.names = FALSE)
  written <- c(scores = scores_path)

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p1 <- file.path(out_dir, "gene_scores.png")
    ggplot2::ggsave(p1, plot_gene_scores(scores, top_n = top_n),
                    width = 7, height = 6, dpi = 120)
    written["gene_scores"] <- p1

    if (length(attr(scores, "metrics")) > 1L) {
      p2 <- file.path(out_dir, "score_contributions.png")
      ggplot2::ggsave(p2, plot_score_contributions(scores, top_n = top_n),
                      width = 7, height = 6, dpi = 120)
      written["score_contributions"] <- p2
    }
  } else {
    message("ggplot2 not available; skipping plots.")
  }

  message("Scored ", nrow(scores), " genes using metric(s): ",
          paste(attr(scores, "metrics"), collapse = ", "))
  message("Wrote: ", paste(written, collapse = ", "))
  invisible(as.list(written))
}

cli_survival <- function(opts) {
  genes_path <- require_opt(opts, "genes", "survival")
  surv_path <- require_opt(opts, "survival", "survival")
  out_dir <- require_opt(opts, "out", "survival")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  genes <- read_gene_list(genes_path)
  cox_df <- read_cohort_csv(surv_path)
  alpha <- if (!is.null(opts[["alpha"]])) as.numeric(opts[["alpha"]]) else 1

  model <- fit_cox_model(cox_df, gene_names = genes, alpha = alpha)
  if (length(model$active_genes) == 0L) {
    stop("The Cox model selected no genes; try a larger gene set or ",
         "different --alpha.")
  }

  risk <- calculate_risk_scores(cox_df, model$active_genes, model$coefficients)

  risk_path <- file.path(out_dir, "risk_scores.csv")
  utils::write.csv(risk, risk_path, row.names = TRUE)
  written <- c(risk_scores = risk_path)

  cindex <- cox_cindex(model$cv_fit)
  if (!is.na(cindex)) message("Cross-validated C-index: ", round(cindex, 4))
  message("Active genes: ", length(model$active_genes))

  if (requireNamespace("survminer", quietly = TRUE)) {
    km_path <- file.path(out_dir, "km_curve.png")
    km <- plot_km_curve(risk)
    grDevices::png(km_path, width = 800, height = 700)
    print(km)
    grDevices::dev.off()
    written["km_curve"] <- km_path
  } else {
    message("survminer not available; skipping Kaplan-Meier plot.")
  }

  message("Wrote: ", paste(written, collapse = ", "))
  invisible(as.list(written))
}

# Extract the cross-validated C-index at lambda.min from a cv.glmnet fit.
cox_cindex <- function(cv_fit) {
  idx <- which(cv_fit$lambda == cv_fit$lambda.min)
  if (length(idx) == 1L && !is.null(cv_fit$cvm)) cv_fit$cvm[idx] else NA_real_
}
