#' Read a gene set from a vector or file
#'
#' Flexible reader that turns a variety of inputs into a plain character
#' vector of gene symbols. Accepts an in-memory character vector, a
#' newline-delimited text file, or the two-column CSV format used by the
#' original paper (a row-index column plus a gene column named `x`).
#'
#' @param x One of:
#'   \describe{
#'     \item{a character vector}{Returned as-is after trimming and
#'       de-duplication.}
#'     \item{a path to a `.csv` file}{If a column named `x` exists (the legacy
#'       format written by `write.csv()`), that column is used; otherwise the
#'       first non-index column is used.}
#'     \item{a path to a `.txt` file}{Read as one gene per line.}
#'   }
#'
#' @return A character vector of unique gene symbols (order preserved).
#' @export
#'
#' @examples
#' read_gene_list(c("TP53", "KRAS", "BRAF"))
#' \dontrun{
#' read_gene_list("inst/extdata/example_genes.csv")
#' }
read_gene_list <- function(x) {
  if (!is.character(x)) {
    stop("`x` must be a character vector of genes or a path to a gene file")
  }
  # A character vector that is empty or has more than one element is treated as
  # gene symbols directly; so is a single value that is not a readable file.
  if (length(x) != 1L || !file.exists(x)) {
    return(clean_gene_vector(x))
  }

  path <- x
  is_csv <- grepl("\\.csv$", path, ignore.case = TRUE)

  if (is_csv) {
    df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    if (ncol(df) == 0L) stop("No columns found in gene file: ", path)
    if ("x" %in% colnames(df)) {
      genes <- df[["x"]]
    } else {
      # Fall back to the first column that is not an unnamed row index.
      first <- if (colnames(df)[1L] %in% c("", "X")) {
        if (ncol(df) >= 2L) df[[2L]] else df[[1L]]
      } else {
        df[[1L]]
      }
      genes <- first
    }
  } else {
    genes <- readLines(path, warn = FALSE)
  }

  clean_gene_vector(as.character(genes))
}

# Trim whitespace/quotes, drop empties, and de-duplicate while preserving order.
clean_gene_vector <- function(genes) {
  genes <- trimws(genes)
  genes <- gsub('^"|"$', "", genes)
  genes <- genes[nzchar(genes) & !is.na(genes)]
  unique(genes)
}
