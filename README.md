# CcSinglecell

An R package for scoring sets of genes from single-cell RNA-seq and miRNA
targeting data, and evaluating how well they predict colorectal cancer
survival. It grew out of the analysis in:

> Andrew Willems, Nicholas Panchy, and Tian Hong. Using Single-Cell RNA
> Sequencing and MicroRNA Targeting Data to Improve Colorectal Cancer Survival
> Prediction. (2023) *Cells* 12(2):228

The package turns that study's scattered scripts into reusable functions, a
high-level `score_gene_set()` wrapper, ggplot2 visualizations, and a `ccscore`
command-line tool.

## Installation

```r
# install.packages("remotes")
remotes::install_github("compbiolover/CcSinglecell")
```

## Score a gene set in R

```r
library(CcSinglecell)

ext   <- function(f) system.file("extdata", f, package = "CcSinglecell")
expr  <- as.matrix(read.csv(ext("example_expression.csv"), row.names = 1, check.names = FALSE))
mirna <- as.matrix(read.csv(ext("example_mirna.csv"),      row.names = 1, check.names = FALSE))
genes <- read_gene_list(ext("example_genes.csv"))

scores <- score_gene_set(genes, expr, mirna_matrix = mirna)
head(scores)

plot_gene_scores(scores, top_n = 15)
plot_score_contributions(scores, top_n = 15)
```

Each metric is optional and additive:

| Metric          | Enabled by         | Captures                                  |
|-----------------|--------------------|-------------------------------------------|
| MAD             | always             | expression variability                    |
| switchDE        | `pseudotime =`     | switch-like dynamics along a trajectory   |
| miRNA           | `mirna_matrix =`   | cancer-miRNA targeting (predicted counts) |
| miRNA activity  | `omics_block(…, "mirna_activity")` | *observed* miRNA repression (anti-correlation) |

`calculate_mirna_activity()` is the expression-anchored upgrade to the simple
target-count metric: given matched miRNA and gene expression plus a predicted
interaction map, it scores each gene by how strongly it is actually
anti-correlated with its targeting cancer miRNAs in your samples (a
ceRNA/repression-evidence view). Unlike the other metrics it is not a
`score_gene_set()` argument; use it via `score_multiomics()` as its own block:

```r
omics_block("mirna_activity",
  list(mirna_expr = mir_expr, gene_expr = gene_expr, target_matrix = targets),
  "mirna_activity")
```

Present metrics are blended (equal weights by default, or pass `weights =`).

## Score across arbitrary omics layers

`score_gene_set()` covers the built-in metrics. For **truly multi-omics**
scoring — any number of layers, any scoring function — use `score_multiomics()`
with one `omics_block()` per layer. A block pairs an omics data object with a
metric (a built-in name or your own `data -> named gene ranking` function):

```r
blocks <- list(
  omics_block("expression",  expr, "mad"),
  omics_block("mirna",        mir,  "mirna"),
  omics_block("methylation",  meth, function(m) sort(apply(m, 1, sd), decreasing = TRUE)),
  omics_block("cnv",          cnv,  function(m) sort(rowMeans(abs(m)), decreasing = TRUE))
)
scores <- score_multiomics(genes, blocks)      # weights default to equal
plot_score_contributions(scores)               # per-layer decomposition
```

`combine_rankings()` and `optimize_weights()` are no longer capped at three
metrics — they generalize to any number of layers. See `ROADMAP.md` for where
this is heading (heavyweight MOFA/DIABLO backends, rigorous survival
evaluation).

### Learn the weights instead of guessing them

Rather than equal weights or a grid search, let the data set the blend:

```r
scores <- score_multiomics(genes, blocks, weights = "learn")
attr(scores, "weights")   # data-driven weights, one per layer
```

`weights = "learn"` (also accepted by `score_gene_set()`) calls
`learn_weights()`, which runs a PCA on the gene-by-metric score matrix and
weights each layer by its loading on the dominant shared axis of variation —
so co-varying, informative metrics count more and flat/idiosyncratic ones count
less. It's dependency-free; the heavier `MOFA2`/`mixOmics` (DIABLO) backends are
on the roadmap.

## Score a gene set from the command line

The package installs an executable `ccscore`:

```sh
# locate it
Rscript -e 'cat(system.file("bin", "ccscore", package = "CcSinglecell"))'

# score genes -> scores.csv + ranking/contribution plots
ccscore score \
  --genes      example_genes.csv \
  --expression example_expression.csv \
  --mirna      example_mirna.csv \
  --out        results/

# survival stratification -> risk_scores.csv + Kaplan-Meier plot
ccscore survival \
  --genes    example_genes.csv \
  --survival example_survival.csv \
  --out      results/

ccscore --help
```

## Key functions

- `score_gene_set()` — score an arbitrary gene set (main entry point)
- `read_gene_list()` — read genes from a vector, CSV, or text file
- `calculate_mad()`, `calculate_switchde()`, `calculate_mirna()` — the metrics
- `combine_rankings()`, `optimize_weights()` — blend / tune metric weights
- `fit_cox_model()`, `calculate_risk_scores()` — survival modeling
- `plot_gene_scores()`, `plot_score_contributions()`, `plot_km_curve()`,
  `plot_weight_optimization()` — visualizations

See `vignette("scoring-genes", package = "CcSinglecell")` for a full walkthrough.

## Project organization

```
CcSinglecell/
├── R/                  # Package functions (scoring, modeling, plotting, CLI)
├── man/                # Documentation (auto-generated by roxygen2)
├── vignettes/          # Usage tutorials
├── tests/testthat/     # Unit tests
├── data-raw/           # Script that generates the bundled example data
├── inst/
│   ├── bin/ccscore     # Command-line tool
│   ├── extdata/        # Small example data (expression, miRNA, survival, genes)
│   └── legacy/         # Original paper code + data for reproducibility
├── DESCRIPTION
└── NAMESPACE
```
