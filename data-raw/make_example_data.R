# Generate small, self-contained example data shipped in inst/extdata/.
#
# The real study data (GSE81861 counts, TCGA survival frames) is large and
# tracked with git-LFS, so it is unsuitable for a package demo. This script
# builds a small, seeded, synthetic-but-realistic data set that exercises
# every scoring metric and the survival path. Re-run with:
#   Rscript data-raw/make_example_data.R

set.seed(42)

out_dir <- file.path("inst", "extdata")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Gene universe ----------------------------------------------------------
# A handful of real colorectal-cancer-associated symbols carry the signal;
# the rest are generic filler so rankings have something to separate.
signal_genes <- c(
  "TP53", "KRAS", "APC", "SMAD4", "PIK3CA", "BRAF", "NRAS", "CTNNB1",
  "MYC", "EGFR", "PTEN", "FBXW7", "TCF7L2", "SOX9", "AXIN2", "LGR5",
  "MKI67", "CDKN2A", "VEGFA", "MMP7"
)
filler_genes <- sprintf("GENE%04d", seq_len(180))
genes <- c(signal_genes, filler_genes)
n_genes <- length(genes)

# --- Single-cell expression matrix (genes x cells) --------------------------
n_cells <- 60L
pseudotime <- sort(runif(n_cells, 0, 1))

expr <- matrix(
  rpois(n_genes * n_cells, lambda = 2),
  nrow = n_genes, ncol = n_cells,
  dimnames = list(genes, sprintf("cell%02d", seq_len(n_cells)))
)

# Highly variable genes (large MAD).
for (g in c("TP53", "KRAS", "MYC", "MKI67", "VEGFA")) {
  expr[g, ] <- rpois(n_cells, lambda = rep(c(1, 12), each = n_cells / 2))
}

# Switch-like genes along pseudotime (sigmoid up / down).
sigmoid <- function(x, k, m) 1 / (1 + exp(-k * (x - m)))
expr["LGR5", ]  <- round(20 * sigmoid(pseudotime, 12, 0.4)) + rpois(n_cells, 1)
expr["SOX9", ]  <- round(18 * sigmoid(pseudotime, 10, 0.6)) + rpois(n_cells, 1)
expr["AXIN2", ] <- round(15 * (1 - sigmoid(pseudotime, 12, 0.5))) + rpois(n_cells, 1)
expr["CDKN2A", ] <- round(16 * sigmoid(pseudotime, 14, 0.7)) + rpois(n_cells, 1)

write.csv(expr, file.path(out_dir, "example_expression.csv"))

# --- Pseudotime per cell ----------------------------------------------------
write.csv(
  data.frame(cell = colnames(expr), pseudotime = pseudotime),
  file.path(out_dir, "example_pseudotime.csv"),
  row.names = FALSE
)

# --- miRNA interaction matrix (genes x miRNAs) ------------------------------
mirnas <- c("miR-21", "miR-34a", "miR-143", "miR-145", "miR-200c")
mirna_mat <- matrix(
  rpois(n_genes * length(mirnas), lambda = 0.3),
  nrow = n_genes, ncol = length(mirnas),
  dimnames = list(genes, mirnas)
)
# Signal genes are more heavily targeted.
mirna_mat[signal_genes, ] <- mirna_mat[signal_genes, ] +
  matrix(rpois(length(signal_genes) * length(mirnas), lambda = 2),
         nrow = length(signal_genes))
write.csv(mirna_mat, file.path(out_dir, "example_mirna.csv"))

# --- Gene set to score (mix of signal + filler) -----------------------------
gene_set <- c(signal_genes, sample(filler_genes, 30))
write.csv(
  data.frame(x = gene_set),
  file.path(out_dir, "example_genes.csv"),
  row.names = TRUE
)

# --- Survival (patients x genes + clinical) ---------------------------------
n_patients <- 90L
surv_genes <- signal_genes
patient_expr <- matrix(
  rnorm(n_patients * length(surv_genes), mean = 6, sd = 2),
  nrow = n_patients, ncol = length(surv_genes),
  dimnames = list(sprintf("TCGA-%03d", seq_len(n_patients)), surv_genes)
)

# Linear predictor: a few genes drive risk.
beta <- stats::setNames(rep(0, length(surv_genes)), surv_genes)
beta[c("TP53", "MYC", "VEGFA")] <- c(0.35, 0.30, 0.25)
beta[c("APC", "SMAD4")] <- c(-0.30, -0.25)
lp <- as.numeric(scale(patient_expr %*% beta))

base_time <- rexp(n_patients, rate = 0.02 * exp(lp))
censor_time <- runif(n_patients, 200, 3000)
time <- pmin(base_time, censor_time)
vital_status <- as.integer(base_time <= censor_time)

surv_df <- data.frame(
  patient_expr,
  days_to_last_follow_up = round(time),
  vital_status = vital_status,
  check.names = FALSE
)
write.csv(surv_df, file.path(out_dir, "example_survival.csv"))

message("Wrote example data to ", normalizePath(out_dir))
message("  genes: ", n_genes, "  cells: ", n_cells,
        "  patients: ", n_patients)
