# CcSinglecell roadmap: toward truly end-to-end multi-omics

This document maps how `CcSinglecell` grows from a single-cell + miRNA gene
scorer into a genuine multi-omics biomarker-discovery toolkit. The guiding
principle is **integration, not enumeration**: the goal is not to mine yet
another cohort-specific signature, but to give users a principled, interpretable
pipeline that takes *any* set of omics layers from raw data to a validated,
explainable prognostic model.

The design is grounded in where the field has moved since the original 2023
paper (see *References* below): from static, prediction-based interaction counts
toward expression-anchored regulatory activity; from hand-tuned weighted sums
toward learned latent-factor and multi-block integration; and from black-box
risk scores back toward per-feature interpretability.

## Design principles

1. **Omics-agnostic core.** Nothing in the integration machinery should hard-code
   "MAD / switchDE / miRNA." Any omics layer (methylation, CNV, proteomics,
   ATAC, …) plugs in as a *block* with a *metric* that turns it into a per-gene
   ranking.
2. **Activity over annotation.** Metrics should reflect what is happening in the
   user's samples (e.g. observed miRNA–mRNA anti-correlation), not just what a
   database predicts could happen.
3. **Learned integration is a first-class option.** Grid-searched weights are the
   naive baseline; data-driven integration (latent factors, sparse multi-block)
   should be selectable, not bolted on.
4. **Interpretability is required, not optional.** However sophisticated the
   backend, the package's output stays a *ranked, attributable gene table*.
5. **Rigor as a supporting layer.** External validation, bootstrap C-index
   intervals, and calibration are how we keep ourselves honest — a quality gate,
   not the headline feature.

## Phases

### Phase 0 — Gene-set scoring, visualization, CLI  *(shipped, v0.2.0)*
`score_gene_set()`, the `ccscore` CLI, ggplot2/survminer visualizations,
bundled example data, and a working survival path.

### Phase 1 — Generalized N-omics integration core  *(this PR, v0.3.0)*
- **`omics_block()`** — describe one omics layer + how to score it (a built-in
  metric name or any user function `data -> named gene ranking`).
- **`score_multiomics()`** — score a gene set across an arbitrary number of
  blocks.
- **`score_rankings()`** — the generic engine: restrict rankings to a gene set,
  renormalize, blend, and emit the tidy per-metric + combined + rank table.
- **N-metric integration** — `combine_rankings()` and `optimize_weights()` are
  generalized from the old 2–3 metric cap to any number of layers (simplex
  weight grid).
- `score_gene_set()` becomes a thin, backward-compatible wrapper over the new
  core.

### Phase 2 — Regulatory-activity metrics  *(next PR, v0.4.0)*
- **`calculate_mirna_activity()`** — replace the static target-count metric with
  an expression-anchored score: for each gene, aggregate the *observed*
  anti-correlation between its expression and that of its targeting cancer
  miRNAs across matched samples (a ceRNA/repression-evidence view).
- Generic adapters so methylation (gene–methylation anti-correlation) and CNV
  (copy-number–expression correlation) enter as blocks with the same interface.

### Phase 3 — Learned integration backends  *(v0.5.0)*
- `score_multiomics(..., method = c("weighted", "mofa", "diablo"))` wrapping
  `MOFA2` (unsupervised latent factors) and `mixOmics::block.splsda` (DIABLO,
  supervised sparse multi-block selection) so integration weights are learned
  from the data rather than grid-searched. Expose the learned loadings as the
  per-gene attribution.

### Phase 4 — Rigorous survival evaluation  *(v0.6.0)*
- Bootstrap C-index confidence intervals, calibration curves, and a
  train-cohort → external-cohort validation harness built into the survival
  path; always report against a clinical-only baseline.

### Phase 5 — Deep integration + interpretability  *(v0.7.0)*
- Optional autoencoder embedding of the joined omics → Cox head (via
  `reticulate`), with SHAP-style per-gene attributions so the model stays
  explainable. Gated behind `Suggests`.

### Phase 6 — Single-cell multi-omics  *(v0.8.0)*
- Extend the trajectory logic (switchDE) to paired single-cell modalities
  (RNA + ATAC / protein via 10x Multiome / CITE-seq), integrating with
  established single-cell methods (WNN, MOFA+, GLUE) for multi-omic pseudotime.

## References (state of the field, from PubMed)

- Deep-learning multi-omics integration with interpretability (Shapley values):
  CustOmics, https://doi.org/10.1371/journal.pcbi.1010921
- Autoencoder multi-omics → survival subtypes across 12 cancers: ProgCAE,
  https://doi.org/10.1093/bib/bbad196
- mRNA+miRNA+methylation autoencoder, C-index 0.748: melanoma risk subtypes,
  https://doi.org/10.1007/s00432-023-05358-x
- Foundational autoencoder multi-omics survival model (Garmire lab), PMID
  29888072 (PMC5961799)
- ceRNA (lncRNA–miRNA–mRNA) prognostic modeling in colorectal cancer:
  https://doi.org/10.3389/fimmu.2023.1279789
- Multi-omics regulatory networks + centrality for biomarker discovery: PLBINs,
  https://doi.org/10.3390/cells12010101
- Bayesian multi-omics integration for prognosis (review):
  https://doi.org/10.21873/cgp.20298

*Article facts above are from PubMed; method toolkits named without a DOI
(MOFA+, DIABLO/mixOmics, SNF, Seurat WNN, GLUE) are established software
referenced from general knowledge.*
