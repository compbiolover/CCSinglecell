# Open-data plan

Addendum to [Multiomics biomarker discovery frontier](Multiomics%20biomarker%20discovery%20frontier.md). Constraint: **every input must be downloadable without an access application**, so anyone can rerun the project from a clean clone. Access status below was checked directly against the source APIs on 2026-11-07 unless marked *verify*.

## Verdict

Yes. Open data loses the randomised-trial and ctDNA questions, but it keeps everything needed for a strong **prognostic** study: four open outcome cohorts with RNA, two of them with matched miRNA (and one with protein), a 371k-cell single-cell atlas, public spatial data, and TCGA H&E slides with Apache-2.0 pathology encoders. Because nothing is gated, the whole pipeline reproduces from one command with no data-use agreements.

What changes from the main plan:

| Main plan | Open-data version | Why |
|---|---|---|
| PETACC-3, S:CORT, TransSCOT for development and confirmation | GSE39582 for development; TCGA, AC-ICAM, curatedCRCData cohorts for leave-one-cohort-out validation | PETACC-3's public deposit (E-MTAB-990) has CEL files, stage and age but **no RFS/OS columns**; the others are consortium-only |
| Treatment-interaction analyses (CALGB 80405, FOCUS) | Dropped; adjuvant-chemo effect modification in GSE39582 is exploratory only (non-randomised) | GSE196576 (CALGB 80405, n=579) is open but carries arm, age, sex and race, **no outcomes** |
| ctDNA-conditional analyses | Dropped | No open individual-level ctDNA with outcomes |
| Baseline: stage + MMR + sidedness (+ ctDNA, Immunoscore) | Stage + MMR/MSI + sidedness + age (+ CMS where called) | All present in the open cohorts |
| H&E with UNI2-h / Virchow2 | H-optimus-0, Virchow (v1), Prov-GigaPath or Hibou-L | These four are **Apache-2.0** (click-through gate); UNI2-h, CONCH, H-optimus-1 are CC-BY-NC-ND with derivative restrictions |

## Open data stack

### Bulk cohorts with outcomes

| Cohort | n | Layers (open) | Outcomes | Source | Role |
|---|---|---|---|---|---|
| GSE39582 (CIT) | 585 samples | Affymetrix mRNA; stage, TNM, MMR, sidedness, KRAS/BRAF/TP53, CIMP/CIN, adjuvant chemo | RFS, OS | GEO (verified fields) | Development |
| TCGA COAD/READ | ~629 | RNA-seq, **miRNA-seq**, methylation, CNV, mutations, RPPA, **H&E slides** | PFI/DSS (TCGA-CDR), DFS, OS | GDC open tier; cBioPortal `coadread_tcga_pan_can_atlas_2018` (verified profiles) | Validation fold; miRNA-activity check; H&E arm |
| AC-ICAM (Sidra-LUMC) | ~348 | RNA-seq normalised counts, CNA, mutations, CMS, MSI score | OS, DSS, PFS | cBioPortal `coad_silu_2022` (verified) | External RNA-seq validation |
| CPTAC-COAD | ~110 | RNA, **miRNA**, proteome, phosphoproteome, CNA | DFS, OS | cBioPortal `coad_cptac_2019` (verified); slides on TCIA (*verify*) | Multi-omics mechanism; small external H&E test |
| curatedCRCData 2.44.0 | many GEO series | Harmonised microarray + curated clinical | Mixed RFS/DFS/OS | Bioconductor (Artistic-2.0, verified in release) | Extra leave-one-cohort-out folds (*verify which series carry RFS*) |
| GSE161158, GSE103479, GSE143985 | 250, 156, 91 | Microarray, stage II/III | DFS per series titles (*verify fields*) | GEO | Extra stage II/III folds |

Pooled, this is roughly 1,500+ patients across 5+ independent cohorts. That is enough for internal-external validation and comfortably clears the Riley budget that sank the TCGA-only model.

### Single-cell, spatial, imaging

| Resource | What | Access |
|---|---|---|
| Pelka 2021, GSE178341 | 371k cells, MMRd vs MMRp, counts h5 + cluster and metadata tables | GEO (verified) |
| Lee 2020, GSE132465 / GSE144735 | SMC and KUL3 cohorts | GEO |
| CELLxGENE Census LTS 2025-11-08 | Harmonised CRC collections | CC BY 4.0 |
| Marteau CRC atlas (4.27M cells) | scANVI-labelled h5ad | GitHub repo has **no licence**; use the CELLxGENE-hosted copy if present (*verify*) |
| Joanito 2022 (iCMS) | 63 patients | Synapse, free registration (*verify*) |
| Visium HD CRC | 10x public datasets; GSE280315 (5 patients) | 10x site / GEO (*verify licence*) |
| TCGA diagnostic H&E | Whole-slide images, COAD/READ | GDC open tier |

## Design decisions

- **Reuse CcSinglecell instead of porting it.** The new repo installs CcSinglecell as an R dependency. `cv_validate_survival()` already fits elastic-net Cox with the clinical variables unpenalised, and `assess_null()`, `power_curve()` and `learning_curve()` exist. The additions belong in CcSinglecell itself: leave-one-cohort-out folds (folds are currently random), time-dependent AUC and IBS, calibration slope, decision curves and a pmsampsize gate.
- **Nextflow DSL2 with nf-test fixtures** for orchestration. Snakemake 9, which the main report suggested, is an equally valid choice; the analysis is file-in, file-out either way.

## Build order

1. **Data spine.** Scaffold the repo (Nextflow, pixi, CI on a small fixture). Pull GSE39582, TCGA (GDC open tier or a pinned cBioPortal datahub commit) and AC-ICAM. Harmonise 3-year RFS/PFI, report clinical-baseline C in each cohort, and commit the pre-registration. Add leave-one-cohort-out folds to the CcSinglecell harness.
2. **Single-cell states.** Run QC on Pelka GSE178341, then scVI/scANVI. Define epithelial states (iCMS2/3, EMP1-high), CAF states and myeloid states. Benchmark scGPT or Geneformer embeddings against PCA and scVI on label transfer.
3. **Into bulk.** Run pseudobulk DE with patients as replicates to get state signatures, then deconvolve every bulk cohort with BayesPrism or InstaPrism. Write the TargetScan 8 loader that replaces hoardeR, run decoupler activity for miRNAs, CollecTRI TFs and PROGENy pathways, and check inferred miRNA activity against measured miRNA in TCGA and CPTAC.
4. **Spatial.** Deconvolve Visium HD CRC spots with the atlas reference (cell2location or RCTD) and map niches with squidpy. Test whether the relapse-linked states sit in specific niches.
5. **Models and locked test.** Compare models against the clinical baseline using leave-one-cohort-out validation, random and proliferation-matched nulls and stability selection. Run one locked test on AC-ICAM and produce the ranked list of states and regulators.
6. **Write-up.** Preprint, a Quarto results site and an interactive explorer of ΔC by cohort and state.

Stretch goals: the H&E arm (attention-MIL survival head on TCGA slides with an Apache-2.0 encoder), the extra curatedCRCData folds, knockoffs and conformal survival bounds.

## Practical constraints

- **Compute.** Everything in the core build runs on a laptop or a modest cloud VM. The stretch H&E arm needs a GPU for slide embedding (single-digit to tens of GPU-hours for ~600 slides with a ViT-L/ViT-g encoder; a judgement, not measured).
- **Storage.** TCGA COAD/READ slides run to hundreds of GB. Stream slide by slide, keep only tile embeddings, delete the SVS.
- **Versioning.** cBioPortal's web API is not versioned; pin downloads to a `cBioPortal/datahub` commit (or checksum the files) and record GDC release numbers.
- **Licences.** Encoders: stay on Apache-2.0 models. Data: CC BY 4.0 (CELLxGENE), GEO/GDC open tier; confirm the Marteau atlas and Visium HD terms before redistributing anything derived.

## Honest limits

- The study can make **prognostic** claims only. With no open randomised cohort carrying outcomes, it cannot claim a biomarker predicts benefit from a specific therapy.
- Microarray (GSE39582, curatedCRCData) and RNA-seq (TCGA, AC-ICAM) cohorts differ by platform, so features must be rank-based or single-sample scores (state fractions, activities, CMS), not raw gene weights.
- No source shows deconvolved cell states add ΔC over stage + MMR in CRC. The most likely molecular result is still a bounded null; the H&E arm is the better bet for a positive ΔC. Either outcome is presentable if pre-registered.
