# Single-cell and spatial omics layer for bulk-applicable CRC prognostic biomarkers (state of the art 2024–2026)

Scope: what should replace each step of the CcSinglecell single-cell layer (GSE81861 FPKM, MAGIC, Monocle3 rooted at low-VIM cells, MAD/switchDE gene ranking, then a TCGA-COAD/READ Cox model whose ΔC over age+stage is about 0), and how single-cell and spatial information is best carried over to bulk survival cohorts. Notes compiled 2026-11-07.

Evidence tags used below:
- **[IND-BENCH]**: independent benchmark, run by people who did not build the methods.
- **[DEV-BENCH]**: the developers comparing their own method with others.
- **[SINGLE]**: a claim from one study.
- **[PREPRINT]**: not peer-reviewed.
- **[ABSTRACT]**: conference abstract only.
- **[PK]**: taken from prior knowledge, with a DOI link I did not re-fetch in this session. These are mostly well-known foundational papers. Check the DOI before citing.
- **[SNIPPET]**: read only from a search-result summary, not from the full text.

## 1. Which large CRC single-cell and spatial atlases now exist (sizes, access, licensing)?

### Takeaway
GSE81861 (about 375 tumour cells, Fluidigm C1, FPKM) is now roughly four orders of magnitude smaller than the largest integrated CRC reference: the Marteau et al. atlas (Cancer Cell 2025), with 4.27 M cells from 650 patients, labelled with scANVI, available as h5ad files and in a cellxgene viewer. The realistic reference stack today is:
- that atlas;
- Pelka 2021 (371k cells, with MMRd/MMRp status);
- Joanito 2022 (the iCMS2/iCMS3 epithelial states, 63 patients);
- Lee 2020 (the SMC and KUL3 cohorts);
- HTAN polyp→CRC data (Chen 2021, Becker 2022);
- small but high-resolution spatial sets: Visium HD n=5, CosMx about 846k cells from n=23 samples, Xenium n=16 patients.

Almost none of these atlases carries survival follow-up, so the outcome signal has to come from bulk cohorts. Processed counts are mostly open. Raw HTAN and Broad reads are behind dbGaP. The licences of most processed objects are not stated clearly.

### Cited Findings
- **Marteau et al. integrated CRC atlas.**
  - Size: 4.27 M cells, 650 patients, 1,670 samples. The preprint says 76 datasets from 48 studies; the GitHub repo says 77 datasets from 49 studies, so the counts conflict slightly. Tissues are normal colon, polyps, primary tumours and metastases. [SNIPPET + fetched repo] — [bioRxiv 10.1101/2024.08.26.609563](https://www.biorxiv.org/content/10.1101/2024.08.26.609563v1); [GitHub icbi-lab/crc-atlas](https://github.com/icbi-lab/crc-atlas)
  - Annotation: 12 coarse cell types and 62 subtypes/states. Cancer cells total 600,740 (primary plus metastatic). Tumours are classified into four immune phenotypes: immune desert, myeloid, B cell and T cell. [SNIPPET] — [bioRxiv 2024.08.26.609563](https://www.biorxiv.org/content/10.1101/2024.08.26.609563v1)
  - Published version (Cancer Cell 2025, a neutrophil-focused paper): adds four cohorts of 266 patients, with spatial single-cell transcriptomics on 3.7 M cells and spatial single-cell protein profiling on 715k cells. [SNIPPET] — [Cancer Cell S1535-6108(25)00537-9](https://www.cell.com/cancer-cell/fulltext/S1535-6108(25)00537-9)
  - Access and licence: cellxgene viewer (via crc.icbi.at) and h5ad objects. Cell-type labels depend on the scANVI embedding. Pre-computed `build_atlas` results are on Zenodo; the input data were "coming soon". The repo has a LICENSE file, but the README does not name the licence type. No survival or bulk-projection analysis is mentioned in the repo. — [GitHub icbi-lab/crc-atlas](https://github.com/icbi-lab/crc-atlas)
- **Pelka et al. 2021 (Human Colon Cancer Atlas, Broad SCP1162).**
  - 371,223 cells from colorectal tumours and adjacent normal tissue. Raw counts are on GEO (GSE178341). Raw sequence data need dbGaP approval. — [Broad Single Cell Portal SCP1162](https://singlecell.broadinstitute.org/single_cell/study/SCP1162)
  - [PK] 62 patients split into MMRd and MMRp. The paper defines spatially organised "multicellular immune hubs", including a CXCL13+ T-cell hub in MMRd tumours. — [Cell 2021, doi:10.1016/j.cell.2021.08.003](https://doi.org/10.1016/j.cell.2021.08.003)
- **Joanito et al. 2022 (Nature Genetics).**
  - Single-cell data from 63 patients, used to define two intrinsic epithelial states, iCMS2 and iCMS3. These refine the CMS classes into the IMF scheme (iCMS, MSI, fibrosis).
  - The authors applied 715 iCMS marker genes to bulk transcriptomes from 15 datasets and ran survival analyses. — [Nat Genet 2022](https://www.nature.com/articles/s41588-022-01100-4); [PMC9279158](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9279158/)
- **Khaliq et al. 2022 (Genome Biology).** 487,829 cells profiled, 49,589 kept after QC. The authors built a reference and deconvolved bulk profiles to pseudo-single-cell resolution, with survival analysis. — [Genome Biol 2022](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02677-z)
- **Lee et al. 2020 (Nature Genetics) [PK].** About 63k cells from 23 Korean patients (SMC, GSE132465) and 6 Belgian patients (KUL3, GSE144735). Describes the SMC/KUL tumour-microenvironment states. — [doi:10.1038/s41588-020-0636-z](https://doi.org/10.1038/s41588-020-0636-z)
- **HTAN Vanderbilt (Chen et al., Cell 2021).** Conventional adenomas and serrated polyps, plus the CRCs that arise from them. The HCA Data Explorer lists CZ CELLxGENE and UCSC Cell Browser as viewers. — [HCA explorer](https://explore.data.humancellatlas.org/projects/50154d1e-2308-44bf-9608-10c7afaa560b)
- **HTAN Stanford (Becker et al., Nature Genetics 2022) [PK].** scRNA-seq plus scATAC-seq along the polyp→CRC continuum, including FAP patients. This is the main CRC chromatin resource. — [doi:10.1038/s41588-022-01088-x](https://doi.org/10.1038/s41588-022-01088-x)
- **HTAN access model.**
  - Single-cell data are harmonised for deposit into CELLxGENE Discover.
  - Controlled Level 1/2 data (raw) need a dbGaP request (phs002371) and are reached through the NCI Cancer Research Data Commons.
  - HTA "CRC Atlas X" (Orion imaging) data were released without restriction under the NCI Moonshot policy. — [HTAN data-sharing preprint, bioRxiv 2024.06.25.598921](https://www.biorxiv.org/content/10.1101/2024.06.25.598921.full.pdf); [HTA CRC Atlas X](https://labsyspharm.github.io/orion-crc)
- **CELLxGENE Census [PK].**
  - A versioned API over all of CELLxGENE Discover, queryable with TileDB-SOMA. It is the practical way to pull every CRC-labelled dataset with harmonised metadata. — [NAR 2025, doi:10.1093/nar/gkae1142](https://doi.org/10.1093/nar/gkae1142)
  - [PK] Discover datasets are generally released under CC BY 4.0. Verify this per dataset.
- **Visium HD (Oliveira, Romero et al., Nature Genetics 2025; 10x Genomics authors).**
  - FFPE CRC samples, n=5, at single-cell-scale resolution.
  - Finds macrophage subpopulations that are spatially distinct and plausibly pro- or anti-tumour. A clonally expanded T-cell population sits near the anti-tumour macrophages.
  - Data: GEO GSE280315. 10x also hosts public CRC Visium HD datasets. — [PubMed 40473992](https://pubmed.ncbi.nlm.nih.gov/40473992/); [10x dataset page](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-human-crc)
- **CosMx atlas (Pattison et al., bioRxiv Dec 2025) [PREPRINT].**
  - About 1,000 genes measured in 846,469 cells from 23 tumour and normal samples. Nine spatial niches.
  - Stromal niches are linked to tumour differentiation, not to outcome. Data: GEO GSE303070. — [bioRxiv 2025.12.09.693138](https://www.biorxiv.org/content/10.64898/2025.12.09.693138v1.full)
- **Xenium (STORM, Yao et al., UT Southwestern, bioRxiv Jul 2026) [PREPRINT].** 32 sections from 16 patients (9 early-onset and 7 average-onset CRC, each with matched normal), on a 375-gene panel. — [PMC13405321](https://pmc.ncbi.nlm.nih.gov/articles/PMC13405321/)
- **Imaging-platform benchmark (Nature Communications 2025).** Compares Xenium, CosMx and other platforms on tissue microarrays for data quality and cell retention, not for clinical outcome. [SNIPPET] — [Nat Commun 2025, s41467-025-64990-y](https://www.nature.com/articles/s41467-025-64990-y.pdf)

### Inferences
- **GSE81861 should be dropped as the primary reference.** Its 11-patient size and C1/FPKM chemistry mean that per-gene, per-cell scores mostly reflect patient and batch structure. Raw counts from Marteau/Pelka/Joanito/Lee allow patient-aware statistics.
- **Suggested role for each dataset:**
  - Marteau atlas: label and reference backbone (scANVI labels, all compartments, metastasis).
  - Pelka: MMR-stratified validation of states.
  - Joanito: the epithelial iCMS axis.
  - Becker/HTAN: the chromatin layer.
  - Visium HD, CosMx, Xenium: niche and co-localisation checks.
- **None of the large atlases carries survival.** The single-cell layer can only define candidate features (states, programs, niches). Their prognostic value must be tested in bulk cohorts with follow-up, and features must never be selected on the same bulk outcome data used to evaluate them.

### Gaps
- Licence type of the Marteau h5ad and Zenodo objects is not stated in the repo. Joanito's data portal and cell counts were not retrieved in this session.
- Which specific CELLxGENE Census collections cover Pelka, Chen or Joanito was not checked. The SCP1162 URL slug ("human-colon-cancer-atlas-c295") indicates that "Human Colon Cancer Atlas" is the Pelka dataset.
- No CRC spatial dataset with more than about 100 patients **and** recurrence follow-up was found, apart from the neutrophil-paper spatial cohorts (266 patients). Their outcome annotations were not verified.
- Availability of raw FASTQ for GSE81861 (needed for any RNA-velocity reanalysis) was not checked.

## 2. What replaced MAGIC-style denoising and Monocle3 + switchDE for trajectory/dynamic-gene detection? What do benchmarks say?

### Takeaway
- **Imputation:** there is no direct successor to MAGIC for gene-level statistics. Imputation inflates gene–gene structure. The field now works on raw counts with count-based latent models (scVI/scANVI), metacells, or patient-level pseudobulk.
- **Fate mapping:** CellRank 2 (several kernels: pseudotime, CytoTRACE, real-time, velocity) and optimal transport (moscot) replace single-rooted pseudotime.
- **Dynamic genes:** GAM-based tradeSeq/condiments, and the multi-sample Lamian framework, replace switchDE. Benchmarks warn that cell-level DE without patient replication gives false discoveries.
- **RNA velocity:** a 15-method, 17-dataset benchmark finds weak agreement between methods, and accuracy falls as dynamics get more complex.
- **The EMT axis:** most important for CcSinglecell, the axis the package rooted on (low VIM inside tumour cells) is biologically doubtful. CRC "EMT/mesenchymal" bulk signal is largely stromal (CAF-derived), a point the GSE81861 authors themselves made.

### Cited Findings
- **MAGIC [PK, foundational 2018].** Graph diffusion that imputes and smooths every gene. — [van Dijk et al., Cell 2018, doi:10.1016/j.cell.2018.05.061](https://doi.org/10.1016/j.cell.2018.05.061)
- **Imputation creates artefacts [PK, foundational].**
  - Andrews & Hemberg showed that imputation methods (including MAGIC) introduce false-positive gene–gene correlations and differential-expression signals. — [F1000Research 2018, doi:10.12688/f1000research.16613.2](https://doi.org/10.12688/f1000research.16613.2)
  - A systematic evaluation of 18 imputation methods found limited and inconsistent benefit for downstream analyses. — [Hou et al., Genome Biol 2020, doi:10.1186/s13059-020-02132-x](https://doi.org/10.1186/s13059-020-02132-x) [IND-BENCH, PK]
- **Trajectory inference baseline [PK, foundational 2019].** Benchmark of 45 trajectory methods. Performance depended on topology, and Slingshot, TSCAN, PAGA and Monocle DDRTree were among the top for various topologies. No newer equally comprehensive trajectory-inference benchmark was found in this session. — [Saelens et al., Nat Biotechnol 2019, doi:10.1038/s41587-019-0071-9](https://doi.org/10.1038/s41587-019-0071-9) [IND-BENCH]
- **switchDE [PK, 2017].** Fits a sigmoid of expression over pseudotime for each gene, at cell level, with no notion of samples. — [Campbell & Yau, Bioinformatics 2017, doi:10.1093/bioinformatics/btw798](https://doi.org/10.1093/bioinformatics/btw798)
- **tradeSeq, condiments, Lamian (dynamic-gene testing).**
  - tradeSeq [PK]: negative-binomial GAMs along one or more lineages, with tests for association, start-vs-end and pattern differences. — [Van den Berge et al., Nat Commun 2020, doi:10.1038/s41467-020-14766-3](https://doi.org/10.1038/s41467-020-14766-3)
  - condiments [PK]: tests whether conditions change topology, progression or differentiation along trajectories. — [Roux de Bézieux et al., Nat Commun 2024, doi:10.1038/s41467-024-44823-0](https://doi.org/10.1038/s41467-024-44823-0)
  - Lamian [PK]: differential pseudotime analysis that treats samples, not cells, as replicates. This matters for patient-level biomarker claims. — [Hou et al., Nat Commun 2023, doi:10.1038/s41467-023-42841-y](https://doi.org/10.1038/s41467-023-42841-y)
- **Pseudoreplication [PK, IND-BENCH].** Single-cell DE methods that treat cells as replicates produce many false discoveries. Pseudobulk methods agree best with ground truth. — [Squair et al., Nat Commun 2021, doi:10.1038/s41467-021-25960-2](https://doi.org/10.1038/s41467-021-25960-2)
- **CellRank 2 (Nature Methods 2024).**
  - A unified fate-mapping framework for multiview data at the scale of millions of cells.
  - Kernels go beyond RNA velocity: PseudotimeKernel, CytoTRACEKernel (developmental potential from gene counts) and real-time (experimental time-point) kernels.
  - A Nature Protocols paper (2025) describes the workflow. — [CellRank 2, Nat Methods 2024](https://www.nature.com/articles/s41592-024-02303-9); [CellRank protocol, Nat Protoc 2025](https://experiments.springernature.com/articles/10.1038/s41596-025-01314-w)
- **RNA velocity benchmark (bioRxiv Aug 2025; Cell Reports Methods 2026) [IND-BENCH].**
  - 15 velocity methods tested on 17 datasets.
  - Agreement exists only between certain method pairs: velocyto with scVelo-stochastic, CellRank with cellDancer, DeepVelo with veloVI.
  - Inference accuracy declines as transcriptional dynamics get more complex.
  - LatentVelo and cell2fate gave the most coherent velocity fields. — [bioRxiv 2025.08.02.668272](https://www.biorxiv.org/content/10.1101/2025.08.02.668272v1.full); [Cell Rep Methods 2026](https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(26)00067-6)
- **cell2fate (Nature Methods 2025).** A velocity model built on modules. It claims better fate resolution than existing velocity methods. [DEV-BENCH] — [Nat Methods 2025](https://www.nature.com/articles/s41592-025-02608-3)
- **Velocity caveats [PK].** Velocity depends on spliced/unspliced counts and on model assumptions (constant kinetic rates, among others) that often fail. — [Bergen et al., Mol Syst Biol 2021, doi:10.15252/msb.202110282](https://doi.org/10.15252/msb.202110282)
- **moscot [PK].** Optimal-transport framework for mapping cells across time points and space, also coupling dissociated cells to spatial data. — [Klein et al., Nature 2025, doi:10.1038/s41586-024-08453-2](https://doi.org/10.1038/s41586-024-08453-2)
- **Palantir [PK, foundational 2019].** Markov-chain fate probabilities and an entropy of differentiation potential. Usable through CellRank's PseudotimeKernel. — [Setty et al., Nat Biotechnol 2019, doi:10.1038/s41587-019-0068-4](https://doi.org/10.1038/s41587-019-0068-4)
- **Latent models and integration [PK].**
  - scVI: a count-based variational autoencoder. — [Lopez et al., Nat Methods 2018, doi:10.1038/s41592-018-0229-2](https://doi.org/10.1038/s41592-018-0229-2)
  - scANVI: the semi-supervised extension of scVI. — [Xu et al., Mol Syst Biol 2021, doi:10.15252/msb.20209620](https://doi.org/10.15252/msb.20209620)
  - scIB benchmark: scANVI, Scanorama, scVI and scGen ranked well on complex integration tasks. — [Luecken et al., Nat Methods 2022, doi:10.1038/s41592-021-01336-8](https://doi.org/10.1038/s41592-021-01336-8) [IND-BENCH]
  - The Marteau CRC atlas itself uses scANVI labels. — [GitHub icbi-lab/crc-atlas](https://github.com/icbi-lab/crc-atlas)
- **Metacells [PK].** SEACells aggregates cells into metacells to cut sparsity without imputing gene–gene structure. — [Persad et al., Nat Biotechnol 2023, doi:10.1038/s41587-023-01716-9](https://doi.org/10.1038/s41587-023-01716-9)
- **Programs instead of trajectories [PK].** Pan-cancer NMF "meta-programs" of intratumour heterogeneity (cell cycle, stress, EMT-like, interferon, and others) built from more than 1,000 tumours. This is a reusable alternative to single-axis pseudotime. — [Gavish et al., Nature 2023, doi:10.1038/s41586-023-06130-4](https://doi.org/10.1038/s41586-023-06130-4)
- **EMT in CRC is stromal [PK, central to the CcSinglecell critique].**
  - The GSE81861 study (Li et al. 2017) reported that EMT-related genes were upregulated in cancer-associated fibroblasts, not in epithelial tumour cells. This explains why bulk tumours were assigned to the mesenchymal CMS4. — [Li et al., Nat Genet 2017, doi:10.1038/ng.3818](https://doi.org/10.1038/ng.3818)
  - Two independent studies showed that the poor-prognosis CRC mesenchymal/stem-like subtype is driven by stromal (CAF) gene expression. — [Isella et al., Nat Genet 2015, doi:10.1038/ng.3225](https://doi.org/10.1038/ng.3225); [Calon et al., Nat Genet 2015, doi:10.1038/ng.3224](https://doi.org/10.1038/ng.3224)
- **A CRC epithelial axis with outcome relevance [PK].** Residual EMP1-high "high-relapse cells" drive metastatic recurrence, and the HRC gene program predicts relapse in human CRC cohorts. — [Cañellas-Socias et al., Nature 2022, doi:10.1038/s41586-022-05402-9](https://doi.org/10.1038/s41586-022-05402-9)

### Inferences
| Current CcSinglecell step | Recommended replacement | Rationale (sources above) |
|---|---|---|
| GSE81861 FPKM, ~375 cells | Raw counts from the Marteau atlas plus Pelka, Joanito and Lee; patient as the unit | Patient-level replication; chemistry that matches modern data |
| MAGIC denoising | No imputation for statistics. Use scVI/scANVI latents for neighbourhoods and labels, and pseudobulk or metacells for expression. MAGIC for visualisation only | Imputation causes false signals; pseudobulk controls false discoveries |
| Monocle3 pseudotime rooted at low-VIM tumour cells | Either (a) drop the trajectory and use discrete states/programs (iCMS2/3, HRC/EMP1, stem/differentiation, Gavish meta-programs), or (b) CellRank 2 with CytoTRACE/Pseudotime kernels on epithelial cells only, with the root chosen on biology (LGR5/stem) rather than VIM | CRC EMT signal is CAF-derived; single-root pseudotime is fragile |
| switchDE + MAD per gene | tradeSeq/condiments at most for description. For biomarker ranking, use Lamian or pseudobulk DE (patients as replicates), or state-abundance features | Cell-level tests ignore patient variance |
| RNA velocity (not used today) | Not recommended as a primary step | Poor cross-method agreement in benchmarks; needs spliced counts |
- **The decisive design change is the unit of analysis.** Move from genes ranked along one per-cell axis to patient-replicated cell states or programs that can be quantified in bulk.

### Gaps
- No independent benchmark was found that compares switchDE, tradeSeq, condiments and Lamian on downstream prognostic usefulness. The benchmarks evaluate topology or DE calibration, not bulk survival transfer.
- No 2024–2026 trajectory-inference benchmark comparable to Saelens 2019 was retrieved. Lineage-tracing approaches (CRISPR barcodes, mtDNA variants) were not researched in this session.
- The CytoTRACE 2 publication venue and its benchmark results were not verified.

## 3. Methods that link single-cell states to bulk patient outcomes: which are validated, which perform best, which are 2024–2026 newcomers?

### Takeaway
**Two families:**
- **(A) Cell-prioritisation by bulk phenotype:** Scissor, scAB, DEGAS, scPAS, SIDISH, PIPET, and PACells for scATAC.
  - Almost all comparisons are run by developers, and each new tool reports beating Scissor.
  - The one independent benchmark found (SigBridgeR, preprint 2026, 8 methods, 4 cancers) rates Scissor and scAB as the most robust and suggests survival endpoints are hard.
  - None of these tools has been shown to add C-index over clinical covariates.
- **(B) Deconvolution and cell-state abundance:** CIBERSORTx, BayesPrism/InstaPrism, EcoTyper, Kassandra, MuSiC2, ReCIDE.
  - A 2026 benchmark on 18 real cohorts (5,891 samples) favours ReCIDE and BayesPrism and warns that pseudobulk-only benchmarks mislead. InstaPrism is fast but weaker on real bulk.

**In CRC,** the single-cell-derived features with the most bulk validation are compositional or state programs: iCMS/IMF, EMP1+ high-relapse cells, FAP+ CAF/SPP1+ macrophage, immune hubs, mCAF. Gene lists ranked along a trajectory are not among them. The one CRC study found that explicitly compared deconvolved immune states against clinical variables reports only modest AUCs, about 0.69.

### Cited Findings
**Family A: cell prioritisation**
- **Scissor [PK, 2022].** Network-regularised regression of bulk phenotype (including Cox) on correlations between cells and bulk samples. It found lung-cancer cell subsets associated with worse survival and with TP53 mutation. — [Sun et al., Nat Biotechnol 2022, doi:10.1038/s41587-021-01091-3](https://doi.org/10.1038/s41587-021-01091-3)
- **scAB (NAR 2022).** Knowledge- and graph-guided matrix factorisation, giving multiresolution phenotype-associated states. It reported the highest C-index against other methods on a liver-cancer test set. [DEV-BENCH] — [NAR 2022](https://academic.oup.com/nar/article/50/21/12112/6849494)
- **DEGAS [PK].** Deep transfer learning that prioritises cells using bulk labels. — [Johnson et al., Brief Bioinform 2022, doi:10.1093/bib/bbab439](https://doi.org/10.1093/bib/bbab439)
- **scPAS (Briefings in Bioinformatics 2025, newcomer).** Developers report that scAB and DEGAS cannot separate phenotype-associated cells from two states of the same cell type, and that Scissor was "noticeably worse" for continuous phenotypes. [DEV-BENCH] — [Brief Bioinform 26(1):bbae655](https://academic.oup.com/bib/article/26/1/bbae655/7925487)
- **SIDISH (Nature Communications 2025, newcomer).** Integrates single-cell and bulk data to find high-risk cells, then runs in-silico perturbation. It claims consistent gains over Scissor, DEGAS and scAB on TCGA-BRCA and TCGA-PDAC. [DEV-BENCH] — [Nat Commun 2025](https://www.nature.com/articles/s41467-025-66162-4)
- **PACells (AJHG 2026, newcomer).** Phenotype-associated cell states from scATAC profiles. — [Am J Hum Genet 2026](https://www.sciencedirect.com/science/article/abs/pii/S0002929726001229)
- **SigBridgeR (Yang et al., bioRxiv 2026-05-12) [PREPRINT, closest to IND-BENCH].** An R framework that unifies 8 phenotype-associated cell screening algorithms, benchmarked across 4 cancer types with binary and survival phenotypes. Scissor and scAB "emerged as the most consistently robust". The snippet suggests several methods struggled with survival phenotypes; the full text could not be retrieved (HTTP 429). — [bioRxiv 10.64898/2026.05.08.723458](https://www.biorxiv.org/content/10.64898/2026.05.08.723458v1.full.pdf); metadata via [Europe PMC](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:10.64898/2026.05.08.723458&resultType=core&format=json)

**Family B: deconvolution and cell-state abundance**
- **CIBERSORTx [PK, foundational].** Signature-matrix deconvolution with batch correction and group-mode or high-resolution purification of cell-type expression. — [Newman et al., Nat Biotechnol 2019, doi:10.1038/s41587-019-0114-2](https://doi.org/10.1038/s41587-019-0114-2)
- **BayesPrism (Nature Cancer 2022).**
  - A Bayesian model that infers cell-type fractions **and** cell-type-specific expression, including the malignant compartment, from an scRNA reference.
  - Its developers claim it outperforms CIBERSORTx on tumour data. [DEV-BENCH] — [Chu et al., Nat Cancer 2022](https://www.nature.com/articles/s43018-022-00356-3)
- **Breast-cancer benchmark on simulated mixtures (Nature Communications 2023).** A later review describes BayesPrism as best overall for 9 major cell types. [IND-BENCH, SNIPPET] — [Nat Commun 2023, s41467-023-41385-5](https://www.nature.com/articles/s41467-023-41385-5)
- **Missing reference cell types (bioRxiv 2024, now in PMC).** NNLS, BayesPrism and CIBERSORTx performed comparably at recovering proportions when a cell type was missing from the reference. BayesPrism lost accuracy when noise was added. The missing type can be detected. [IND-BENCH] — [PMC11974051](https://pmc.ncbi.nlm.nih.gov/articles/PMC11974051/)
- **Reference–mixture matching (Briefings in Bioinformatics 2026).** Seven methods tested. The choice of algorithm mattered much more than matching the reference to the mixture. InstaPrism was comparable on simulated bulk but worse on real bulk, and fastest (about 2–2.5 min). [IND-BENCH, SNIPPET] — [Brief Bioinform 27(2):bbag184](https://academic.oup.com/bib/article/27/2/bbag184/8659159)
- **Real-bulk prognostic benchmark (Li, Su, Tang, Lee, Tian; Genome Biology, published 2026-01-21).**
  - 18 real bulk cohorts (5,891 samples) across 9 cancer types, 5 methods.
  - ReCIDE and BayesPrism were the most robust. Matrix CAFs (mCAF) were prognostic across cancers.
  - Traditional pseudobulk benchmarks "may not always be reliable in real-world settings". [IND-BENCH] — [Genome Biol 2026, doi:10.1186/s13059-026-03942-1](https://link.springer.com/article/10.1186/s13059-026-03942-1); abstract via [Europe PMC](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:10.1186/s13059-026-03942-1&resultType=core&format=json)
- **Single-cell plus single-nucleus references (Cell Reports Methods 2026).** Combining scRNA and snRNA references improves bulk deconvolution, studied with BayesPrism and InstaPrism. [SINGLE, SNIPPET] — [Cell Rep Methods 2026](https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(26)00046-9)
- **InstaPrism [PK].** A fast re-implementation of the BayesPrism model. — [Hu & Chikina, Bioinformatics 2024, doi:10.1093/bioinformatics/btae440](https://doi.org/10.1093/bioinformatics/btae440)
- **MuSiC2 [PK].** Deconvolution across conditions, accounting for genes that differ between conditions. — [Fan et al., Brief Bioinform 2022, doi:10.1093/bib/bbac430](https://doi.org/10.1093/bib/bbac430)
- **Kassandra [PK].** A tree-based deconvolution model trained on large artificial mixtures, validated against cytometry and IHC. — [Zaitsev et al., Cancer Cell 2022, doi:10.1016/j.ccell.2022.07.006](https://doi.org/10.1016/j.ccell.2022.07.006)
- **EcoTyper [PK].** CIBERSORTx-based cell states and multicellular "carcinoma ecotypes" (CE1–CE10) across 16 carcinomas including CRC. Ecotypes were associated with survival and immunotherapy response. — [Luca et al., Cell 2021, doi:10.1016/j.cell.2021.09.014](https://doi.org/10.1016/j.cell.2021.09.014)
- **Older pipeline benchmark [PK, foundational].** Data transformation and normalisation choices strongly affect deconvolution accuracy. — [Avila Cobos et al., Nat Commun 2020, doi:10.1038/s41467-020-19015-1](https://doi.org/10.1038/s41467-020-19015-1) [IND-BENCH]

**CRC-specific links from single cells to bulk outcomes**
- **iCMS/IMF.**
  - 715 iCMS markers were applied to bulk from 15 datasets, and survival analysis was run. — [Joanito et al., Nat Genet 2022](https://www.nature.com/articles/s41588-022-01100-4)
  - A stage IV study used iCMS2/iCMS3 scores from single-cell deconvolution with overall survival as the primary outcome. [SINGLE; venue unverified] — [ResearchGate](https://www.researchgate.net/publication/376893955_Single-cell_deconvolution_reveals_high_lineage-and_location-_dependent_heterogeneity_in_mesenchymal_multivisceral_stage_IV_colorectal_cancer)
- **Deconvolved immune states versus clinical variables for CRC recurrence (Kim et al., iScience 2022).**
  - Method: MethylCIBERSORT with a sorted-cell methylation reference (not single-cell).
  - Cohorts: internal SMC (114 patients) and TCGA (106 patients).
  - Immune-only AUC: about 0.67–0.69 in SMC and 0.69 in TCGA.
  - SMC: combined model 0.74 against clinical-only 0.61. TCGA: clinical-only 0.50, and TNM added little.
  - Small and retrospective. — [PMC9636036](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9636036/)
- **FAP+ CAF–SPP1+ macrophage axis [PK].** Single-cell and spatial data show the two co-localise. A high FAP/SPP1 score was associated with poor prognosis in bulk CRC cohorts. — [Qi et al., Nat Commun 2022, doi:10.1038/s41467-022-29366-6](https://doi.org/10.1038/s41467-022-29366-6)
- **EMP1+ high-relapse-cell program [PK].** Predicts relapse in human CRC cohorts. — [Cañellas-Socias et al., Nature 2022](https://doi.org/10.1038/s41586-022-05402-9)
- **Typical 2025 papers that pair single-cell and bulk data with machine learning.**
  - Fibroblast prognostic markers validated on GSE17536 with LOOCV. — [PMC12499553](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12499553/)
  - SPINK1/TIMP1 epithelial markers from TISCH linked to survival. — [PMC12733225](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12733225/)
  - Both are [SINGLE]: low evidence, with no ΔC over clinical reported in the summaries.
- **Random-signature baseline [PK, foundational].** Most random gene signatures are significantly associated with outcome in breast cancer. Any single-cell-derived signature therefore needs a comparison against random or proliferation baselines. — [Venet et al., PLoS Comput Biol 2011, doi:10.1371/journal.pcbi.1002240](https://doi.org/10.1371/journal.pcbi.1002240)
- **TCGA endpoint guidance [PK].** The TCGA Clinical Data Resource cautions about OS for some cancer types and recommends PFI/DFI where follow-up is short. This is relevant to TCGA-COAD/READ, which has few OS events. — [Liu et al., Cell 2018, doi:10.1016/j.cell.2018.02.052](https://doi.org/10.1016/j.cell.2018.02.052)

### Inferences
- **Recommended transfer path:**
  1. Define CRC cell states and programs on the integrated atlas (patient-replicated).
  2. Deconvolve bulk cohorts with BayesPrism (InstaPrism for speed and sensitivity runs; ReCIDE worth testing given the 2026 benchmark). Use a reference at state resolution: iCMS2/3 epithelium, HRC/EMP1, CAF subtypes such as mCAF/FAP+, SPP1+/C1QC+ macrophages, CXCL13+ T cells.
  3. Enter the fractions, or state-specific expression scores (BayesPrism's malignant-compartment expression), as covariates in a penalised Cox model alongside age, stage, MSI/MMR and possibly CMS.
  4. Report ΔC and time-dependent AUC under nested CV, plus external validation (cohorts with RFS such as GSE39582, PETACC-3, CIT).
- **Use Scissor and scAB for hypothesis generation only.** Their cell selection is tuned on the same bulk phenotype, so evaluating on that cohort is circular.
- **Expected effect size is small.** The CRC MethylCIBERSORT study and the 2026 real-bulk benchmark suggest modest discrimination (AUC about 0.7). A ΔC of about 0 over age+stage is a realistic outcome even for a well-built pipeline. The successor should pre-register ΔC as the primary metric and include negative controls: random signatures, a proliferation score, and CMS4/stromal content.
- **Malignant-state expression may beat gene-level rankings.** BayesPrism's estimate of malignant-compartment expression is the most direct way to carry tumour-intrinsic single-cell signal into bulk without stromal confounding, and stromal confounding is exactly the EMT/VIM issue.

### Gaps
- PIPET could not be verified (venue, year, benchmark); it was not retrieved this session.
- No head-to-head comparison of Scissor, scAB, scPAS, SIDISH and DEGAS was found that reports incremental C-index over clinical covariates in CRC.
- The full text of SigBridgeR (rate-limited) and the per-cancer results of the Genome Biology 2026 benchmark (paywall redirect) were not read. Whether CRC is among the 9 cancer types is unconfirmed.
- No independent benchmark of EcoTyper or Kassandra for CRC survival was found.

## 4. Do single-cell foundation models actually help biomarker discovery or prognosis?

### Takeaway
No demonstrated gain. Independent benchmarks from 2025–2026 agree on three points:
- Zero-shot embeddings from scGPT, Geneformer, UCE and TranscriptFormer are at best on par with HVG/PCA, scVI or Harmony, and are often worse.
- Perturbation-prediction models do not beat additive or linear baselines (Ahlmann-Eltze et al., Nature Methods 2025). A 2026 Nature Biotechnology rebuttal argues that they do beat uninformative baselines under calibrated metrics.
- Performance degrades under dataset or platform shift.

No study was found in which foundation-model features improved bulk survival prediction over clinical covariates. Nicheformer (spatial) is the most relevant new model, but it is validated on spatial-label tasks, not outcomes. Treat foundation models as an optional ablation arm with mandatory simple baselines.

### Cited Findings
- **Kedzierska et al., Genome Biology 2025 [IND-BENCH].** Zero-shot scGPT and Geneformer embeddings performed inconsistently against HVG selection, scVI and Harmony, and trailed scVI and Harmony on batch integration. — [Genome Biol 2025](https://link.springer.com/article/10.1186/s13059-025-03574-x)
- **Ahlmann-Eltze, Huber & Anders, Nature Methods 2025 (doi:10.1038/s41592-025-02772-6) [IND-BENCH].**
  - Five foundation models and two other deep-learning models were tested on single and double perturbations. None beat a "no change" or additive baseline.
  - A simple linear model pretrained on a different perturbation dataset beat all of them.
  - The authors suggest observational pretraining data may be one cause. — [PMC12328236](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12328236/); [code/data Zenodo](https://zenodo.org/records/16092690); [GitHub](https://github.com/const-ae/linear_perturbation_prediction-Paper)
- **Rebuttal: deep perturbation models beat uninformative baselines on calibrated metrics.** Eight models across four datasets, in a bioRxiv preprint (Oct 2025) and a Nature Biotechnology version (2026). This is a contested point, not settled. — [bioRxiv 2025.10.20.683304](https://www.biorxiv.org/content/10.1101/2025.10.20.683304.full.pdf); [Nat Biotechnol 2026, s41587-026-03307-w](https://www.nature.com/articles/s41587-026-03307-w)
- **"One PCA still rules them all" (arXiv 2410.13956) [IND-BENCH, PREPRINT].** Benchmarked scGPT, Geneformer, CellPLM and UCE against PCA and scVI for perturbation analysis. — [arXiv 2410.13956](https://arxiv.org/pdf/2410.13956)
- **Souza & Mehta, Feb 2026 [PREPRINT].** Careful normalisation plus linear methods reaches state-of-the-art or near it on several benchmarks, beating foundation models on some out-of-distribution tasks. — [PMC12918925](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12918925/)
- **Zero-shot benchmark, bioRxiv Aug 2026 [IND-BENCH, PREPRINT].** Covers scGPT, SCimilarity, UCE and TranscriptFormer. No embedding consistently beats the baseline for supervised cell-type annotation. TranscriptFormer showed near-zero cell-type silhouette. — [bioRxiv 2026.08.03.739553](https://www.biorxiv.org/content/10.64898/2026.08.03.739553v1.full)
- **CellBench-LS (bioRxiv Apr 2026) [PREPRINT].** Seven foundation models against PCA, UMAP and scVI across five tasks, zero-shot and few-shot. — [bioRxiv 2026.04.01.714123](https://www.biorxiv.org/content/10.64898/2026.04.01.714123v1.full)
- **Harmonised single-cell and spatial foundation-model benchmark (arXiv 2607.17227) [PREPRINT].** Widely used models can underperform simpler baselines and degrade under dataset or platform shift. [SNIPPET] — [arXiv 2607.17227](https://arxiv.org/pdf/2607.17227)
- **scFME in-silico perturbation evaluation (bioRxiv Sep 2025) [PREPRINT].** A second-generation Geneformer trained with extra cancer data beat Geneformer v1 and two scGPT variants. This is a relative ranking, not a comparison with simple baselines. — [bioRxiv 2025.09.22.677811](https://www.biorxiv.org/content/10.1101/2025.09.22.677811.full.pdf)
- **Nicheformer (Nature Methods 2025).** A transformer pretrained on SpatialCorpus-110M (more than 110 M dissociated and spatial cells, human and mouse). Evaluated on spatial density, niche and region label prediction, and on transferring spatial context onto dissociated data. [DEV-BENCH] — [Nat Methods 2025](https://www.nature.com/articles/s41592-025-02814-z); [PMC12695652](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12695652/)
- **TranscriptFormer (CZI, bioRxiv Apr 2025).** Claims F1 above 0.65 across distant species and better transfer than UCE. [DEV-BENCH] Contradicted on structure preservation by the Aug 2026 zero-shot benchmark above. — [bioRxiv 2025.04.25.650731](https://www.biorxiv.org/content/10.1101/2025.04.25.650731v1.full.pdf)
- **Original foundation-model papers [PK].**
  - Geneformer. — [Theodoris et al., Nature 2023, doi:10.1038/s41586-023-06139-9](https://doi.org/10.1038/s41586-023-06139-9)
  - scGPT. — [Cui et al., Nat Methods 2024, doi:10.1038/s41592-024-02201-0](https://doi.org/10.1038/s41592-024-02201-0)
  - scFoundation. — [Hao et al., Nat Methods 2024, doi:10.1038/s41592-024-02305-7](https://doi.org/10.1038/s41592-024-02305-7)
- **Commentary [opinion, not evidence].** Single-cell foundation models have underdelivered for drug discovery. — [Drug Target Review](https://www.drugtargetreview.com/home/why-single-cell-foundation-models-have-underdelivered-and-what-drug-discovery-needs-instead/2136187.article)

### Inferences
- **Hype versus demonstrated gains.** The claimed gains come from developer benchmarks. Independent benchmarks across 2025–2026 consistently show PCA, HVG, scVI and linear baselines matching or beating zero-shot foundation-model representations. Survival-relevant biomarker discovery has not been benchmarked at all.
- **What would justify including a foundation model:**
  - it is cheap to run, for example frozen embeddings used as an alternative cell-state space;
  - it is compared against scVI/PCA using the **same downstream ΔC** metric;
  - it is reported as an ablation, never as the core method.
- **Spatial models may matter later.** Nicheformer-style models that carry spatial context onto dissociated cells could eventually let niche labels be inferred for dissociated atlases. There is no outcome evidence yet.

### Gaps
- STATE (Arc Institute, 2025), CellFM, UCE (peer-reviewed version) and the Virtual Cell Challenge results were not retrieved or verified in this session.
- No paper was found that tests foundation-model-derived features for bulk CRC (or any cancer) survival prediction against clinical covariates.

## 5. Spatial biomarkers in CRC: niches, cell–cell communication, recurrence-predictive domains, image-derived features

### Takeaway
- **Validated so far:** the only clinically validated spatial CRC biomarkers remain immune-density scores (Immunoscore: CD3/CD8 in tumour core and invasive margin) and H&E deep-learning scores (for example DoMore-v1-CRC).
- **2024–2026 sequencing-based spatial studies are small and mostly descriptive:** Visium HD n=5, CosMx 23–25 patients, Xenium 16 patients.
- **Prognostic translation goes through H&E classifiers or bulk projection.** STORM's "desmoplastic fibrotic barrier" gives HR 1.50 in stage III/IV TCGA, stratified rather than adjusted for stage, and is a preprint. A CosMx-derived COL1A1/αSMA stromal signature reached 5-year RFS 0.88 vs 0.49 in a GeoMx cohort, from a conference abstract.
- **Cell–cell communication inference is useful for hypotheses, not as a biomarker.** Methods disagree substantially, so consensus frameworks (LIANA+) are recommended.

### Cited Findings
- **Immunoscore [PK, foundational, IND-validated].**
  - International validation in stage I–III colon cancer: CD3/CD8 density in tumour core and invasive margin was prognostic for recurrence, independent of TNM and other factors.
  - This is the benchmark any new spatial CRC biomarker has to beat. — [Pagès et al., Lancet 2018, doi:10.1016/S0140-6736(18)30789-X](https://doi.org/10.1016/S0140-6736(18)30789-X)
- **H&E deep learning [PK].** DoMore-v1-CRC, an H&E deep-learning prognostic marker, was validated in independent cohorts. — [Skrede et al., Lancet 2020, doi:10.1016/S0140-6736(19)32998-8](https://doi.org/10.1016/S0140-6736(19)32998-8)
- **CODEX neighbourhoods [PK, foundational].** Coupling of cellular neighbourhoods at the CRC invasive front differs between Crohn's-like-reaction and diffuse-inflammatory-infiltration patients, which differ in survival. — [Schürch et al., Cell 2020, doi:10.1016/j.cell.2020.07.005](https://doi.org/10.1016/j.cell.2020.07.005)
- **Orion [PK].** Same-section high-plex immunofluorescence plus H&E in CRC. An immune-profile score was associated with progression-free survival. — [Lin et al., Nat Cancer 2023, doi:10.1038/s43018-023-00576-1](https://doi.org/10.1038/s43018-023-00576-1); data: [HTA CRC Atlas X](https://labsyspharm.github.io/orion-crc)
- **STORM (Xenium, bioRxiv Jul 2026) [PREPRINT].**
  - Method: a graph-attention VAE that finds 10 spatial motifs without labels. The desmoplastic fibrotic barrier (DFB) is associated with restricted CD8+ access to tumour cores, independent of CD8+ abundance.
  - H&E transfer: a ViT-B/16 trained on about 300k H&E tiles reaches tile accuracy 0.82 and was applied to TCGA CRC slides.
  - Full TCGA cohort (n=412, 90 events): not significant after FDR.
  - Stage III/IV subset (n=177, 58 events): DFB-barrier composite HR 1.50 (95% CI 1.14–1.97, q=0.032). VBMS fraction HR 0.71.
  - No C-index reported. Not adjusted for stage, MSI or CMS. — [PMC13405321](https://pmc.ncbi.nlm.nih.gov/articles/PMC13405321/)
- **CosMx collagen signature (JCO 2025 supplement) [ABSTRACT].** CosMx 1,000-plex on 25 resectable CRC patients. In the GeoMx validation cohort, a COL1A1 signature in the αSMA compartment gave 5-year RFS 0.88 vs 0.49 (p<0.005). [SNIPPET] — [JCO 2025;43(4_suppl):227](https://ascopubs.org/doi/10.1200/JCO.2025.43.4_suppl.227)
- **CosMx stromal niches (Pattison et al.) [PREPRINT].** Nine niches. Stromal niches influence tumour differentiation; no outcome link. — [bioRxiv 2025.12.09.693138](https://www.biorxiv.org/content/10.64898/2025.12.09.693138v1.full)
- **Visium HD macrophage niches.** n=5, descriptive. — [Oliveira et al., Nat Genet 2025](https://pubmed.ncbi.nlm.nih.gov/40473992/)
- **FAP+ CAF/SPP1+ TAM spatial co-localisation with bulk prognostic association [PK].** — [Qi et al., Nat Commun 2022](https://doi.org/10.1038/s41467-022-29366-6)
- **Pelka multicellular immune hubs [PK].** Spatially organised immune hubs identified by co-variation of cell programs and confirmed with multiplex imaging. — [Cell 2021](https://doi.org/10.1016/j.cell.2021.08.003)
- **Cell–cell communication benchmark [PK, IND-BENCH].** 16 resources and 7 methods compared: low agreement between methods and between resources. This motivated the LIANA consensus framework. — [Dimitrov et al., Nat Commun 2022, doi:10.1038/s41467-022-30755-0](https://doi.org/10.1038/s41467-022-30755-0)
- **LIANA+ [PK].** A scverse-compatible framework that unifies cell–cell communication methods and resources, including spatially informed and multi-view inference. — [Dimitrov et al., Nat Cell Biol 2024, doi:10.1038/s41556-024-01469-w](https://doi.org/10.1038/s41556-024-01469-w)
- **NicheNet [PK, foundational].** Infers ligands that regulate target genes, using prior signalling and gene-regulatory knowledge. — [Browaeys et al., Nat Methods 2020, doi:10.1038/s41592-019-0667-5](https://doi.org/10.1038/s41592-019-0667-5)
- **CellChat v2 [PK].** Adds spatial-proximity-constrained inference and an expanded database. — [Jin et al., Nat Protoc 2025, doi:10.1038/s41596-024-01045-4](https://doi.org/10.1038/s41596-024-01045-4)
- **Spatial deconvolution benchmark [PK, IND-BENCH].** cell2location, SpatialDWLS and RCTD were among the top methods for deconvolving spot-based spatial data. — [Li et al., Nat Methods 2022, doi:10.1038/s41592-022-01480-9](https://doi.org/10.1038/s41592-022-01480-9)
- **Spatial domain benchmark [PK, IND-BENCH].** Benchmark of spatial clustering and domain methods: no universal winner, and performance is dataset-dependent. — [Yuan et al., Nat Methods 2024, doi:10.1038/s41592-024-02215-8](https://doi.org/10.1038/s41592-024-02215-8)

### Inferences
**For a bulk-applicable successor, spatial data have three roles:**
1. **Validation.** Check that the candidate prognostic states co-localise in plausible niches, for example a FAP+ CAF/SPP1+ macrophage barrier or DFB excluding CD8 cells. This gives biological credibility to features that bulk deconvolution can estimate.
2. **Niche signatures for bulk.** Derive niche-level gene signatures or niche-state references and project them onto bulk, using them as Cox covariates.
3. **An H&E route.** TCGA COAD/READ has diagnostic slides, so ST-trained H&E classifiers are a parallel route.

**Bar to clear:** every spatial feature should be tested against stage, MSI and an Immunoscore-like CD8/CD3 proxy, plus stromal content (CMS4), not against stage alone. STORM's stratified-only analysis is exactly the gap a rigorous ΔC framework would expose.

**Cell–cell communication** should be used to annotate mechanisms (consensus LIANA+), not as a source of prognostic features. Its poor inter-method agreement means any ligand–receptor score needs robustness checks across methods.

### Gaps
- No CRC study was found that combines a sequencing-based spatial niche feature with recurrence in more than about 100 patients and reports ΔC over stage plus MSI.
- The outcome analyses of the neutrophil Cancer Cell 2025 spatial cohorts (266 patients) were not checked.
- Recent 2024–2026 H&E-to-transcriptome models (HEST-1k, iStar and others) and their CRC prognostic validation were not researched in this session.

## 6. Single-cell multi-omics (10x Multiome, CITE-seq, scNMT) and integration tools: relevance to regulatory or miRNA-adjacent signals

### Takeaway
**Integration tools:** no universal winner. The 2025 Nature Methods multitask benchmark favours:
- Seurat WNN, sciPENN and Multigrate for vertical integration of RNA+ADT;
- GLUE and Seurat v5 for diagonal integration;
- MOFA+, Matilda and scMoMaT as the few methods that also select features.

A 2024 RNA+ATAC benchmark ranks GLUE first, then MultiVI.

**miRNA:** standard 3′ droplet scRNA-seq (and the C1 data in GSE81861) does not measure mature miRNAs. Single-cell miRNA-adjacent signal is reachable only indirectly:
- pri-miRNA/host-gene expression or chromatin accessibility (Multiome);
- TF regulons (SCENIC+);
- miRNA target-set activity scores (decoupler-style).

None of these was found to be validated for CRC prognosis. For CRC, the main chromatin resource is the HTAN Stanford polyp→CRC scATAC/scRNA dataset.

### Cited Findings
- **Multitask benchmark of single-cell multimodal integration (Nature Methods 2025) [IND-BENCH].**
  - Vertical integration (RNA+ADT): Seurat WNN, sciPENN and Multigrate were top for preserving cell types.
  - Feature selection: only Matilda, scMoMaT and MOFA+ support selecting molecular markers.
  - Diagonal integration: GLUE and Seurat v5 were top performers; Seurat used much more time and memory.
  - Code: scMultiBench. — [Nat Methods 2025, s41592-025-02856-3](https://www.nature.com/articles/s41592-025-02856-3); [GitHub PYangLab/scMultiBench](https://github.com/PYangLab/scMultiBench)
- **RNA+ATAC integration benchmark (Briefings in Bioinformatics 2024) [IND-BENCH].** GLUE was best for omics mixing, cell-type conservation and alignment, followed by MultiVI. MOFA+ was good at cell-type conservation. — [Brief Bioinform 25(2):bbae095](https://academic.oup.com/bib/article/25/2/bbae095/7630478); [PMC10944570](https://pmc.ncbi.nlm.nih.gov/articles/PMC10944570/)
- **Mosaic/bridge integration benchmark [PREPRINT, SNIPPET].** Seurat v5 bridge and MultiVI beat GLUE on the fraction of samples closer than threshold; GLUE had the highest cell-match score. — [ResearchGate](https://www.researchgate.net/publication/393572896_Benchmarking_single-cell_multi-modal_data_integrations)
- **Integration tool papers [PK].**
  - GLUE: graph-linked unified embedding for unpaired multi-omics. — [Cao & Gao, Nat Biotechnol 2022, doi:10.1038/s41587-022-01284-4](https://doi.org/10.1038/s41587-022-01284-4)
  - MultiVI. — [Ashuach et al., Nat Methods 2023, doi:10.1038/s41592-023-01909-9](https://doi.org/10.1038/s41592-023-01909-9)
  - totalVI (CITE-seq). — [Gayoso et al., Nat Methods 2021, doi:10.1038/s41592-020-01050-x](https://doi.org/10.1038/s41592-020-01050-x)
  - MOFA+. — [Argelaguet et al., Genome Biol 2020, doi:10.1186/s13059-020-02015-1](https://doi.org/10.1186/s13059-020-02015-1)
  - Seurat v5 bridge integration. — [Hao et al., Nat Biotechnol 2024, doi:10.1038/s41587-023-01767-y](https://doi.org/10.1038/s41587-023-01767-y)
- **SCENIC+ [PK].** Enhancer-driven gene-regulatory networks from paired or unpaired scRNA+scATAC. This is the main route to regulatory programs. — [Bravo González-Blas et al., Nat Methods 2023, doi:10.1038/s41592-023-01938-4](https://doi.org/10.1038/s41592-023-01938-4)
- **decoupler [PK].** Scores pathway, TF or any gene-set activity per cell or sample, so it could score miRNA target sets the way CcSinglecell's targeting data would. — [Badia-i-Mompel et al., Bioinform Adv 2022, doi:10.1093/bioadv/vbac016](https://doi.org/10.1093/bioadv/vbac016)
- **Single-cell miRNA measurement [PK].**
  - Dedicated single-cell small-RNA sequencing is low-throughput (foundational). — [Faridani et al., Nat Biotechnol 2016, doi:10.1038/nbt.3701](https://doi.org/10.1038/nbt.3701)
  - Smart-seq-total captures small and long RNAs from the same cell. — [Isakova et al., PNAS 2021, doi:10.1073/pnas.2113568118](https://doi.org/10.1073/pnas.2113568118)
- **CRC chromatin atlas [PK].** HTAN Stanford scATAC plus scRNA along the polyp→CRC continuum. — [Becker et al., Nat Genet 2022](https://doi.org/10.1038/s41588-022-01088-x)
- **Phenotype-associated cells from scATAC (PACells, AJHG 2026)** extend family-A bulk linking to chromatin data. — [AJHG 2026](https://www.sciencedirect.com/science/article/abs/pii/S0002929726001229)

### Inferences
- **The miRNA layer.** CcSinglecell's miRNA-targeting idea is best re-expressed at single-cell level as **target-set activity scores** (decoupler with a curated miRNA–target network) computed on pseudobulk per patient × state.
  - These scores can be checked against pri-miRNA host-gene expression or Multiome accessibility at miRNA loci (for example MIR17HG, MIR21) in Becker/HTAN data.
  - In bulk, miRNA expression itself is available (TCGA miRNA-seq), so the single-cell layer only needs to say *which cell compartment* a miRNA program acts in.
  - This is a design inference. No paper was found demonstrating that this improves CRC prognosis.
- **Tool choice for the successor:**
  - Multiome (RNA+ATAC, paired): MultiVI, or MOFA+ when interpretable factors and feature selection are wanted.
  - Unpaired RNA with ATAC from another study: GLUE.
  - CITE-seq: totalVI or WNN.
  - Multi-omics is a secondary, mechanism-annotation layer. It is not needed to build bulk survival features.

### Gaps
- No CRC single-cell or spatial miRNA dataset was found. Spatial miRNA methods were not researched.
- No scNMT-seq (methylation plus accessibility plus RNA) CRC dataset was found in this session.
- No evidence was found that single-cell multi-omic or regulatory features add prognostic value over clinical covariates in CRC bulk cohorts.
