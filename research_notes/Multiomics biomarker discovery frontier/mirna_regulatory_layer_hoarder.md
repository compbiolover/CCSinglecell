# miRNA / ncRNA regulatory layer for CcSinglecell successor: hoardeR status and the 2024–2026 state of the art

Notes as of 2026-11-07. "Verified live" means I checked it directly this session (curl or CRAN/Bioconductor page). Citations marked "(prior knowledge)" are well-known papers whose DOI I did not re-fetch this session. The report writer should treat them as solid but not re-checked here.

## 1. hoardeR: maintainer, CRAN status, what targetScan() does, does it still work?

### Takeaway
Do not port hoardeR. It was archived from CRAN on 2026-04-13. Its `targetScan()` is a line-by-line HTML scraper of the TargetScan **7.1** web CGI, and that endpoint now returns HTTP 403. The function also never looks at release versions: `release` is accepted but ignored. Everything it gave CcSinglecell (the top-N targets per miRNA by cumulative weighted context++) is available as a single bulk download from TargetScan 8.0, or offline through Bioconductor.

### Cited Findings
- CRAN status: "Package 'hoardeR' was removed from the CRAN repository … Archived on 2026-04-13 as issues were not corrected despite reminders." — [CRAN hoardeR page](https://cran.r-project.org/package=hoardeR)
- Final CRAN checks (version 0.10, 2026-04-13): "Installation failed" on seven flavors (r-devel Debian clang/gcc, Fedora clang/gcc, macOS arm64, Windows x86_64; r-patched Linux). The vignette rebuild failed on r-oldrel-windows (missing PNG, LaTeX failure). r-release passed on Linux, macOS and Windows. — [CRAN check results archive 2026-04-13](https://cran-archive.r-project.org/web/checks/2026/2026-04-13_check_results_hoardeR.html)
- Release history: 0.0-1 (2014-02-20), 0.0-2 (2014-06), 0.1 (2015-06), 0.9 (2016-09), 0.9.1 and 0.9.2 (2016-10), 0.9.4 and 0.9.4-2 (2019-02), **0.10 (2024-04-08, last)**. — [CRAN archive](https://cran.r-project.org/src/contrib/Archive/hoardeR/)
- Maintainer: Daniel Fischer (Luke / Natural Resources Institute Finland), sole aut/cre, `daniel.fischer@luke.fi`. License GPL (>= 2). Title: "Collect and Retrieve Annotation Data for Various Genomic Data". Description: "Cross-species identification of novel gene candidates using the NCBI web service is provided. Further, sets of miRNA target genes can be identified by using the targetscan.org API." Imports include bamsignals and Biostrings. — DESCRIPTION file in [hoardeR_0.10.tar.gz](https://cran.r-project.org/src/contrib/Archive/hoardeR/hoardeR_0.10.tar.gz) (verified live)
- Other functions: `blastSeq`, `getBlastResult`, `findSpecies`, `getEnsgInfo`, `getGeneLocation`, `getGeneSeq`, `getSequenceFromNCBI`, `getFastaFromBed`, `getAnnotation`, `intersectXMLAnnot`, `coverageDensity`, `plotCoverage`, `plotHit`, `createSpecies`, `tableSpecies`, plus summary and print methods for `fa`, `fq` and `blastRes`. Most of these do cross-species NCBI BLAST mining. Only `targetScan()` is relevant to CcSinglecell. — R/ directory and NAMESPACE of [hoardeR_0.10.tar.gz](https://cran.r-project.org/src/contrib/Archive/hoardeR/hoardeR_0.10.tar.gz) (verified live)
- What `targetScan()` does:
  - Signature: `targetScan(mirna=NULL, species=NULL, release="7.1", maxOut=NULL)`.
  - It builds `http://www.targetscan.org/cgi-bin/targetscan/vert_71/targetscan.cgi?species=<Species>&mirg=<miRNA>`. **vert_71 is hard-coded and `release` is not used in the URL.**
  - It reads the page with `scan()` and splits it with `grepl`/`strsplit` on fixed HTML markers (`<th>total</th>`, `target=new>`, `<td nowrap>`, `<B>`). It returns a data.frame of Ortholog, geneName, consSites and poorlySites. `maxOut` keeps the first N rows.
  - It guesses species from the hsa/mmu prefix. If the name matches several families, it retries alternative names containing "/".
  - — [GitHub fischuu/hoardeR R/targetScan.R](https://raw.githubusercontent.com/fischuu/hoardeR/master/R/targetScan.R); same lines (56, 86) in CRAN 0.10 source (verified live)
- Live endpoint test (2026-11-07):
  - `https://www.targetscan.org/cgi-bin/targetscan/vert_71/targetscan.cgi?species=Human&mirg=hsa-miR-21-5p` returns **HTTP 403 "403 forbidden"**.
  - The same query on vert_72 and vert_80 returns 200 with a full table.
  - hoardeR's plain-http URL failed in this sandbox with a proxy DNS error, so that path is inconclusive.
  - The vert_72 and vert_80 result pages say "Table sorted by cumulative weighted context++ score". The default view lists targets with conserved sites, and a separate link shows targets "irrespective of site conservation".
  - — [TargetScan 8.0 CGI](https://www.targetscan.org/cgi-bin/targetscan/vert_80/targetscan.cgi?species=Human&mirg=hsa-miR-21-5p) (verified live)
- TargetScan 8.0 bulk files are downloadable: `vert_80_data_download/Summary_Counts.default_predictions.txt.zip` returned HTTP 200 (verified live). Other files in use include Conserved_Site_Context_Scores, Nonconserved_Site_Context_Scores and miR_Family_Info. — [enrichMiR DownloadFunctions.R](https://rdrr.io/github/ETHZ-INS/enrichMiR/src/R/DownloadFunctions.R); [TargetScan 8.0](https://www.targetscan.org/vert_80/)
- CcSinglecell legacy usage: `inst/legacy/Code/mirna_calculator_original.R` does `require(hoardeR)` and then calls `targetScan(mirna=common_mirnas[...], ...)` in a loop over the dbDEMC-high ∩ miRmap miRNAs. `inst/legacy/paper.Rmd:1664` does the same. — local repo /home/user/CCSinglecell (verified)

### Inferences
- The Cells 2023 paper says "TargetScan 7.2", but hoardeR has always queried vert_71. The published miRNA layer was very likely built on **TargetScan 7.1 web output**, not 7.2. This is worth stating in the successor's methods and changelog.
- With `maxOut = 10`, CcSinglecell took the top 10 targets per miRNA from the conserved-site table, sorted by cumulative weighted context++ score. The same data comes from TargetScan 8.0's `Summary_Counts.default_predictions.txt` (or the context-score files): filter by miR family and slice the 10 rows with the lowest cumulative weighted context++. That is one `data.table` read and offline, with no scraping, rate limits or HTML fragility. "Porting and improving hoardeR" means about 10 lines of R that read a versioned bulk file, not a web client.
- Even a fixed scraper would be the wrong design: one HTTP request per miRNA, results tied to the web UI's default (conserved-only) view, and no record of the release.
- The remaining hoardeR functions (BLAST and NCBI mining) are irrelevant to the successor.

### Gaps
- I did not find a public statement from TargetScan on when or why vert_71 CGI started returning 403. A hoardeR GitHub issue about the breakage was not searched.
- I did not test whether hoardeR's marker-based parser still works on vert_80 HTML. vert_80 pages contain `target=new>` but did not show `<th>total</th>` in a quick grep, so the parser would probably need changes. This is moot if bulk files are used.

## 2. Current miRNA target resources: versions, bulk access, licensing, update cadence

### Takeaway
Pick resources by evidence class.
- **Prediction** (TargetScan 8.0 context++; scanMiR biochemical KdModels): free bulk download or a Bioconductor package.
- **Experimentally supported:** TarBase v9 (2024, ~6M entries, cell-type-annotated, includes AGO-CLIP) and miRTarBase 2025 (NAR 2025, LLM-assisted curation).
- **CLIP-supported:** ENCORI (still online).
- **Naming:** standardize on miRBase v22 accessions, with MirGeneDB 3.0 (NAR 2025) family and seed-curated names.

multiMiR 1.34 bundles a 2024-08-28 DB snapshot with TargetScan 8, TarBase 9 and miRTarBase 9. That is convenient but lags miRTarBase 2025. targetscan.Hs.eg.db is no longer in Bioconductor release.

### Cited Findings
- **TargetScan 8.0:** web CGI and bulk files live (verified live). Release 8.0 human file: 264,563 conserved-site interactions, 405 miRNAs, 13,077 target genes, as packaged by CyTargetLinker. — [CyTargetLinker TargetScan linksets](https://cytargetlinker.github.io/pages/linksets/targetscan). TargetScan 8 was released with McGeary et al. 2022, eLife 11:e69803 — [doi:10.7554/eLife.69803](https://doi.org/10.7554/eLife.69803) (prior knowledge). I found no newer release than 8.0 (2021).
- **miRTarBase 2025:**
  - Citation: Li Z. et al., NAR 53(D1):D147–D156, published 2024-11-23, DOI 10.1093/nar/gkae1072, PMID 39578692.
  - Content: more than 3,817,550 validated MTIs from 13,690 articles. First release to put an LLM (Meta Llama 3) into its NLP curation pipeline. Adds miRNA–drug interactions, biomarker annotations and oxidized miRNAs.
  - Hosted at mirtarbase.cuhk.edu.cn.
  - "Version 10" appears only on a Wikipedia-derived page and is unverified.
  - — [PMC11701613](https://pmc.ncbi.nlm.nih.gov/articles/11701613); [Wikipedia](https://en.wikipedia.org/wiki/MiRTarBase)
- **DIANA-TarBase v9.0:**
  - Paper: NAR 2024, 52(D1):D304, gkad1071.
  - Size: about 6 million entries, described as roughly 3× miRTarBase and 6× TarBase v8. Adds ~34,000 viral-miRNA/host pairs, TDMD events, and 172 tissues/cell types across 37 protocols.
  - Free, no registration. RNAcentral release 25 (June 2025) moved from TarBase v7 to v9.
  - — [NAR](https://academic.oup.com/nar/article/52/D1/D304/7438911); [PMC10767993](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10767993/); [TarBase site](https://dianalab.e-ce.uth.gr/tarbasev9); [RNAcentral blog](https://blog.rnacentral.org/2025/06/rnacentral-release-25.html)
- **ENCORI (starBase):** online at rnasysu.com/encori (HTTP 200, verified live). It integrates CLIP-seq–supported miRNA–mRNA, ceRNA, RBP–RNA and other interaction types. — [ENCORI](https://rnasysu.com/encori/) (verified live). I found no 2024–2026 version or update paper.
- **miRDB:** online (verified live), with predictions from MirTarget. multiMiR lists the source version as miRDB **6 (June 2019)**. — [miRDB](https://mirdb.org/); [multiMiR vignette](https://bioconductor.org/packages/release/bioc/vignettes/multiMiR/inst/doc/multiMiR.html)
- **mirDIP:** online at ophid.utoronto.ca/mirDIP (HTTP 200). The page refers to "Version 5" and has a license statement in which use "of the database is explicitly granted to non-profit organi[zations]", implying commercial restrictions. — [mirDIP](https://ophid.utoronto.ca/mirDIP/) (verified live). Exact sub-version and date were not extracted.
- **multiMiR (Bioconductor 3.23):**
  - Package v1.34.0, maintainer Spencer Mahaffey (CU Anschutz).
  - Bundled DB v2.4.0, updated 2024-08-28.
  - Sources: TargetScan 8 (Sept 2021), TarBase 9 (2023), miRTarBase 9 (Sept 2021), miRDB 6 (June 2019), DIANA-microT 5 (2013). Legacy sources are frozen: PITA 2008, MicroCosm 2009, miRanda 2010, ElMMo 2011, PicTar 2012, miR2Disease, PhenomiR. The vignette warns that old databases "carry over" old target IDs.
  - It is queried through a remote server (multimir.org), not fully offline.
  - — [Bioconductor multiMiR](https://bioconductor.org/packages/release/bioc/html/multiMiR.html); [vignette](https://bioconductor.org/packages/release/bioc/vignettes/multiMiR/inst/doc/multiMiR.html)
- **scanMiR (Bioconductor 3.23, v1.18.0)** and **scanMiRData (v1.18.0):**
  - Maintainer Pierre-Luc Germain (ETH Zürich).
  - KdModels (12-mer binding affinity models) for all human, mouse and rat miRBase miRNAs. Scans by seed, full sequence (3′ alignment) or KdModel. Predicts repression, slicing sites and TDMD sites.
  - — [scanMiR](https://bioconductor.org/packages/release/bioc/html/scanMiR.html); [scanMiRData](https://bioconductor.org/packages/release/data/experiment/html/scanMiRData.html) (verified live)
  - Paper: Soutschek et al., Bioinformatics 2022, [doi:10.1093/bioinformatics/btac110](https://doi.org/10.1093/bioinformatics/btac110) (prior knowledge). It builds on McGeary et al., Science 2019, [doi:10.1126/science.aav1741](https://doi.org/10.1126/science.aav1741) (prior knowledge), which reported that the biochemical affinity model outperformed TargetScan 7 context++ at predicting repression.
- **targetscan.Hs.eg.db:** returned 404 under Bioconductor release data/annotation, data/experiment and bioc. It is no longer in the current release (3.23). — [Bioconductor URL tested](https://bioconductor.org/packages/release/data/annotation/html/targetscan.Hs.eg.db.html) (verified live). The removal date was not found.
- **miRBaseConverter (Bioconductor 3.23, v1.36.0):** converts names and accessions across miRBase versions. Maintainer Taosheng Xu. — [Bioconductor](https://bioconductor.org/packages/release/bioc/html/miRBaseConverter.html) (verified live)
- **MirGeneDB 3.0:**
  - Paper: Clarke AW et al., NAR 2025, 53(D1):D116–D128, DOI 10.1093/nar/gkae1094, PMID 39673268.
  - Content: 21,822 miRNA genes in 1,743 families. Adds 33 invertebrate and 6 mammal species. Takes over naming of conserved novel metazoan miRNAs. Retrained covariance models usable with Infernal `cmsearch`. Part of RNAcentral and ELIXIR Norway.
  - — [PubMed 39673268](https://pubmed.ncbi.nlm.nih.gov/39673268/); [MirGeneDB](https://mirgenedb.org/covariance)
  - bayesReact (NAR 2026) uses both miRBase v22 and MirGeneDB 3.0 seeds, which shows the two are used together. — [PMC12884093](https://pmc.ncbi.nlm.nih.gov/articles/PMC12884093/)

### Inferences
- Recommended successor stack. Every source is bulk, offline, versioned and checksummed in a data manifest:
  1. TargetScan 8.0 bulk context++ files for the "predicted" layer, a drop-in replacement for hoardeR.
  2. scanMiR KdModels as a biochemical alternative or ensemble member, and for isomiR/non-canonical sites.
  3. TarBase v9, filtered to AGO-CLIP and to colon/CRC cell types where annotated, plus miRTarBase 2025 strong-evidence entries for the "experimentally supported" layer.
  4. miRBase v22 accession or MirGeneDB 3.0 family as the join key, with miRBaseConverter for legacy dbDEMC names.
- multiMiR is useful for quick cross-database lookups. Its miRTarBase is v9 (2021), not 2025, and it depends on a remote server, so it is weaker for reproducibility than pinned bulk files.
- Licensing differs by resource and is a practical constraint for a redistributable package. Do not ship third-party tables inside the package. Download at runtime with recorded versions.

### Gaps
- Exact current license terms for TargetScan, miRTarBase 2025 and TarBase v9 bulk downloads were not verified. mirDIP's non-profit clause was seen but not read in full.
- No 2024–2026 update papers were found for ENCORI, miRDB or miRWalk. miRWalk's current version was not checked.
- dbDEMC (used for CcSinglecell's cancer-miRNA list): the latest version (dbDEMC 3.0, 2022?) and its status were not checked.
- Whether miRBase has released anything after v22.1 (2018) was not verified this session.

## 3. Deep-learning miRNA target prediction (2023–2026): does it beat TargetScan context++ or biochemical models on held-out data?

### Takeaway
On site-level CLIP/chimeric benchmarks (miRBench 2025), deep models clearly beat seed rules and older CNNs. Gains in the newest models are small, e.g. AP 0.85 vs 0.84. On the question that matters for gene ranking (predicting **repression** of a gene), the best 2026 evidence (miRBind2) shows only modest gains over TargetScan 8 context++: Pearson 0.30 vs 0.24, ROC-AUC 0.60 vs 0.56. That comparison is on 7 miRNAs from the same transfection source used for training. No source I found shows a deep model decisively beating the McGeary biochemical model on held-out repression. RNA foundation models (RNA-FM, RiNALMo) have barely been applied to this task, and only self-reported results exist.

### Cited Findings
- **miRBench:**
  - Paper: Bioinformatics 2025, DOI 10.1093/bioinformatics/btaf233, published about 2025-07-15; bioRxiv 2024.11.14.623628.
  - Finding: a "miRNA frequency class bias" in negative-sample generation inflates reported performance of chimeric/CLIP-trained models. The authors propose bias-mitigated negatives, building more than 2.5M class-balanced interactions including a new AGO2 chimeric dataset.
  - Datasets: AGO2_CLASH_Hejret2023, AGO2_eCLIP_Klimentova2022, AGO2_eCLIP_Manakov2022. Python package `miRBench` (v1.0.1), bundling TargetScanCnn, miRBind and seed baselines.
  - — [PMC12261448](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12261448/); [bioRxiv](https://www.biorxiv.org/content/10.1101/2024.11.14.623628.full.pdf); [PyPI](https://pypi.org/project/mirbench)
- **miRBind2:**
  - Paper: Čechák … Alexiou, Bioinformatics 42 Suppl 2, btag410, published 2026-08-21.
  - Architecture: sequence-only CNN with about 92% fewer parameters than the previous baseline.
  - Site-level vs miRBenchCNN: AP 0.85 vs 0.84 and ROC-AUC 0.82 vs 0.81 on Manakov test; AP 0.83 vs 0.81 on left-out miRNAs.
  - Gene-level repression (miRBind2-3UTR vs TargetScan 8.0 weighted context++), on 50,549 miRNA–gene pairs (7 miRNAs × 7,486 genes; Agarwal 2015 HeLa transfections): Pearson 0.30 vs 0.24; Spearman 0.19 vs 0.15; ROC-AUC 0.60 vs 0.56; AP 0.47 vs 0.41.
  - Caveats: TargetScan assigns no score to 79.5% of pairs. Gene-level training and test data both come from Agarwal 2015.
  - — [Bioinformatics](https://academic.oup.com/bioinformatics/article/42/Supplement_2/btag410/8767265)
- **Expression-anchored benchmark:**
  - Paper: Jin … Qin, Ran (Emory), arXiv 2610.33960, 2026-10-25.
  - Design: TargetScan, miRDB, miRmap, DIANA-microT, RNA22, miRWalk and mirDIP were scored against paired miRNA/mRNA-seq plus miRTarBase/TarBase, on iPSC 22q11.2 deletion lines (n = 10).
  - Youden's J: mirDIP 0.566, DIANA-microT 0.481, miRmap 0.457, miRDB 0.123, TargetScan 0.040. TargetScan had the lowest Type I error (0.001) but power of only 0.042. "No single tool is optimal across all metrics."
  - Deep-learning tools were **excluded**.
  - — [arXiv 2610.33960](https://arxiv.org/html/2610.33960)
- **Older deep models (pre-2024):** TargetNet (Bioinformatics 2022) reported F1 0.774 vs TargetScan 0.471 on functional vs non-functional classification, with an author-chosen test set. miRAW (PLoS CB 2018) and chimiRic (PLoS CB 2016) also reported beating TargetScan. — [TargetNet](https://academic.oup.com/bioinformatics/article/38/3/671/6408435); [miRAW](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1006185); [chimiRic](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1005026)
- **DeepExoMir** (IJMS, 2026): an RNA language-model framework with mean AU-PRC 0.855 on miRBench CLIP test sets, beating 8 retrained miRBench baselines. Its comparison with TargetScan 8.0 is retrospective only. — [doi:10.3390/ijms27146184](https://doi.org/10.3390/ijms27146184)
- **RNA foundation models:**
  - RiNALMo (650M parameters; Nature Communications 16:5671, 2025) was evaluated on structure, splicing and translation tasks, not miRNA targeting.
  - DeepMiRT (GitHub) uses an RNA-FM encoder with miRNA–target cross-attention and claims first place among 12 methods on eCLIP benchmarks. This is self-reported, with no paper found.
  - — [RiNALMo arXiv](https://arxiv.org/html/2403.00043v2); [DeepMiRT](https://github.com/zichengll/DeepMiRT)
- **Biochemical model:** McGeary et al. 2019 reported that the CNN-derived affinity (Kd) model predicted repression better than TargetScan 7 context++. — [doi:10.1126/science.aav1741](https://doi.org/10.1126/science.aav1741) (prior knowledge)

### Inferences
- For a prognostic gene-ranking pipeline, the predictor matters less than (a) which miRNAs are included and (b) whether the miRNA is actually active in the tissue. Gene-level repression predictability tops out at about r = 0.3 for every method. Spend effort on activity evidence (Q4), not on swapping in a deep target predictor.
- A sensible successor design is pluggable target sets (TargetScan 8 context++, scanMiR Kd, TarBase-CLIP) plus a sensitivity analysis showing that the gene ranking is robust to the choice. Deep models (miRBind2, which is sequence-only and light) can be an optional backend. They should not be the default until an independent repression benchmark exists.

### Gaps
- I found no independent head-to-head of deep models (TargetNet, miTAR, Mimosa, TEC-miTarget, miRBind2, DeepMiRT) against scanMiR/McGeary Kd and TargetScan 8 on held-out **repression** data (AGO-CLIP + transfection + CRISPR knockouts). Mimosa, miTAR and TEC-miTarget were not specifically checked.
- I found no published application of RiNALMo to miRNA targeting.

## 4. Inferring miRNA activity from expression; is bulk gene–miRNA anti-correlation a valid repression readout?

### Takeaway
The field has moved from pairwise miRNA–mRNA correlation to **motif- or target-set enrichment over the transcriptome**, which needs only mRNA data:
- miReact (Sci Rep 2021)
- miTEA-HiRes (Communications Biology 2025, single-cell and spatial)
- bayesReact (NAR 2026; R package; uses miRBase v22 and MirGeneDB 3.0 seeds; better on sparse and single-cell data, while miReact is slightly better on large bulk)

These scores can be computed per cell in scRNA-seq, which suits CcSinglecell. Bulk anti-correlation in TCGA is weak in tumors and heavily confounded by tumor purity and composition. It should not be used as a repression readout without purity adjustment.

### Cited Findings
- **bayesReact:**
  - Paper: Rasmussen, Bouchard-Côté, Pedersen, NAR 2026 54(4):gkag072, DOI 10.1093/nar/gkag072, PMID 41657247; CC BY.
  - Method: unsupervised Bayesian model (Stan) of motif-position skew along 3′UTRs ranked by fold-change. Target sets are 7-mer seed motifs (positions 2–8) from miRBase v22 and MirGeneDB 3.0, not TargetScan predictions.
  - Benchmark vs miReact, Sylamer, cWords and miTEA-HiRes: "outperforms existing methods on sparse bulk data". miReact has higher expression–activity correlation on full bulk TCGA. bayesReact holds up better under downsampling and in single cells, and "the activity score does not exhibit a similar zero inflation" as expression.
  - The authors caution that "activity inference methods remain challenging to benchmark", and that activity "does not directly correspond to its expression level".
  - Example: miR-124-3p activity is strongly anti-correlated with REST in low-grade glioma.
  - Software: github.com/JakobSkouPedersenLab/bayesReact.
  - — [PMC12884093](https://pmc.ncbi.nlm.nih.gov/articles/PMC12884093/); [NAR](https://academic.oup.com/nar/article/54/4/gkag072/8468457)
- **miTEA-HiRes:** "Inferring single-cell and spatial microRNA activity from transcriptomics data", Communications Biology 2025. — [Commun Biol s42003-025-07454-9](https://www.nature.com/articles/s42003-025-07454-9)
- **Paired-expression methods** named in bayesReact's introduction: BIRTA, biRte, ActMiR. decoupleR is not discussed. — [PMC12884093](https://pmc.ncbi.nlm.nih.gov/articles/PMC12884093/)
- **decoupleR (Bioconductor 3.23, v2.17.0):** a generic framework (ulm, mlm, wsum, viper, etc.) for activity inference from any prior-knowledge network. It can take a signed, weighted miRNA→target network (e.g. weight = −context++), though that is not its documented default use. — [Bioconductor decoupleR](https://bioconductor.org/packages/release/bioc/html/decoupleR.html) (verified live)
- **Tumor-vs-normal correlation:** across 11 TCGA cancers, average miRNA–target correlation was −0.38 to −0.61 in normal tissue but only −0.04 to −0.26 in tumors. — [PMC6304515](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6304515/)
- **Purity confounding:** tumor purity confounds correlation and clustering across 21 TCGA cancer types and "should be accounted for" in all omics analyses (Aran et al., Nat Commun 2015). Cellular composition "has strong effects on the correlation analysis" of miRNA regulatory landscapes, and purity adjustment changed the results (F1000Research 10:633). DNA- and RNA-derived purity estimates agree poorly. — [Nat Commun ncomms9971](https://www.nature.com/articles/ncomms9971); [F1000Research](https://f1000research.com/articles/10-633); [JCO PO purity assessment](https://compbio.hms.harvard.edu/publications/pdf/po.20.00016.pdf)
- **miReact:** Nielsen MM & Pedersen JS, "miRNA activity inferred from single cell mRNA expression", Sci Rep 2021, [doi:10.1038/s41598-021-88480-5](https://doi.org/10.1038/s41598-021-88480-5) (prior knowledge). Downstream users treat its output as "a relative measure for the level of target repression" and select targets by correlation below −0.5 with activity, e.g. in early-stage bladder cancer. — [Nat Genet 2025, s41588-024-02030-z](https://www.nature.com/articles/s41588-024-02030-z)

### Inferences
- CcSinglecell's "miRNA activity" metric (gene–miRNA anti-correlation in matched TCGA samples) is the weakest-evidenced design among current options. Three confounders apply:
  1. Purity and stroma: many CRC miRNAs are stromal, such as miR-143/145 and miR-199/214.
  2. CMS subtype structure.
  3. Copy-number co-variation.
  
  At minimum, use partial correlation conditioned on purity (ESTIMATE/ABSOLUTE) plus CMS. Better, replace it with seed-enrichment activity (bayesReact/miReact) computed **per cell type in the scRNA-seq data**. That needs no matched miRNA-seq and works in the epithelial compartment directly.
- A clean decoupleR route: build a signed, weighted network miRNA → targets (weight = −context++ or −log Kd), run `ulm`/`mlm` per cell or pseudobulk, and treat the scores like TF activities. This puts the miRNA, TF and pathway layers in one framework (see Q6).
- isomiR-level analysis matters for activity because a 5′ isomiR shifts the seed. scanMiR and bayesReact (7-mer motifs) can handle custom seeds. TargetScan bulk files are canonical-family only.

### Gaps
- I found no study that directly quantifies how much TCGA miRNA–mRNA anti-correlation remains after purity adjustment in COAD/READ.
- CoMeTa, Sylamer (current status), ceRNA-network methods and isomiR-specific activity tools were not checked individually this session.
- I found no benchmark of decoupleR with miRNA target sets against bayesReact or miReact.

## 5. miRNAs (tissue; circulating/exosomal/cfRNA) as prognostic biomarkers in CRC, 2024–2026

### Takeaway
Evidence is strongest for **diagnostic** circulating miRNA panels (pooled AUC about 0.90 in 2025 meta-analyses), and even these mostly lack external validation. **Prognostic** miRNA evidence remains single-marker (miR-21, miR-20a, miR-675-5p), small and heterogeneous. I found no externally validated multi-miRNA prognostic panel in clinical use. The best-replicated CRC miRNA is **predictive**, not prognostic: low miR-31-3p predicts anti-EGFR benefit in RAS-wild-type mCRC, across several cohorts that share authors.

### Cited Findings
- **Cross-specimen diagnostic panel meta-analysis** (2025): pooled sensitivity 0.85, specificity 0.84, AUC 0.90. Most panels lack external validation, designs are mostly single-center, and 26 studies (76%) had problematic patient-selection applicability. — [PMC12578759](https://pmc.ncbi.nlm.nih.gov/articles/PMC12578759/)
- **Circulating miRNA diagnostic meta-analysis** (Frontiers in Genetics, Feb 2025): 37 studies, 2,775 patients. — [Frontiers 10.3389/fgene.2025.1574586](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2025.1574586/full)
- **Small 2025 diagnostic meta-analysis** (MDPI Applied Sciences): 7 studies. — [doi:10.3390/app15042111](https://doi.org/10.3390/app15042111)
- **miR-21 pooled analysis** (2017): high miR-21 is associated with shorter DFS. Tissue miR-21 is likely the better prognostic marker, circulating miR-21 the better diagnostic one. — [PMC5546529](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5546529/)
- **miR-20a meta-analysis:** pooled OS HR 2.02 with 95% CI 0.90–3.14, which crosses 1. — [PMC6721456](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6721456/)
- **miR-675-5p:** single-cohort study; high expression is independently associated with shorter DFS and OS after TNM adjustment. — [PMC10298325](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10298325/)
- **miR-31-3p (predictive for anti-EGFR):**
  - Manceau 2014 discovery: validation HR 4.1 (frozen) and 2.44 (FFPE).
  - PROSPECT-C 2019 (ISH, 45 patients): low miR-31-3p means better ORR, PFS and OS.
  - New EPOC post hoc: patients with low miR-31-3p were not harmed by adding cetuximab.
  - Cohorts are small and the groups overlap.
  - — [ASCO Post 2014](https://ascopost.com/issues/september-15-2014/potential-marker-for-response-to-anti-egfr-therapy-in-kras-wild-type-metastatic-colorectal-cancer/); [ICR repository (PROSPECT-C)](https://repository.icr.ac.uk/handle/internal/3193); [Oncotarget (New EPOC)](https://www.oncotarget.com/article/21291/text/)
- **Reviews:** a 2025 BMC Cancer mini-systematic review covers studies through Nov 2023 but does not pool signature validation. — [BMC Cancer s12885-025-14343-1](https://link.springer.com/article/10.1186/s12885-025-14343-1); [updated review PMC12429582](https://pmc.ncbi.nlm.nih.gov/articles/PMC12429582/)
- **Methylation contrast:** a methylation-marker external-validation study in DACHS notes that prognostic markers "are often not translated", hampered by a lack of validation in large independent cohorts. This is the same failure mode as miRNA signatures. — [medRxiv 2022.11.03.22281595](https://www.medrxiv.org/content/10.1101/2022.11.03.22281595.full.pdf)

### Inferences
- CcSinglecell's miRNA layer selects miRNAs from dbDEMC (differential expression, largely diagnostic evidence) and uses them to weight genes for **survival** prediction. The literature does not support assuming that DE-listed cancer miRNAs carry prognostic information. The successor should (a) select miRNAs by CRC/colon activity or expression in the actual data, not dbDEMC lists, and (b) validate any miRNA-derived feature on held-out cohorts against a TNM + CMS baseline.
- Circulating/exosomal miRNAs are a separate use case (detection, not survival). They are probably out of scope unless the successor adds liquid-biopsy data.

### Gaps
- I found no 2024–2026 meta-analysis specifically of CRC **prognostic** multi-miRNA signatures with PROBAST/REMARK assessment.
- I found no data on whether any miRNA prognostic assay (e.g. a commercial panel) is in CRC guidelines. My belief is that none is in NCCN/ESMO, but I did not verify it.
- Exosome- or cfRNA-specific prognostic evidence for 2024–2026 was not found in these searches.

## 6. Beyond miRNA: are lncRNA/circRNA/RBP/TF regulatory layers better-evidenced for prognosis?

### Takeaway
For a modern regulatory layer, **TF activity (CollecTRI/DoRothEA via decoupleR) and pathway activity (PROGENy)** are the most mature, benchmarked and maintained options. They work in bulk and single-cell data, use footprint scoring like miRNA activity, and share one API (decoupleR). SCENIC+ adds cell-type-resolved enhancer-GRNs but needs multiome data. lncRNA and circRNA prognostic signatures share or exceed miRNA's reproducibility problems. I found no 2024–2026 evidence that they are better prognostic signals.

### Cited Findings
- **decoupleR:** Bioconductor 3.23 v2.17.0, Saez-Rodriguez lab (Heidelberg). Explicitly built to infer activities from any omics with any prior-knowledge resource, including signed and weighted networks. — [Bioconductor decoupleR](https://bioconductor.org/packages/release/bioc/html/decoupleR.html) (verified live). Paper: Badia-i-Mompel et al., Bioinformatics Advances 2022, [doi:10.1093/bioadv/vbac016](https://doi.org/10.1093/bioadv/vbac016) (prior knowledge).
- **CollecTRI:** a curated signed TF–gene regulon collection that outperformed DoRothEA in TF-activity benchmarks on perturbation data. Müller-Dott et al., NAR 2023, [doi:10.1093/nar/gkad841](https://doi.org/10.1093/nar/gkad841) (prior knowledge).
- **PROGENy:** perturbation-response footprints for 14 pathways. Schubert et al., Nat Commun 2018, [doi:10.1038/s41467-017-02391-6](https://doi.org/10.1038/s41467-017-02391-6) (prior knowledge).
- **SCENIC+:** enhancer-driven GRNs from single-cell multiome. Bravo González-Blas et al., Nat Methods 2023, [doi:10.1038/s41592-023-01938-4](https://doi.org/10.1038/s41592-023-01938-4) (prior knowledge).
- **ENCORI** also hosts RBP–RNA (CLIP) and ceRNA interaction tables, so it can supply an RBP layer from the same resource. — [ENCORI](https://rnasysu.com/encori/) (verified live)

### Inferences
- A good successor design is a single "regulator activity" module on decoupleR. Its networks would be CollecTRI (TF), PROGENy (pathway), and miRNA → target from TargetScan 8 / scanMiR / TarBase-CLIP, weighted and signed negative. Run it per cell type on the scRNA-seq data and on TCGA/GEO bulk. Feed the activity scores, not raw target counts, into the survival model, and report incremental C-index over TNM + CMS. This replaces both hoardeR and the row-sum target-count ranking. Row-sum counts mostly reward genes with long, site-rich 3′UTRs that are targeted by many miRNAs, regardless of whether those miRNAs are active.
- lncRNA/circRNA layers are poorly suited to a scRNA-seq-first package. Coverage is low in 3′ 10x data (circRNAs need total-RNA/back-splice reads), and prognostic signatures from them are typically TCGA-only. Treat them as optional at most.

### Gaps
- I did not search this session for a 2024–2026 head-to-head study in CRC comparing the prognostic value of TF/pathway activity with miRNA activity or lncRNA signatures. The claim that TF/pathway layers are "better-evidenced for prognosis" rests on tool maturity and benchmarking, not on a CRC survival comparison.
- CMS/CRIS subtype prognostic baselines were not re-sourced here.
- I did not check whether decoupleR or OmniPath currently ship a ready-made miRNA–target resource (OmniPath has historically included miRNA–mRNA interactions).
