# Multi-omics and multimodal integration models for cancer prognosis (2024–2026 frontier): what adds value beyond clinical baselines

Scope note: research run 2026-11-07. Items marked **[older]** are from before 2024. Items marked **[snippet-only]** were seen only in search-result excerpts, without reading the full text, so their numbers should be re-checked before anyone quotes them. "Demonstrated" means the result came from a held-out or external cohort with a clinical or simple baseline. "Claimed" means a self-reported gain over weak or untuned baselines, or a gain measured only in internal cross-validation.

## 1. Benchmarks of multi-omics survival prediction against clinical baselines (Herrmann 2021 and follow-ups)

### Takeaway
Every large benchmark finds the same thing. Clinical variables plus (at most) gene expression carry almost all the signal. Group-structured statistical models (block forest, clinical-favouring CoxBoost or priority-Lasso) match or barely beat a clinical-only Cox model, adding more omics layers usually hurts, and deep models do not beat statistical ones (they calibrate worse). The CcSinglecell result of ΔC ≈ 0 for mRNA+miRNA over age+stage is the expected outcome, not a failure of that particular pipeline.

### Cited Findings
- **[older] Herrmann et al. 2021, Brief Bioinform 22(3):bbaa167 (DOI 10.1093/bib/bbaa167).**
  - Design: 11 methods (penalised regression, boosting, random forests) on 18 TCGA multi-omics datasets, repeated 5-fold CV, with Kaplan–Meier and a clinical-only Cox model as references.
  - Clinical-only Cox: mean Uno C = 0.618, integrated Brier score (IBS) 0.175.
  - Block forest was the only method to beat it on average, by a tiny margin: C = 0.620, IBS 0.174. CoxBoost favouring clinical features got C = 0.618, IBS 0.174. IPF-Lasso got IBS 0.176 but C = 0.578. Priority-Lasso variants got IBS 0.180 and 0.181.
  - Paired t-tests found no significant difference between block forest, clinical-favouring CoxBoost and the clinical-only Cox model.
  - In 6 of 18 datasets no omics method beat clinical Cox. Multi-omics "helped substantially" for some datasets, including LAML and COAD.
  - Respecting the group structure and favouring clinical variables both helped on average. The authors judge the predictive utility of multi-omics to be "limited".
  - Source: [Herrmann et al. 2021](https://academic.oup.com/bib/article/22/3/bbaa167/5895463)
- **Li et al. 2024, BMC Med Inform Decis Mak (DOI 10.1186/s12911-024-02642-9), from the same research line.** They tested all 31 combinations of 5 omics types across 14 TCGA datasets, with clinical data favoured in every model. Adding more omics types usually made performance worse, so the practice of including as many data types as possible deserves "a second look". — [Springer](https://link.springer.com/article/10.1186/s12911-024-02642-9); [PubMed](https://pubmed.ncbi.nlm.nih.gov/39223659/) **[snippet-only; full text blocked]**
- **SurvBoard (Wissel/Boeva lab), Brief Bioinform 2025, 26:bbaf521 (DOI 10.1093/bib/bbaf521).**
  - Design: 28 datasets (TCGA 21, ICGC, TARGET, METABRIC) and 12 models: PriorityLasso (L1+L2), BlockForest and 10 deep-learning models.
  - "Regularized linear models and random-forest-based methods outperformed deep learning methods" when trained on clinical plus gene expression.
  - Clinical and gene expression were the most predictive modalities, and many models got worse as more modalities were added.
  - BlockForest had the best IBS rank and PriorityLasso was second. The deep model with the best Antolini-C rank was among the worst calibrated (IBS and D-calibration), and C-index differences between methods were small.
  - Pan-cancer training (clinical + expression) improved C and IBS for most methods, deep models most and BlockForest least.
  - Adding samples with missing modalities as extra training data improved C for the deep models. Only PriorityLasso gained in IBS and D-calibration.
  - The authors say multi-omics survival analysis "currently has limited clinical applications".
  - The main text gives no separate clinical-only baseline; the numbers are only in the figures and supplements.
  - Sources: [SurvBoard PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12486238/); code and web service at [GitHub BoevaLab/survboard](https://github.com/BoevaLab/survboard)
- **[older] Yuan et al. 2014, Nat Biotechnol 32:644 (DOI 10.1038/nbt.2940).** Covered 953 TCGA samples from 4 cancer types. Adding molecular data to clinical variables (age, stage, grade, gender) gave statistically significant but modest gains in only 3 cancers, and molecular data alone showed "little predictive power" with one exception. — [eScholarship](https://escholarship.org/uc/item/0th055mz)
- **Genetic-algorithm feature-selection study (arXiv 2604.00065, 2026).** Used clinical + mRNA + miRNA in TCGA-KIRC, LGG and SARC (not COAD/READ). In KIRC, multi-omic panels looked better internally but did "not exceed the clinical baseline" externally. LGG was the only cohort with a clear gain from clinical + mRNA. Overestimation (internal vs nested-external) was a recurring problem. — [arXiv 2604.00065](https://arxiv.org/html/2604.00065)
- **TCGA colorectal cohort (COADREAD/"CRC") in the MMP benchmark (ICML 2024, 5-fold site-stratified CV, disease-specific survival).**
  - Clinical-only C = 0.655 ± 0.119.
  - Gene expression MLP 0.588 ± 0.105; pathway SNN 0.590 ± 0.104.
  - Multimodal (histology + genomics): MCAT 0.578 ± 0.136, SurvPath 0.539 ± 0.150, MOTCat 0.590 ± 0.130, CMTA 0.559 ± 0.195, MMP 0.630–0.636 ± 0.12.
  - So in TCGA CRC the clinical baseline beat every genomic and multimodal model, and fold SDs were about 0.10–0.20.
  - Source: [MMP, Song et al., arXiv 2407.00224](https://arxiv.org/html/2407.00224v1)

### Inferences
- The CcSinglecell age+stage C ≈ 0.69 sits at or above the clinical baseline published for TCGA COADREAD (0.655 in MMP). With 78 events and fold SDs of about 0.1, a molecular increment below about 0.03–0.05 C is statistically undetectable in TCGA COAD/READ alone. A successor should treat "ΔC over clinical, with a CI" as the primary endpoint and power it with external or pooled cohorts.
- The most defensible ways to "integrate" for prediction are clinical-favouring structured penalisation (priority-Lasso with an unpenalised clinical block, CoxBoost with mandatory clinical covariates, IPF-Lasso) or block forests, judged on IBS and calibration as well as C. Blending per-gene rankings across omics before fitting the model has no benchmark support in this literature.
- Pan-cancer pretraining or transfer, and training on samples with missing modalities, are the two "tricks" with benchmark-level support (SurvBoard) for raising performance.

### Gaps
- I found no 2024–2026 benchmark that reports ΔC over a clinical-only model specifically for TCGA COAD/READ with mRNA+miRNA. SurvBoard's per-dataset numbers are in its supplements and Zenodo outputs, which I did not read.
- Li et al. 2024 effect sizes (C and IBS per combination) could not be extracted because the full text was blocked.

## 2. Classical and statistical integration (MOFA+/MEFISTO, DIABLO, SNF, iCluster, intNMF, priority-Lasso, IPF-Lasso, block forests, cooperative learning, late fusion)

### Takeaway
For supervised survival prediction, block forests and priority-Lasso (with clinical favoured) are the methods with benchmark support. Unsupervised factor and clustering methods (MOFA+, iCluster, SNF, intNMF) are for exploration and subtyping, not prediction, and have not beaten PCA in head-to-head tests. DIABLO is a classification method with no survival-benchmark support. Cooperative learning is principled but has only thin survival evidence, none of it in cancer benchmarks.

### Cited Findings
- **[older] Herrmann 2021, method by method.**
  - Block forest: best average C (0.620) and IBS (0.174).
  - IPF-Lasso: best IBS among penalised methods that do not favour clinical variables (0.176), but low C (0.578).
  - Priority-Lasso: worse on average (IBS 0.180–0.181).
  - Group-structure-aware methods beat naive concatenation, but the authors advise against recommending a single method because rankings depend on dataset and metric.
  - Source: [Herrmann 2021](https://academic.oup.com/bib/article/22/3/bbaa167/5895463)
- **SurvBoard 2025.** BlockForest was best and PriorityLasso second on calibration-type metrics (IBS, D-calibration) across 28 datasets. Both beat all 10 deep models. — [SurvBoard](https://pmc.ncbi.nlm.nih.gov/articles/PMC12486238/)
- **[older] CCLE benchmark 2022 (MOFA2, DIABLO, iClusterPlus, iClusterBayes, moCluster).** "No single best method" emerged, the iCluster family and MOFA2 "did not perform better than PCA", and DIABLO was among the fastest. — [Frontiers Genet 2022](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2022.854752/full)
- **Single-cell multimodal integration benchmark, Nat Methods 2025.** MOFA+ had the highest reproducibility of feature selection. This is single-cell data, not bulk survival. — [Nat Methods 2025](https://www.nature.com/articles/s41592-025-02856-3)
- **STAR Protocols 2026, practical guidance.** When missingness exceeds about 30% across layers, use MOFA rather than forcing DIABLO. — [STAR Protocols](https://www.cell.com/star-protocols/fulltext/S2666-1667(26)00187-5)
- **[older] Cooperative learning (Ding, Li, Narasimhan, Tibshirani, PNAS 2022).** An agreement penalty between view-specific predictions gives a continuum from early to late fusion and pairs with the lasso for sparsity. — [arXiv 2112.12337](https://arxiv.org/pdf/2112.12337)
- **Survival extensions of cooperative learning.**
  - Hahn et al. 2024, Brief Bioinform: cooperative learning with Cox polygenic hazard models, AUC about 0.7. This is Alzheimer's disease-free survival, not cancer. — [PMC11151121](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11151121/)
  - A 2026 conference abstract proposes a cooperative accelerated-failure-time model ("CoopAFT"). — [SMOD 2026 abstract](https://www.dipmat2.unisa.it/iniziative/smod2026/Abstracts/Iuliano.pdf)
- **[older] Integration that "blends rankings" before modelling** (the CcSinglecell approach) is not among the method families any benchmark evaluated. Every benchmarked integrator either uses group-wise penalties or splits (priority-Lasso, IPF-Lasso, block forest) or fuses learned representations. — [Herrmann 2021](https://academic.oup.com/bib/article/22/3/bbaa167/5895463); [SurvBoard](https://pmc.ncbi.nlm.nih.gov/articles/PMC12486238/)

### Inferences
- Recommended stack for a successor's statistical arm:
  1. Clinical-only Cox as the reference.
  2. Priority-Lasso (prioritylasso R package) or CoxBoost with the clinical block unpenalised.
  3. IPF-Lasso.
  4. Block forest (blockForest R package).
  5. Optionally, cooperative-learning Cox (multiview) as a principled early/late-fusion dial.
  Judge each on C, IBS and D-calibration under nested CV.
- Keep MOFA+ and DIABLO only as exploratory or feature-generation layers (factors used as candidate covariates). They should not be the main route to prognostic signatures.

### Gaps
- I found no 2024–2026 bulk-cancer survival benchmark that pits SNF, iCluster, intNMF or MEFISTO against clinical baselines. The most recent head-to-heads I found are from 2022 or are single-cell only.
- I found no published head-to-head of cooperative learning against priority-Lasso or block forest on TCGA survival.

## 3. Deep learning for omics survival (autoencoders, GNNs, biologically informed networks, transformers, DeepSurv/Cox-Time/DeepHit heads)

### Takeaway
On bulk multi-omics survival, deep models do not beat regularised Cox or random forests, and they calibrate worse (SurvBoard 2025). Biologically informed networks' gains come from sparsity, not from the pathway wiring: randomly rewired networks of equal sparsity do as well or better. The case for P-NET-style models now rests on interpretability, which is itself contested.

### Cited Findings
- **SurvBoard 2025.** Six state-of-the-art deep multi-omics survival methods, plus other neural models, were compared with PriorityLasso and BlockForest. "Statistical models tend to outperform deep learning methods, especially for metrics measuring survival function calibration." — [SurvBoard](https://pmc.ncbi.nlm.nih.gov/articles/PMC12486238/)
- **"Sparsity is all you need" (Brief Bioinform 2026, bbag425; arXiv 2505.04300).**
  - Tested five pathway-informed models: BINN, DeepKEGG, PASNet, PathCNN and PINNet.
  - Randomly rewired networks with the same sparsity performed as well. Statistically significant differences *favouring random sparsity* appeared for BINN and DeepKEGG.
  - Randomised P-NET reached AUC 0.896 against 0.899 for the original.
  - Conclusion: the predictive gains attributed to pathways "do not require biologically meaningful connectivity".
  - Sources: [Brief Bioinform 2026](https://academic.oup.com/bib/article/27/4/bbag425/8753615); [arXiv](https://arxiv.org/html/2505.04300v1)
- **Simulation and empirical evaluation of BiNNs (2025).** Small samples, weak signal and extreme feature sparsity limit BiNN performance. P-NET beat random forest in less sparse settings but was slightly worse than logistic regression. — [PMC12642320](https://pmc.ncbi.nlm.nih.gov/articles/PMC12642320/)
- **P-NET expansion preprint (bioRxiv 2026).** P-NET importance concentrates on meaningful pathways while random architectures spread it evenly. However, the importance stability of sparse networks is "still far off from traditional models". — [bioRxiv 2026.04.19.719454](https://www.biorxiv.org/content/10.64898/2026.04.19.719454v1.full)
- **Genomics-only deep baselines in multimodal papers are often weak or untuned.** The MIST paper reports its genomics-only baseline near chance (C ≈ 0.50) in every cohort. — [MIST arXiv 2609.21811](https://arxiv.org/html/2609.21811v1)
- **SurvPath paper (CVPR 2024).** "A simple MLP trained on our set of transcriptomics constitutes a strong baseline that outperforms several multimodal methods." — [arXiv 2304.06819](https://arxiv.org/html/2304.06819v2)

### Inferences
- A successor gains little in prediction from a DeepSurv, Cox-Time or DeepHit head on bulk mRNA/miRNA at n ≈ 350 and 78 events. If a neural arm is included, it should be benchmarked against elastic-net Cox and block forest with a matched tuning budget, and on IBS and calibration, not C alone.
- If pathway structure is used, include a "random sparse wiring" control, as the 2026 Brief Bioinform paper did. Otherwise biological claims are unfalsifiable.

### Gaps
- I found no 2024–2026 head-to-head of DeepSurv, Cox-Time or DeepHit against elastic-net Cox on high-dimensional cancer omics beyond SurvBoard. I found no rigorous evaluation of GNNs over PPI graphs or omics transformers with clinical baselines.
- I did not verify per-model deep-learning results inside SurvBoard (they are in the figures).

## 4. Multimodal pathology + omics, pathology foundation models, and H&E biomarkers in CRC

### Takeaway
- **Pan-cancer multimodal (WSI+genomics) models.** Papers report average C-index gains of about 0.02–0.05 over unimodal models, but always in internal CV with fold SDs of about 0.05–0.15. In TCGA CRC, a clinical-only model beat every multimodal model, and these models are badly calibrated.
- **Pathology foundation models.** They give survival C of roughly 0.66–0.77 internally, with small, non-significant differences between models.
- **H&E in CRC.** The demonstrated value is H&E-derived prognostic scores such as DoMore-v1/CAPAI. They are externally validated and *independent of stage* (adjusted HR about 3), and they reclassify stage II/III patients. The evidence is retrospective; prospective trials are still lacking.

### Cited Findings
**Multimodal survival models (MCAT, PORPOISE, SurvPath, MMP and successors)**
- **MMP (ICML 2024), using a UNI encoder and 5-fold site-stratified CV.**
  - Average C over 6 TCGA cohorts: clinical 0.585, gene expression 0.612, best unimodal histology about 0.605–0.611, SurvPath 0.629, MOTCat 0.631, MMP 0.665.
  - In CRC, "unimodal histology baselines outperform multimodal baselines": AttnMISL got 0.728 ± 0.110 against MMP 0.636 ± 0.120, and clinical got 0.655.
  - Source: [arXiv 2407.00224](https://arxiv.org/html/2407.00224v1)
- **Structural prognostic event modelling (arXiv 2512.01116, 2025).** PORPOISE, PIBD and MMP did not beat unimodal models in their benchmark. — [arXiv 2512.01116](https://arxiv.org/html/2512.01116)
- **PS3 (ICCV 2025).** Claims a 6.72% improvement over the next-best multimodal method. Gene-based methods were the strongest unimodal baselines. — [arXiv 2509.20022](https://arxiv.org/html/2509.20022) **[snippet-only for numbers]**
- **Calibration audit (arXiv 2604.04239, Apr 2026).**
  - On TCGA-BRCA: SurvPath C 0.669 ± 0.063, MCAT 0.645 ± 0.057, MMP 0.598 ± 0.081.
  - 1-calibration failed in 4/5, 5/5 and 3/5 folds respectively. Across 290 tests, 166 rejected calibration. MCAT reached C = 0.817 on GBMLGG but failed calibration on all five folds.
  - A genomics-only SNN failed in only 1/25 folds, and Cox-PH passed.
  - Gating fusion calibrated better than bilinear or concatenation fusion. Platt scaling fixed most of the miscalibration without changing C.
  - The SurvPath–MCAT gap (0.024) is well inside the fold SD (about 0.06).
  - Source: [arXiv 2604.04239](https://arxiv.org/html/2604.04239)

**Pathology foundation models (UNI2, CONCH 1.5, Virchow2, H-optimus-1, Prov-GigaPath, etc.)**
- **PathBench (arXiv 2505.20202).** 19 models; survival tasks used 3×5-fold CV in single-institution cohorts.
  - CRC overall survival best: CTransPath 0.757. CRC DFS best: H-optimus-1 0.738. CRC disease-specific survival best: Virchow2 0.768.
  - No model led all CRC survival tasks. Breast CIs were about 0.35 wide, so rankings are not separable.
  - I found no external survival cohort in the text.
  - Source: [PathBench](https://arxiv.org/html/2505.20202)
- **Breast survival FM benchmark (arXiv 2604.24679, 2026).** 13 models. H-optimus-1 was best overall, but "absolute differences between top-performing models are relatively small". Second-generation models beat first-generation ones (Virchow2 over Virchow, UNI2-h over UNI). — [arXiv 2604.24679](https://arxiv.org/html/2604.24679v1)
- **Nat Biomed Eng 2025 weakly supervised benchmark.** CONCH had the highest average AUROC, followed by Virchow2 and Prov-GigaPath. These are mostly classification and biomarker tasks, not survival. — [Neidlinger et al.](https://www.nature.com/articles/s41551-025-01516-3)
- **31-model medRxiv benchmark (2025).** Virchow2 performed best across TCGA, CPTAC and external tasks. — [medRxiv 2025.05.08.25327250](https://www.medrxiv.org/content/10.1101/2025.05.08.25327250v1.full)

**CRC H&E prognostic biomarkers**
- **[older] DoMore-v1-CRC, Skrede et al. Lancet 2020 (DOI 10.1016/S0140-6736(19)32998-8).**
  - External validation in 1,122 stage II/III patients from the QUASAR 2 trial (capecitabine). Poor vs good prognosis: HR 3.84 (95% CI 2.72–5.43) unadjusted, and HR 3.04 (2.07–4.47) after adjusting for established markers including pN stage.
  - Endpoint: cancer-specific survival.
  - Source: [PubMed 32007170](https://pubmed.ncbi.nlm.nih.gov/32007170/)
- **[older] Kleppe et al. Lancet Oncol 2022, 23(9):1221–1232.**
  - Combined DoMore-v1 with pT, pN and lymph-node count into a decision-support system. External validation in 1,075 QUASAR 2 patients.
  - Risk groups: 41% low, 32% intermediate, 27% high, against 13%, 57% and 30% under guidelines. The new low-risk group had prognosis similar to the guideline low-risk group.
  - DoMore-v1 gave prognostic information beyond guideline risk.
  - Sources: [PubMed 35964620](https://pubmed.ncbi.nlm.nih.gov/35964620/); [ASCO Post 2022](https://ascopost.com/news/september-2022/new-risk-stratification-system-for-completely-resected-stage-ii-or-iii-colorectal-cancer/)
- **CAPAI Dutch nationwide validation (ESMO 2025, abstract 726O, interim).**
  - 453 untreated high-risk stage II/III colon cancer patients from 65 hospitals.
  - 3-year cancer-specific survival: low risk 93.7%, intermediate 87.5%, high 60.4%. HR high vs low 7.9 (4.1–15.2); intermediate vs low 2.1 (1.0–4.3).
  - No head-to-head with stage alone was reported, and it is a congress abstract with no DOI.
  - Source: [ASCO Post Dec 2025](https://ascopost.com/issues/december-10-2025/validation-of-pathology-ai-integration-tool-for-risk-stratification-in-locally-advanced-colon-cancer-without-chemotherapy)
- **Jiang et al. Lancet Digit Health (online 2023, issue 2024), end-to-end attention-MIL.**
  - 4,428 patients from Australia, Germany and the USA.
  - Disease-specific survival HR in external tests: 3.08 (2.44–3.89; n = 1,395), 2.23 (1.23–4.04) and 3.07 (1.78–5.3). Reported as "independent of established clinical risk factors".
  - Sources: [White Rose eprint](https://eprints.whiterose.ac.uk/203369/10/End-to-end%20prognostication%20in%20colorectal%20cancer%20by%20deep.pdf); [Lancet DH](https://www.thelancet.com/journals/landig/article/PIIS2589-7500(23)00208-X/fulltext) **[abstract-level]**
- **HistoNav (Cancers 2026, DOI 10.3390/cancers18172863), a DINO-ViT plus graph network.**
  - Trained on QUASAR 2: CV C-index only 0.592 ± 0.039.
  - External HRs: TCGA-COAD 5.46; SCOT 2.35 (95% CI 0.83–6.61, p = 0.148, *not significant*).
  - No head-to-head with TNM was reported, and no treatment-interaction analysis.
  - Source: [PMC13564889](https://pmc.ncbi.nlm.nih.gov/articles/PMC13564889/)
- **Systematic review of DL H&E prognosis in CRC (Research Square preprint, Nov 2026).** Pooled high- vs low-risk HR was 2.36 in 17 externally validated studies against 2.74 in 10 internally validated studies, so effects shrink on external validation. — [Research Square rs-11297991](https://www.researchsquare.com/article/rs-11297991/v1) **[snippet-only]**
- **Multiview DL for stage II CRC (2025).** Two external Chinese cohorts (n ≈ 438 and 405). — [PMC12801286](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12801286/) **[snippet-only; no effect sizes extracted]**
- **MSI-H from H&E.** A systematic review and meta-analysis exists. — [PMC12274608](https://pmc.ncbi.nlm.nih.gov/articles/PMC12274608/) **[not read; no numbers extracted]**

### Inferences
- In CRC, H&E carries the strongest *validated* non-staging prognostic signal among the modalities in this note: an adjusted HR of about 3 against pN stage in a trial cohort. That is stronger than any published mRNA/miRNA increment.
- Adding H&E to TCGA COAD/READ (about 600 slides) is feasible with frozen FM features (UNI2, Virchow2, H-optimus-1) and ABMIL. Expect C in the 0.6–0.7 range with SD of about 0.1 in TCGA. An offset model (clinical linear predictor plus H&E score) and an external cohort are needed to show any ΔC.
- Fancy cross-attention fusion (MCAT, SurvPath) has no demonstrated advantage in CRC. Late fusion with recalibration (Platt) is the safer default.

### Gaps
- I found no study reporting ΔC (stage + H&E against stage alone) with CIs on an external CRC cohort. Most CRC studies report HRs of dichotomised risk groups.
- I found no prospective randomised trial of H&E-guided adjuvant decisions in CRC.
- TITAN, THREADS and other 2025–2026 slide-level FMs were not covered for survival in the sources I read.

## 5. Omics and bulk-transcriptome foundation models; LLM and agentic systems for biomarker discovery

### Takeaway
Bulk-transcriptome foundation models (BulkFormer and others) report prognosis gains only against other foundation models. I found no comparison against raw-expression elastic-net Cox, and the one independent benchmark has internally inconsistent numbers. Agentic "AI scientist" systems (Biomni, STELLA, etc.) score about 40–65% on general bio-analysis benchmarks and are judged "more reliable as assistants than fully autonomous agents". There is no rigorous biomarker-discovery evaluation.

### Cited Findings
- **BulkFormer (bioRxiv 2025.06.11.659222).** About 150M parameters covering about 20k protein-coding genes, pretrained on more than 500k bulk RNA-seq profiles. Prognosis is one of six downstream tasks, and the paper claims it beats single-cell LLMs on all of them. — [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.06.11.659222.full.pdf) **[snippet-only; full text rate-limited]**
- **Independent benchmark (Kim et al., LG AI Research and Vanderbilt, AACR 2026).**
  - Linear probing on frozen embeddings of 6 single-cell FMs plus BulkFormer, survival across 14 TCGA cancers.
  - Reported C: scFoundation 0.672 ± 0.081, CellFM 0.665 ± 0.078, BulkFormer 0.839 ± 0.086, scBERT 0.599, scLong 0.589.
  - Yet the abstract calls scFoundation and CellFM "comparable" to BulkFormer, and its scGPT text contradicts its own table.
  - No raw-expression baseline was reported.
  - Sources: [AACR 2026 abstract mirror](https://aacr2026-wiki.vercel.app/references/5261); [LG AI Research](https://www.lgresearch.ai/publication/view?seq=157)
- **Downstream use.** A preprint built a deep Cox model on the 37M BulkFormer variant and got mean C 0.686 ± 0.002. — [search result, arXiv 2605.16571](https://arxiv.org/pdf/2605.16571) **[snippet-only; context unclear]**
- **BAISBench (Luo et al., Bioinformatics 2026, DOI 10.1093/bioinformatics/btag227).**
  - Tasks: scRNA-seq only. Cell annotation on 15 datasets, and 193 discovery multiple-choice questions from 41 studies.
  - Ranking: STELLA best (about human-expert level on discovery), then Pantheon, then Biomni. Biomni with Claude Opus 4.5 approached human level.
  - Failure modes: brittle execution, error accumulation and answering from prior knowledge. Performance dropped on shuffled data.
  - Conclusion: AI scientists are "more reliable as assistants than fully autonomous agents".
  - Source: [Bioinformatics 2026](https://academic.oup.com/bioinformatics/article/42/Supplement_1/btag227/8726318)
- **HeurekaBench (arXiv 2601.01678).** On single-cell tasks, Biomni scored about 50%, CellVoyager 27.8% and BixBench-Agent 44.4%. — [arXiv 2601.01678](https://arxiv.org/pdf/2601.01678) **[snippet-only]**
- **BiomniBench (bioRxiv May 2026).** A process-level evaluation of agent analyses. A secondary write-up says agents fail consistently at *method selection*. — [bioRxiv](https://www.biorxiv.org/content/10.64898/2026.05.12.724604v1.full.pdf); [SEK write-up](https://www.sekbio.com/blog-9-bioomnibench-ai-evaluation.html) **[snippet-only]**
- **Vendor-adjacent claims.** Biomni Lab scored 52.2% on BixBench against 39.5% for Claude Code (Opus 4.6), and 65.0 against 68.5 for senior pharma scientists on a rubric. These come from a developer-connected source. — [Phylo blog](https://phylo.bio/blog/evaluating-ai-agents-in-biology)
- **Biomni original preprint (2025).** — [ResearchGate](https://www.researchgate.net/publication/392392722_Biomni_A_General-Purpose_Biomedical_AI_Agent)

### Inferences
- Treat the BulkFormer survival claims as unproven for this use until it beats elastic-net Cox on raw expression *with clinical covariates* in held-out data. A successor could test this cheaply by using BulkFormer embeddings as one more block in priority-Lasso or block forest.
- Agentic systems fit as *workflow assistants* (literature triage, pipeline scaffolding). They should not decide biomarkers. Their documented weakness, method selection, is exactly the failure that produces optimistic survival signatures.

### Gaps
- I found no rigorous evaluation of any LLM agent on survival-biomarker discovery with held-out validation.
- I found no bulk FM evaluated with clinical baselines in CRC. I also found no multi-omics (beyond RNA) foundation model with a survival benchmark.

## 6. Proteomics, proteogenomics, methylation and cfDNA as layers with stronger prognostic evidence

### Takeaway
In CRC, the layer with by far the strongest validated prognostic evidence is post-operative ctDNA molecular residual disease (MRD): multivariable DFS HR of about 12, dominating stage, nodal status and RAS/BRAF. Tissue or plasma proteomics shows promise only in small single-cohort studies. I found no evidence that CPTAC proteomics beats mRNA for CRC survival in matched patients.

### Cited Findings
- **GALAXY / CIRCULATE-Japan, Nat Med 2024 (DOI 10.1038/s41591-024-03254-6).**
  - Cohort: 2,240 patients with stage II–IV CRC, median follow-up 23 months.
  - ctDNA positivity in the MRD window: DFS HR 11.99 (10.02–14.35), OS HR 9.68.
  - In multivariable analysis, ctDNA was "the single most significant prognostic factor": DFS HR 12.08 against lymph-node positivity 1.56, BRAF V600E 2.03, RAS 1.43 and pT 1.67.
  - Adjuvant chemotherapy benefit was confined to MRD-positive patients (adjusted HR 0.23) and was not significant in MRD-negative patients (HR 0.70, p = 0.091).
  - Source: [PMC11564113](https://pmc.ncbi.nlm.nih.gov/articles/PMC11564113/)
- **Olink Explore plasma proteomics (West China Hospital, Carcinogenesis 2026).** 150 CRC patients and 367 neuro-related proteins. LASSO-Cox gave 5-protein OS and 11-protein DFS signatures, with the OS signature validated in UK Biobank (n = 1,133). — [Olink publication page](https://olink.com/publication/plasma-proteomic-profiling-identified-prognostic-indicators-with-therapeutic-potential-for-colorectal-cancer)
- **UK Biobank pre-diagnosis plasma proteomics and colon cancer survival (BMC Cancer 2025).** Samples were collected about 7.9 years before diagnosis. A 10-protein profile plus age and sex gave AUCs of 0.85, 0.82 and 0.89 at 1, 5 and 10 years for early stage, and 0.71–0.79 for late stage. Incremental value over age and sex alone was not stated in the excerpt. — [UK Biobank publications](https://www.ukbiobank.ac.uk/publications/association-of-pre-diagnosis-plasma-proteomic-contexture-with-overall-survival-of-early-and-late-stage-colon-cancer-patients/) **[snippet-only]**
- **[older] Olink Target 96 (2018).** In 261 stage II–IV patients, only osteoprotegerin out of 92 proteins was associated with OS, and none with recurrence. — [Olink](https://olink.com/publication/plasma-protein-profiling-reveal-osteoprotegerin-as-a-marker-of-prognostic-impact-for-colorectal-cancer)
- **[older] TCPA reverse-phase protein array (RPPA) protein signature (2021).** A 4-protein LASSO-Cox model in 315 colon cancer patients had the best 1-year time-dependent ROC among the indicators compared. Not externally validated. — [PMC7974900](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7974900/)
- **[older] CPTAC proteogenomic recurrence signature (Front Oncol 2022).** Only 34 CPTAC COAD protein samples, 6 of them recurrent, combined with TCGA RNA-seq. The design does not isolate protein against RNA. — [Frontiers 2022](https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2022.871568/full)

### Inferences
- If the successor wants a layer that demonstrably moves prognosis beyond stage in CRC, ctDNA MRD (not public in TCGA) and H&E (public in TCGA) are the evidence-backed choices. Proteomics is a hypothesis-generating layer; CPTAC COAD is too small for survival modelling on its own.
- UK Biobank Pharma Proteomics Project data address pre-diagnostic *risk* and survival, not post-resection prognosis. This is a different clinical question from CcSinglecell's.

### Gaps
- I found no 2024–2026 source on tissue methylation classifiers or cfDNA fragmentomics (e.g. DELFI-style) as validated *prognostic* (rather than detection) layers in CRC. This was not searched in depth.
- A source for the 2026 ALTAIR randomised trial (ctDNA-guided therapy, primary endpoint not met, HR 0.79) appeared only in a search summary with no verifiable URL, so it is omitted from the findings.
- I found no matched CPTAC–TCGA comparison of protein against mRNA prognostic models.

## 7. Interpretability that holds up (SurvSHAP(t), knockoffs, stability selection, attention caveats)

### Takeaway
SurvSHAP(t) is the established model-agnostic, time-dependent attribution for survival models. Attention heatmaps in MIL pathology are not reliably faithful and should be checked with perturbation tests (patch flipping, confounder insertion). For biomarker *selection*, stability under resampling is the real test. Sparse and pathway networks remain far less stable than classical models, and I found no 2024–2026 survival-specific knockoff or stability-selection benchmark.

### Cited Findings
- **SurvSHAP(t) (Krzyziński et al., Knowledge-Based Systems 262, 2023) [older].** The first time-dependent SHAP explanation for black-box survival models. Values sum to the predicted survival function, it detects time-varying effects, and it beat SurvLIME in the authors' tests. A Python implementation is available. — [arXiv 2208.11080](https://arxiv.org/abs/2208.11080)
- **Attention maps in digital pathology (Albuquerque et al., Roche, 2024).** There is no standard evaluation framework. When artificial confounders were inserted, attention generally highlighted relevant regions, but robustness depended on the type and number of confounders. — [arXiv 2407.02484](https://arxiv.org/abs/2407.02484)
- **xMIL (2024).** Proposes an explanation method with better faithfulness scores than attention on biomarker tasks. — [arXiv 2406.04280](https://arxiv.org/pdf/2406.04280v3.pdf)
- **Patch-flipping faithfulness benchmark for MIL heatmaps (2026).** — [arXiv 2603.08328](https://arxiv.org/pdf/2603.08328)
- **Stability of biologically informed and sparse networks.** Feature-importance stability is "still far off from traditional models". — [bioRxiv 2026 P-NET expansion](https://www.biorxiv.org/content/10.64898/2026.04.19.719454v1.full)
- **Feature selection in the 2024–2025 single-cell integration benchmark.** MOFA+ had the most reproducible feature selection. — [Nat Methods 2025](https://www.nature.com/articles/s41592-025-02856-3)

### Inferences
- For a successor, the "biomarker" claim should be tied to:
  - selection frequency across bootstrap or CV resamples (stability selection with a per-family error bound, run on the elastic-net Cox path);
  - a permutation or random-wiring null;
  - ideally replication of the direction of effect in an external cohort.
- SHAP rankings without a stability check repeat the same optimism problem.
- Attention weights from MIL models should not be reported as biology without perturbation-based faithfulness checks.

### Gaps
- Searches turned up no 2024–2026 papers specifically on model-X knockoffs or stability selection for Cox models in multi-omics, so I cannot cite current benchmarks for them. The canonical references (Meinshausen & Bühlmann 2010; Candès et al. 2018) were not fetched this session.
- I found no study testing whether SurvSHAP(t) rankings are stable across resamples in high-dimensional omics.
