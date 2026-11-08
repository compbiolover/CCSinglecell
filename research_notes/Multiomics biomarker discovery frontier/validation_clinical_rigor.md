# Validation standards, clinical endpoints and cohorts for multi-omics CRC biomarkers (2024–2026)

Scope note: written 2026-11-07 for planning the successor to CcSinglecell (TCGA-COAD/READ OS model; ΔC ≈ 0 over age+stage at n=350, 78 events). Items marked **[memory, not re-fetched]** come from the researcher's prior knowledge with a DOI, and were not re-checked in this session. Treat them as lower confidence than items backed by a link fetched or searched here.

## Q1. Where is the unmet clinical need in CRC prognosis/prediction, and which questions remain open for multi-omics?

### Takeaway
Post-surgical ctDNA has become the dominant prognostic biomarker in stage II/III CRC. Randomized trials (DYNAMIC, DYNAMIC-III, COBRA) show that it does **not** yet work as a reliable *predictive* tool for de-escalation or escalation, and NCCN calls it "prognostic, not predictive". The open questions where tissue or multi-omics could add value are: (1) relapse in ctDNA-negative patients (about 12–15% at 3 years in stage III); (2) which ctDNA-positive patients benefit from which regimen, since escalation failed; (3) predicting response to chemotherapy duration or regimen (IDEA/Immunoscore, still unvalidated); (4) anti-EGFR selection beyond RAS/sidedness; (5) pMMR/MSS immunotherapy subsets. In dMMR localized colon cancer, neoadjuvant/adjuvant ICI response is now so high that a prognostic biomarker has little room to matter there.

### Cited Findings
**ctDNA / MRD (stage II–III)**
- DYNAMIC (stage II, NEJM 2022): ctDNA-guided management was non-inferior to standard care for recurrence-free survival and reduced adjuvant chemo use. **[memory: 455 randomized; chemo 15% vs 28%; 2-y RFS 93.5% vs 92.4%]** — [Tie et al. NEJM 2022, DOI 10.1056/NEJMoa2200075](https://www.nejm.org/doi/full/10.1056/NEJMoa2200075)
- DYNAMIC-III (stage III, phase 2/3; ESMO 2025 LBA9; *Nature Medicine* 31(12):4291–4300, Dec 2025, DOI 10.1038/s41591-025-04030-w):
  - 1,002 randomized at 66 sites (AUS/NZ/Canada); 968 evaluable; 702 (72.5%) ctDNA-negative.
  - ctDNA-negative de-escalation cut oxaliplatin use (34.8% vs 88.6%), grade ≥3 AEs (6.2% vs 10.6%) and hospitalisations (8.5% vs 13.2%).
  - 3-y RFS was 85.3% vs 88.1%, so the **non-inferiority margin was not met**.
  - In the pre-planned clinically low-risk (T1-3 N1) subgroup, 3-y RFS was 91% vs 93.2%.
  - In ctDNA-positive patients, escalation did not improve outcomes (2-y RFS 51% vs 61%).
  - Median follow-up is reported as 47 months in the paper and 45 months in news coverage.
  - Sources: [ASCO Post Oct 2025](https://ascopost.com/news/october-2025/postsurgical-ctdna-testing-in-stage-iii-colon-cancer-for-treatment-de-escalation/); [CCTG release](https://www.ctg.queensu.ca/cctg_media_releases/dynamic-iii-trial-results-ctdna-blood-test-powerful-tool-guide-treatment); [themednet record of Nat Med paper](https://www.themednet.org/publications/circulating-tumor-dna-guided-adjuvant-therapy-in-locally-advanced-colon-cancer-the-randomized-phase-2-3-dynamic-iii-trial); [ESMO commentary "not ready for prime time"](https://www.esmo.org/insights-by-esmo/all-news/ctdna-guided-adjuvant-chemotherapy-in-colon-cancer-not-ready-for-prime-time).
- Interpretation conflict: sponsors framed DYNAMIC-III as supporting ctDNA-guided de-escalation, while ESMO's reporter headlined it "not ready for prime time". The same ESMO coverage notes that **PEGASUS** also missed its primary endpoint, partly because of reduced power — [CCTG](https://www.ctg.queensu.ca/cctg_news/trial-results-show-blood-test-can-guide-treatment-decisions-advanced-colon-cancer); [ESMO](https://www.esmo.org/insights-by-esmo/all-news/ctdna-guided-adjuvant-chemotherapy-in-colon-cancer-not-ready-for-prime-time)
- COBRA (NRG-GI005, low-risk stage IIA; ASCO GI 2024): 635 randomized. In the phase II analysis of the first 16 ctDNA-positive patients, ctDNA clearance was 3/7 (43%) with observation vs 1/9 (11%) with chemotherapy (P=0.98). The trial closed per protocol. Commentators note that ctDNA clearance is not a validated surrogate — [ASCO Post Jan 2024](https://ascopost.com/news/january-2024/role-of-postoperative-ctdna-in-patients-with-low-risk-stage-iia-colon-cancer/); [ASCO Daily News](https://dailynews.ascopubs.org/do/closure-cobra-trial-creates-uncertainty-ctdna-guided-adjuvant-therapy-stage-ii-colon)
- CIRCULATE-Japan GALAXY (observational, Signatera tumor-informed assay):
  - 2023 analysis: 1,039 patients, 18% ctDNA-positive at 4 weeks. Post-surgical ctDNA positivity identified the patients most likely to benefit from adjuvant chemotherapy (HR 6.59 reported) — [Kotani et al. Nat Med 2023, PMC](https://pmc.ncbi.nlm.nih.gov/articles/9873552/)
  - 2024 analysis (Nat Med, 16 Sep 2024, DOI 10.1038/s41591-024-03254-6): 2,240 patients. The NCC release reports 2-y recurrence-free proportions of 20.6% (ctDNA+) vs 85.1% (ctDNA−). ctDNA predicted recurrence independently of tumor genomics — [NCC press release](https://www.ncc.go.jp/en/information/press_release/2024/0917_1/index.html)
  - In ctDNA-negative patients, adjuvant chemo showed no significant benefit (18-mo DFS 94.9% vs 91.5%; adjusted HR 1.71, 0.80–3.7). This is non-randomized — [OncLive, CIRCULATE-Japan update](https://www.onclive.com/view/updated-results-of-circulate-japan-in-colorectal-cancer)
  - Company-reported (Natera, ESMO 2024) 24-mo OS: 100% with sustained clearance, 82% with transient clearance, 61% with no clearance — [Clinical Trials Arena](https://www.clinicaltrialsarena.com/?p=216609)
- VEGA (CIRCULATE-Japan; ctDNA-negative high-risk stage II / low-risk stage III; surgery alone vs 3-mo CAPOX, phase III non-inferiority): no readout found as of these searches — [CIRCULATE-Japan design, Cancer Sci 2021](https://onlinelibrary.wiley.com/doi/10.1111/cas.14926)
- CIRCULATE-US (NRG-GI008, NCT05174169): recruiting as of an Aug 2026 registry snapshot; target 1,912 participants; estimated primary completion 2029-03. ctDNA-negative patients are randomized to immediate chemo or serial ctDNA monitoring; ctDNA-positive patients receive intensification — [trial record summary](https://dashboard.datalookout.com/trial/NCT05174169); [NRG Oncology](https://www.nrgoncology.org/Home/News/Post/testing-ctdna-as-a-biomarker-for-adjuvant-therapy-intensity-in-stage-3-colon-cancer-nrg-gi008-the-circulate-us-study/)
- NCCN Colon v1.2025: ctDNA is "prognostic, but not predictive" in the adjuvant MSS/pMMR setting. Evidence is "insufficient" for routine use outside trials. De-escalation based on ctDNA is **not recommended**, and ctDNA is "not recommended for surveillance". Sources disagree on the exact version/date — [OncLive](https://www.onclive.com/view/nccn-updates-ctdna-stance-in-colon-rectal-and-mcc-guidelines); [GenomeWeb](https://www.genomeweb.com/node/619036)
- The ESMO localised colon cancer guideline is still the 2020 version (Argilés et al., Ann Oncol, DOI 10.1016/j.annonc.2020.06.022; marked "Current", last updated 1 Oct 2020). It includes Immunoscore as a prognostic refinement and restricts ctDNA to studies. No 2024–2026 ESMO update was found — [ESMO guideline page](https://www.esmo.org/Guidelines/Gastrointestinal-Cancers/Early-Colon-Cancer); [UCL copy](https://discovery.ucl.ac.uk/id/eprint/10106865/)

**Adjuvant duration (IDEA) and Immunoscore**
- IDEA final OS (pooled 6 RCTs, Lancet Oncol 2020): 5-y OS was 82.4% (3 mo) vs 82.8% (6 mo), HR 1.02 (0.95–1.11); non-inferiority was not formally met overall (FDR-adj p=0.058), but it was met for CAPOX (HR 0.96). The absolute difference was 0.4%. In high-risk patients (T4 and/or N2, n=5,273), 5-y OS was 72.0% vs 74.1% (HR 1.08) — [André et al., PMC7786835](https://pmc.ncbi.nlm.nih.gov/articles/PMC7786835); [ASCO Post](https://ascopost.com/news/january-2021/final-overall-survival-results-of-the-idea-collaboration-3-vs-6-months-of-adjuvant-chemotherapy-for-stage-iii-colon-cancer)
- Immunoscore in IDEA-France (n=1,322 stage III): low IS gave HR 1.54 for relapse/death (3-y DFS 66.8% vs 77.1%). A predictive interaction for 3 vs 6 months of FOLFOX was reported (P=0.0015 in T/N subgroups), but the authors say it "needs to be confirmed in an external validation cohort" — [ASCO abstract](https://ace.asco.org/content/258839/abstract/171171); [cancer.fr summary](https://www.cancer.fr/professionnels-de-sante/veille/nota-bene-cancer/bulletin-n-447/prognostic-and-predictive-value-of-the-immunoscore-in-stage-iii-colon-cancer-patients-treated-with-oxaliplatin-in-the-prospective-idea-france-prodi)
- The pooled SCOT + IDEA-HORG analysis (JCO 2024) confirmed Immunoscore's prognostic value in CAPOX-treated stage III patients, but was **not powered** for the predictive (duration) claim — [Christie repository record](https://christie.openrepository.com/entities/publication/14c50cd9-2df3-4e98-9cf6-1e5e88e79cf7); [TransSCOT publications](https://www.gla.ac.uk/research/az/transscot/publications/)
- Original Immunoscore international validation: **[memory, not re-fetched: Pagès et al. Lancet 2018, n≈2,681 stage I–III, DOI 10.1016/S0140-6736(18)30789-X]**

**Immunotherapy (dMMR/MSI-H)**
- NICHE-2 (single-arm phase 2, locally advanced dMMR colon; neoadjuvant nivolumab + ipilimumab):
  - 115 enrolled; 113 (98%) had timely surgery; pCR 68%. NEJM 2024 — [Oncodaily, NEJM publication note](https://oncodaily.com/blog/79062); [2 Minute Medicine](https://2minutemedicine.com/neoadjuvant-immunotherapy-in-mismatch-repair-deficient-colon-cancer)
  - 3-y DFS 100% (ESMO 2024 LBA24); ctDNA negative at the MRD timepoint in all patients — [ESMO OncologyPRO](https://oncologypro.esmo.org/meeting-resources/esmo-congress-2024/neoadjuvant-immunotherapy-in-locally-advanced-mmr-deficient-colon-cancer-3-year-disease-free-survival-from-niche-2)
  - Lancet long-term update: data cutoff 19 Jun 2026; median follow-up 58.7 mo; all 111 per-protocol patients alive and disease-free; exploratory 5-y OS 100%; ITT 3-y DFS 99%. Seen only via secondary sources — [Oncodaily NICHE-2](https://oncodaily.com/gastrointestinal-oncology/niche-2)
- ATOMIC (Alliance A021502; stage III dMMR; adjuvant FOLFOX ± atezolizumab; n=712): 3-y DFS 86.3% vs 76.2%, HR 0.50. Grade ≥3 TRAEs were 71.7% vs 62.1%. NEJM 2026;394(12):1155–1166, DOI 10.1056/NEJMoa2507874. NCCN has incorporated it and extends it to T4bN0 — [Targeted Oncology Mar 2026](https://www.targetedonc.com/view/atezolizumab-plus-chemotherapy-improves-dfs-in-dmmr-colon-cancer); [ecancer](https://ecancer.org/en/news/27961-atezolizumab-plus-chemo-after-surgery-shows-benefit-for-stage-iii-dmmr-colon-cancer)

**Anti-EGFR and subtypes (metastatic)**
- PARADIGM ctDNA analysis (Nat Med 2024, DOI 10.1038/s41591-023-02791-w):
  - 733 RAS-WT mCRC patients, with a baseline plasma panel covering KRAS/NRAS/PTEN/EGFR-ECD mutations, HER2/MET amplifications and ALK/RET/NTRK1 fusions.
  - In "negative hyperselected" patients, median OS was 40.7 vs 34.4 mo with panitumumab vs bevacizumab. Patients with any alteration did not benefit (HR 1.13).
  - Left-sided hyperselected patients: 42.1 vs 35.5 mo (HR 0.76), interaction P=0.171 (hypothesis-generating).
  - 79.4% of left-sided vs 50.3% of right-sided tumors were hyperselected.
  - Sources: [Nature Medicine](https://www.nature.com/articles/s41591-023-02791-w); [PMC10957476](https://pmc.ncbi.nlm.nih.gov/articles/PMC10957476)
- CALGB/SWOG 80405 CMS analysis (JCO 2019, DOI 10.1200/JCO.18.02258):
  - NanoString CMS on 581 primary tumors; CMS was prognostic.
  - Strong treatment × CMS interaction for OS (P<0.001): CMS1 favoured bevacizumab, CMS2 favoured cetuximab (P=0.0046).
  - Unstratified arms did not differ (~30 vs 29 mo; HR 0.88, P=0.08).
  - A CMS classifier had to be retrained because of panel mismatch.
  - Sources: [KU Leuven repository](https://lirias.kuleuven.be/2798750); [ASCO](https://ace.asco.org/content/10.1200JCO.18.02258/journal/JCO)
- iCMS (Joanito et al. Nat Genet 2022, DOI 10.1038/s41588-022-01100-4):
  - Based on 373,058 single-cell transcriptomes from 63 patients, the study identified two epithelial states, iCMS2 and iCMS3.
  - iCMS3 includes MSI-H tumors and about one-third of MSS tumors. CMS4 with iCMS3 epithelium had the worst prognosis.
  - The proposed 5-class IMF scheme combines iCMS, MSI and fibrosis, validated on bulk data from 3,614 patients.
  - Sources: [SKKU record](https://pure.skku.edu/en/publications/single-cell-and-bulk-transcriptome-sequencing-identifies-two-epit/); [Clinical Research News](https://www.clinicalresearchnewsonline.com/news/2022/08/11/classification-system-for-colorectal-cancer-needs-an-update)
- Existing gene-expression risk scores add little over clinical factors: in PETACC-3 (688 FFPE), four tested gene-expression risk scores "provide prognostic information but contribute only marginally to improving models based on established risk factors" — [Di Narzo et al. JNCI 2014, dju247](https://doc.rero.ch/record/299791/files/dju247.pdf) (older; same pattern as CcSinglecell's ΔC≈0)
- Real-world outcome prediction: MSK-CHORD found that NLP-derived clinical features (e.g., sites of disease) outperformed genomic data or stage alone for OS prediction, both in cross-validation and in an external multi-institution dataset — [Jee et al. Nature 2024, PubMed 39506116](https://pubmed.ncbi.nlm.nih.gov/39506116/)

### Inferences
- Endpoints where a successor could matter, in priority order:
  1. **Relapse among ctDNA-negative stage III patients.** In DYNAMIC-III's ctDNA-negative de-escalated arm, about 15% relapsed by 3 years. A tissue or multi-omic classifier that finds the ctDNA-negative-but-high-risk group would hit a clear decision gap (the "residual risk after negative MRD" problem). The required ΔC must be judged *conditional on ctDNA status*, not on stage alone.
  2. **Regimen choice for ctDNA-positive patients.** Escalation to FOLFOXIRI-type therapy did not help (51% vs 61% 2-y RFS), so a predictive marker of chemo-refractoriness would be useful. This needs treatment-randomized cohorts (IDEA, SCOT, FOCUS, PETACC-3, CALGB 80405).
  3. **Treatment-by-biomarker interaction** (CMS × anti-EGFR/bevacizumab, Immunoscore × duration). These are the only designs that answer *predictive* questions. TCGA OS cannot do this.
- The baseline to beat is **stage + ctDNA (+ MMR + Immunoscore where available)**, not age + stage. A biomarker that gives ΔC≈0 over age + stage is very unlikely to add over ctDNA.
- Localized dMMR colon cancer is close to "solved" for neoadjuvant ICI (NICHE-2 has no recurrences). The biomarker need moves to the MSS/pMMR majority and to the small group of ICI-refractory dMMR patients.

### Gaps
- Exact DYNAMIC-III non-inferiority margin and HR with CI: full text not fetched.
- GALAXY 2024 multivariable HR and ACT-benefit HRs by ctDNA status: Nature paywall; only the press release was read.
- No VEGA, ALTAIR or PEGASUS final numbers were verified. ALTAIR appears as an accepted Nat Med 2026 manuscript ([Oxford ORA](https://ora.ox.ac.uk/objects/uuid:b00c9f1b-f05d-4de5-bda1-5bbd61a0d8df/files/rpn89d9157)), with results unconfirmed.
- No 2025–2026 data were found on ColoPrint (PARSC) or Oncotype DX Colon. Their current guideline status beyond the 2020 ESMO text was not established.
- CheckMate 8HW / KEYNOTE-177 (metastatic dMMR) and primary ICI resistance rates were not researched here.

## Q2. Better endpoints than TCGA overall survival

### Takeaway
TCGA-CDR (Liu et al., Cell 2018) calls **PFI** the most reliable TCGA endpoint and recommends it without reservation for COAD and READ. TCGA OS is weak because follow-up is short: median 22.0 months (COAD) and 20.0 months (READ). READ OS is only "use with caution", READ DFI is not recommended, and OS includes non-cancer deaths. Switching to PFI adds only about 20% more events, so TCGA is still underpowered for discovery and should be used, at most, as one site in an internal-external validation.

### Cited Findings
All from [Liu et al., Cell 2018;173(2):400–416.e11, DOI 10.1016/j.cell.2018.02.052](https://escholarship.org/content/qt9ng938g3/qt9ng938g3.pdf) (pre-2024; still the reference standard):
- Table 3, COAD (N=459): OS ✓ (102 events / 357 censored); PFI ✓ (123/336); DFI ✓ (24/166); DSS ✓ approximate (64/379).
- Table 3, READ (N=170): OS ✓* (26 events), noted "need a longer follow-up for OS, DSS, and DFI"; PFI ✓ (39 events); DFI × not recommended (7 events, "too small"); DSS ✓ approximate* (15 events).
- Table 2: median follow-up is 22.0 mo (COAD) and 20.0 mo (READ). For READ, OS passed the step-1 tests but "may need a longer follow-up".
- "Among the 4 endpoints, PFI is generally considered a more informative endpoint for TCGA pan-cancer studies." PFI is recommended without reservation in all but 4 of 33 cancer types.
- Rationale: OS "may weaken a clinical study as deaths because of non-cancer causes do not necessarily reflect tumor biology". OS and DSS need longer follow-up than DFI or PFI.
- Caveats for molecular analyses:
  - Molecular subgrouping "partition[s] the sample sets, potentially compromising the statistical significance".
  - Conclusions "should be further validated on independent tumor datasets".
  - Use molecular data from primary tumors only.
  - DSS is only "approximately defined" (derived from cause-of-death fields).
- Comparable cohort endpoints: GENIE BPC CRC defines OS as death censored at last known alive contact, and ascertainment of death varies by institution — [AACR GENIE BPC CRC public release](https://www.aacr.org/professionals/research/aacr-project-genie/bpc/crc-public-release/)

### Inferences
- With TCGA COAD+READ PFI (about 162 events in 629 patients before QC), a modest-effect model is still undersized by Riley criteria (see Q4 table). TCGA should be a validation or "internal-external" fold, not the discovery engine.
- For stage II/III questions, the clinically standard endpoints are DFS/RFS (IDEA, DYNAMIC-III, ATOMIC, NICHE-2) and TTR (time to recurrence). The successor should harmonize to **RFS or TTR at 3 years** across trial cohorts (GSE39582 RFS, PETACC-3 RFS/OS, SCOT/IDEA DFS) and handle competing non-cancer death explicitly (Fine–Gray or cause-specific hazards).

### Gaps
- No 2024–2026 re-curation of TCGA-COAD/READ follow-up (e.g., updated GDC clinical data) was checked. Whether newer GDC releases extend follow-up is unknown.

## Q3. External validation cohorts for CRC with omics and outcomes

### Takeaway
The most useful cohorts for a successor are those with **randomized treatment and RFS/DFS**: PETACC-3 (open array data), FOCUS/S:CORT, SCOT/TransSCOT, QUASAR2, CALGB/SWOG 80405 and IDEA sub-studies. Large clinicogenomic real-world resources add scale: MSK-CHORD (5,543 CRC, open via cBioPortal, non-commercial), GENIE BPC CRC (1,485), 100kGP (2,023 WGS), Hartwig (metastatic WGS). Deep multi-omics cohorts (AC-ICAM, CPTAC-COAD) are small but rich. Most trial cohorts are controlled-access through consortium agreements.

### Cited Findings
| Cohort | n / setting | Omics | Endpoints | Access | Source |
|---|---|---|---|---|---|
| GSE39582 (CIT/Marisa) | 566 stage I–IV colon | Affymetrix U133Plus2 | RFS, OS | Open (GEO) | [PMC5111107](https://pmc.ncbi.nlm.nih.gov/articles/PMC5111107); **[memory: Marisa et al. PLoS Med 2013, DOI 10.1371/journal.pmed.1001453]** |
| GSE17536 / GSE14333 | **[memory: 177 / 290]** | Affymetrix | DFS/OS (17536), DFS (14333) | Open (GEO) | [PMC5891048](https://pmc.ncbi.nlm.nih.gov/articles/PMC5891048); **[memory: Smith 2010 Gastroenterology; Jorissen 2009 Clin Cancer Res]** |
| PETACC-3 | 688 FFPE stage II–III, RCT (5FU/LV ± irinotecan) | Almac ADXCRC array; E-MTAB-990 (+64 normals E-MTAB-1026) | RFS, OS | ArrayExpress accession cited; open/controlled status not confirmed | [Di Narzo JNCI 2014](https://doc.rero.ch/record/299791/files/dju247.pdf); [PMC3840702](https://pmc.ncbi.nlm.nih.gov/articles/PMC3840702) |
| FOCUS / S:CORT | ~2,000 patients' transcriptomes in S:CORT; FOCUS mCRC subset n=360 (GSE156915) | RNA (array/seq), DNA | PFS/OS, treatment-randomized | Consortium; contact Prof Tim Maughan | [Oxford S:CORT](https://www.cancer.ox.ac.uk/research/networks/scort); [PMC10937375](https://pmc.ncbi.nlm.nih.gov/articles/PMC10937375) |
| SCOT / TransSCOT | >3,000 patients' tissue + blood; stage III/high-risk II, 3 vs 6 mo | IHC/multiplex, Immunoscore, MMR AI, genomics (in progress) | DFS | Consortium (Glasgow) | [TransSCOT](https://gla.ac.uk/research/az/transscot); [projects](https://www.gla.ac.uk/research/az/transscot/projects/) |
| QUASAR2 | 1,069 in Lancet Oncol 2024 validation | Multiplex immune (CD8/CD20/FoxP3/CD68) | DFS | Consortium (Oxford) | [TransSCOT publications](https://www.gla.ac.uk/research/az/transscot/publications/) |
| CALGB/SWOG 80405 | 581 (NanoString CMS); 554 (GSE196576 immune); ~1,284 RNA-seq in a 2024 preprint | NanoString, RNA-seq, germline | OS, PFS, randomized cetuximab vs bevacizumab | GSE196576 on GEO; rest via Alliance/NCTN (dbGaP ID not found) | [KU Leuven](https://lirias.kuleuven.be/2798750); [OmicsDI GSE196576](https://www.omicsdi.org/dataset/geo/GSE196576); [Research Square preprint](https://www.researchsquare.com/article/rs-3673588/v1) |
| AC-ICAM (Sidra-LUMC) | 348 treatment-naive colon | RNA-seq, WES, TCR-seq, 16S microbiome; tumor + matched normal | OS/DSS-type outcomes | EGA/dbGaP study phs002978 (controlled); cBioPortal coad_silu_2022 | [Roelands Nat Med 2023, DOI 10.1038/s41591-023-02324-5](https://pmc.ncbi.nlm.nih.gov/articles/PMC10202816); [EGA](https://ega-archive.org/studies/phs002978); [cBioPortal KB](https://jim-bo.github.io/cbio-kb/datasets/coad_silu_2022.html) |
| CPTAC-COAD | 110 tumors + normals | WES, RNA, proteome, phosphoproteome | Limited outcomes | Public, under a CPTAC data use agreement (PDC); earlier TCGA-CPTAC set: 64 COAD + 31 READ | [Vasaikar Cell 2019, DOI 10.1016/j.cell.2019.03.030](https://pmc.ncbi.nlm.nih.gov/articles/PMC6768830); [NCI blog](https://www.cancer.gov/news-events/cancer-currents-blog/2019/proteogenomics-cptac-colon-cancer) |
| MSK-CHORD | 5,543 CRC (of 24,950 MSK-IMPACT patients) | MSK-IMPACT panel + NLP-derived clinical features | OS, treatment lines, metastatic sites | Open on cBioPortal (msk_chord_2024), CC BY-NC-ND 4.0; commercial use via datarequests@mskcc.org | [Jee et al. Nature 2024, DOI 10.1038/s41586-024-08167-5](https://pubmed.ncbi.nlm.nih.gov/39506116/); [MSK data catalog](https://datacatalog.mskcc.org/dataset/11458) |
| AACR GENIE BPC CRC v2.0-public | 1,485 (MSK, DFCI, VICC; sequenced 2014–2018) | Panel NGS | OS, PFS (imaging/notes-based), treatment histories | Synapse (terms of use); subset on cBioPortal | [AACR BPC CRC release](https://www.aacr.org/professionals/research/aacr-project-genie/bpc/crc-public-release/) |
| UK 100,000 Genomes CRC | 2,023 WGS (1,898 primary, 122 metastases, 3 recurrences); 100× tumor / 33× blood | WGS | Prognostic associations of 4 new MSS subgroups | Controlled; Genomics England Research Environment **[memory]** | [Cornish et al. Nature 2024;633:127–136, DOI 10.1038/s41586-024-07747-9](https://publications.mpi-cbg.de/Cornish_2024_8887.pdf) |
| Hartwig Medical Foundation | 429 mCRC biopsies (CPCT-02); pan-cancer 2,520 genomes / 2,399 patients (2019) | WGS (+RNA for subsets) | Treatment response/OS for some cohorts | Data Access Board review + licence, academic use; use limited to the approved project | [Utrecht dspace mCRC paper](https://dspace.library.uu.nl/handle/1874/442157); [Priestley Nature 2019](https://link.springer.com/10.1038/s41586-019-1689-y); [Erasmus repository](https://repub.eur.nl/pub/134166/Repub_134166_O-A.pdf) |

- Data-access caveat for AC-ICAM: the paper says the data were deposited in dbGaP and cBioPortal, but the exact application route was not confirmed — [PMC10202816](https://pmc.ncbi.nlm.nih.gov/articles/PMC10202816)
- Note: the CALGB 80405 preprint's 79 "OS-associated" genes use p<0.1. Treat them as hypothesis-generating — [Research Square](https://www.researchsquare.com/article/rs-3673588/v1)

### Inferences
- **Recommended validation stack for a stage II/III RFS model:**
  - Discovery: GSE39582 + PETACC-3 (both array-based with RFS; PETACC-3 is an RCT).
  - Internal-external folds: TCGA PFI, GSE17536/14333, AC-ICAM.
  - Confirmatory: a consortium cohort with ctDNA or Immunoscore (TransSCOT/IDEA-HORG, S:CORT).
  - Platform heterogeneity (Affymetrix vs Almac vs RNA-seq vs NanoString) means any signature must be built on platform-robust features, e.g., rank-based or single-sample scores, or CMS/iCMS calls.
- For metastatic or predictive work: CALGB 80405 and FOCUS (randomized), with MSK-CHORD/GENIE BPC for real-world treatment-line outcomes. MSK-CHORD's own result (clinical NLP > genomics) is a warning that panel genomics alone is unlikely to beat rich clinical features.

### Gaps
- No usable search hit was found for an **HTAN** CRC atlas with outcomes; HTAN's coverage and endpoints are unverified.
- PETACC-3 access conditions (open vs controlled on ArrayExpress), the CALGB 80405 dbGaP accession, and exact current Hartwig CRC counts were not confirmed.
- No CRC ctDNA trial was found that releases patient-level ctDNA + tissue omics openly. DYNAMIC, GALAXY and DYNAMIC-III data-sharing routes were not researched.

## Q4. Reporting and methodology standards that define "rigorous" now

### Takeaway
The 2024–2025 standard set is: **TRIPOD+AI** (BMJ 2024, 27 items) for reporting; **PROBAST+AI** (BMJ 2025) for risk of bias; the **BMJ 2024 "Evaluation of clinical prediction models" series** (parts 1–3) for validation design and external-validation sample size; and **Riley pmsampsize** criteria for development sample size. REMARK, STARD and BRISQ are older, biomarker-specific complements. Performance must be reported as discrimination (time-dependent C/AUC), calibration (slope, calibration-in-the-large, curves at a fixed horizon) and clinical utility (decision-curve net benefit). Validation should use internal-external CV across sites, not a single random split.

### Cited Findings
- **TRIPOD+AI** (Collins et al., BMJ 2024;385:e078378; erratum q902) has a 27-item checklist plus a 13-item abstract checklist. It covers regression and ML and **supersedes TRIPOD 2015** — [LMU PDF](https://epub.ub.uni-muenchen.de/116670/1/bmj-2023-078378.full.pdf); [Oxford ORA](https://ora.ox.ac.uk/objects/uuid:414e5fb8-5575-4b67-b94d-a0833a7cc58a)
- **PROBAST+AI** (Moons et al., BMJ 2025;388:e082505, 24 Mar 2025, DOI 10.1136/bmj-2024-082505) has separate parts for development (16 signalling questions on quality and applicability) and evaluation (18 questions on risk of bias and applicability). It uses four domains (participants/data sources, predictors, outcome, analysis) and may replace PROBAST 2019 — [Oxford ORA](https://ora.ox.ac.uk/objects/uuid:8b15c087-459b-43fa-ad03-1f66ed6ba370); [Maastricht](https://cris.maastrichtuniversity.nl/en/publications/probastai-an-updated-quality-risk-of-bias-and-applicability-asses/)
- **BMJ 2024 evaluation series:**
  - Part 1 (Collins, e074819): internal vs external validation rationale.
  - Part 2 (e074820): how to run an external validation, including clinical usefulness.
  - Part 3 (Riley, e074821): sample size for external validation. It notes many published validations are too small and advocates tailored calculations over rules of thumb.
  - Sources: [Birmingham PDF part 1](https://research.birmingham.ac.uk/files/223789914/bmj-2023-074819.full.pdf); [part 3 record](https://research.birmingham.ac.uk/en/publications/evaluation-of-clinical-prediction-models-part-3-calculating-the-s)
- **Development sample size (pmsampsize, survival option):**
  - Inputs: overall event rate, mean follow-up, prediction horizon, anticipated Cox–Snell R² (or optimism-adjusted C to approximate it), and target shrinkage of 0.9.
  - It returns the minimum n and number of events across three criteria.
  - Sources: Riley et al. Stat Med 2019 (Part II, DOI 10.1002/sim.7992) and BMJ 2020 — [pmsampsize docs](https://rdrr.io/cran/pmsampsize/man/pmsampsize.html); [CRAN manual](https://cran.r-universe.dev/pmsampsize/doc/manual.html)
- A newer Riley extension targets precise *individual-level* risk at specified time points, implemented in **pmstabilityss** — [Oxford ORA](https://ora.ox.ac.uk/objects/uuid:7df66d0c-a278-4e43-8223-75685a5b3bad/files/rfj2364154)
- External validation for time-to-event outcomes has its own sample-size method, which needs assumptions about the linear-predictor distribution and the event and censoring distributions — **[memory: Riley et al. Stat Med 2022, DOI 10.1002/sim.9275]**
- **Illustrative computation (this session).** This applies the closed-form shrinkage criterion n = p / ((S−1)·ln(1 − R²cs/S)) with S=0.9 from the Riley criteria above, at CcSinglecell's event fraction (78/350 = 22%):
  - R²cs=0.05: p=5 needs n=875 (~195 events); p=10 needs 1,750 (~390); p=20 needs 3,500 (~780).
  - R²cs=0.10: p=5 needs 425 (~95 events); p=10 needs 850 (~190); p=20 needs 1,699 (~379).
  - R²cs=0.15: p=5 needs 275 (~62 events); p=10 needs 549 (~123); p=20 needs 1,097 (~245).
  - This is criterion 1 only; pmsampsize takes the maximum over three criteria, so real minima are greater than or equal to these.
- **Decision-curve analysis** for survival: net benefit must account for censoring. **dcurves** (R) implements time-to-event DCA at a fixed horizon, following Vickers et al. 2008. **bayesDCA** offers Bayesian DCA with censored outcomes — [dcurves](https://cloud.r-project.org/web/packages/dcurves/refman/dcurves.html); [bayesDCA arXiv 2308.02067](https://arxiv.org/pdf/2308.02067)
- Older or still-current biomarker standards (pre-2024):
  - REMARK (McShane et al. 2005; explanation and elaboration Altman et al. BMC Med 2012), cited as the tumour-marker reporting standard in TCGA-CDR's references — [TCGA-CDR PDF](https://escholarship.org/content/qt9ng938g3/qt9ng938g3.pdf)
  - **[memory, not re-fetched]** STARD 2015 (Bossuyt et al. BMJ 2015, DOI 10.1136/bmj.h5527; diagnostic accuracy)
  - **[memory, not re-fetched]** BRISQ (Moore et al. Cancer Cytopathol 2011, DOI 10.1002/cncy.20147; pre-analytical biospecimen reporting)
- **[memory, not re-fetched]** Survival calibration: graphical calibration curves and ICI/E50/E90 for Cox models at a horizon (Austin, Harrell & van Klaveren, Stat Med 2020, DOI 10.1002/sim.8570). Internal-external cross-validation across sites or studies (Steyerberg & Harrell, J Clin Epidemiol 2016, DOI 10.1016/j.jclinepi.2015.04.005).
- **[memory, not re-fetched]** Benchmarks show molecular data add little over clinical data:
  - Yuan et al. (Nat Biotechnol 2014, DOI 10.1038/nbt.2940) found modest added value of TCGA molecular data over clinical variables.
  - Herrmann et al. (Brief Bioinform 2021, DOI 10.1093/bib/bbaa167) benchmarked multi-omics survival methods on TCGA and found that clinical-only models are often hard to beat.

### Inferences
- CcSinglecell's harness already matches much of this (bootstrap ΔC, repeated CV, calibration at a horizon, power curve). Missing pieces to add:
  - a Riley sample-size gate *before* fitting;
  - internal-external CV by cohort or site rather than random folds;
  - time-dependent AUC and integrated Brier score alongside C;
  - calibration slope and ICI;
  - decision-curve net benefit at a clinically defined threshold (e.g., the 3-y recurrence risk at which adjuvant oxaliplatin would be given or withheld);
  - a TRIPOD+AI checklist and PROBAST+AI self-assessment shipped with each report.
- The table above makes the core issue explicit: with 78 events, even a model with R²cs=0.15 supports about 6 parameters before shrinkage exceeds 10%. Penalized Cox on thousands of genes is in heavy-shrinkage territory, so ΔC≈0 is the expected outcome, not a surprise.

### Gaps
- No 2024–2026 dedicated guidance on DCA for time-to-event models (e.g., a Vickers/Van Calster consensus paper) was located.
- The criteria 2 and 3 formulas for survival were not reproduced. The table above is criterion 1 only.
- No CRC-specific consensus (e.g., ESMO Precision Medicine Working Group) on validation standards for multi-omics prognostic signatures was found.

## Q5. Modern statistical tools for biomarker selection with error control

### Takeaway
For high-dimensional survival selection with formal error control, the current options are: model-X knockoffs adapted to Cox models (CoxKnockoff, multiple-knockoff CoxMKF, derandomized knockoffs with k-FWER control, 2025); stability selection with error bounds; and conformalized survival analysis for calibrated, distribution-free lower predictive bounds under censoring. Venet et al. (2011) is the reference showing that random signatures are often "significant". Any new signature must therefore beat random and meta-gene (e.g., proliferation) null distributions, not just a p-value.

### Cited Findings
- **Conformalized survival analysis** (Candès, Lei & Ren, arXiv 2103.09763; JRSSB) wraps any survival model to give covariate-dependent lower predictive bounds.
  - Coverage is guaranteed in finite samples under Type I censoring.
  - Under conditionally independent censoring it is doubly robust: coverage holds approximately if either the censoring model or the survival model is well estimated.
  - It was demonstrated on UK Biobank COVID-19 data.
  - Sources: [arXiv](https://arxiv.org/abs/2103.09763v3); [UK Biobank publication page](https://www.ukbiobank.ac.uk/publications/conformalized-survival-analysis/)
- **Knockoffs for Cox regression:**
  - CoxKnockoff (Li, Yu & Zhao, arXiv 2308.00269) gives finite-sample FDR control — [arXiv summary](https://arxiv.deeppaper.ai/papers/2308.00269v1)
  - Derandomized knockoffs with Cox regression give **k-FWER** control for "confirmatory biomarker identification" (arXiv 2504.03907, 2025) — [arXiv v2](https://arxiv.org/html/2504.03907v2)
  - Multi-knockoff filter for high-dimensional survival omics (CoxMKF) — [Columbia Tech Ventures](https://inventions.techventures.columbia.edu/technologies/multi-knockoff-filter-for-high--CU25376)
  - Robust knockoffs for Cox models with an unknown precision matrix (JUSTC 2024) — [JUSTC](https://justc.ustc.edu.cn/en/article/pdf/preview/10.52396/JUSTC-2024-0075.pdf)
- **[memory, not re-fetched]** Foundational papers:
  - Model-X knockoffs: Candès et al. JRSSB 2018, DOI 10.1111/rssb.12265.
  - Stability selection: Meinshausen & Bühlmann, JRSSB 2010, DOI 10.1111/j.1467-9868.2010.00740.x; complementary-pairs bounds by Shah & Samworth, JRSSB 2013, DOI 10.1111/j.1467-9868.2011.01034.x.
  - Random-signature problem: Venet, Dumont & Detours, PLoS Comput Biol 2011;7(10):e1002240, DOI 10.1371/journal.pcbi.1002240. Most random gene sets were significantly associated with breast-cancer outcome, largely through correlation with proliferation.
- Bayesian decision analysis is available for net benefit under censoring (bayesDCA) — [arXiv 2308.02067](https://arxiv.org/pdf/2308.02067)

### Inferences
- Concrete guards for the successor:
  1. **Random-signature null:** compare the candidate's ΔC against ΔC from size-matched random gene sets, and from sets matched on correlation with a proliferation/meta-PCNA score, drawn from the same expression matrix. CcSinglecell's block-permutation p-value tests "signal vs none". This tests "signal vs generic transcriptome signal".
  2. **Selection with error control:** use knockoff-Cox (FDR or k-FWER) or stability selection for feature lists. Then refit a low-dimensional model so the parameter count respects the Riley sample-size calculation.
  3. **Conformal lower bounds on time-to-recurrence** as a patient-level, distribution-free uncertainty output. Report coverage per external cohort, which is also a transportability check.
  4. **Pre-registration of a small number of candidate signatures**, with multiplicity control (e.g., Holm or FDR across signatures × cohorts), to avoid the forking-paths version of Venet's problem.

### Gaps
- No CRC application of knockoff-Cox or conformal survival was found in these searches.
- Newer conformal survival variants (adaptive cut-offs, doubly-robust under general censoring, 2024–2025) were not verified this session.
- No head-to-head benchmark of these methods on omics survival data was found.

## Q6. Causal and genetic validation of candidate biomarkers

### Takeaway
Large pQTL resources make proteome-wide Mendelian randomization with colocalization routine: UKB-PPP has 54,219 participants and 2,923 Olink proteins; deCODE adds SomaScan proteins. CRC proteome-MR has prioritized proteins such as GREM1, CLSTN3, CSF2RA, CD86 and ANGPTL4. These are almost entirely studies of **CRC risk (incidence)**, not prognosis. For prognostic or predictive biomarkers, MR evidence is weak (index-event/collider bias, small survival GWAS), so perturbation screens and treatment-randomized cohorts are the more relevant orthogonal evidence.

### Cited Findings
- **UKB-PPP** (Sun et al., Nature 2023;622:329–338, DOI 10.1038/s41586-023-06592-6): 54,219 UKB participants, Olink Explore (~3,000 proteins), pQTLs for 2,923 proteins. About 14,287 primary associations (81% novel per Olink). Includes non-European pQTL mapping; data available to UKB researchers — [Nature PDF](https://www.nature.com/articles/s41586-023-06592-6.pdf); [UK Biobank](https://www.ukbiobank.ac.uk/publications/plasma-proteomic-associations-with-genetics-and-health-in-the-uk-biobank/); [Olink](https://olink.com/application/population-scale-proteogenomics/new-truths-discovered-in-worlds-largest-proteogenomic-study)
- **CRC proteome-wide MR** (Sun et al., Genome Med 2023):
  - Used pQTLs from 7 proteome GWAS (4,853 proteins); outcome data from 16,871 cases / 26,328 controls plus FinnGen and UKB (9,276 cases).
  - 13 proteins were associated with CRC risk.
  - After colocalization and SMR/HEIDI, GREM1, CLSTN3, CSF2RA and CD86 were prioritized.
  - Sources: [PMC10508028](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10508028/); [PubMed 37726845](https://pubmed.ncbi.nlm.nih.gov/37726845/)
- JNCI 2025: cis-MR + colocalization suggested that genetically proxied ANGPTL4 inhibition lowers CRC risk, replicated in EPIC and UKB. Seen via a search summary; the linked record is Oxford NDPH — [NDPH publication](https://www.ndph.ox.ac.uk/publications/publication_modal/2130523)
- MR applied to *outcome within a trial*: a CALGB/SWOG 80405 preprint ran MR on RNA-seq and germline data (n≈1,284) for OS. It is not peer reviewed and uses a liberal p<0.1 — [Research Square](https://www.researchsquare.com/article/rs-3673588/v1)
- **[memory, not re-fetched]**
  - deCODE plasma proteome pQTLs: Ferkingstad et al., Nat Genet 2021, DOI 10.1038/s41588-021-00978-w (~35,559 Icelanders; ~4,907 SomaScan aptamers).
  - Index-event bias in GWAS of disease progression/survival: Dudbridge et al., Nat Commun 2019, DOI 10.1038/s41467-019-09381-w.
  - Genome-scale Perturb-seq: Replogle et al., Cell 2022, DOI 10.1016/j.cell.2022.05.013.

### Inferences
- Use MR and colocalization as a **filter for risk-biology biomarkers** (e.g., plasma proteins for screening or early detection), not as validation of a prognostic tissue signature. For prognosis, germline instruments rarely explain tumor-intrinsic expression, and survival GWAS among cases suffer index-event bias.
- Better orthogonal evidence for a prognostic or predictive gene module:
  - (a) CRISPR/DepMap dependency or drug-sensitivity association in CRC lines and organoids;
  - (b) Perturb-seq showing the module responds to a putative driver;
  - (c) a **treatment × biomarker interaction** in randomized cohorts (CALGB 80405 CMS, PARADIGM ctDNA, IDEA/Immunoscore).
  - For the successor this could be an optional "evidence ladder" score, but the primary claim should rest on external validation in Q3 cohorts.

### Gaps
- No 2024–2026 MR study of CRC *survival or recurrence* using UKB-PPP was found.
- No CRC-specific Perturb-seq or CRISPR study linking a prognostic signature to function was retrieved.
- Whether UKB-PPP was used in any CRC prognosis MR remains unverified.
