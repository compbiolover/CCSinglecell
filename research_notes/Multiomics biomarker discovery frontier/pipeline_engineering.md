# Pipeline engineering for a 2026 multi-omics biomarker-discovery successor to CcSinglecell

Scope: software, data access, workflow and QA ecosystem as of 2026-11-07. Research method: ~40 web searches/fetches, preferring official docs, release notes, GitHub and peer-reviewed papers. Where search snippets were the only evidence, or where sources disagree, that is flagged inline. Version numbers move fast; every pin below should be re-checked at repo-creation time.

## 1. R vs Python vs hybrid in 2026 (scverse vs Bioconductor, interop, what consortia standardize on)

### Takeaway
The centre of gravity for single-cell, spatial and imaging is Python/scverse on AnnData/MuData/SpatialData files (h5ad and Zarr), and the consortia that matter (HTAN, CELLxGENE/CZI) standardize on h5ad that follows the CELLxGENE schema. Bioconductor is still strong and current (3.23, R 4.6, April 2026) for bulk RNA statistics, and since 2025–2026 R can read and write AnnData natively through anndataR without embedding Python. So a hybrid is practical: Python-primary, with R steps that exchange data through files rather than through in-process reticulate bridges.

### Cited Findings
**scverse core (Python)**
- anndata 0.12.0 was released 2025-07-16 and announced 2025-07-29. Its focus was Zarr v3 support, full lazy loading (`read_lazy` / `read_elem_lazy` replace `read_elem_as_dask`), API extensibility and obstore for remote access, while staying fully compatible with Zarr v2. — [scverse blog](https://scverse.org/blog/2025-anndata-012/)
- The 0.12 line was still getting patches in 2026 (0.12.8 on 2026-01-27; 0.12.19 on 2026-07-02). anndata 0.13 followed, moving X into layers and making Zarr v3 the default. The stable docs are titled 0.13.3.post0. — [anndata release notes](https://anndata.scverse.org/en/stable/release-notes/)
- Latest release dates on scverse's PyPI organisation page: spatialdata and spatialdata-io on 2026-07-02, mudata on 2026-08-17, scvi-tools on 2026-07-09, rapids-singlecell on 2026-07-27. — [PyPI scverse org](https://pypi.org/org/scverse/)
- scvi-tools changelog: 1.4.3 (2026-05-12) added `scvi.external.DIAGVI`; 1.5.0 (2026-07-08) added rapids-singlecell support; 1.5.1 (2026-09-10) fixed MuData saving with AutotuneExperiment. — [scvi-tools changelog](https://docs.scvi-tools.org/en/latest/changelog.html). Conflict: the GitHub Releases page reportedly showed 1.4.3 as the latest tag. — [scvi-tools releases](https://github.com/scverse/scvi-tools/releases)
- decoupler is now an scverse package, with 2.1.1 the latest build seen. — [PyPI decoupler](https://pypi.org/p/decoupler). LIANA+ 1.7.3 was released 2026-05-26 and requires `decoupler>=2.0.0`. — [LIANA+ changelog](https://liana-py.readthedocs.io/en/latest/changelog.html)
- MOFA in Python: mofapy2 is at 0.7.2. — [PyPI mofapy2](https://pypi.org/project/mofapy2/). muon provides a `muon.tl.mofa()` interface on MuData. — [muon changelog](https://muon.readthedocs.io/en/stable/changelog.html). The R package MOFA2 is version 1.22.1 in Bioconductor 3.23. — [Bioconductor MOFA2](https://www.bioconductor.org/packages/release/bioc/html/MOFA2.html)

**Bioconductor (R)**
- Bioconductor 3.22 was released 2025-10-30 with 2,361 software packages, 435 experiment-data, 926 annotation, 29 workflows and 6 books, for R 4.5, including Linux arm64 and macOS arm64 builds. — [Bioc 3.22 announcement](https://bioconductor.org/news/bioc_3_22_release/)
- Bioconductor 3.23 was released 2026-04-29 with 2,418 software packages on R 4.6. It is the current release; on the usual six-month schedule 3.24 would follow in late October 2026, but I did not verify that it happened. — [Bioconductor release announcements](https://bioconductor.org/about/release-announcements)

**Interop**
- anndataR is a pure-R implementation that reads and writes H5AD and Zarr natively with no Python dependency, converts to and from SingleCellExperiment and Seurat, and is on Bioconductor (R 4.5+). Published in *Bioinformatics* 42(6), btag288, 2026. — [Bioinformatics paper](https://academic.oup.com/bioinformatics/article/42/6/btag288/8675377); [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC13242924/)
- anndataR 0.2.0 (2025-08-01) supports Python anndata 0.12.0 and has an option to stop writing NULLs so that files stay readable by anndata < 0.12. — [Zenodo anndataR 0.2.0](https://zenodo.org/records/16678013)
- zellkonverter, the reticulate+basilisk route, now defaults to an anndata 0.12.3 environment and prints a message recommending anndataR. — [zellkonverter releases](https://github.com/theislab/zellkonverter/releases)
- The anndata docs list the older R `anndata` (reticulate wrapper) as superseded by anndataR. — [anndata interoperability](https://anndata.readthedocs.io/en/latest/interoperability.html)
- Bioconductor's OSTA (spatial transcriptomics) book has a Python-interoperability chapter covering reticulate, basilisk and anndataR. — [OSTA ch. 8](https://bioconductor.org/books/release/OSTA/pages/bkg-python-interoperability.html)

**Spatial and imaging formats**
- OME-Zarr 0.5.0 (2024-11-21) moved OME-NGFF onto the Zarr v3 spec; patches 0.5.1 and 0.5.2 came in January 2025. — [OME-NGFF 0.5 spec](https://ngff.openmicroscopy.org/0.5/)
- SpatialData's on-disk format is a Zarr container that follows OME-NGFF where possible. By design it is not compatible with OME-Zarr 0.4 and targets what 0.5 should be (coordinate transformations). Its roadmap still lists Zarr v3 sharding and moving to ome-zarr-models-py as open items. — [SpatialData design doc](https://spatialdata.scverse.org/en/latest/design_doc.html)
- spatialdata 0.7.3 writes Zarr v3 stores (container format 0.2, images in an "OME 0.5-dev" layout). Stores at format 0.1 are Zarr v2, so readers must check the version. An R reader/writer (spatialdataR) is in development. — [spatialdataR issue #115](https://github.com/HelenaLC/spatialdataR/issues/115)
- The vizarr viewer cannot open SpatialData stores (open issue). — [vizarr #169](https://github.com/BioNGFF/vizarr/issues/169)

**What consortia standardize on**
- HTAN (Nature Methods 2025) puts open-access processed Level 3 and 4 data on Synapse, controlled lower-level data in the NCI CRDC, and imaging in original formats such as OME-TIFF through the CDS. The data model is published as CSV, JSON-LD, JSONSchema and YAML. — [HTAN Nat Methods 2025](https://www.nature.com/articles/s41592-025-02643-0); [bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2024.06.25.598921v1.full)
- HTAN Phase 2 sc/snRNA-seq Level 3/4 files are h5ad modelled on CELLxGENE's schema requirements, in AnnData 0.10 format, and are checked with the HTAN-h5ad-validator PyPI package. — [HTAN sc submission docs](https://docs.humantumoratlas.org/data_submission/scrnaseq_data_submission/); [HTAN portal docs](https://docs.humantumoratlas.org/data_access/portal/)
- CZI's CELLxGENE Census is distributed as TileDB-SOMA with Python and R APIs (`cellxgene_census` / `cellxgene.census`), using dataset schema 7.0.0 and Census schema 2.4.0 for the 2025-11-08 LTS. — [Census data releases](https://chanzuckerberg.github.io/cellxgene-census/cellxgene_census_docsite_data_release_info.html)

### Inferences
- Recommended stack: a Python-primary analysis layer (scanpy + anndata ≥0.12 + MuData + SpatialData + scvi-tools + decoupler 2.x + LIANA+), with R kept for what it does best (DESeq2/edgeR/limma for bulk; MOFA2 or mixOmics if needed). Hand data across the language boundary as h5ad or Zarr files that anndataR and zellkonverter can read. Avoid CcSinglecell's current pattern of driving mofapy2 from R through reticulate: on-disk exchange removes a whole class of environment-coupling bugs.
- Write h5ad to the CELLxGENE/HTAN schema, validated with HTAN-h5ad-validator, so that outputs are deposit-ready and Census/HTAN inputs need no conversion.
- Pin anndata to one minor line and test that both directions (anndataR read/write) work in CI. 0.12 vs 0.13 changes (X moved into layers, Zarr v3 default) and spatialdata format 0.1 vs 0.2 are the most likely silent breakages.
- MOFA can now run entirely in Python on MuData via muon + mofapy2, which removes the R→Python bridge CcSinglecell needs today.

### Gaps
- I found no recent primary source on CellRank 2 release status, nor on MuData/muon release cadence beyond PyPI dates.
- I did not verify whether Bioconductor 3.24 (expected ~late October 2026) has shipped.
- Whether TileDB-SOMA itself adopts Zarr v3 or OME-Zarr was not answered by any source found.
- HCA's current processed-matrix format standard was not researched directly. It is reasonable to assume it aligns with CELLxGENE h5ad, but that is unsourced.

## 2. Data access: which sources give reproducible, versioned snapshots

### Takeaway
Only some sources give citable, immutable versions: CELLxGENE Census LTS builds (kept at least 5 years), ISB-CGC BigQuery `TCGA_versioned` tables keyed by GDC release number, and Bioconductor-versioned packages such as recount3. Raw GDC API pulls through TCGAbiolinks are not snapshots. The GDC changes data in place across releases (DR42–DR46 all touched TCGA clinical, slide or COAD data), and TCGAbiolinks keeps breaking when the GDC API changes. The successor should therefore pin a GDC data-release number and record file UUIDs and md5s, or use ISB-CGC versioned tables. The CPTAC Python package looks stale, and CcSinglecell's hoardeR dependency (used for TargetScan) is failing CRAN checks.

### Cited Findings
**GDC / TCGA**
- Recent GDC data releases: v42.0 (2025-01-30), v43.0 (2025-05-07), v44.0 (2025-10-29), v45.0 (2025-12-04), v46.0 (2026-08-10). — [GDC Data Release Notes](https://docs.gdc.cancer.gov/Data/Release_Notes/Data_Release_Notes/)
- Changes in those releases that matter for a TCGA-COAD/READ survival model:
  - DR42 made TCGA clinical fields that were previously only in supplements queryable.
  - DR43 migrated survival-plot data back to `diagnoses.days_to_last_follow_up`, fixing a DR42 issue.
  - DR46 removed ensemble MAFs from 397 contaminated TCGA alignments and noted one TCGA-COAD methylation aliquot still pending.
  - DR44 says some replaced MuTect2 VCFs "will be versioned in a future release".
  - Source: [GDC Data Release Notes](https://docs.gdc.cancer.gov/Data/Release_Notes/Data_Release_Notes/)
- The GDC API release notes record that data migrated into `follow_ups.days_to_follow_up` for some TCGA projects in DR42, which can change clinical parsing. — [GDC API release notes](https://docs.gdc.cancer.gov/API/Release_Notes/API_Release_Notes)
- DR43 (May 2025) reportedly introduced new CNV category labels ("gain/amplification", "heterozygous/homozygous deletion"), which can break code that filters on the old labels. — [CRDC updates June 2025](https://datacommons.cancer.gov/news/crdc-components-updates-june-2025) (from a search snippet)
- DR45 added an RC-PTCL project, Burkitt lymphoma slides and a Mutect2 fix for chromosomes 10 and 20. — [GDC DR45 news](https://gdc.cancer.gov/gdc-releases-new-data-dr45). A search snippet said GDC is planning DR47 for January 2027; I did not verify this on a primary page. — [GDC access-data tag page](https://gdc.cancer.gov/subject-tag/access-data)
- TCGAbiolinks remains on Bioconductor. — [Bioconductor TCGAbiolinks](https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html). Its issue tracker has 2025 reports including "Issues connecting to GDC API" (opened 2025-11-07), clinical-download errors, and a missing `disease_response` column. — [TCGAbiolinks issues](https://github.com/BioinformaticsFMRP/TCGAbiolinks/issues). An older documented failure mode is GDCquery parsing errors when the GDC `projects` TSV endpoint changed. — [Bioc support 102119](https://support.bioconductor.org/p/102119); [issue #512](https://github.com/BioinformaticsFMRP/TCGAbiolinks/issues/512)

**ISB-CGC (BigQuery)**
- `isb-cgc-bq.TCGA` holds the latest tables; `isb-cgc-bq.TCGA_versioned` holds all prior releases plus the current one, with the GDC release number in the table name (e.g. `clinical_gdc_r24`, `RNAseq_hg38_gdc_r40`). — [ISB-CGC TCGA docs](https://isb-cancer-genomics-cloud.readthedocs.io/en/latest/sections/data/TCGA_top.html); [ISB-CGC data release notes](https://isb-cancer-genomics-cloud.readthedocs.io/en/latest/sections/ReleaseNotes/ISB-CGCDataReleases.html)
- ISB-CGC loaded GDC r42-derived data in May 2025 and r43 in June 2025. The legacy `isb-cgc.TCGA_*_v0` datasets are deprecated but still exist. — [ISB-CGC data release notes](https://isb-cancer-genomics-cloud.readthedocs.io/en/latest/sections/ReleaseNotes/ISB-CGCDataReleases.html)

**cBioPortal / UCSC Xena**
- cBioPortal lists both TCGA PanCancer Atlas studies and newer "TCGA GDC, 2025" study versions. — [cBioPortal](https://www.cbioportal.org/)
- UCSC Xena hosts GDC-derived cohorts and a curated Pan-Cancer Atlas hub, accessible from R through UCSCXenaTools. — [Xena TCGA hub docs](https://ucsc-xena.gitbook.io/project/public-data-we-host/tcga); [UCSCXenaTools](https://shixiangwang.github.io/home/en/tools/ucscxenatools-intro/)
- PanCancer Atlas data are "similar but not the same" as the original TCGA marker-paper data, which had extra filtering, so results depend on which source is used. — reported in search results alongside [NCI BTEP TCGA download guide](https://bioinformatics.ccr.cancer.gov/docs/btep-coding-club/CC2024/TCGA/TCGA_download/) (attribution from snippet; verify)

**CELLxGENE Census**
- LTS builds: 2023-05-15, 2023-07-25, 2023-12-15, 2024-07-01, 2025-01-30 and 2025-11-08. The 2025-11-08 LTS has 1,845 datasets, 162,025,130 human cells (99.6M unique) and 46.3M mouse cells.
- LTS builds stay public for at least 5 years; weekly builds for 1 month.
- Pin with `open_soma(census_version="2025-11-08")`. `"stable"` resolves to the latest LTS, which is a moving target.
- The 2025-11-08 LTS is compatible with package 1.17.x.
- Source: [Census data releases](https://chanzuckerberg.github.io/cellxgene-census/cellxgene_census_docsite_data_release_info.html)
- GitHub shows cellxgene-census v1.18.0 as the latest package release (dated 24 Jun, presumably 2026), which conflicts with the docs' "1.17.x" statement. — [cellxgene-census releases](https://github.com/chanzuckerberg/cellxgene-census/releases)
- The 2023-05-15 LTS has an erratum: 243,569 cells are duplicated as `is_primary_data=True`. — [Census data releases](https://chanzuckerberg.github.io/cellxgene-census/cellxgene_census_docsite_data_release_info.html)

**HTAN**
- Open-access Level 3/4 data on Synapse (Level 4 = h5ad); controlled data in CRDC; imaging (OME-TIFF) through CDS. — [HTAN Nat Methods 2025](https://www.nature.com/articles/s41592-025-02643-0); [HTAN portal docs](https://docs.humantumoratlas.org/data_access/portal/). An unofficial HTAN metadata dataset exists on HuggingFace. — [seandavis/htan](https://huggingface.co/datasets/seandavis/htan)

**CPTAC / PDC (proteomics)**
- The PayneLab `cptac` Python package exposes per-cancer dataset objects such as `cptac.Colon()`. Its latest PyPI release is 1.5.14 (2024-06-06), with nothing newer found. — [PyPI cptac](https://pypi.org/project/cptac/); [GitHub PayneLab/cptac](https://github.com/PayneLab/cptac)
- CPTAC proteomic data itself is served from the Proteomic Data Commons (PDC). — [TCIA CPTAC wiki](https://wiki.cancerimagingarchive.net/display/Public/CPTAC+Imaging+Proteomics)

**Uniformly processed bulk RNA compendia**
- recount3 has over 750,000 uniformly processed human and mouse RNA-seq samples. — [recount3 Genome Biology 2021 (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8628444). The Bioconductor package is at 1.20.0 (2025). — [recount3 manual](https://free.nchc.org.tw/bioconductor/packages/release/bioc/manuals/recount3/man/recount3.pdf)
- ARCHS4 has continuing NIH funding (2024 project record), but I found no 2025–2026 version release. — [NIH grant record](https://api.ai-analytics.org/grant/nih/10909127.md). DEE2 claims more than 1.3M datasets as of January 2025, a self-reported comparison against ARCHS4 and recount3. — [DEE2 blog](https://genomespot.blogspot.com/2025/01/dee2-2025-updates-growth-and-development.html)

**miRNA target scraping (current CcSinglecell dependency)**
- hoardeR 0.10 (2024-04-08) queries TargetScan through targetscan.org. CRAN mirrors show "issues need fixing before 2026-03-30 / 2026-04-13" notices, meaning it was failing CRAN checks; archival was not confirmed. — [CRAN hoardeR](https://cran.r-project.org/web/packages/hoardeR/); [Duke CRAN mirror](https://archive.linux.duke.edu/cran/web/packages/hoardeR/index.html)

### Inferences
- Reproducibility recipe for TCGA-COAD/READ:
  - Record the GDC data-release number, the file UUIDs and the md5s in a manifest committed to the repo.
  - Download with gdc-client or the GDC API against that manifest.
  - Optionally cross-check against ISB-CGC `TCGA_versioned.*_gdc_rNN` tables, which are the only truly frozen, queryable TCGA snapshots found.
  - Treat TCGAbiolinks as a convenience wrapper, not a reproducibility boundary.
- For clinical survival endpoints, prefer the curated TCGA Pan-Cancer Clinical Data Resource (via Xena or cBioPortal PanCancer Atlas) over raw GDC fields, given the DR42/DR43 back-and-forth on follow-up fields. The specific CDR resource was not separately sourced here.
- Pin the Census by explicit LTS date (2025-11-08), never `"stable"`.
- For proteomics, the stale `cptac` package is a risk. Plan to pull from PDC directly and freeze files with checksums.
- Replace hoardeR web-scraping with a downloaded, versioned TargetScan release file, or a Bioconductor miRNA-target resource, checksummed in the data manifest. Live-website scraping is inherently irreproducible.

### Gaps
- No current sources gathered on curatedTCGAData versioning, GEOquery status, the gdc-client version, or the GDC's policy on accessing superseded file versions.
- I did not confirm the date of DR47 or whether a DR47 exists before 2027.
- I found no ARCHS4 v2.x release information.
- The HTAN Phase 2 data-release cadence and versioning were not established.
- I did not determine whether `cptac` is formally unmaintained or merely quiet.

## 3. Workflow managers, pipelines, containers and environment pinning

### Takeaway
Nextflow plus nf-core remains the default for heavy upstream processing, and 2026 brought big changes: strict syntax is the default from Nextflow 26.04, nf-core/tools 4.0.0 (April 2026) requires Nextflow ≥25.10.4, and modules are migrating from Biocontainers to Seqera Containers. Snakemake 9 is actively maintained (9.27.0, September 2026) and fits a Python-centric downstream analysis layer. Flyte, Metaflow and Dagster are ML/data-engineering orchestrators, not genomics-native. Pixi is emerging as the lockfile manager that spans conda-forge, bioconda and PyPI (through uv).

### Cited Findings
**Nextflow / nf-core**
- Nextflow's strict-syntax parser is opt-in in 25.04/25.10 and default from 26.04, with `NXF_SYNTAX_PARSER=v1` to opt out. — [Nextflow strict syntax](https://nextflow.io/docs/latest/strict-syntax.html); [Migrating to 25.10](https://nextflow.io/docs/latest/migrations/25-10.html)
- nf-core/tools 4.0.0 (2026-04-28) makes the template strict-syntax compliant and sets `nextflowVersion = '!>=25.10.4'`. It also:
  - adds a `nextflow lint` pre-commit hook;
  - makes version topic channels mandatory for module linting;
  - auto-generates Docker, Singularity/Apptainer and Conda container configs for amd64 and arm64 from module `meta.yml`;
  - removes prefix-free commands and webhook notifications (now Nextflow plugins);
  - prepares for the Seqera-containers transition.
  - Source: [nf-core/tools 4.0.0 blog](https://nf-co.re/blog/2026/tools-4_0_0)
- The nf-core roadmap makes strict syntax mandatory in Q2 2026 and adds static types and records to the template in Q4 2026. — [nf-core syntax roadmap](https://nf-co.re/blog/2025/nextflow_syntax_nf-core_roadmap)
- nf-core/rnaseq PR #1680 raised the floor to Nextflow 25.10.4 with nf-schema 2.6.1, validated against 26.04.0. — [rnaseq #1680](https://github.com/nf-core/rnaseq/pull/1680). A known friction: boolean CLI params can be parsed as strings under the strict parser after moving to template 4.0.2. — [nf-core Slack archive](https://www.linen.dev/s/nf-core/t/33296025/nextflow-v26-04-0-uses-the-strict-syntax-parser-by-default-a)
- Latest Nextflow stable is reported as 26.10.3 (2026-10-20). — [Wikipedia](https://en.wikipedia.org/wiki/Nextflow) (not verified on nextflow.io). This is consistent with nf-prov 1.8.0 (2026-10-13) requiring Nextflow ≥26.10.0. — [nf-prov registry](https://registry.nextflow.io/plugins/nf-prov)
- nf-core lists 160 pipelines. — [nf-co.re/pipelines](https://nf-co.re/pipelines)
- Pipeline status for this project:
  - differentialabundance: v2.0.0 (2026-06-23). — [releases](https://github.com/nf-core/differentialabundance/releases)
  - scrnaseq: 4.2.0 per a third-party catalogue. — [Cirro catalogue](https://docs.cirro.bio/pipelines/nf-core-community-curated-nextflow-pipelines/)
  - methylseq: sources conflict; GitHub lists 4.2.0 among releases, while third-party catalogues lag. — [methylseq releases](https://github.com/nf-core/methylseq/releases)
  - smrnaseq: only 2.4.1 (about two years old) surfaced. — [smrnaseq release stats](https://nf-co.re/smrnaseq/releases_stats)
  - spatialvi: "not yet released" (dev only). — [spatialvi dev](https://nf-co.re/spatialvi/dev/)
  - nf-core/sopa is listed for spatial omics. — [Cirro catalogue](https://docs.cirro.bio/pipelines/nf-core-community-curated-nextflow-pipelines/)
- Containers: nf-core/modules issue #13143 tracks migrating the remaining Biocontainers modules to Seqera Containers (`community.wave.seqera.io`), tied to the October 2026 hackathon. Conversions use `nf-core modules container create` and were merging in early November 2026. Several PRs report Wave linux/arm64 build failures (amd64-only images), plus a missing conda lock file in one amd64 build. — [modules #13143](https://github.com/nf-core/modules/issues/13143); [admixture PR #13529](https://github.com/nf-core/modules/pull/13529); [seqtk PR #13210](https://github.com/nf-core/modules/pull/13210)

**Snakemake**
- Snakemake 9.0.0 (2025-03-14) was a major release with breaking changes. It added `--replace-workflow-config`, dynamic module names and IOCache save/reload. — [Snakemake changelog v9.2.1](https://snakemake.readthedocs.io/en/v9.2.1/project_info/history.html)
- The stable changelog is at 9.27.0 (2026-09-11), including diffs shown in generated unit tests and dry-run validation fixes. — [Snakemake stable changelog](https://snakemake.readthedocs.io/en/stable/project_info/history.html)

**Flyte / Metaflow / Dagster (opinion sources)**
- A common pattern is NGS preprocessing in Nextflow or Snakemake with downstream modelling in Flyte, Dagster or Ray. — [Tracer review (vendor marketing)](https://tracer-marketing-website-v3.vercel.app/resources/a-comparative-review-of-bioinformatics)
- Flyte is recommended for GPU-heavy ML on Kubernetes and Metaflow for notebook-to-cloud data-science teams. — [agentpedia 2026](https://agentpedia.codes/blog/ai-ml-orchestration-tools-2026)
- Dagster is Python-centric with weak support for other languages. — [Deep Origin (2022)](https://deeporigin.com/blog/why-are-bioinformatics-workflows-different). It is good for chaining pipelines but a poor choice for CLI-heavy bioinformatics. — [Seqera podcast ep. 56](https://seqera.io/podcasts/episode-56-pipeline-chaining-meta-pipelines-part-1/)
- No head-to-head benchmarks were found.

**Environment pinning**
- Pixi is a lockfile-first manager on the conda ecosystem (`pixi.toml` + `pixi.lock` committed to git, multi-environment), presented at SciPy 2025. — [SciPy 2025 talk](https://cfp.scipy.org/scipy2025/talk/GJACX9/)
- Pixi resolves PyPI dependencies through uv, and maps `uv lock` to `pixi lock`. — [Pixi: switching from uv](https://pixi.prefix.dev/latest/switching_from/uv/)
- Bioinformatics tutorials use conda-forge first (including R) and bioconda second. — [corebio pixi tutorial](https://corebio.info/tutorials/pixi/); [NGS101 pixi scRNA-seq setup](https://ngs101.com/setting-up-single-cell-rna-seq-analysis-environment-with-pixi-10x-faster-setup-zero-version-conflicts/)

### Inferences
- Recommended layering:
  1. **Upstream (FASTQ → counts):** nf-core pipelines run unmodified with pinned `-r` revisions: rnaseq, scrnaseq, methylseq; smrnaseq with caution because of its stale release.
  2. **Data acquisition:** a small, explicit downloader driven by a manifest.
  3. **Analysis/modelling layer:** Snakemake 9 or Nextflow, calling Python/R scripts in one pixi environment. Pixi can hold both R (conda-forge r-base plus bioconda Bioconductor packages) and Python, with uv handling PyPI.
- If the team adopts Nextflow for its own steps, start directly on strict syntax with Nextflow ≥25.10.4 (ideally 26.10.x) and the nf-core/tools 4.x template, to avoid a migration in 2027.
- The successor does not need Flyte, Metaflow or Dagster. Their value (Kubernetes ML orchestration, data assets) outweighs their cost only at platform scale.
- For containers, use the Seqera Containers or Wave image of each nf-core module, plus one project image built from `pixi.lock`. Expect amd64-only images for some tools (arm64 Wave failures), which matters on Apple Silicon and Graviton.
- renv vs pixi: pixi can replace renv when all R packages come from conda-forge/bioconda. Keep renv only if CRAN-only or GitHub-only R packages are needed. This is an inference, not sourced.

### Gaps
- No primary sources gathered on conda-lock's status, Apptainer versions, or Seqera's Wave pricing/rate limits for public use.
- Current release versions of nf-core/rnaseq, scrnaseq, methylseq and smrnaseq were not verified on nf-co.re: the pipelines page fetch was truncated, and third-party catalogues conflict.
- I could not pin a URL for informal reports that installing Bioconductor packages through pixi sometimes fails to resolve R dependencies. Test it early.
- Nextflow 26.10.3 is from Wikipedia only.

## 4. Experiment tracking, data versioning and leakage-safe model selection for survival

### Takeaway
MLflow 3 is the mature, self-hostable default for tracking (local file store, no SaaS). DVC still works but changed stewardship in November 2025 (lakeFS acquired it from Iterative.ai), which is a mild continuity risk. For leakage-free nested cross-validation of survival models, scikit-survival (sklearn-compatible, so Pipeline + GridSearchCV inside an outer CV) and mlr3proba (native nested resampling) are the two most composable frameworks. TorchSurv (Novartis/FDA, JOSS 2024, MIT) is the lightweight option for deep survival models.

### Cited Findings
- MLflow 3 positions itself around experiment tracking, observability and evaluation, adding GenAI features (prompt optimisation, feedback tracking). Models are first-class entities and `log_model` no longer needs an active run. — [MLflow 3 docs](https://mlflow.org/docs/latest/genai/mlflow-3); [MLflow 3.1.3 docs](https://mlflow.org/docs/3.1.3/genai/mlflow-3)
- A third-party source mentions MLflow 3.10 adding multi-workspace support (unverified on mlflow.org). — [letsdatascience](https://letsdatascience.com/blog/mlflow-experiment-tracking-and-ml-lifecycle-management)
- Managed MLflow 3 is available on Amazon SageMaker AI and Databricks. — [AWS blog](https://aws.amazon.com/blogs/machine-learning/accelerating-generative-ai-development-with-fully-managed-mlflow-3-0-on-amazon-sagemaker-ai); [Databricks MLflow](https://docs.databricks.com/en/mlflow)
- lakeFS acquired the DVC open-source project from Iterative.ai (announced November 2025). DVC remains open source at dvc.org, positioned for smaller datasets, while lakeFS targets enterprise scale. — [DVC: joins lakeFS FAQ](https://dvc.org/blog/dvc-joins-lakefs-your-questions-answered/); [lakeFS blog](https://lakefs.io/blog/celebration-shared-vision-lakefs-dvc/). This followed a $20M lakeFS round in July 2025. — [Dealroom](https://app.dealroom.co/news/feed/lakefs-acquires-dvc-raises-20m)
- scikit-survival follows the scikit-learn API, so standard sklearn CV and model-selection tooling (Pipelines, GridSearchCV, nested CV) applies directly. — [scikit-survival GitHub](https://github.com/sebp/scikit-survival); [JMLR paper](http://jmlr.org/papers/volume21/20-729/20-729.pdf)
- mlr3proba supports nested resampling and tuning on any survival measure. A tutorial recommends 3-fold inner CV with 60-iteration random search for real analyses. — [mlr3proba paper (arXiv)](https://arxiv.org/pdf/2008.08080); [mlr3gallery survival networks](https://mlr3gallery.mlr-org.com/posts/2021-04-20-survival-networks/)
- The R package survcompare 0.3.0 (2025-06-25) runs repeated nested CV comparing Cox-PH, Cox-Lasso, survival random forests and their ensemble. — [survcompare manual](https://archive.linux.duke.edu/cran/web/packages/survcompare/refman/survcompare.html)
- TorchSurv was developed by Novartis with the US FDA. It allows custom PyTorch survival models without imposing a parametric form, was published in JOSS in 2024, is MIT-licensed, and has latest release 0.2.0 per ecosyste.ms. — [JOSS paper](https://theoj.org/joss-papers/joss.07341/10.21105.joss.07341.pdf); [arXiv 2404.10761](https://arxiv.org/pdf/2404.10761); [ecosyste.ms](https://science.ecosyste.ms/projects/338)

### Inferences
- CcSinglecell's original design selected features (DE genes, pseudotime-associated genes, miRNA targets) on the full cohort before survival modelling. The successor should move all data-driven feature selection (DE filtering, MOFA factor fitting, glmnet λ) inside the outer CV loop. With sklearn this means an sksurv `Pipeline`, possibly with custom transformers, wrapped in an outer `cross_validate`. With mlr3 it means `AutoTuner` inside `resample`.
- Simplest tracking setup: MLflow with a local `mlruns/` file store or SQLite backend, logging params, fold IDs, metrics (Harrell's/Uno's C, integrated Brier score) and the data-manifest hash. Use DVC only if large intermediate artifacts must be versioned in git. Otherwise a checksummed manifest plus Zenodo deposit of frozen inputs is lighter and has no vendor risk.
- An external validation cohort (e.g. a GEO CRC cohort, or CPTAC-COAD for proteomics) should be kept completely out of model selection. Nested CV alone is not a substitute.

### Gaps
- No sources gathered on Weights & Biases' 2025–2026 ownership or licensing status, pycox maintenance status, lifelines' latest release, or tidymodels `censored` versions. These need verification before being listed as options.
- The exact current MLflow version was not confirmed.

## 5. GPU/compute: single-cell and pathology foundation models, RAPIDS, cloud

### Takeaway
The leading pathology foundation models (UNI2-h, Virchow2, CONCH) are gated on HuggingFace under CC-BY-NC-ND 4.0: academic and non-commercial use only, with restrictions that extend to derived models. This is fine for a university repo but must be documented, and weights cannot be redistributed. Single-cell foundation models (Geneformer, scGPT, and in 2026 benchmarks also UCE, TranscriptFormer and Arc State) repeatedly fail to beat simple baselines zero-shot, so they belong in an optional, benchmarked branch. rapids-singlecell 0.15–0.16 (2026) is the mature GPU path for scanpy-style workflows, and scvi-tools 1.5 integrates with it.

### Cited Findings
**Pathology**
- UNI2-h (MahmoodLab) is licensed CC-BY-NC-ND 4.0, non-commercial academic research only. Commercial use, including of models trained on its outputs, needs prior approval. Access is gated and requests from personal email domains (e.g. Gmail) are denied. — [UNI2-h model card](https://huggingface.co/MahmoodLab/UNI2-h)
- Virchow2 (Paige) is CC-BY-NC-ND 4.0 and gated, and commercial-entity requests are rejected by default. — [Virchow2 model card](https://huggingface.co/paige-ai/Virchow2/blob/main/README.md). A third-party page claims Virchow2 permits commercial use; this contradicts the model card and should be ignored. — [aimodels.fyi (unreliable)](https://www.aimodels.fyi/models/huggingFace/uni2-h-mahmoodlab)
- TRIDENT, MahmoodLab's toolkit for large-scale whole-slide-image processing, lists CONCH under CC-BY-NC-ND-4.0 and notes that several supported encoders need HuggingFace gated-access approval. — [TRIDENT GitHub](https://github.com/mahmoodlab/TRIDENT)

**Single-cell foundation models**
- Genome Biology 2025: zero-shot Geneformer and scGPT underperform simpler methods on cell-type clustering and batch integration. They were chosen because pretrained weights were available. — [Genome Biology 2025](https://link.springer.com/content/pdf/10.1186/s13059-025-03574-x.pdf)
- VCBench (bioRxiv, June 2026) benchmarks Geneformer, scGPT, UCE, TranscriptFormer and Arc State against pre-registered linear and nearest-neighbour baselines. Excerpts show Geneformer beating TranscriptFormer on at least one zero-shot metric. — [VCBench](https://www.biorxiv.org/content/10.64898/2026.06.18.733146v1.full.pdf). A separate July 2026 harmonised benchmark across single-cell and spatial models reports context-dependent generalisation. — [arXiv 2607.17227](https://arxiv.org/pdf/2607.17227)
- scPRINT (Nature Communications 2025) was pre-trained on 50M cells for gene-network inference. — [scPRINT](https://www.nature.com/articles/s41467-025-58699-1)
- CZI's TranscriptFormer is publicly accessible on CZI's Virtual Cells Platform, with code on GitHub. — [Biohub TranscriptFormer overview](https://biohub.org/blog/transcriptformer-model-overview/). CZI and NVIDIA announced a collaboration on virtual-cell models in October 2025. — [News-Medical](https://www.news-medical.net/news/20251028/CZI-and-NVIDIA-collaborate-to-accelerate-life-science-research-through-virtual-cell-models.aspx)

**GPU single-cell**
- rapids-singlecell 0.15.0 (2026-04-30) moved to prebuilt CUDA wheels and compiled (nanobind) kernels. Install now uses `-cu12` or `-cu13` package variants instead of plain `pip install rapids-singlecell`. 0.16.0 (2026-07-16) and 0.16.1 (2026-07-27) followed. — [rapids-singlecell 0.15 blog](https://scverse.org/blog/2026-rsc-goes-nanobind/); [release notes](https://rapids-singlecell.readthedocs.io/en/latest/release-notes/index.html)
- scvi-tools 1.5.0 (2026-07-08) added rapids-singlecell support. — [scvi-tools changelog](https://docs.scvi-tools.org/en/latest/changelog.html)

**Cloud**
- ISB-CGC exposes TCGA/GDC data in Google BigQuery. — [ISB-CGC TCGA docs](https://isb-cancer-genomics-cloud.readthedocs.io/en/latest/sections/data/TCGA_top.html)
- Managed MLflow is available on SageMaker and Databricks (see Q4). — [AWS blog](https://aws.amazon.com/blogs/machine-learning/accelerating-generative-ai-development-with-fully-managed-mlflow-3-0-on-amazon-sagemaker-ai)

### Inferences
- Make the H&E branch "bring your own token": a config key for the HuggingFace token, cached embeddings (not weights) stored as Zarr or Parquet, and a LICENSE note stating that downstream models inherit CC-BY-NC-ND restrictions. Provide a permissively licensed fallback encoder so CI and outside users can run end-to-end. Which open encoder to use is not researched here.
- Treat scFM embeddings as one more candidate feature block inside nested CV, always reported against a PCA/scVI baseline. Current benchmarks do not justify making them the default.
- GPU is needed for WSI feature extraction and optionally for scVI and rapids-singlecell. Everything in the core CRC survival path (bulk RNA, miRNA, methylation, Cox/glmnet) runs on CPU. Keep GPU steps in separate, optional workflow targets with their own pixi environment (CUDA variants).

### Gaps
- I did not obtain licence terms for TranscriptFormer weights, scGPT, Geneformer or UCE, nor the current status of CONCH/UNI successors (e.g. TITAN, CONCH v1.5).
- No sources on permissively licensed pathology encoders, on HuggingFace gating turnaround, or on cloud GPU pricing.
- Arc Institute's State model licence is not established.

## 6. Packaging and QA best practice (monorepo, CI, testing, docs, DOIs, FAIR/RO-Crate)

### Takeaway
The tooling for a FAIR, citable pipeline repo is mature:
- pixi or uv lockfiles;
- nf-test plus `nextflow lint` for Nextflow, and pytest/testthat for code;
- CITATION.cff plus GitHub→Zenodo release integration for versioned DOIs, with a concept DOI for "all versions";
- nf-prov (v1.8.0, October 2026) emitting Workflow Run RO-Crate provenance with no pipeline changes.

anndataR (R) and anndata (Python) show the cross-language pattern: each language has a native implementation of a shared on-disk spec, rather than one language wrapping the other.

### Cited Findings
- The nf-core template ships nf-test (bumped from 0.9.4), a `nextflow lint` pre-commit hook and strict meta/ext key linting. — [nf-core/tools 4.0.0](https://nf-co.re/blog/2026/tools-4_0_0)
- Snakemake can generate unit tests per rule, and 9.x shows diffs on failure. — [Snakemake stable changelog](https://snakemake.readthedocs.io/en/stable/project_info/history.html)
- nf-prov 1.8.0 (2026-10-13, requires Nextflow ≥26.10.0) writes Workflow Run RO-Crate (`ro-crate-metadata.json`, covering the Process, Workflow and Provenance profiles) and BioCompute Objects, configured in `nextflow.config` without changing the pipeline script. — [nf-prov registry](https://registry.nextflow.io/plugins/nf-prov). WRROC has been supported since nf-prov 1.4.0. — [RO-Crate nf-prov page](https://www.researchobject.org/ro-crate/nf-prov). WorkflowHub supports the Workflow Crate profile. — [Workflow Run Crate](https://researchobject.org/workflow-run-crate/). The nf-core/provenancereport pipeline produces RO-Crates with nf-prov. — [provenancereport output](https://nf-co.re/provenancereport/output)
- GitHub→Zenodo: enable the repo in Zenodo and each GitHub release gets an archived DOI. A `CITATION.cff` on the default branch gives GitHub's "Cite this repository" button and supplies Zenodo metadata (otherwise Zenodo uses GitHub API metadata). Zenodo also mints a concept DOI that resolves to the latest release. — [spec-kit PR #2291](https://github.com/github/spec-kit/pull/2291); [CASRAI](https://casrai.org/wp/tag/github/)
- Pixi plus uv covers mixed conda/PyPI lockfiles. — [Pixi: switching from uv](https://pixi.prefix.dev/latest/switching_from/uv/)
- Small-fixture testing of AnnData I/O across languages is exactly the compatibility surface anndataR's paper addresses, with round-trips between H5AD/Zarr and SingleCellExperiment/Seurat. — [anndataR paper](https://academic.oup.com/bioinformatics/article/42/6/btag288/8675377)

### Inferences
- Suggested repo shape (monorepo):
  - `pyproject.toml` package (src layout) for the Python core;
  - an optional thin R package (`DESCRIPTION`) only for R-only statistics;
  - `pixi.toml`/`pixi.lock` at the root defining `default`, `r` and `gpu` environments plus task aliases;
  - `workflow/` (Snakefile or `main.nf`);
  - `tests/fixtures/` containing tiny synthetic h5ad/MuData and TCGA-shaped TSVs, under 1 MB each, with no controlled data.
- GitHub Actions jobs: lint (ruff, lintr, `nextflow lint`); unit tests (pytest, testthat); one end-to-end "stub" workflow run on the fixtures; the anndataR↔anndata round-trip test; and a release job that relies on the Zenodo integration plus CITATION.cff validation.
- Docs: a single Quarto site can render both R and Python notebooks, which avoids running both pkgdown and mkdocs. This is an inference; the Quarto, pkgdown and mkdocs specifics were not sourced.
- Provenance: if Nextflow is used, enable nf-prov `wrroc` and archive the crate with each results release. If Snakemake is used, an equivalent RO-Crate export path needs checking (see Gaps).

### Gaps
- No sources gathered on Quarto, pkgdown or mkdocs-material versions or recommended practice, nor on Snakemake's RO-Crate/provenance support.
- An informal claim that an invalid CITATION.cff makes the Zenodo release fail could not be pinned to a primary URL. Verify against Zenodo's help pages.
- No sources on R+Python monorepo exemplars beyond anndata/anndataR.

## 7. LLM-agent tooling for bioinformatics: mature enough as an optional layer?

### Takeaway
Agents exist and are active but immature. Biomni (Stanford) is at PyPI 0.0.8 (October 2025), supports MCP, and warns that generated code runs with full system privileges. BioMCP (MIT licence, about 13 tools) covers literature, trials and variants rather than TCGA or single-cell. Agent "skills" for CELLxGENE Census and pixi are appearing. These are reasonable as an optional, sandboxed exploration layer (e.g. literature or annotation lookup), but should never run inside the reproducible pipeline path.

### Cited Findings
- Biomni is a general-purpose biomedical agent (bioRxiv, May 2025) with benchmarks on causal gene prioritisation, drug repurposing and rare-disease diagnosis. — [bioRxiv Biomni](https://www.biorxiv.org/content/10.1101/2025.05.30.656746v1). Latest PyPI is 0.0.8 (2025-10-27). — [PyPI biomni](https://pypi.org/project/biomni/)
- The Biomni repo documents loading external MCP servers through config. It warns that executed code runs with full system privileges and that the frozen April 2025 release differs from the web platform. — [Biomni GitHub](https://github.com/snap-stanford/Biomni)
- AWS published a walkthrough exposing Biomni's database tools through Bedrock AgentCore Gateway (MCP). — [AWS blog](https://aws.amazon.com/blogs/machine-learning/build-a-biomedical-research-agent-with-biomni-tools-and-amazon-bedrock-agentcore-gateway/)
- BioMCP (GenomOncology) is an MCP server over ClinicalTrials.gov, PubMed and MyVariant.info, listed from 2025-04-05, MIT licence, about 13 tools. — [PulseMCP listing](https://www.pulsemcp.com/servers/genomoncology-biomcp); [Augment MCP registry](https://www.augmentcode.com/mcp/biomcp)
- Census-specific agent tooling exists: Harvard Zitnik lab's ToolUniverse has CELLxGENE Census tools, and K-Dense publishes a cellxgene-census agent skill (targeting 1.17.x and the 2025-11-08 LTS). — [ToolUniverse Census tools](https://zitniklab.hms.harvard.edu/ToolUniverse/en/tools/cellxgene_census_tools.html); [K-Dense skill](https://github.com/K-Dense-AI/scientific-agent-skills/blob/main/skills/cellxgene-census/SKILL.md). There is also a pixi-based "agent-safe execution" skill for bioinformatics environments. — [pixi-skill](https://github.com/cheahhl814/pixi-skill)
- Evaluation infrastructure is appearing, e.g. BioMedArena (May 2026), an open toolkit for building and evaluating biomedical deep-research agents. — [arXiv 2605.06177](https://arxiv.org/pdf/2605.06177)
- A hosted, separate "Biomni Lab" MCP connector exists and is a different product from the open-source repo. — [Gumloop listing](https://www.gumloop.com/connectors/biomni-lab-mcp)

### Inferences
- Include agents as an optional `agents/` or docs-level integration, not a dependency:
  - an MCP config pointing at BioMCP for literature and variant lookups during biomarker interpretation;
  - a README section on running Biomni in a container with no credentials and read-only data mounts.
- Agents' value in this project is interpretation and curation (literature support for candidate biomarkers, annotation lookups), not executing the frozen pipeline. Reproducibility requires deterministic code paths, which agents don't provide.
- A well-structured repo (pixi tasks, manifest-driven data, CLI entry points) is what makes an agent useful later. That is a better investment than wiring a specific agent now.

### Gaps
- No primary sources found on MCP servers for GDC, cBioPortal or CELLxGENE that are officially maintained by those projects. The Census tools found are third-party.
- No independent evaluation of Biomni's error rate on TCGA-style analyses was found.
- Biomni's 2026 release status beyond 0.0.8 is unknown.
