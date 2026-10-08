# Dog-human gene network colocalization links synaptic signalling and ubiquitination to compulsive disorders

> [!NOTE]
> **All code used to generate the results in the manuscript is included in this repository**, apart
> from SNP-heritability estimation and the regional GWAS plots, which live in the external repositories
> linked in the figure table below. We are
> additionally connecting these analysis scripts into end-to-end Nextflow pipelines, so that the
> full analysis can be re-run from raw data with a single set of commands. This pipeline layer is
> still being tested stage by stage; in the meantime, the original analysis scripts (in each stage's
> `code/` folder and in `archive/`) remain the reference record of how each result was produced. Each
> stage's `README.md` under `code/` describes what it does and how it has been tested so far.

> Katarina Tengvall, Vista Sohrab, Matthew J. Christmas, Gustaf Brander, Eric Pederson, Raphaela Pensch, Brittney Kenney, Chengcheng Song, Ola Wallerman, Elisabeth Sundström, Chao Wang, Åsa Karlsson, Linnéa L. Monsén, Nordic OCD and Related Disorders Consortium (NORDiC), PGC TS/OCD working group, Kathryn A. Lord, Nora I. Strom, Christian Rück, David Mataix-Cols, Manuel Mattheisen, James J. Crowley, Jennifer R. S. Meadows, Elinor K. Karlsson, Kerstin Lindblad-Toh
>
> **Code:** Katarina Tengvall, Vista Sohrab, Elinor K. Karlsson, Mahesh Binzer-Panchal
>
> Manuscript submitted to *Nature*. This repository is tagged at submission and will be tagged again at publication, to track versions against any later updates or resubmission.

This repository contains the computational pipeline for genomic analyses of compulsive disorder in dogs (dogCD), and cross-species network studies between humans and dogs. dogCD is phenotyped using exploratory factor analysis of 14 compulsive-behaviour items from the survey on the Darwin's Ark platform. The pipeline processes sequencing data from 3,032 dogs sequenced at low coverage and 11 dogs sequenced at high coverage, and integrates canine genotype array data (Axiom Canine Genotyping Array Sets A and B, n = 411) lifted from CanFam3.1 to CanFam4. Low-coverage sequencing and array genotypes are then imputed using the Dog10K reference panel. The pipeline then performs GWAS (genome-wide association studies) for dogCD (n = 2,584). GWAS regions are functionally annotated against a dog brain cis-regulatory element (cCRE) atlas derived from ATAC-seq, CUT&RUN, and ChIP-seq chromatin profiling across eight brain regions (wet-lab assays not covered by this repository). Network propagation (NetColoc) then colocalizes dogCD with human OCD and related psychiatric traits to identify shared gene networks. Human and dog brain expression data are used for subsequent enrichment analysis of network genes and dogCD GWAS regions.

------------------------------------------------------------------------

## Repository Contents

This repository contains code and instructions necessary to reproduce the analysis, figures, and tables in the paper.

| Folder | Description |
|----------------------------|--------------------------------------------|
| `code/` | Analysis logic, one numbered stage per pipeline step (`00`–`07`) |
| `analyses/` | Per-stage Nextflow run configuration (`params.yml`, `run_nextflow.sh`, `nextflow.config`) |
| `data/` | Inputs, intermediate/derived files, and outputs, mirroring `code/`'s stage numbering |
| `env/` | Conda environment specs behind each stage's frozen Seqera Wave container |
| `archive/` | Superseded/original-author material and pre-conversion working directories, kept for provenance — not live code |

The same three-way split (`code/` = logic, `analyses/` = per-stage run config, `data/` = inputs/outputs) repeats identically across all numbered stages.

Some READMEs and code comments refer to the authors' internal working notes (`REPRODUCIBILITY_AUDIT.md`, `REVIEW_AND_SUGGESTIONS.md`, `AGENTS.md`) or to the manuscript draft folder (`data/Manuscript_draft/`). These are not part of this public copy; the corresponding author can provide details on request.

------------------------------------------------------------------------

## Pipeline Stages

A folder-by-folder map of what each stage does. Each stage's own `README.md` (in both its `code/` and `analyses/` directories) has full detail — inputs/outputs, process list, how to run it standalone.

### 00 — Reference data & full-cohort mapping (`analyses/00_fetch-raw-data`)

No `code/` counterpart — plain fetch scripts, not a Nextflow pipeline. `00a-fetch-reference-data` downloads reference files (genome, BQSR/depth sites, Dog10K GLIMPSE panel, Dog10K gene models, from which it also derives the GWAS clumping gene-range file) with checksum verification. `00c-fetch-scilifelab-deposit` downloads the pipeline inputs we deposited on the SciLifeLab Data Repository (Axiom genotypes, survey responses, breed/sex metadata), also checksum-verified. `00b-mapping-batched` fetches the full sequencing cohort from SRA (3,043 runs, \~3.5 TB total) and maps it, batched via `run_batched.sh` (which self-submits as a SLURM job) rather than a single `nf-core/fetchngs` run, specifically to stay within HPC disk quota. Cohort genotyping is a deliberately separate, manually-run step (`01b-cohort-genotyping`) once every batch is confirmed done — not auto-chained after `00b`.

``` bash
pixi run 00a-fetch-reference-data  # reference files only — small, run this first
pixi run 00c-fetch-scilifelab-deposit  # deposited Axiom genotypes, survey responses, breed/sex metadata
pixi run 00b-mapping-batched    # fetches + maps the full cohort, in disk-bounded batches
```

------------------------------------------------------------------------

## Correspondence Between Paper and Code

| Figure | Source | Pixi Task |
|-----------------------|-------------------------|------------------------|
| Fig 1a | Illustrations — no code | — |
| Fig 1b (heritability) | Plotting: `code/07_additional_plots_and_analyses/1_heritability_plot`. **Calculation code lives in a different repository**: `github.com/VistaSohrab/dog-gwas-heritability-nextflow` | `07a-heritability-plot` |
| Fig 1c (Manhattan/QQ) | `code/04_gwas/1_mlma-loco`, `2_polmm`, `5_plotting` | `04a-mlma-loco`, `04b-polmm`, `04e-gwas-plotting` |
| Fig 1d (dogCD vs. human OCD GWAS power) | `code/07_additional_plots_and_analyses/4_power_analysis` | `07e-power-analysis` |
| Fig 1e (seed-gene Venn) | No code needed — raw seed-gene lists, no overlap by construction | — |
| Fig 1f (NPS-threshold gene sets) | `code/05_network/2_CrossSpeciesBMI` (`CONSERVED_NETWORK` computes the gene counts + hypergeometric overlap p-value, `PLOT_CONSERVED_NETWORK_VENN` renders the panel) | `05b-cross-species-bmi` |
| Fig 1g (dogCD-vs-OCD z-score scatter) | `code/05_network/1_NetColoc` (`PLOT_NETCOLOC_SCATTER`) | `05a-netcoloc` |
| Fig 1h-i (colocalization permutation test, control comparisons) | `code/05_network/2_CrossSpeciesBMI` | `05b-cross-species-bmi` |
| Fig 2a-b (hierarchy-capture bar chart, gene-overlap Venn) | `code/05_network/4_Plotting` (`BUILD_NETWORK_GENE_COUNTS` + `PLOT_FIG2A` for 2a, `PLOT_VENN_HIERARCHY` for 2b) — underlying systems-map gene lists trace back to `2_CrossSpeciesBMI`'s merged hierarchy output | `05d-network-overlap` |
| Fig 2c-e (GO semantic clustering) | `code/05_network/5_GO_semantic_clustering` | `05c-go-semantic-clustering` |
| Fig 3a (hierarchy diagram) | `code/05_network/2_CrossSpeciesBMI` | `05b-cross-species-bmi` |
| Fig 3b (GO fold-enrichment) | `code/05_network/5_GO_semantic_clustering` | `05c-go-semantic-clustering` |
| Fig 3c-e (MP enrichment, community detail panels) | Cytoscape/Illustrator — no code, manual. HiDeF community network + cosine-similarity subgraph inputs come from `code/05_network/1_NetColoc`'s own source notebook — see `code/05_network/3_Cytoscape/README.md`. The Results sentence quantifying Fig 3c's term counts is checked (not reproduced exactly — see `2_CrossSpeciesBMI/README.md`'s Notes) by `MGD_NETWORK_ENRICHED_TERMS` | — (panel itself manual; the check runs under `05b-cross-species-bmi`) |
| Fig 4a-b (cross-species snRNA-seq) | `code/07_additional_plots_and_analyses/3_cross_species_expression`. Input is a single deposited Figshare object, `Siletti_human_subset_normalised.rds` (a subset of the Siletti et al. 2023 human brain snRNA-seq atlas, fetched by this substage's own `fetch-siletti-subset.sh`) | `07d-cross-species-expression` (needs `07c-fetch-siletti-subset`) |
| Fig 4c (emission-probability heatmap) | `code/06_cCRE/emission_state_UU.R` | `06a-emission-state-plot` |
| Fig 4d (UU vs. EpicDog) | `code/06_cCRE/compare_UU_epicdog` (`BUILD_ALL_TISSUES_FILTERED`, `COMPARE_UU_EPICDOG_ELEMENTS`) — element-count comparison; EpicDog side fetched live by `06d-fetch-epicdog-bed`, UU side (`all_tissues_filtered.bed`) built live from `06e-fetch-chromatin-states-bed`'s 8 per-region BEDs | `06f-compare-uu-epicdog` (needs `06d-fetch-epicdog-bed` + `06e-fetch-chromatin-states-bed`) |
| Fig 4e (cCRE odds ratios per UU brain region, dogCD vs. dog size); EpicDog brain vs. other-tissue odds ratios in the Results text and Supplementary Table 13; per-region overlaps in Supplementary Table 14 | `code/06_cCRE/cre_overlap` (`PLOT_REGION_FOREST` for Fig 4e; `PLOT_TISSUE_SPECIFICITY` for the EpicDog odds ratios). Uses GWAS regions defined for dogCD and size. | `06c-cre-overlap` (needs `06b-clump-regions-100kb`, `06d-fetch-epicdog-bed`, `06e-fetch-chromatin-states-bed`) |
| Fig 5 | Illustration — no code | — |

### Extended Data Figures

| Figure | Source | Pixi Task |
|-----------------------|-------------------------|------------------------|
| ED Fig 1 (per-region SuSiE lead-SNP effects/implicated genes) | **Outsourced to a different repository**: detailed analyses and plotting of dog GWAS regions were performed using scripts from [`github.com/JenniferMeadowsLab/Plotting_GWAS_results`](https://github.com/JenniferMeadowsLab/Plotting_GWAS_results) ([doi.org/10.5281/zenodo.23054418](https://doi.org/10.5281/zenodo.23054418)) — no code in this repo | — |
| ED Fig 2 (GO Biological Process enrichment across OCD/dogCD network communities) | `code/05_network/5_GO_semantic_clustering` | `05c-go-semantic-clustering` |
| ED Fig 3 (cross-species transcriptomic concordance heatmap) | `code/07_additional_plots_and_analyses/5_cross_species_transcriptomic_heatmap` — verified 2026-09-30, exact match against the published figure | `07g-cross-species-transcriptomic-heatmap` (input `human_dog_5regions.rds` is deposited at [doi.org/10.6084/m9.figshare.33787129](https://doi.org/10.6084/m9.figshare.33787129); `07f-fetch-human-dog-5regions` does not download it yet — see that stage's README) |
| ED Fig 4 (example regional plots for dogCD GWAS hits) | **Outsourced to a different repository**, same as ED Fig 1: [`github.com/JenniferMeadowsLab/Plotting_GWAS_results`](https://github.com/JenniferMeadowsLab/Plotting_GWAS_results) ([doi.org/10.5281/zenodo.23054418](https://doi.org/10.5281/zenodo.23054418)) — no code in this repo | — |

------------------------------------------------------------------------

## Reproducing the Analysis

### 1. Requirements

- [**pixi**](https://pixi.sh) — installs Nextflow (pinned version, see `pixi.toml`) and is the entry point for running anything: `pixi run <task>`, `pixi task list`.
- **Docker** — every substage's `stub_test` profile hardcodes `docker.enabled = true`, so Docker specifically is what the everyday stub-test/dev loop needs installed and running locally. Singularity/Apptainer profiles also exist for HPC systems without Docker (UPPMAX/Dardel profiles are predefined too).
- **Nothing else needs installing locally** — every tool an actual pipeline step needs (BWA, GATK, samtools, bcftools, PLINK/PLINK2, GLIMPSE, GCTA, PolyFun/SuSiE, bedtools, liftOver, R, Python, …) is containerized per-process via Seqera Wave; pixi/Nextflow pull the right container automatically.
- [**Seqera Wave CLI**](https://github.com/seqeralabs/wave-cli) — only needed if you're changing a container spec or adding a new one; not needed to just run the existing pipeline.

### 2. Configuration

The pipeline is built around the **Dog10K** reference panel and the UU_Cfam_GSD_1.0/CanFam4 assembly. No `.env` file or manual path configuration is needed — every stage's own `analyses/<stage>_*/params.yml` already points at the correct deposited/fetched input paths (relative to that `analyses/` subdirectory), and `pixi.toml`'s task `cwd` handles running each one from the right place. To point a stage at different data, edit that stage's `params.yml` directly — see its own `README.md` for what each parameter expects.

One exception: several `params.yml` files include a `project: uppmax2025-2-420` line — this is the manuscript authors' own UPPMAX/SLURM billing allocation, not a generic default. Anyone else running this on UPPMAX needs to swap in their own project ID; it's ignored on non-SLURM infrastructure (e.g. plain Docker/local execution).

### 3. Execution Order

Every step below is a `pixi run <task>` command — run `pixi task list` to see the full set with descriptions. Task names are prefixed `<stage><letter>-`, matching the numbered `analyses/<stage>_*` folder each one runs from (00-07, same numbering `code/`, `analyses/`, and `data/` all share) plus a letter ordering tasks within that stage — so `pixi task list`'s alphabetical output already reads top-to-bottom as a run order. Each subfolder's own README has more detail on that step's inputs/outputs. Steps on the same numbered line below have no dependency on each other and can run in parallel; later steps depend on everything above them unless noted.

Two tasks are deliberately **not** numbered and don't appear below because they aren't pipeline steps: `clean-nextflow` (removes stale Nextflow work directories/logs/caches — run it anytime) and `manual-fetch-sequence-data` (small-scale/manual sequence-data fetches only — never use it for the full cohort, see `analyses/00_fetch-raw-data/README.md`).

``` bash
# 00 + 01 — reference data, then fetch+map the full sequencing cohort, then cohort-genotype it
pixi run 00a-fetch-reference-data
pixi run 00c-fetch-scilifelab-deposit
pixi run 00b-mapping-batched      # full-cohort fetch+map only, in disk-bounded batches; 01a-fetch-map
                                   # is the same workflow exposed as a task for manual/small-scale runs
pixi run 01b-cohort-genotyping    # run yourself once every batch above is confirmed done -- not
                                   # auto-chained after 00b (see run_batched.sh's own header for why)

# 02 — survey factor analysis (independent of mapping), Axiom imputation, then merge + plot
pixi run 02a-fetch-axiom-jars && pixi run 02b-axiom-imputation
pixi run 02c-survey-data
pixi run 02d-merging-filtering    # needs 01b-cohort-genotyping + 02b-axiom-imputation
pixi run 02e-plotting             # needs 02d-merging-filtering

# 03 — genotype-concordance QC
pixi run 03a-qc                   # needs 01b-cohort-genotyping + 02b-axiom-imputation

# 04 — GWAS, fine-mapping, and downstream region work
pixi run 04a-mlma-loco            # needs 02d-merging-filtering
pixi run 04b-polmm                # needs 02d-merging-filtering
pixi run 04c-finemap-susie        # needs 02d-merging-filtering + 04a-mlma-loco + 04b-polmm
pixi run 04d-qc-signals-igv       # needs 02d-merging-filtering + 04c-finemap-susie (manual IGV review)
pixi run 04e-gwas-plotting        # needs 04a-mlma-loco + 04b-polmm; Fig 1c Manhattan/QQ plots

# 05 — network propagation, cross-species colocalization, and GO semantic clustering
pixi run 05a-netcoloc                # also produces Fig 1g; NPS-threshold sensitivity analysis (ST7)
                                      # is not produced by this pipeline — deposited/frozen only,
                                      # see 1_NetColoc/README.md
pixi run 05b-cross-species-bmi       # needs 05a-netcoloc; control-trait comparison (Fig 1i, ST10)
                                      # reloads a frozen deposited result, not re-executed live
                                      # (unseeded null permutation) — see 2_CrossSpeciesBMI/README.md
pixi run 05c-go-semantic-clustering  # reads deposited data/05_network/*.tsv directly, no pipeline
                                      # dependency; the GO BP enrichment behind those TSVs isn't
                                      # re-run live either — deterministic tests, but against
                                      # unpinned external reference databases — see 2_CrossSpeciesBMI/README.md
pixi run 05d-network-overlap         # Fig 2a-b; reads deposited data/05_network/*_systemsmap_genes.txt
                                      # directly, no pipeline dependency

# 06 — chromatin-state annotation (independent of the GWAS/network chain above)
pixi run 06a-emission-state-plot  # Fig 4c
pixi run 06b-clump-regions-100kb  # needs 04a-mlma-loco + 04b-polmm
pixi run 06c-cre-overlap          # Fig 4e (+ EpicDog tissue ORs); needs 06b-clump-regions-100kb, 06d-fetch-epicdog-bed
                                   # + 06e-fetch-chromatin-states-bed — all pulled in automatically
                                   # via depends-on
pixi run 06f-compare-uu-epicdog   # Fig 4d; needs 06d-fetch-epicdog-bed + 06e-fetch-chromatin-
                                   # states-bed — both pulled in automatically via depends-on

# 07 — additional plots and analyses (each independent, deposited/archived inputs only)
pixi run 07a-heritability-plot        # Fig 1b
pixi run 07b-gwas-catalog-overlap     # "12 of 27" GWAS-catalog claim, mechanical join
pixi run 07d-cross-species-expression # Fig 4a-b; needs 07c-fetch-siletti-subset (2.58 GB,
                                       # run once) — pulled in automatically via depends-on
pixi run 07e-power-analysis           # Fig 1d; dogCD vs. human OCD GWAS power comparison
```

Estimated total runtime: <!-- FILL IN: e.g. "48 hours on 32-core node with 128GB RAM" -->

------------------------------------------------------------------------

## Data

Both — most raw sequencing data is external (SRA), but a substantial amount of intermediate/ deposited data (seed gene lists, frozen network-analysis results, chromatin-state BEDs, etc.) is committed directly under `data/`, split by stage to mirror `code/`'s numbering. See each stage's own `data/<stage>/README.md` for what's committed vs. fetched vs. pipeline-generated.

### Genomic Data Sources

| Data Type | Source | Access |
|----------------------------|----------------------|----------------------|
| Low-coverage dog sequencing (3,032 dogs) | NCBI SRA | [PRJNA675863](https://www.ncbi.nlm.nih.gov/seqread/PRJNA675863) |
| High-coverage dog sequencing (11 dogs) | NCBI SRA | [PRJNA683923](https://www.ncbi.nlm.nih.gov/seqread/PRJNA683923) |
| Axiom canine genotype array (411 dogs, CanFam3.1), lifted/merged/imputed genotypes (CanFam4), phenotype data, GWAS sumstats, Darwin's Ark breed/sex metadata for the genotyped dogs (`DarwinsArk_20220715_dogs_genotyped_breed_sex.csv`), recoded dogCD survey answers for 25,302 dogs (`response_df_dogCDitems.txt`) | SciLifeLab Data Repository (Figshare) | [doi.org/10.17044/scilifelab.33339309](https://doi.org/10.17044/scilifelab.33339309) |
| ATAC-seq/CUT&RUN/ChIP-seq raw reads (cCRE atlas) | ENA | [PRJEB124064](https://www.ebi.ac.uk/ena/browser/view/PRJEB124064) |
| Processed cCRE peak files + ChromHMM state assignments (UU and CanFam4-remapped EpicDog) | Figshare | [doi.org/10.6084/m9.figshare.31971747](https://doi.org/10.6084/m9.figshare.31971747) |
| EpicDog cCRE raw data (pre-CanFam4 remapping) | NCBI GEO | [GSE203104](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203104) |
| PacBio Iso-seq hypothalamic expression | NCBI SRA | [PRJNA657719](https://www.ncbi.nlm.nih.gov/seqread/PRJNA657719) |
| Dog and human brain single nuclei RNA-seq datasets, as Seurat objects (incl. the Siletti et al. 2023 human subset) | Figshare | [doi.org/10.6084/m9.figshare.33787129](https://doi.org/10.6084/m9.figshare.33787129) |
| Dog brain snRNA-seq Expression dataset (DBEx) — explorable via web app | SciLifeLab | [dog-brain-single-cell.serve.scilifelab.se](https://dog-brain-single-cell.serve.scilifelab.se) |
| Survey data (Darwin's Ark) | Darwin's Ark platform (datafreeze 2022-07-15) | [darwinsark.org](https://darwinsark.org) — phenotype data also bundled in the SciLifeLab Axiom deposit above |

### Reference Panels

| Resource | Version | Usage |
|----------------------------|-------------------------|--------------------|
| Dog10K panel | 1,929-canine phased imputation panel ([EVA PRJEB62420](https://www.ebi.ac.uk/ena/browser/view/PRJEB62420), [kiddlabshare.med.umich.edu/dog10K](https://kiddlabshare.med.umich.edu/dog10K/)) | GLIMPSE reference for imputation |
| CanFam4 genome | UU_Cfam_GSD_1.0 autosomes/X (GCF_011100685.1) + ROS_Cfam_1.0 Y chromosome (GCF_014441545.1) | Reference assembly for mapping, LiftOver target |
| TOGA | Human-dog 1:1 orthologues, `orthologsClassification.tsv.gz` ([genome.senckenberg.de/download/TOGA](http://genome.senckenberg.de/download/TOGA/)) | Human orthologue mapping (`1_NetColoc`) |
| Zoonomia phastCons elements | 241-way mammalian alignment, constrained elements ([doi.org/10.5281/zenodo.20136916](https://doi.org/10.5281/zenodo.20136916)) | Evolutionary-conservation overlap for cCREs and GWAS regions |
| Zoonomia phyloP constraint scores | 241-way mammalian alignment ([doi.org/10.5281/zenodo.20137136](https://doi.org/10.5281/zenodo.20137136)) | GWAS functional annotation — constrained sites (phyloP \> 2.27) |
| NCBI gene models | *Canis lupus familiaris* annotation release 106 ([ftp.ncbi.nlm.nih.gov](https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9615/106/)) | Gene models for GWAS/cCRE functional annotation |

------------------------------------------------------------------------

## Repository DOI / Archival Version

A static snapshot of this repository at the time of publication is archived at: \> <!-- FILL IN: Zenodo / Figshare / OSF DOI here, or note as missing --> \> \[DOI not yet assigned\]

------------------------------------------------------------------------

## License

Code: Apache License 2.0 — see [`LICENSE`](LICENSE).\
Data: data files in this repository and our deposited datasets (SciLifeLab Data Repository / Figshare, DOIs in the Data section) are licensed under CC BY 4.0. Third-party data (e.g. NCBI SRA, Dog10K, EpicDog, the Gene Ontology annotation, the GWAS Catalog, published GWAS summary statistics and the gene lists derived from them) remain under their original providers' terms; see the sources listed above and in each `data/<stage>/README.md`.

------------------------------------------------------------------------

## Citation

If you use this code, please cite the article:

> Tengvall K, Sohrab V, Christmas MJ, *et al.* Dog-human gene network colocalization links synaptic signalling and ubiquitination to compulsive disorders. *Submitted* (2026).

The full author list and citation metadata are in [`CITATION.cff`](CITATION.cff); on GitHub, use **Cite this repository** in the sidebar.

------------------------------------------------------------------------

## AI Assistance

The conversion of the original analysis scripts into Nextflow/pixi pipelines (containers, workflow wiring, stub tests and documentation) was carried out with the help of [Claude Code](https://claude.com/claude-code) (Anthropic; Claude Opus models), under the direction of the authors. All changes were reviewed by the authors, and results were checked against the original analyses where data allowed. The authors take full responsibility for the code.

------------------------------------------------------------------------

> \[!TIP\] The dog genome mapping pipeline is optimized for the **Dog10K** reference panel. Each stage's own `analyses/<stage>_*/params.yml` is where its input paths are set — no environment variables needed.