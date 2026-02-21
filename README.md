# skin-aging

Code and reproducible analysis workflow for the manuscript:

**“A single-cell atlas of age-associated changes in sun-protected human penile skin reveals stromal depletion and inflammatory fibroblast remodeling”**  
(submitted to *GeroScience*)

This repository contains scripts for reprocessing, integration, and downstream analyses of human penile skin scRNA-seq data, including donor-aware pseudo-bulk modeling, continuous-age trend modeling, temporal module analysis, and senescence signature anchoring (SenSkin™).

---

## Data availability

### Controlled-access raw data
Raw sequencing data (FASTQ) were obtained from the **Genome Sequence Archive (GSA)** under accession **HRA007611** via controlled access. This repository does **not** contain controlled-access data.

### Processed data
Processed data products generated in this work (expression matrices, integrated object, and donor metadata) are available on Figshare:  
- **DOI:** 10.6084/m9.figshare.30498725

---

## Repository structure
skin-aging/
├─ README.md
├─ LICENSE
├─ config/
├─ scripts/
│ ├─ 00_cellranger.sh
│ ├─ 01_Data_import_and_Processing.py
│ ├─ 02_decontx.R
│ ├─ 03_Basic_Process_for_the_remove_ambientRNA_count.py
│ ├─ 04_integration_scvi.py
│ ├─ 05_composition_analysis.py
│ ├─ 06_pseudobulk_DE.R
│ ├─ 07_GAM_age_model.R
│ ├─ 08_modules_clusterGvis.R
│ └─ 09_external_benchmarking.R
└─ .gitignore

> Note: large intermediate results (e.g., `*.h5ad`, `*.rds`) should not be committed to GitHub. See `.gitignore`.

---

## Computational environment

- Cell Ranger: v6.1.2  
- Python: 3.9 (Scanpy v1.9.3; scvi-tools v0.19.0; Scrublet v0.2.3)  
- R: 4.2.1 (DaMiRseq; limma; mgcv; clusterProfiler; clusterGvis)
