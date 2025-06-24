# TNBC scRNA-seq: Metastatic Subpopulations Project

This project analyzes single-cell RNA sequencing data from triple-negative breast cancer (TNBC) to identify metastatic subpopulations and tumor heterogeneity.

## 📁 Project Structure

- `data/` – Contains placeholders and instructions to download raw data
- `scripts/` – Analysis scripts
- `output/` –  Generated plots and results
- `TNBC_scRNAseq_Analysis.Rmd` – Main RMarkdown analysis 
- `PowerPoint/` – Final presentation

## 🧬 Dataset

- [GEO: GSE118389](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118389)
- Original paper: [Nature Communications, 2018](https://www.nature.com/articles/s41467-018-06052-0)
- GitHub repo: [Michorlab/tnbc_scrnaseq](https://github.com/Michorlab/tnbc_scrnaseq)

## 🔬 Analysis Goals

1. Perform QC and filtering
2. Identify cell types and tumor cells
3. Discover tumor subpopulations
4. Detect metastatic signatures
5. Correlate metastasis scores with other genes/processes

## 📦 Requirements

- R (>= 4.2)
- Seurat
- dplyr, ggplot2, patchwork
- clusterProfiler, Matrix, readr

