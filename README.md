# End-to-End Single-Cell RNA-Seq Analysis Pipeline (Seurat v5)
**Profiling 3k Peripheral Blood Mononuclear Cells (PBMCs) with Unsupervised Clustering and Cell Type Annotation**

![R](https://img.shields.io/badge/R-4.3+-blue.svg)
![Seurat](https://img.shields.io/badge/Seurat-v5-green.svg)
![License](https://img.shields.io/badge/License-MIT-yellow.svg)

---

## Overview
This repository contains a reproducible, end-to-end single-cell RNA-sequencing (scRNA-seq) workflow implemented in **R** using **Seurat v5**. The pipeline processes the standard 10x Genomics **PBMC 3k** dataset (Peripheral Blood Mononuclear Cells from a healthy donor), performing quality control, normalization, feature selection, dimensionality reduction, unsupervised graph-based clustering, differential expression analysis, and biological cell type annotation.

---

## Biological Insights & Results

The analysis identified **9 distinct immune cell populations** based on well-established canonical markers:

| Cluster | Identified Cell Type | Canonical / Top Markers | Cell Count | Percentage (%) |
| :---: | :--- | :--- | :---: | :---: |
| **0** | Naive CD4+ T cells | `CCR7`, `LEF1`, `MAL`, `PIK3IP1` | 684 | 25.9% |
| **1** | CD14+ Monocytes | `CD14`, `S100A8`, `S100A9`, `FOLR3` | 481 | 18.2% |
| **2** | Memory CD4+ T cells | `AQP3`, `CD40LG`, `CD2`, `TRAT1` | 476 | 18.0% |
| **3** | B cells | `CD79A`, `VPREB3`, `TCL1A`, `LINC00926` | 344 | 13.0% |
| **4** | CD8+ T cells | `CD8A`, `GZMK`, `GZMH`, `CCL5` | 291 | 11.0% |
| **5** | FCGR3A+ (CD16+) Monocytes | `CKB`, `CDKN1C`, `MS4A4A`, `HES4` | 162 | 6.1% |
| **6** | Natural Killer (NK) cells | `GNLY`, `GZMB`, `SPON2`, `AKR1C3` | 155 | 5.9% |
| **7** | Dendritic Cells (DCs) | `FCER1A`, `CLEC10A`, `SERPINF1` | 32 | 1.2% |
| **8** | Platelets | `PPBP`, `ITGA2B (CD41)`, `GP9`, `PF4` | 13 | 0.5% |

<p align="center">
  <img src="figures/08_umap_annotated_celltypes.png" width="65%" alt="UMAP Cell Annotation">
  <img src="figures/09_cell_type_proportions.png" width="65%" alt="Cell Type Proportions">
</p>

---

## Pipeline Workflow

1. **Quality Control & Filtering (`01_qc_filter.R`)**:
   - Filtered low-quality cells and potential doublets/empty droplets:
     - `nFeature_RNA`: 200 to 2,500 genes
     - Mitochondrial read percentage: < 5%
   - Retained **2,638 high-quality cells** out of 2,700 raw cells.

2. **Normalization & Feature Selection (`02_normalize_pca.R`)**:
   - Applied global scaling normalization (`LogNormalize`, scale factor = 10,000).
   - Identified **top 2,000 highly variable genes (HVGs)** using variance-stabilizing transformation (`vst`).
   - Scaled data and performed linear dimensionality reduction via **PCA**.

3. **Graph-based Clustering & UMAP (`03_clustering_umap.R`)**:
   - Constructed Shared Nearest Neighbor (SNN) graph using top 10 principal components.
   - Clustered cells using the Louvain algorithm (`resolution = 0.5`).
   - Projected cells onto 2D space using **UMAP**.

4. **Biomarker Identification (`04_find_markers.R`)**:
   - Wilcoxon Rank Sum test with Bonferroni correction (`FindAllMarkers`).
   - Extracted positive cluster markers (`log2FC >= 0.25`, `min.pct = 0.25`).
   - Generated marker expression Heatmaps and DotPlots.

5. **Cell Type Annotation (`05_cell_annotation.R`)**:
   - Curated and assigned immunological identities based on known markers.
   - Generated final publication-ready figures.

---

## Repository Structure

```text
├── data/
│   └── filtered_gene_bc_matrices/hg19/
├── figures/
├── results/
│   ├── qc/
│   ├── clustering/
│   ├── markers/
│   └── final_annotated_pbmc.rds
├── scripts/
│   ├── 01_qc_filter.R
│   ├── 02_normalize_pca.R
│   ├── 03_clustering_umap.R
│   ├── 04_find_markers.R
│   └── 05_cell_annotation.R
├── environment.yml
├── .gitignore
└── README.md
```

---

## Installation & Usage

### 1. Clone the repository
```bash
git clone https://github.com/shayesteh68/scrnaseq-seurat-pbmc3k.git
cd scrnaseq-seurat-pbmc3k
```

### 2. Set up environment
```bash
conda env create -f environment.yml
conda activate r_seurat_env
```

### 3. Run the pipeline
```bash
Rscript scripts/01_qc_filter.R
Rscript scripts/02_normalize_pca.R
Rscript scripts/03_clustering_umap.R
Rscript scripts/04_find_markers.R
Rscript scripts/05_cell_annotation.R
```

---

## Author
**Narges Shayesteh**
- LinkedIn: [linkedin.com/in/narges-shayesteh](https://linkedin.com/in/narges-shayesteh)
- GitHub: [@shayesteh68](https://github.com/shayesteh68)
