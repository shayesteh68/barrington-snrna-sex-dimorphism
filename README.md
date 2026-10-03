# Single-Nucleus RNA-Seq Analysis of Barrington's Nucleus (PMC)

## Analysis of GSE226809 — Neuronal Subclustering with Seurat v5

This repository contains a reproducible **single-nucleus RNA-seq (snRNA-seq) analysis** of the **Barrington's nucleus** (also referred to as the pontine micturition center, PMC), based on the publicly available GEO dataset **GSE226809**.

**Author:** Narges Shayesteh — Bioinformatics Researcher
**LinkedIn:** https://www.linkedin.com/in/narges-shayesteh

---

## 1. Biological Background

The Barrington's nucleus is a compact group of neurons in the dorsomedial pons that acts as the **pontine micturition center (PMC)**. It coordinates the bladder–urethra reflex circuit and communicates with the spinal cord, periaqueductal gray (PAG), and hypothalamic and cortical regions.

Single-nucleus profiling of this small, deep brain structure resolves its molecularly defined neuronal populations (glutamatergic, GABAergic, cholinergic, and neuromodulatory neurons) as well as resident glia, which cannot be captured reliably in dissociated bulk tissue.

## 2. Dataset

| Property | Value |
|---|---|
| GEO accession | GSE226809 |
| Assay | Single-nucleus RNA-seq (snRNA-seq) |
| Tissue / region | Brainstem — Barrington's nucleus (pontine micturition center) |
| Technology | 10x Genomics Chromium, nuclear transcriptomes |
| Framework | Seurat v5 (R) |
| Focus | Unsupervised clustering and neuronal/glial subclustering |

Raw count matrices and cell-level metadata were downloaded directly from GEO and parsed with `Read10X_h5` / `Read10X` depending on the deposited format.

## 3. Analysis Pipeline

The workflow was implemented in **R 4.x with Seurat v5**:

1. **Loading** — construction of a Seurat object from the deposited count matrices (`CreateSeuratObject`).
2. **Quality control** — per-nucleus filtering on:
   - `nFeature_RNA` (detected genes)
   - `nCount_RNA` (UMI depth)
   - percent mitochondrial reads (`percent.mt`)
3. **Normalization** — log-normalization (`NormalizeData`, `LogNormalize`, scale factor 10,000) and detection of highly variable genes (`FindVariableFeatures`).
4. **Scaling and dimensionality reduction** — `ScaleData`, PCA (`RunPCA`), and an elbow/intrinsic-dimensionality survey of the top PCs.
5. **Neighbor graph and clustering** — `FindNeighbors` on significant PCs followed by multi-resolution Leiden/Louvain clustering (`FindClusters`).
6. **Visualization** — UMAP (`RunUMAP`) plus QC violin plots, feature plots, and cluster marker dot plots.
7. **Cell-type annotation** — canonical marker-based assignment of neuronal and non-neuronal populations:
   - **Glutamatergic neurons:** *Slc17a6 (Vglut2), Slc17a7*
   - **GABAergic neurons:** *Gad1, Gad2, Slc32a1*
   - **Cholinergic neurons:** *Chat, Slc18a3*
   - **Neuropeptidergic / PMC-enriched:** *Adcyap1, Crh, Nts, Pnoc*
   - **Astrocytes:** *Aqp4, Gfap*; **Oligodendrocytes:** *Mbp, Mog, Plp1*
   - **Microglia:** *C1qa, Cx3cr1*; **Endothelial:** *Pecam1*
8. **Differential expression** — cluster and population markers via `FindAllMarkers` (Wilcoxon, `min.pct = 0.25`), filtered to biologically meaningful markers.

## 4. Repository Structure

```text
GSE226809/
├── data/
│   ├── raw/                      # GEO-deposited count matrices (h5 / mtx)
│   └── metadata/                 # sample / cell metadata
├── scripts/
│   ├── 01_qc_filtering.R
│   ├── 02_integration_clustering.R
│   ├── 03_annotation_subclustering.R
│   └── 04_markers_plots.R
├── results/
│   ├── seurat_obj/               # saved .rds objects
│   ├── markers/                  # cluster and cell-type marker tables (csv)
│   └── plots/
│       ├── qc_violins.png
│       ├── elbow_plot.png
│       ├── umap_clusters.png
│       ├── umap_celltypes.png
│       └── marker_dotplot.png
├── README.md
└── environment.yml
```

## 5. Reproducibility

The full analysis runs from the scripts in `scripts/` in numerical order against `data/raw/`. Key reproducibility anchors:

- Seurat **v5** object model (layered counts) with `JoinLayers` before downstream plotting.
- Fixed random seed for PCA/UMAP/clustering reproducibility.
- QC thresholds documented per sample in `scripts/01_qc_filtering.R`.
- Cell-type calls are marker-based and should be revisited if a curated reference atlas for the Barrington's nucleus becomes available.

### Environment

```r
library(Seurat)      # v5
library(ggplot2)
library(dplyr)
library(patchwork)
```

## 6. Caveats

- Nuclei-only protocols under-sample cytoplasmic transcripts; absolute expression levels differ from whole-cell scRNA-seq.
- Deep brainstem dissection may carry admixed neighboring nuclei; subclustering and marker checks guard against this.
- Marker-based annotation is an inference, not ground truth; functional claims require electrophysiology or retrograde-tracing validation.

## 7. Author

**Narges Shayesteh**
Bioinformatics Researcher — single-cell and single-nucleus genomics
LinkedIn: https://www.linkedin.com/in/narges-shayesteh
