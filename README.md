# Sexual Dimorphism in Barrington’s Nucleus (`Penk+`) Neurons: snRNA-seq Pipeline

[![R](https://img.shields.io/badge/R-4.x-blue.svg)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-v5.0-green.svg)](https://satijalab.org/seurat/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Author](https://img.shields.io/badge/Author-Narges%20Shayesteh-orange.svg)](https://github.com/shayesteh68)

An end-to-end, memory-optimized single-nucleus RNA sequencing (snRNA-seq) pipeline analyzing pontine micturition center (PMC / Barrington's nucleus) cell types, with a primary focus on sexually dimorphic gene expression within **`Penk+` Barrington projection neurons**.

Based on public snRNA-seq cohort data from **Dr. Anne M.J. Verstegen's laboratory** ([GEO: GSE226809](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226809)).

---

## 🔬 Scientific Background & Rationale

Barrington’s nucleus (pontine micturition center; PMC) provides primary supraspinal control over the lower urinary tract, orchestrating detrusor contraction and urinary sphincter coordination. Proenkephalin-expressing (`Penk+`) and corticotropin-releasing hormone-expressing (`Crh+`) neurons within Barrington’s nucleus show significant clinical and translational relevance regarding sex-biased bladder dysfunction, interstitial cystitis, and stress-induced urinary incontinence.

This project delivers:
1. **Large-scale QC & Integration** across 30 snRNA-seq samples (66,994 high-quality nuclei).
2. **High-resolution cluster annotation** resolving Barrington excitatory, inhibitory, and glial populations.
3. **Targeted Sex-Differential Expression Analysis** within confirmed `Penk+` Barrington neurons ($N = 1,632$ nuclei).

---

## 📊 Key Results & Biological Findings

### 1. Robust Cell-Type Deconvolution
Across 23 resolved clusters, canonical cell-type identities were clearly established:
- **Cluster 8 (`Penk+` Barrington Neurons):** High-level marker expression (`Penk` = 8.36, `Slc17a6` = 1.22, `Gad1` = 1.49).
- **Clusters 13 & 21 (Glutamatergic Neurons):** Enriched in `Slc17a6` (2.58 - 2.77) with baseline `Gad1` expression.
- **Clusters 7, 12 & 22 (GABAergic Interneurons):** Marked by robust `Gad1` expression (2.27 - 2.56).
- **Clusters 2, 19 & 20 (Astrocytes):** Confirmed by high `Aqp4` levels (6.44 - 11.03).
- **Clusters 3, 4 & 18 (Oligodendrocytes):** Confirmed by pan-myelin `Mbp` expression (59.3 - 161.7).

### 2. Validated Sex-Dimorphic Profiles in `Penk+` Neurons (Wilcoxon Test)
Differential expression between female ($n = 949$) and male ($n = 683$) nuclei in Cluster 8 revealed high internal validation:
* **Female-Enriched:** 
  * `Xist` ($\log_2\text{FC} = +8.49$, $p_{adj} = 5.58 \times 10^{-198}$, $81.6\%$ vs $0.1\%$)
  * `Tsix` ($\log_2\text{FC} = +5.69$, $p_{adj} = 7.71 \times 10^{-37}$, $23.8\%$ vs $0.3\%$)
* **Male-Enriched:** 
  * `Eif2s3y` ($\log_2\text{FC} = -11.26$, $p_{adj} = 1.50 \times 10^{-56}$)
  * `Kdm5d` ($\log_2\text{FC} = -10.22$, $p_{adj} = 2.83 \times 10^{-28}$)
  * `Uty` ($\log_2\text{FC} = -10.03$, $p_{adj} = 7.36 \times 10^{-23}$)
  * `Ddx3y` ($\log_2\text{FC} = -9.73$, $p_{adj} = 7.46 \times 10^{-21}$)
* **Physiological / Autosomal Candidates:**
  * `mt-Nd1` ($\log_2\text{FC} = -0.58$, $p_{adj} = 0.027$): Male-biased baseline mitochondrial respiratory complex I subunit expression.
  * `Nr3c2` ($p = 1.82 \times 10^{-4}$): Male-biased mineralocorticoid receptor expression implicated in osmoregulation and volume-dependent micturition reflexes.

---

## 📁 Repository Structure
```text
barrington-snrna-sex-dimorphism/
├── data/
│   ├── raw/               # GSE226809 DGE matrices (ignored in git)
│   └── processed/         # Joined Seurat RDS objects (ignored in git)
├── scripts/
│   ├── 01_load_qc_integrate.R       # QC filtering, Seurat v5 Layer joining, PCA/UMAP
│   ├── 02_annotate_celltypes.R      # Marker dot plots, feature plots & cluster profiling
│   └── 03_penk_sex_differential.R   # Cluster 8 isolation, Wilcoxon test, Volcano plot
├── results/
│   ├── figures/
│   │   ├── barrington_umap_clusters_sex.png
│   │   ├── barrington_canonical_markers_dotplot.png
│   │   └── penk_neurons_sex_volcano.png
│   └── tables/
│       ├── cluster_canonical_markers_summary.tsv
│       └── penk_neurons_sex_degs.tsv
├── .gitignore
└── README.md
⚙️ Installation & Reproducibility
Conda Environment
bash
conda create -n r_seurat_env -c conda-forge -c bioconda r-base=4.3 r-seurat=5.0 r-ggplot2 r-dplyr r-ggrepel
conda activate r_seurat_env
Execution
bash
# 1. Quality Control, Normalization, Integration & UMAP
Rscript scripts/01_load_qc_integrate.R

# 2. Cell-type Profiling & DotPlot Generation
Rscript scripts/02_annotate_celltypes.R

# 3. Cluster 8 Penk+ Sex-Differential Analysis & Volcano Plot
Rscript scripts/03_penk_sex_differential.R
👩‍💻 Author
Narges Shayesteh

Bioinformatics & Computational Biology Specialist

GitHub: @shayesteh68
LinkedIn: Narges Shayesteh
