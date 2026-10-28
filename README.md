# Barrington's nucleus single-nucleus RNA-seq

This repository analyzes the Barrington's nucleus (PMC) subset of GEO series GSE226809, a single-cell atlas of the mouse dorsal pons. The workflow loads GEO library count matrices, performs QC and clustering, then estimates a sex contrast for the selected Penk-high cluster.

## Data and environment

Run Rscript scripts/00_download_geo.R to download and unpack GEO supplementary files into data/raw/. The analysis expects files named GSM..._Mm_PMC_..._dge.tsv.gz.

Create and activate the environment with conda env create -f environment.yml and conda activate barrington-snrna.

## Pipeline

Run these commands from the repository root:

1. Rscript scripts/00_download_geo.R
2. Rscript scripts/01_load_qc_integrate.R
3. Rscript scripts/02_annotate_celltypes.R
4. Rscript scripts/03_penk_sex_differential.R
5. Rscript scripts/04_pathway_enrichment.R

The previous cell-level Wilcoxon p-values treated nuclei as independent biological replicates. The updated differential analysis uses one pseudobulk count column per GEO batch and requires at least two batches per sex. Old DEG and enrichment outputs are retained for provenance but are not validated by the corrected model; regenerate them before citing.

Cluster 8 is the default Penk-positive population carried over from the existing annotation. Confirm the cluster-to-cell-type mapping after rerunning clustering because cluster numbers can change.

The GEO records identify PMC libraries with a `Bar_*` batch and an F/M sex suffix. The loader records every GSM library and batch in `results/tables/GSE226809_PMC_library_manifest.tsv`, then aggregates libraries from the same GEO batch for pseudobulk. The GEO series describes each batch as a pool of multiple mice, so the batch is the replicate unit; cells and GSM libraries are not counted as independent replicates. Cluster 8 is used only if it remains the Penk-high cluster in the current object. The MIT license applies to repository code; GEO data remain subject to the source terms.
