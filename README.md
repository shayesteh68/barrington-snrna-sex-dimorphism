# Barrington's nucleus snRNA-seq analysis

This repository contains a focused, independent reanalysis of the 30 Barrington's nucleus (Bar/PMC) count matrices in GEO series GSE226809. It is not a reproduction of either the full dorsal pons atlas or the later Bar-Penk paper.

This is an independent portfolio case study using public data. It was not commissioned by, or conducted in collaboration with, the source-study authors or their institutions.

## Sources and scope

The source atlas contains 72 DroNc-seq DGE matrices: 42 from parabrachial nucleus-centered dissections and 30 from Bar-centered dissections. This workflow analyzes the 30 Bar matrices only. The accession, GEO title, region, `Bar_*` group label, and sex for each library are recorded in [`metadata/GSE226809_PMC_GEO_sample_manifest.tsv`](metadata/GSE226809_PMC_GEO_sample_manifest.tsv); the loader checks each file against this manifest.

Data source and related study:

- Source atlas: Nardone et al., “A spatially-resolved transcriptional atlas of the murine dorsal pons at single-cell resolution,” *Nature Communications* (2024), [doi:10.1038/s41467-024-45907-7](https://doi.org/10.1038/s41467-024-45907-7).
- Related Bar-Penk study: Klymko et al., “Descending Enkephalinergic neurons in Barrington’s nucleus gate sex-biased control of micturition,” *Nature Communications* (2026), [doi:10.1038/s41467-026-77949-4](https://doi.org/10.1038/s41467-026-77949-4). The authors' processed objects and Seurat workflow are at [Zenodo record 15225406](https://doi.org/10.5281/zenodo.15225406).

The later Bar-Penk study uses a staged analysis to identify an excitatory Penk population. Its filters, doublet handling, SCTransform workflow, cell exclusions, and reclustering differ from this repository. The analysis here therefore does not claim to reproduce that population or its published figures.

## Methods and interpretation

The workflow applies the source atlas's mouse QC limits of 400-10,000 detected genes, at most 10% mitochondrial counts, and at most 5% hemoglobin counts. It uses LogNormalize, 2,000 variable genes, 20 PCs, resolution 0.4, fixed random seeds, and no batch correction. In the 2026-10-29 audit run, 53,035 cells passed these filters and formed 23 clusters. These settings intentionally define a lower-memory independent analysis; results are expected to differ from the later paper.

The marker workflow calculates canonical-marker detection and mean expression from all QC-passing cells. The genome-wide Wilcoxon marker scan uses a reproducible sample of at most 250 cells per cluster to bound memory. Its cell-level p-values are descriptive and are not biological-replicate tests or automatic cell-type assignments.

For the secondary sex contrast, the pipeline selects the cluster with the highest mean log-normalized Penk expression in the current run, then sums raw counts by the nine GEO `Bar_*` group labels (five female and four male). The selected cluster is a candidate, not a validated match to the published excitatory Bar-Penk population. In the audit run it was cluster 7: 38.4% Penk-positive, 8.7% Slc17a6-positive, 47.2% Gad2-positive, and 44.5% Mbp-positive. These mixed markers make its cell identity provisional.

The GEO record [GSE226809](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226809) describes the full experiment as 9 male and 8 female batches, with 3–5 mice per batch. This Bar-focused subset contains five female and four male batch pools. The DESeq2 contrast therefore uses the pool/batch as its unit; it does not recover individual-mouse measurements or estimate within-pool mouse-to-mouse variation. Nuclei and GSM libraries within a pool are aggregated rather than treated as independent replicates. The unpaired pool-level contrast is exploratory. In the audit run, six genes passed the formal absolute log2 fold-change >0.5 test at BH-adjusted p <0.05; all six were sex-chromosome genes. Directional GO enrichment is skipped when fewer than three significant mapped genes are available, to avoid interpreting pathways driven by only one or two genes.

## Workflow

Run commands from the repository root. Create the Conda environment with `conda env create -f environment.yml`. For the pinned Linux package solution, use `environment-linux-64.lock.txt`.

1. `Rscript scripts/00_download_geo.R` downloads and validates the GEO supplementary archive, then extracts the matrices to `data/raw/`.
2. `Rscript scripts/01_load_qc_integrate.R` validates integer counts, saves sparse per-library ingest checkpoints, checks GEO metadata, applies QC, clusters the cells, and saves the processed object and UMAP.
3. `Rscript scripts/02_annotate_celltypes.R` writes canonical-marker evidence from all cells, cluster sizes, a memory-bounded descriptive marker scan, and marker plots.
4. `Rscript scripts/03_penk_sex_differential.R` records the current highest-Penk cluster, aggregates raw counts by GEO group, and runs the exploratory batch-group pseudobulk DESeq2 contrast.
5. `Rscript scripts/04_pathway_enrichment.R` uses the tested, Entrez-mapped genes as the GO background and writes a summary even when no direction has enough genes for enrichment.

## Outputs and data handling

The workflow writes sample and QC manifests, cluster sizes, marker evidence, cluster-selection evidence, pseudobulk sample QC, DESeq2 results, enrichment summaries, figures, and R session information under `results/` and `run_logs/`. Raw GEO files, extracted count matrices, sparse ingest checkpoints, and Seurat objects stay under the ignored `data/` directory and are not committed. The MIT license applies to repository code; source data remain subject to GEO and source-study terms.
