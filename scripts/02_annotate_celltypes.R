suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(Matrix)
})
set.seed(42)

cat("=== Step 1: Loading Processed Seurat Object ===\n")
pmc <- readRDS("data/processed/pmc_barrington_seurat.rds")
DefaultAssay(pmc) <- "RNA"

cat("Number of cells loaded:", ncol(pmc), "\n")
cat("Identified Clusters:", length(unique(Idents(pmc))), "\n")

cat("=== Step 2: Defining Canonical Cell Markers ===\n")
# Canonical markers relevant to Barrington's PMC & LPBN neuroanatomy
markers_to_check <- list(
  "Barrington_PMC"  = c("Crh", "Penk", "Tac1", "Gal"),
  "Glutamatergic"   = c("Slc17a6"),
  "GABAergic"       = c("Gad1", "Gad2", "Slc32a1"),
  "Astrocytes"      = c("Aqp4", "Gja1"),
  "Oligodendrocytes"= c("Mbp", "Mog"),
  "OPC"             = c("Pdgfra", "Cspg4"),
  "Microglia"       = c("Cx3cr1", "Tmem119"),
  "Endothelial"     = c("Flt1", "Cldn5")
)

# Flatten list and verify presence in dataset
all_markers <- unlist(markers_to_check, use.names = FALSE)
existing_markers <- all_markers[all_markers %in% rownames(pmc)]
if (!length(existing_markers)) stop("None of the configured marker genes are present in the RNA assay.")
cat("Markers detected in RNA assay:", paste(existing_markers, collapse = ", "), "\n")
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("run_logs", recursive = TRUE, showWarnings = FALSE)

cluster_id <- as.character(Idents(pmc))
cluster_sizes <- as.data.frame(table(cluster = cluster_id), stringsAsFactors = FALSE)
names(cluster_sizes)[2L] <- "cells"
write.table(cluster_sizes, "results/tables/cluster_sizes.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)

cat("=== Step 3: Generating Marker DotPlot ===\n")
p_dot <- DotPlot(pmc, features = existing_markers, cols = c("lightgrey", "darkblue"), dot.scale = 6) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5)
  ) +
  ggtitle("Canonical Cell Type Marker Expression Across Barrington PMC Clusters")

ggsave("results/figures/barrington_canonical_markers_dotplot.png", p_dot,
       width = 12, height = 7, dpi = 300, bg = "white")
cat("DotPlot saved to: results/figures/barrington_canonical_markers_dotplot.png\n")

cat("=== Step 4: Generating PMC Key Markers FeaturePlot ===\n")
# Highlighting PMC specific markers: Crh, Penk, Slc17a6, Gad1
key_pmc_genes <- c("Crh", "Penk", "Slc17a6", "Gad1")
key_pmc_genes <- key_pmc_genes[key_pmc_genes %in% rownames(pmc)]
if (!length(key_pmc_genes)) stop("None of the configured key PMC genes are present in the RNA assay.")

p_feat <- FeaturePlot(pmc, features = key_pmc_genes, ncol = 2, pt.size = 0.3, cols = c("lightgrey", "#b2182b"))
ggsave("results/figures/barrington_key_markers_featureplot.png", p_feat,
       width = 11, height = 9, dpi = 300, bg = "white")
cat("FeaturePlot saved to: results/figures/barrington_key_markers_featureplot.png\n")

cat("=== Step 5: Exporting Cluster Average Expression Table ===\n")
avg_exp <- AverageExpression(pmc, features = existing_markers, assays = "RNA")$RNA
avg_exp_df <- data.frame(Gene = rownames(avg_exp), avg_exp, check.names = FALSE)

write.table(avg_exp_df, file = "results/tables/canonical_markers_avg_expression.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("=== Step 6: Exporting Cluster Marker Evidence ===\n")
normalized <- GetAssayData(pmc, assay = "RNA", layer = "data")
marker_sets <- data.frame(
  gene = unlist(markers_to_check, use.names = FALSE),
  marker_set = rep(names(markers_to_check), lengths(markers_to_check)),
  stringsAsFactors = FALSE
)
marker_sets <- marker_sets[marker_sets$gene %in% rownames(normalized), , drop = FALSE]
cluster_levels <- sort(unique(cluster_id))
evidence <- do.call(rbind, lapply(seq_len(nrow(marker_sets)), function(i) {
  gene <- marker_sets$gene[i]
  do.call(rbind, lapply(cluster_levels, function(cluster) {
    cells <- which(cluster_id == cluster)
    values <- as.numeric(normalized[gene, cells, drop = TRUE])
    data.frame(
      marker_set = marker_sets$marker_set[i], gene = gene, cluster = cluster,
      cells = length(cells), detected_cells = sum(values > 0),
      percent_detected = 100 * mean(values > 0),
      mean_log_normalized_expression = mean(values),
      stringsAsFactors = FALSE
    )
  }))
}))
write.table(evidence, "results/tables/cluster_canonical_marker_evidence.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
penk_cluster_means <- evidence[evidence$gene == "Penk", c("cluster", "mean_log_normalized_expression")]
penk_cluster_means <- penk_cluster_means[order(-penk_cluster_means$mean_log_normalized_expression), ]

cat("=== Step 7: Finding Descriptive Cluster Markers ===\n")
marker_cell_cap <- 250L
cells_by_cluster <- split(colnames(pmc), cluster_id)
marker_cells_by_cluster <- lapply(cells_by_cluster, function(cells) {
  if (length(cells) > marker_cell_cap) sample(cells, marker_cell_cap) else cells
})
marker_scan_qc <- data.frame(
  cluster = names(cells_by_cluster),
  cells_available = lengths(cells_by_cluster),
  cells_tested = lengths(marker_cells_by_cluster),
  stringsAsFactors = FALSE
)
write.table(marker_scan_qc, "results/tables/marker_scan_sample_counts.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
marker_cells <- unlist(marker_cells_by_cluster, use.names = FALSE)
pmc_marker_scan <- subset(pmc, cells = marker_cells)
all_markers <- FindAllMarkers(
  pmc_marker_scan, assay = "RNA", only.pos = TRUE, test.use = "wilcox",
  min.pct = 0.25, logfc.threshold = 0.25,
  max.cells.per.ident = marker_cell_cap, random.seed = 42
)
write.table(all_markers, "results/tables/cluster_markers_wilcoxon.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(
  "Canonical marker detection and mean log-normalized expression are descriptive evidence, not automatic cell-type calls.",
  paste("FindAllMarkers uses a fixed random seed and tests at most", marker_cell_cap,
        "cells per identity; cells were sampled from the full object before testing to bound memory."),
  paste("Marker scan cells:", length(marker_cells), "of", ncol(pmc), "QC-passing cells."),
  paste("Highest mean Penk cluster:", penk_cluster_means$cluster[1]),
  paste("Top mean log-normalized Penk expression:", signif(penk_cluster_means$mean_log_normalized_expression[1], 4)),
  "FindAllMarkers uses a cell-level Wilcoxon rank-sum test only to describe cluster markers; its p-values are not biological-replicate inference."
), "results/tables/cluster_annotation_review.txt")
writeLines(capture.output(sessionInfo()), "run_logs/annotation_sessionInfo.txt")
rm(pmc_marker_scan, marker_cells, marker_cells_by_cluster, cells_by_cluster)
cat("Canonical marker evidence and descriptive cluster markers saved.\n")

cat("=== Annotation profiling completed. Review evidence before naming clusters. ===\n")
