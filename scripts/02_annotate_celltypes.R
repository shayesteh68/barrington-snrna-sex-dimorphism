suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(Matrix)
})

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

cat("=== Step 3: Generating Marker DotPlot ===\n")
p_dot <- DotPlot(pmc, features = existing_markers, cols = c("lightgrey", "darkblue"), dot.scale = 6) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5)
  ) +
  ggtitle("Canonical Cell Type Marker Expression Across Barrington PMC Clusters")

ggsave("results/figures/barrington_canonical_markers_dotplot.png", p_dot, width = 12, height = 7, dpi = 300)
cat("DotPlot saved to: results/figures/barrington_canonical_markers_dotplot.png\n")

cat("=== Step 4: Generating PMC Key Markers FeaturePlot ===\n")
# Highlighting PMC specific markers: Crh, Penk, Slc17a6, Gad1
key_pmc_genes <- c("Crh", "Penk", "Slc17a6", "Gad1")
key_pmc_genes <- key_pmc_genes[key_pmc_genes %in% rownames(pmc)]
if (!length(key_pmc_genes)) stop("None of the configured key PMC genes are present in the RNA assay.")

p_feat <- FeaturePlot(pmc, features = key_pmc_genes, ncol = 2, pt.size = 0.3, cols = c("lightgrey", "#b2182b"))
ggsave("results/figures/barrington_key_markers_featureplot.png", p_feat, width = 11, height = 9, dpi = 300)
cat("FeaturePlot saved to: results/figures/barrington_key_markers_featureplot.png\n")

cat("=== Step 5: Exporting Cluster Average Expression Table ===\n")
avg_exp <- AverageExpression(pmc, features = existing_markers, assays = "RNA")$RNA
avg_exp_df <- data.frame(Gene = rownames(avg_exp), avg_exp, check.names = FALSE)

write.table(avg_exp_df, file = "results/tables/canonical_markers_avg_expression.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Average expression table saved to: results/tables/canonical_markers_avg_expression.tsv\n")

cat("=== Step 2 Annotation Profiling Completed Successfully! ===\n")
