suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
})

cat("=== Step 1: Loading Processed Seurat Object ===\n")
pmc <- readRDS("data/processed/pmc_barrington_seurat.rds")
DefaultAssay(pmc) <- "RNA"

cat("=== Step 2: Isolating Penk+ Barrington Nucleus Neurons (Cluster 8) ===\n")
# Subset Cluster 8
penk_neurons <- subset(pmc, idents = 8)

cat("Cluster 8 total nuclei count:", ncol(penk_neurons), "\n")
sex_counts <- table(penk_neurons$sex)
print(sex_counts)

cat("=== Step 3: Running Sex-Differential Expression Analysis (Female vs Male) ===\n")
Idents(penk_neurons) <- "sex"

# Perform differential expression testing using Wilcoxon rank-sum test
# comparing Female vs Male within Penk+ neurons
deg_sex <- FindMarkers(
  object = penk_neurons,
  ident.1 = "Female",
  ident.2 = "Male",
  min.pct = 0.1,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

# Format DEG table
deg_sex$gene <- rownames(deg_sex)
deg_sex$significant <- ifelse(deg_sex$p_val_adj < 0.05 & abs(deg_sex$avg_log2FC) >= 0.5, "Significant", "Not Significant")

cat("Significant Sex-DEGs detected (FDR < 0.05, |log2FC| >= 0.5):", sum(deg_sex$significant == "Significant"), "\n")

# Save table
write.table(deg_sex, file = "results/tables/penk_neurons_sex_degs.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("DEG results saved to: results/tables/penk_neurons_sex_degs.tsv\n")

cat("=== Step 4: Generating Publication-Grade Volcano Plot ===\n")
# Select top genes for labeling
top_labeled_genes <- deg_sex %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC))) %>%
  head(15)

p_volcano <- ggplot(deg_sex, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = significant), alpha = 0.7, size = 1.8) +
  scale_color_manual(values = c("Not Significant" = "#bdbdbd", "Significant" = "#d95f02")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40", size = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", size = 0.5) +
  geom_text_repel(
    data = top_labeled_genes,
    aes(label = gene),
    size = 4,
    fontface = "bold.italic",
    box.padding = 0.4,
    max.overlaps = 20
  ) +
  labs(
    title = "Sex-Biased Gene Expression in Penk+ Barrington Nucleus Neurons",
    subtitle = "Single-Nucleus RNA-seq Analysis (PMC Cohort: Female vs Male)",
    x = "log2 Fold-Change (Positive = Female Enriched, Negative = Male Enriched)",
    y = "-log10 (FDR-adjusted p-value)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey30"),
    axis.title = element_text(face = "bold")
  )

ggsave("results/figures/penk_neurons_sex_volcano.png", p_volcano, width = 9, height = 7, dpi = 300)
cat("Volcano plot saved to: results/figures/penk_neurons_sex_volcano.png\n")

cat("=== Penk Sex DEG Pipeline Complete! ===\n")
