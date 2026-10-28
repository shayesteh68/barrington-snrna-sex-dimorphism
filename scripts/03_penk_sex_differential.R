#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
})

input_rds <- "data/processed/pmc_barrington_seurat.rds"
if (!file.exists(input_rds)) stop("Processed object missing: ", input_rds)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

pmc <- readRDS(input_rds)
DefaultAssay(pmc) <- "RNA"
cluster_id <- "8"
if (!cluster_id %in% levels(Idents(pmc))) stop("Cluster 8 is absent; review the cluster annotation before DE.")
if (!"Penk" %in% rownames(pmc)) stop("Penk is absent from the current RNA assay.")
penk_expression <- FetchData(pmc, vars = "Penk", layer = "data")[["Penk"]]
penk_by_cluster <- tapply(penk_expression, as.character(Idents(pmc)), mean, na.rm = TRUE)
top_penk_cluster <- names(penk_by_cluster)[which.max(penk_by_cluster)]
if (!identical(top_penk_cluster, cluster_id)) {
  stop("Cluster 8 is not the current Penk-high cluster; review markers and update the selected cluster before DE.")
}
penk <- subset(pmc, idents = cluster_id)
meta <- penk[[]][colnames(penk), , drop = FALSE]
if (!all(c("sample_id", "sex") %in% names(meta))) stop("Object metadata must contain sample_id and sex.")
if (anyNA(meta$sample_id) || anyNA(meta$sex)) stop("Missing sample_id or sex metadata.")

sample_meta <- unique(meta[, c("sample_id", "sex"), drop = FALSE])
if (anyDuplicated(sample_meta$sample_id)) stop("Each sample must map to exactly one sex.")
sample_meta$sex <- factor(sample_meta$sex, levels = c("Male", "Female"))
sample_meta <- sample_meta[!is.na(sample_meta$sex), , drop = FALSE]
if (any(table(sample_meta$sex) < 2L)) stop("Pseudobulk DE requires at least two biological samples per sex.")

counts <- GetAssayData(penk, assay = "RNA", layer = "counts")
sample_factor <- factor(meta$sample_id, levels = sample_meta$sample_id)
aggregation <- Matrix::sparse.model.matrix(~ 0 + sample_factor)
colnames(aggregation) <- levels(sample_factor)
pseudobulk <- as.matrix(counts %*% aggregation)
storage.mode(pseudobulk) <- "integer"
sample_meta <- sample_meta[match(colnames(pseudobulk), sample_meta$sample_id), , drop = FALSE]
rownames(sample_meta) <- sample_meta$sample_id

dds <- DESeqDataSetFromMatrix(countData = pseudobulk, colData = sample_meta, design = ~ sex)
keep <- rowSums(counts(dds) >= 10) >= 2L
dds <- dds[keep, ]
if (nrow(dds) < 2L) stop("Too few genes remain after pseudobulk filtering.")
dds <- DESeq(dds)
res <- results(dds, contrast = c("sex", "Female", "Male"))
deg <- as.data.frame(res)
deg$gene <- rownames(deg)
deg$significant <- !is.na(deg$padj) & deg$padj < 0.05 & abs(deg$log2FoldChange) >= 0.5
deg <- deg[order(deg$padj, na.last = TRUE), , drop = FALSE]

write.table(deg, "results/tables/penk_neurons_sex_degs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(
  "GEO-batch-level pseudobulk DESeq2",
  paste("Biological GEO batches:", nrow(sample_meta)),
  paste("Male GEO batches:", sum(sample_meta$sex == "Male")),
  paste("Female GEO batches:", sum(sample_meta$sex == "Female")),
  "GSM libraries within the same Bar_* batch are aggregated before differential expression.",
  "Each GEO batch pools multiple mice; cells and GSM libraries are not treated as independent replicates.",
  "Contrast: Female versus Male",
  "Threshold: padj < 0.05 and |log2FoldChange| >= 0.5",
  "Cluster identity: cluster 8; current Penk expression check passed. Review after reclustering."
), "results/tables/penk_neurons_sex_degs_summary.txt")

plot_data <- deg[!is.na(deg$padj), , drop = FALSE]
plot_data$neg_log10_padj <- -log10(pmax(plot_data$padj, .Machine$double.xmin))
top_labels <- plot_data %>% filter(significant) %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(15)
volcano <- ggplot(plot_data, aes(log2FoldChange, neg_log10_padj, color = significant)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#d95f02")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(data = top_labels, aes(label = gene), max.overlaps = 20) +
  labs(title = "Penk-high cluster: GEO-batch-level sex contrast",
       x = "DESeq2 log2 fold change (Female versus Male)",
       y = "-log10 adjusted p-value", color = "Significant") +
  theme_bw()
ggsave("results/figures/penk_neurons_sex_volcano.png", volcano, width = 9, height = 7, dpi = 300)
message("GEO-batch-level pseudobulk results saved.")
