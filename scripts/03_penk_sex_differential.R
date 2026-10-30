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
dir.create("run_logs", recursive = TRUE, showWarnings = FALSE)

pmc <- readRDS(input_rds)
DefaultAssay(pmc) <- "RNA"
if (!"Penk" %in% rownames(pmc)) stop("Penk is absent from the current RNA assay.")
penk_expression <- FetchData(pmc, vars = "Penk", layer = "data")[["Penk"]]
penk_by_cluster <- tapply(penk_expression, as.character(Idents(pmc)), mean, na.rm = TRUE)
if (!length(penk_by_cluster) || any(!is.finite(penk_by_cluster))) {
  stop("Penk expression is missing or non-finite across clusters.")
}
top_penk_clusters <- names(penk_by_cluster)[penk_by_cluster == max(penk_by_cluster)]
if (length(top_penk_clusters) != 1L) stop("The highest-Penk cluster is tied; review cluster evidence before DE.")
cluster_id <- top_penk_clusters[[1L]]

marker_genes <- intersect(c("Penk", "Crh", "Slc17a6", "Slc32a1", "Gad2", "Mbp", "Mog"), rownames(pmc))
marker_expression <- FetchData(pmc, vars = marker_genes, layer = "data")
cluster_labels <- as.character(Idents(pmc))
cluster_selection <- do.call(rbind, lapply(sort(unique(cluster_labels)), function(cluster) {
  cells <- which(cluster_labels == cluster)
  values <- marker_expression[cells, , drop = FALSE]
  row <- data.frame(cluster = cluster, cells = length(cells), stringsAsFactors = FALSE)
  for (gene in marker_genes) {
    row[[paste0(gene, "_percent_detected")]] <- 100 * mean(values[[gene]] > 0)
    row[[paste0(gene, "_mean_log_normalized_expression")]] <- mean(values[[gene]])
  }
  row
}))
cluster_selection$selected_as_highest_Penk <- cluster_selection$cluster == cluster_id
write.table(cluster_selection, "results/tables/penk_enriched_cluster_selection.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Selected candidate cluster ", cluster_id, " by highest mean log-normalized Penk expression.")

penk_cluster <- subset(pmc, idents = cluster_id)
meta <- penk_cluster[[]][colnames(penk_cluster), , drop = FALSE]
if (!all(c("sample_id", "sex") %in% names(meta))) stop("Object metadata must contain sample_id and sex.")
if (anyNA(meta$sample_id) || anyNA(meta$sex)) stop("Missing sample_id or sex metadata.")

sample_meta <- unique(meta[, c("sample_id", "sex"), drop = FALSE])
if (anyDuplicated(sample_meta$sample_id)) stop("Each sample must map to exactly one sex.")
sample_meta$sex <- factor(sample_meta$sex, levels = c("Male", "Female"))
sample_meta <- sample_meta[!is.na(sample_meta$sex), , drop = FALSE]
if (any(table(sample_meta$sex) < 2L)) stop("Pseudobulk DE requires at least two GEO batch-pool units per sex.")

counts <- GetAssayData(penk_cluster, assay = "RNA", layer = "counts")
sample_factor <- factor(meta$sample_id, levels = sample_meta$sample_id)
aggregation <- Matrix::sparse.model.matrix(~ 0 + sample_factor)
colnames(aggregation) <- levels(sample_factor)
pseudobulk <- as.matrix(counts %*% aggregation)
if (anyNA(pseudobulk) || any(!is.finite(pseudobulk)) || any(pseudobulk < 0) ||
    any(pseudobulk != floor(pseudobulk))) stop("Pseudobulk counts must be finite non-negative integers.")
storage.mode(pseudobulk) <- "integer"
sample_meta <- sample_meta[match(colnames(pseudobulk), sample_meta$sample_id), , drop = FALSE]
rownames(sample_meta) <- sample_meta$sample_id
sample_meta$cluster_cells <- as.integer(table(factor(meta$sample_id, levels = sample_meta$sample_id)))
sample_meta$total_UMIs <- as.integer(colSums(pseudobulk))
write.table(sample_meta, "results/tables/penk_enriched_cluster_pseudobulk_sample_qc.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

dds <- DESeqDataSetFromMatrix(countData = pseudobulk, colData = sample_meta, design = ~ sex)
keep <- rowSums(counts(dds) >= 10) >= 2L
dds <- dds[keep, ]
if (nrow(dds) < 2L) stop("Too few genes remain after pseudobulk filtering.")
dds <- DESeq(dds)
contrast <- c("sex", "Female", "Male")
res <- results(dds, contrast = contrast)
res_threshold <- results(dds, contrast = contrast, lfcThreshold = 0.5,
                         altHypothesis = "greaterAbs")
deg <- as.data.frame(res)
deg$gene <- rownames(deg)
deg$pvalue_lfc_threshold <- res_threshold$pvalue
deg$padj_lfc_threshold <- res_threshold$padj
deg$significant <- !is.na(deg$padj_lfc_threshold) & deg$padj_lfc_threshold < 0.05 &
  abs(deg$log2FoldChange) > 0.5
deg <- deg[order(deg$padj, na.last = TRUE), , drop = FALSE]

write.table(deg, "results/tables/penk_enriched_cluster_sex_degs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(
  "GEO Bar batch-pool pseudobulk DESeq2",
  paste("GEO Bar batch-pool units:", nrow(sample_meta)),
  paste("Male batch-pool units:", sum(sample_meta$sex == "Male")),
  paste("Female batch-pool units:", sum(sample_meta$sex == "Female")),
  "GSM libraries sharing a Bar_* GEO label are aggregated before differential expression.",
  "One pseudobulk column represents one GEO Bar batch pool; nuclei and GSM libraries within that pool are not treated as independent replicates.",
  paste("Unpaired comparison:", sum(sample_meta$sex == "Female"), "female and",
        sum(sample_meta$sex == "Male"), "male GEO batch pools after QC."),
  "The GEO DGE matrices expose Bar_* batch labels but no individual-animal IDs or paired animal-level measurements. This unpaired batch-group contrast is exploratory and cannot recover within-group animal variation.",
  "Contrast: Female versus Male",
  "Candidate threshold: DESeq2 Wald test of |log2FoldChange| > 0.5, BH-adjusted p < 0.05.",
  "The standard zero-effect Wald p-value and the threshold-test p-value are both reported.",
  paste("Cluster identity:", cluster_id, "selected by highest mean Penk expression in this workflow; it is a Penk-enriched candidate, not a validated match to the published excitatory Bar-Penk cluster."),
  "Canonical marker evidence shows mixed inhibitory and oligodendrocyte-associated signal in the selected cluster; interpret its cell identity cautiously."
), "results/tables/penk_enriched_cluster_sex_degs_summary.txt")

plot_data <- deg[!is.na(deg$padj_lfc_threshold), , drop = FALSE]
plot_data$neg_log10_padj <- -log10(pmax(plot_data$padj_lfc_threshold, .Machine$double.xmin))
top_labels <- plot_data %>% filter(significant) %>%
  arrange(padj_lfc_threshold, desc(abs(log2FoldChange))) %>% head(15)
volcano <- ggplot(plot_data, aes(log2FoldChange, neg_log10_padj, color = significant)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#d95f02")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(data = top_labels, aes(label = gene), max.overlaps = 20) +
labs(title = paste0("Penk-enriched candidate cluster ", cluster_id, ": batch-pool sex contrast"),
       x = "DESeq2 log2 fold change (Female versus Male)",
       y = "-log10 BH-adjusted p-value (|LFC| > 0.5 test)", color = "Significant") +
  theme_bw()
ggsave("results/figures/penk_enriched_cluster_sex_volcano.png", volcano,
       width = 9, height = 7, dpi = 300, bg = "white")
writeLines(capture.output(sessionInfo()), "run_logs/differential_sessionInfo.txt")
message("GEO batch-pool pseudobulk results saved.")
