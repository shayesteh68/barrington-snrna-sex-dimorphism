#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggplot2)
})

deg_file <- "results/tables/penk_neurons_sex_degs.tsv"
if (!file.exists(deg_file)) stop("Pseudobulk DEG table not found: ", deg_file)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

deg <- read.delim(deg_file, stringsAsFactors = FALSE)
required <- c("gene", "log2FoldChange", "padj")
if (!all(required %in% names(deg))) stop("DEG file must contain: ", paste(required, collapse = ", "))
deg <- deg[!is.na(deg$gene) & !is.na(deg$log2FoldChange), ]
all_map <- bitr(unique(deg$gene), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
universe <- unique(all_map$ENTREZID)
if (!length(universe)) stop("No gene symbols mapped to Entrez IDs.")

for (direction in c("Female", "Male")) {
  if (direction == "Female") {
    selected <- deg$gene[!is.na(deg$padj) & deg$padj < 0.05 & deg$log2FoldChange > 0]
  } else {
    selected <- deg$gene[!is.na(deg$padj) & deg$padj < 0.05 & deg$log2FoldChange < 0]
  }
  if (!length(selected)) {
    message(direction, ": no significant genes; skipping.")
    next
  }
  mapped <- bitr(unique(selected), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  genes <- unique(mapped$ENTREZID)
  if (length(genes) < 2L) {
    message(direction, ": fewer than two significant mapped genes; skipping.")
    next
  }

  go <- enrichGO(gene = genes, universe = universe, OrgDb = org.Mm.eg.db,
                 keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
  if (!is.null(go) && nrow(as.data.frame(go))) {
    write.table(as.data.frame(go), paste0("results/tables/penk_", tolower(direction), "_go.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    ggsave(paste0("results/figures/penk_", tolower(direction), "_go.png"),
           dotplot(go, showCategory = 15), width = 8, height = 6, dpi = 300)
  }
}
message("Direction-specific pseudobulk enrichment completed.")
