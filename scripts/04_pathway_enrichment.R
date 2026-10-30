#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggplot2)
})

deg_file <- "results/tables/penk_enriched_cluster_sex_degs.tsv"
if (!file.exists(deg_file)) stop("Pseudobulk DEG table not found: ", deg_file)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("run_logs", recursive = TRUE, showWarnings = FALSE)

deg <- read.delim(deg_file, stringsAsFactors = FALSE)
required <- c("gene", "log2FoldChange", "pvalue_lfc_threshold", "padj_lfc_threshold")
if (!all(required %in% names(deg))) stop("DEG file must contain: ", paste(required, collapse = ", "))
deg <- deg[!is.na(deg$gene) & !is.na(deg$log2FoldChange) & !is.na(deg$pvalue_lfc_threshold), ]
all_map <- bitr(unique(deg$gene), fromType = "SYMBOL", toType = "ENTREZID",
                OrgDb = org.Mm.eg.db, drop = TRUE)
universe <- unique(all_map$ENTREZID)
if (!length(universe)) stop("No gene symbols mapped to Entrez IDs.")

output_files <- c(
  "results/tables/penk_enriched_cluster_female_go.tsv", "results/tables/penk_enriched_cluster_male_go.tsv",
  "results/figures/penk_enriched_cluster_female_go.png", "results/figures/penk_enriched_cluster_male_go.png"
)
unlink(output_files[file.exists(output_files)])
summary_rows <- list()

for (direction in c("Female", "Male")) {
  if (direction == "Female") {
    selected <- deg$gene[!is.na(deg$padj_lfc_threshold) & deg$padj_lfc_threshold < 0.05 & deg$log2FoldChange > 0.5]
  } else {
    selected <- deg$gene[!is.na(deg$padj_lfc_threshold) & deg$padj_lfc_threshold < 0.05 & deg$log2FoldChange < -0.5]
  }
  mapped <- if (length(selected)) bitr(unique(selected), fromType = "SYMBOL", toType = "ENTREZID",
                                      OrgDb = org.Mm.eg.db, drop = TRUE) else data.frame()
  genes <- if (nrow(mapped)) unique(mapped$ENTREZID) else character()
  status <- "no_significant_directional_genes"
  enriched_terms <- 0L
  if (!length(selected)) {
    message(direction, ": no significant genes; skipping.")
    summary_rows[[direction]] <- data.frame(
      direction = direction, tested_symbols = length(unique(deg$gene)),
      mapped_background_entrez = length(universe), significant_symbols = 0L,
      mapped_significant_entrez = 0L, enriched_terms = 0L, status = status
    )
    next
  }
  if (length(genes) < 3L) {
    message(direction, ": fewer than three significant mapped genes; skipping enrichment.")
    status <- "fewer_than_three_mapped_genes"
    summary_rows[[direction]] <- data.frame(
      direction = direction, tested_symbols = length(unique(deg$gene)),
      mapped_background_entrez = length(universe), significant_symbols = length(unique(selected)),
      mapped_significant_entrez = length(genes), enriched_terms = 0L, status = status
    )
    next
  }

  go <- enrichGO(gene = genes, universe = universe, OrgDb = org.Mm.eg.db,
                 keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
  go_df <- if (!is.null(go)) as.data.frame(go) else data.frame()
  if (nrow(go_df)) {
    enriched_terms <- nrow(go_df)
    status <- "enrichment_terms_found"
    ggsave(paste0("results/figures/penk_enriched_cluster_", tolower(direction), "_go.png"),
           dotplot(go, showCategory = 15), width = 8, height = 6, dpi = 300, bg = "white")
  } else {
    status <- "no_enriched_terms"
  }
  write.table(go_df, paste0("results/tables/penk_enriched_cluster_", tolower(direction), "_go.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  summary_rows[[direction]] <- data.frame(
    direction = direction, tested_symbols = length(unique(deg$gene)),
    mapped_background_entrez = length(universe), significant_symbols = length(unique(selected)),
    mapped_significant_entrez = length(genes), enriched_terms = enriched_terms, status = status
  )
}
summary <- do.call(rbind, summary_rows)
summary$mapping_fraction_background <- length(unique(all_map$ENTREZID)) / length(unique(deg$gene))
write.table(summary, "results/tables/penk_enriched_cluster_pathway_enrichment_summary.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), "run_logs/enrichment_sessionInfo.txt")
message("Direction-specific pseudobulk enrichment completed using the tested mapped-gene universe.")
