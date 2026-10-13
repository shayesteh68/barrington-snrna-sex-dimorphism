#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggplot2)
})

deg_file <- "results/tables/penk_neurons_sex_degs.tsv"
if (!file.exists(deg_file)) {
  stop("DEG file not found: ", deg_file)
}

deg_df <- read.delim(deg_file, stringsAsFactors = FALSE)
cat(sprintf("Loaded %d genes from %s\n", nrow(deg_df), deg_file))

# Separate significant genes (p_val_adj < 0.05)
sig_genes <- deg_df[deg_df$p_val_adj < 0.05 & !is.na(deg_df$p_val_adj), ]
cat(sprintf("Identified %d significant DEGs (adj.p < 0.05)\n", nrow(sig_genes)))

# Map gene symbols to Entrez IDs
sig_mapped <- bitr(
  sig_genes$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

sig_merged <- merge(sig_genes, sig_mapped, by.x = "gene", by.y = "SYMBOL")

female_entrez <- sig_merged$ENTREZID[sig_merged$avg_log2FC > 0]
male_entrez <- sig_merged$ENTREZID[sig_merged$avg_log2FC < 0]

cat(sprintf("Female-biased mapped: %d | Male-biased mapped: %d\n", length(female_entrez), length(male_entrez)))

# 1. Over-Representation Analysis (ORA) for GO Biological Process
if (length(sig_merged$ENTREZID) > 0) {
  go_ora <- enrichGO(
    gene = sig_merged$ENTREZID,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  if (!is.null(go_ora) && nrow(as.data.frame(go_ora)) > 0) {
    write.table(
      as.data.frame(go_ora),
      file = "results/tables/penk_sex_go_enrichment.tsv",
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    p_go <- dotplot(go_ora, showCategory = 15, title = "GO BP Enrichment: Penk+ Sex DEGs")
    ggsave("results/figures/penk_sex_go_dotplot.png", plot = p_go, width = 8, height = 6, dpi = 300)
    cat("GO ORA enrichment completed and saved.\n")
  } else {
    cat("No significant terms found via strict ORA cutoffs (expected with sparse DEG list).\n")
  }
}

# 2. Gene Set Enrichment Analysis (GSEA) on full ranked gene list
all_mapped <- bitr(
  deg_df$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)
deg_all_merged <- merge(deg_df, all_mapped, by.x = "gene", by.y = "SYMBOL")
deg_all_merged <- deg_all_merged[!duplicated(deg_all_merged$ENTREZID), ]

gene_list <- deg_all_merged$avg_log2FC
names(gene_list) <- deg_all_merged$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

cat("Running GSEA (Gene Ontology: BP) on ranked gene list...\n")
set.seed(42)
gsea_go <- gseGO(
  geneList = gene_list,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)

if (!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
  write.table(
    as.data.frame(gsea_go),
    file = "results/tables/penk_sex_gsea_go.tsv",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  p_gsea <- dotplot(gsea_go, showCategory = 15, title = "GSEA GO BP: Sex-Biased Profiles in Penk+ Neurons")
  ggsave("results/figures/penk_sex_gsea_dotplot.png", plot = p_gsea, width = 9, height = 7, dpi = 300)
  cat("GSEA GO analysis completed and saved.\n")
} else {
  cat("GSEA completed; no pathways passed cutoff.\n")
}

cat("Pipeline step 04 execution finished.\n")
