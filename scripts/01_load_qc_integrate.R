#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(data.table)
    library(ggplot2)
    library(dplyr)
    library(patchwork)
})

cat("=== Step 1: Scanning and Loading PMC Samples ===\n")
raw_dir <- "/home/digital_home/barrington_snrna_verstegen/data/raw"
pmc_files <- list.files(raw_dir, pattern = "Mm_PMC_.*_dge\\.tsv\\.gz$", full.names = TRUE)

if(length(pmc_files) == 0) {
    stop("No PMC dge files found in ", raw_dir)
}

cat("Found", length(pmc_files), "PMC sample files.\n")

seurat_list <- list()

for (i in seq_along(pmc_files)) {
    f <- pmc_files[i]
    base_name <- basename(f)
    sample_id <- gsub("^GSM[0-9]+_Mm_PMC_", "", base_name)
    sample_id <- gsub("_dge\\.tsv\\.gz$", "", sample_id)
    
    sex <- ifelse(grepl("[0-9]F-", sample_id), "Female",
           ifelse(grepl("[0-9]M-", sample_id), "Male", "Unknown"))
    
    cat(sprintf("[%d/%d] Reading %s (%s)...\n", i, length(pmc_files), sample_id, sex))
    
    dt <- fread(f, data.table = FALSE)
    gene_names <- make.unique(gsub("_", "-", dt[[1]]))
    dt[[1]] <- NULL
    mat <- as(as.matrix(dt), "dgCMatrix")
    rownames(mat) <- gene_names
    rm(dt)
    
    s_obj <- CreateSeuratObject(counts = mat, project = "Barrington_PMC", min.cells = 3, min.features = 200)
    rm(mat)
    
    s_obj$sample_id <- sample_id
    s_obj$sex <- sex
    s_obj[["percent.mt"]] <- PercentageFeatureSet(s_obj, pattern = "^mt-")
    
    seurat_list[[sample_id]] <- s_obj
}

cat("=== Step 2: Merging Samples and Consolidating Layers ===\n")
merged_pmc <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list), project = "Barrington_PMC")
rm(seurat_list)
invisible(gc())

# In Seurat v5, join layers into a single sparse matrix to avoid multi-layer memory overhead
merged_pmc[["RNA"]] <- JoinLayers(merged_pmc[["RNA"]])

cat(sprintf("Total raw nuclei: %d | Total genes: %d\n", ncol(merged_pmc), nrow(merged_pmc)))
print(table(merged_pmc$sex))

cat("=== Step 3: Nuclei QC & Filtering ===\n")
qc_plot <- VlnPlot(merged_pmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                   group.by = "sex", ncol = 3, pt.size = 0)
ggsave("/home/digital_home/barrington_snrna_verstegen/results/figures/qc_vlnplot_raw.png", 
       plot = qc_plot, width = 11, height = 4.5, dpi = 300)

pmc_filtered <- subset(merged_pmc, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 5)
rm(merged_pmc)
invisible(gc())

cat(sprintf("Nuclei passing QC: %d\n", ncol(pmc_filtered)))
print(table(pmc_filtered$sex))

cat("=== Step 4: Normalization & Variable Feature Selection ===\n")
pmc_filtered <- NormalizeData(pmc_filtered, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
pmc_filtered <- FindVariableFeatures(pmc_filtered, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

cat("=== Step 5: Scaling Only Highly Variable Genes (Memory Efficient) ===\n")
top_features <- VariableFeatures(pmc_filtered)
pmc_filtered <- ScaleData(pmc_filtered, features = top_features, verbose = FALSE)
invisible(gc())

cat("=== Step 6: PCA, UMAP, and Clustering ===\n")
pmc_filtered <- RunPCA(pmc_filtered, features = top_features, npcs = 30, verbose = FALSE)
pmc_filtered <- FindNeighbors(pmc_filtered, dims = 1:20, verbose = FALSE)
pmc_filtered <- FindClusters(pmc_filtered, resolution = 0.4, verbose = FALSE)
pmc_filtered <- RunUMAP(pmc_filtered, dims = 1:20, verbose = FALSE)

cat("=== Step 7: Saving Processed Seurat Object and Figures ===\n")
saveRDS(pmc_filtered, file = "/home/digital_home/barrington_snrna_verstegen/data/processed/pmc_barrington_seurat.rds")

p1 <- DimPlot(pmc_filtered, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 4) + 
      ggtitle("Barrington's Nucleus (PMC) - Clusters") +
      theme_minimal()

p2 <- DimPlot(pmc_filtered, reduction = "umap", group.by = "sex", cols = c("Female" = "#E64B35", "Male" = "#4DBBD5"), alpha = 0.5) + 
      ggtitle("Nuclei by Sex") +
      theme_minimal()

p_umap <- p1 + p2
ggsave("/home/digital_home/barrington_snrna_verstegen/results/figures/barrington_umap_clusters_sex.png", 
       plot = p_umap, width = 13, height = 5.5, dpi = 300)

cat("=== Pipeline Execution Complete Successfully! ===\n")
