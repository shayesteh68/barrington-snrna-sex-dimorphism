suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(data.table)
    library(ggplot2)
    library(dplyr)
    library(patchwork)
})

raw_dir <- "data/raw"
pmc_files <- sort(list.files(raw_dir, pattern = "Mm_PMC_.*_dge\\.tsv\\.gz$", full.names = TRUE))

if(length(pmc_files) == 0) {
    stop("No PMC dge files found in ", raw_dir)
}

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

cat("=== Step 1: Scanning and Loading PMC Samples (Sparse-Safe Loader) ===\n")
cat("Found", length(pmc_files), "PMC sample files.\n")

load_dge_as_sparse <- function(file_path, sample_id) {
    # Read header to count cells
    header <- strsplit(readLines(gzfile(file_path), n = 1), "\t")[[1]]
    barcodes <- header[-1]
    cell_names <- paste0(sample_id, "_", barcodes)
    
    # Read entire table via fread as pure data.table
    dt <- fread(file_path, header = TRUE, sep = "\t", showProgress = FALSE)
    gene_names <- make.unique(gsub("_", "-", dt[[1]]))
    
    # Calculate column sums / detected genes per cell to discard empty droplets immediately
    # Using matrix extraction only for non-empty columns
    dt_vals <- dt[, -1, with = FALSE]
    rm(dt)
    
    # Identify non-empty cells (at least 200 expressed genes)
    col_nz <- vapply(dt_vals, function(col) sum(col > 0), integer(1))
    valid_cols <- which(col_nz >= 200)
    
    if (length(valid_cols) == 0) {
        rm(dt_vals)
        gc(verbose = FALSE)
        return(NULL)
    }
    
    # Subset to keep only valid cells BEFORE dense conversion
    dt_valid <- dt_vals[, valid_cols, with = FALSE]
    rm(dt_vals)
    
    mat <- as(as.matrix(dt_valid), "dgCMatrix")
    rm(dt_valid)
    
    rownames(mat) <- gene_names
    colnames(mat) <- cell_names[valid_cols]
    
    gc(verbose = FALSE)
    return(mat)
}

merged_pmc <- NULL

for (i in seq_along(pmc_files)) {
    f <- pmc_files[i]
    base_name <- basename(f)
    sample_id <- gsub("^GSM[0-9]+_Mm_PMC_", "", base_name)
    sample_id <- gsub("_dge\\.tsv\\.gz$", "", sample_id)

    sex <- ifelse(grepl("[0-9]F-", sample_id), "Female",
           ifelse(grepl("[0-9]M-", sample_id), "Male", "Unknown"))

    cat(sprintf("[%d/%d] Reading %s (%s)...\n", i, length(pmc_files), sample_id, sex))

    mat <- load_dge_as_sparse(f, sample_id)
    if (is.null(mat) || ncol(mat) == 0) {
        cat(sprintf("  Skipping %s (no cells passed 200-feature threshold)\n", sample_id))
        next
    }

    s_obj <- CreateSeuratObject(counts = mat, project = "Barrington_PMC", min.cells = 3)
    rm(mat)

    s_obj$sample_id <- sample_id
    s_obj$sex <- sex
    s_obj[["percent.mt"]] <- PercentageFeatureSet(s_obj, pattern = "^mt-")

    if (is.null(merged_pmc)) {
        merged_pmc <- s_obj
    } else {
        merged_pmc <- merge(merged_pmc, y = s_obj, project = "Barrington_PMC")
    }
    rm(s_obj)
    gc(verbose = FALSE)
}

cat("=== Step 2: Consolidating Layers ===\n")
merged_pmc[["RNA"]] <- JoinLayers(merged_pmc[["RNA"]])
gc(verbose = FALSE)

cat(sprintf("Total raw nuclei: %d | Total genes: %d\n", ncol(merged_pmc), nrow(merged_pmc)))
print(table(merged_pmc$sex))

cat("=== Step 3: Nuclei QC & Filtering ===\n")
qc_plot <- VlnPlot(merged_pmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                   group.by = "sex", ncol = 3, pt.size = 0)
ggsave("results/figures/qc_vlnplot_raw.png",
       plot = qc_plot, width = 11, height = 4.5, dpi = 300)

pmc_filtered <- subset(merged_pmc, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 5)
rm(merged_pmc)
gc(verbose = FALSE)

cat(sprintf("Nuclei passing QC: %d\n", ncol(pmc_filtered)))
print(table(pmc_filtered$sex))

cat("=== Step 4: Normalization & Variable Feature Selection ===\n")
pmc_filtered <- NormalizeData(pmc_filtered, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
pmc_filtered <- FindVariableFeatures(pmc_filtered, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

cat("=== Step 5: Scaling Only Highly Variable Genes ===\n")
top_features <- VariableFeatures(pmc_filtered)
pmc_filtered <- ScaleData(pmc_filtered, features = top_features, verbose = FALSE)
gc(verbose = FALSE)

cat("=== Step 6: PCA, UMAP, and Clustering ===\n")
pmc_filtered <- RunPCA(pmc_filtered, features = top_features, npcs = 30, verbose = FALSE)
pmc_filtered <- FindNeighbors(pmc_filtered, dims = 1:20, verbose = FALSE)
pmc_filtered <- FindClusters(pmc_filtered, resolution = 0.4, verbose = FALSE)
pmc_filtered <- RunUMAP(pmc_filtered, dims = 1:20, verbose = FALSE)

cat("=== Step 7: Saving Processed Seurat Object and Figures ===\n")
saveRDS(pmc_filtered, file = "data/processed/pmc_barrington_seurat.rds")

p1 <- DimPlot(pmc_filtered, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 4) +
      ggtitle("Barrington's Nucleus (PMC) - Clusters") +
      theme_minimal()

p2 <- DimPlot(pmc_filtered, reduction = "umap", group.by = "sex") +
      ggtitle("Barrington's Nucleus (PMC) - Sex") +
      theme_minimal()

combined_plot <- p1 + p2
ggsave("results/figures/umap_clusters_sex.png",
       plot = combined_plot, width = 14, height = 6, dpi = 300)

cat("=== Step 1 Pipeline Completed Successfully! ===\n")
