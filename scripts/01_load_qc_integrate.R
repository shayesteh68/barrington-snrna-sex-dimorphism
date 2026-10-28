#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

raw_dir <- "data/raw"
pmc_files <- sort(list.files(raw_dir, pattern = "^GSM[0-9]+_Mm_PMC_.*_dge\\.tsv\\.gz$", full.names = TRUE))
if (!length(pmc_files)) stop("No PMC DGE files found in data/raw/. Run scripts/00_download_geo.R first.")

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

load_dge_as_sparse <- function(file_path, sample_id, min_features = 200L) {
  con <- gzfile(file_path, open = "rt")
  on.exit(close(con), add = TRUE)
  header <- strsplit(readLines(con, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1]]
  if (length(header) < 2L) stop("Invalid DGE header in ", file_path)
  barcodes <- header[-1]
  n_cells <- length(barcodes)
  if (!n_cells || any(!nzchar(barcodes)) || anyDuplicated(barcodes)) {
    stop("DGE header contains no cells or duplicated/empty barcodes in ", file_path)
  }
  genes <- character()
  feature_count <- integer(n_cells)
  i_chunks <- list()
  j_chunks <- list()
  x_chunks <- list()
  row_index <- 0L

  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) break
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) != n_cells + 1L) stop("Row width mismatch in ", file_path)
    row_index <- row_index + 1L
    genes[row_index] <- fields[1]
    values <- suppressWarnings(as.numeric(fields[-1]))
    if (length(values) != n_cells || anyNA(values) || any(!is.finite(values)) ||
        any(values < 0) || any(values != floor(values))) {
      stop("Invalid count values in ", file_path, " at row ", row_index + 1L)
    }
    nz <- which(values != 0)
    if (length(nz)) {
      i_chunks[[length(i_chunks) + 1L]] <- rep.int(row_index, length(nz))
      j_chunks[[length(j_chunks) + 1L]] <- nz
      x_chunks[[length(x_chunks) + 1L]] <- values[nz]
      feature_count <- feature_count + as.integer(values > 0)
    }
  }
  if (!row_index) stop("No genes found in ", file_path)
  keep_cells <- which(feature_count >= min_features)
  if (!length(keep_cells)) return(NULL)

  i <- unlist(i_chunks, use.names = FALSE)
  j_original <- unlist(j_chunks, use.names = FALSE)
  x <- unlist(x_chunks, use.names = FALSE)
  keep_nz <- j_original %in% keep_cells
  j <- match(j_original[keep_nz], keep_cells)
  Matrix::sparseMatrix(
    i = i[keep_nz],
    j = j,
    x = x[keep_nz],
    dims = c(row_index, length(keep_cells)),
    dimnames = list(make.unique(gsub("_", "-", genes)), paste0(sample_id, "_", barcodes[keep_cells]))
  )
}

merged_pmc <- NULL
sample_manifest <- data.frame(
  library_id = character(), batch_id = character(), sex = character(),
  retained_cells = integer(), stringsAsFactors = FALSE
)

for (i in seq_along(pmc_files)) {
  file_path <- pmc_files[i]
  file_name <- basename(file_path)
  library_id <- sub("_.*$", "", file_name)
  batch_code <- sub("^GSM[0-9]+_Mm_PMC_([0-9]+[FM])-.*_dge\\.tsv\\.gz$", "\\1", file_name)
  if (identical(batch_code, file_name)) stop("Could not parse the GEO batch from sample title: ", file_name)
  batch_id <- paste0("Bar_", batch_code)
  sex_code <- substr(batch_code, nchar(batch_code), nchar(batch_code))
  sex <- if (sex_code == "F") "Female" else if (sex_code == "M") "Male" else "Unknown"
  if (sex == "Unknown") stop("Could not infer sex from the GEO batch title: ", file_name)

  message(sprintf("[%d/%d] Loading %s as %s (%s)", i, length(pmc_files), library_id, batch_id, sex))
  mat <- load_dge_as_sparse(file_path, library_id)
  if (is.null(mat) || !ncol(mat)) {
    message("Skipping ", library_id, ": no cells passed the feature threshold.")
    sample_manifest <- rbind(sample_manifest, data.frame(
      library_id = library_id, batch_id = batch_id, sex = sex,
      retained_cells = 0L, stringsAsFactors = FALSE
    ))
    next
  }

  obj <- CreateSeuratObject(counts = mat, project = "Barrington_PMC", min.cells = 3)
  obj$sample_id <- batch_id
  obj$batch_id <- batch_id
  obj$library_id <- library_id
  obj$sex <- sex
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  merged_pmc <- if (is.null(merged_pmc)) obj else merge(merged_pmc, y = obj, project = "Barrington_PMC")
  sample_manifest <- rbind(sample_manifest, data.frame(
    library_id = library_id, batch_id = batch_id, sex = sex,
    retained_cells = ncol(obj), stringsAsFactors = FALSE
  ))
  rm(mat, obj)
  gc(verbose = FALSE)
}
if (is.null(merged_pmc)) stop("No sample passed the initial cell filter.")
write.table(sample_manifest, "results/tables/GSE226809_PMC_library_manifest.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
retained_libraries <- sample_manifest[sample_manifest$retained_cells > 0L, , drop = FALSE]
batch_manifest <- unique(retained_libraries[, c("batch_id", "sex"), drop = FALSE])
if (anyDuplicated(batch_manifest$batch_id)) stop("A GEO batch maps to more than one sex.")
if (any(table(batch_manifest$sex) < 2L)) stop("Need at least two GEO batches per sex for pseudobulk DE.")
message("Using ", nrow(batch_manifest), " GEO batches as biological replicates; GSM libraries within a batch will be aggregated.")

merged_pmc[["RNA"]] <- JoinLayers(merged_pmc[["RNA"]])
qc_plot <- VlnPlot(merged_pmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "sex", ncol = 3, pt.size = 0)
ggsave("results/figures/qc_vlnplot_raw.png", qc_plot, width = 11, height = 4.5, dpi = 300)

pmc_filtered <- subset(merged_pmc, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 5)
pmc_filtered <- NormalizeData(pmc_filtered, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
pmc_filtered <- FindVariableFeatures(pmc_filtered, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
features <- VariableFeatures(pmc_filtered)
pmc_filtered <- ScaleData(pmc_filtered, features = features, verbose = FALSE)
pmc_filtered <- RunPCA(pmc_filtered, features = features, npcs = 30, verbose = FALSE)
pmc_filtered <- FindNeighbors(pmc_filtered, dims = 1:20, verbose = FALSE)
pmc_filtered <- FindClusters(pmc_filtered, resolution = 0.4, verbose = FALSE)
pmc_filtered <- RunUMAP(pmc_filtered, dims = 1:20, verbose = FALSE)

saveRDS(pmc_filtered, "data/processed/pmc_barrington_seurat.rds")
p1 <- DimPlot(pmc_filtered, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle("Barrington's nucleus clusters")
p2 <- DimPlot(pmc_filtered, reduction = "umap", group.by = "sex") + ggtitle("Sample sex")
ggsave("results/figures/umap_clusters_sex.png", p1 + p2, width = 14, height = 6, dpi = 300)
message("QC, clustering, and sample-level metadata saved.")
