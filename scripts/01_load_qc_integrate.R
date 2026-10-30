#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

raw_dir <- "data/raw"
pmc_files <- sort(list.files(raw_dir, pattern = "^GSM[0-9]+_Mm_PMC_.*_dge\\.tsv\\.gz$", full.names = TRUE))
if (!length(pmc_files)) stop("No PMC DGE files found in data/raw/. Run scripts/00_download_geo.R first.")
metadata_file <- "metadata/GSE226809_PMC_GEO_sample_manifest.tsv"
if (!file.exists(metadata_file)) stop("GEO sample metadata manifest is missing: ", metadata_file)
geo_metadata <- read.delim(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(geo_metadata) != 30L || anyDuplicated(geo_metadata$geo_accession) ||
    !all(c("geo_accession", "title", "region", "batch", "sex") %in% names(geo_metadata))) {
  stop("GEO PMC metadata manifest must contain 30 unique records with accession, title, region, batch, and sex.")
}
file_accessions <- sub("_.*$", "", basename(pmc_files))
if (!setequal(file_accessions, geo_metadata$geo_accession)) {
  stop("The 30 PMC count files do not match the accession list in the GEO metadata manifest.")
}

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
cache_dir <- file.path("data/processed", "ingest_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

load_dge_as_sparse <- function(file_path, sample_id, min_features = 200L) {
  cache_path <- file.path(cache_dir, paste0(sample_id, "_min_features_", min_features, ".rds"))
  if (file.exists(cache_path)) {
    cached <- readRDS(cache_path)
    if (inherits(cached, "sparseMatrix") && ncol(cached) > 0L && nrow(cached) > 0L) {
      message("Using validated sparse ingest cache for ", sample_id)
      return(cached)
    }
    stop("Invalid sparse ingest cache: ", cache_path)
  }

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
  chunk_size <- 256L

  repeat {
    lines <- readLines(con, n = chunk_size, warn = FALSE)
    if (!length(lines)) break
    fields <- strsplit(lines, "\t", fixed = TRUE)
    if (any(lengths(fields) != n_cells + 1L)) stop("Row width mismatch in ", file_path)
    chunk_genes <- vapply(fields, `[[`, character(1), 1L)
    if (anyNA(chunk_genes) || any(!nzchar(chunk_genes))) {
      stop("Empty or missing feature IDs in ", file_path)
    }
    values <- suppressWarnings(as.numeric(unlist(lapply(fields, `[`, -1L), use.names = FALSE)))
    if (length(values) != length(lines) * n_cells || anyNA(values) ||
        any(!is.finite(values)) || any(values < 0) || any(values != floor(values))) {
      stop("Invalid count values in ", file_path, " near feature row ", row_index + 1L)
    }
    chunk_matrix <- matrix(values, nrow = length(lines), ncol = n_cells, byrow = TRUE)
    nonzero <- which(chunk_matrix > 0, arr.ind = TRUE)
    if (nrow(nonzero)) {
      i_chunks[[length(i_chunks) + 1L]] <- as.integer(nonzero[, 1L] + row_index)
      j_chunks[[length(j_chunks) + 1L]] <- as.integer(nonzero[, 2L])
      x_chunks[[length(x_chunks) + 1L]] <- chunk_matrix[nonzero]
    }
    feature_count <- feature_count + as.integer(colSums(chunk_matrix > 0))
    genes <- c(genes, chunk_genes)
    row_index <- row_index + length(lines)
    rm(fields, values, chunk_matrix, nonzero, chunk_genes)
  }
  if (!row_index) stop("No genes found in ", file_path)
  keep_cells <- which(feature_count >= min_features)
  if (!length(keep_cells)) return(NULL)

  i <- unlist(i_chunks, use.names = FALSE)
  j_original <- unlist(j_chunks, use.names = FALSE)
  x <- unlist(x_chunks, use.names = FALSE)
  cell_map <- integer(n_cells)
  cell_map[keep_cells] <- seq_along(keep_cells)
  mapped_j <- cell_map[j_original]
  keep_nz <- mapped_j > 0L
  sparse_counts <- Matrix::sparseMatrix(
    i = i[keep_nz], j = mapped_j[keep_nz], x = x[keep_nz],
    dims = c(row_index, length(keep_cells)),
    dimnames = list(make.unique(gsub("_", "-", genes)),
                    paste0(sample_id, "_", barcodes[keep_cells]))
  )
  saveRDS(sparse_counts, cache_path, compress = FALSE)
  sparse_counts
}

library_cache_paths <- rep(NA_character_, length(pmc_files))
library_cell_metadata <- vector("list", length(pmc_files))
sample_manifest <- data.frame(
  library_id = character(), batch_id = character(), sex = character(),
  retained_cells = integer(), stringsAsFactors = FALSE
)

for (i in seq_along(pmc_files)) {
  file_path <- pmc_files[i]
  file_name <- basename(file_path)
  library_id <- sub("_.*$", "", file_name)
  metadata_row <- match(library_id, geo_metadata$geo_accession)
  metadata <- geo_metadata[metadata_row, , drop = FALSE]
  file_title <- sub("_dge\\.tsv\\.gz$", "", sub("^GSM[0-9]+_", "", file_name))
  if (is.na(metadata_row) || metadata$region != "Bar" || metadata$title != file_title) {
    stop("GEO metadata does not match this PMC count file: ", file_name)
  }
  batch_id <- as.character(metadata$batch)
  sex <- if (tolower(metadata$sex) == "female") "Female" else if (tolower(metadata$sex) == "male") "Male" else "Unknown"
  if (!grepl("^Bar_[0-9]+[FM]$", batch_id) || sex == "Unknown") {
    stop("Invalid batch or sex metadata for ", library_id)
  }

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

  library_cache_paths[i] <- file.path(cache_dir, paste0(library_id, "_min_features_200.rds"))
  library_cell_metadata[[i]] <- data.frame(
    sample_id = rep(batch_id, ncol(mat)), batch_id = rep(batch_id, ncol(mat)),
    library_id = rep(library_id, ncol(mat)), sex = rep(sex, ncol(mat)),
    row.names = colnames(mat), stringsAsFactors = FALSE
  )
  sample_manifest <- rbind(sample_manifest, data.frame(
    library_id = library_id, batch_id = batch_id, sex = sex,
    retained_cells = ncol(mat), stringsAsFactors = FALSE
  ))
  rm(mat)
}
active_libraries <- which(!is.na(library_cache_paths))
if (!length(active_libraries)) stop("No sample passed the initial cell filter.")
write.table(sample_manifest, "results/tables/GSE226809_PMC_library_manifest.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
retained_libraries <- sample_manifest[sample_manifest$retained_cells > 0L, , drop = FALSE]
batch_manifest <- unique(retained_libraries[, c("batch_id", "sex"), drop = FALSE])
if (anyDuplicated(batch_manifest$batch_id)) stop("A GEO batch maps to more than one sex.")
if (any(table(batch_manifest$sex) < 2L)) stop("Need at least two GEO batches per sex for pseudobulk DE.")
message("Using ", nrow(batch_manifest), " GEO batch units; GSM libraries within a batch will be aggregated.")

active_cache_paths <- library_cache_paths[active_libraries]
all_features <- character()
all_cell_names <- unlist(lapply(library_cell_metadata[active_libraries], rownames), use.names = FALSE)
cell_metadata <- do.call(rbind, library_cell_metadata[active_libraries])
if (anyDuplicated(all_cell_names) || !identical(rownames(cell_metadata), all_cell_names)) {
  stop("Cell identifiers are duplicated or misaligned before sparse matrix assembly.")
}

library_nnz <- integer(length(active_libraries))
for (k in seq_along(active_libraries)) {
  mat <- readRDS(active_cache_paths[k])
  if (!inherits(mat, "sparseMatrix") || !ncol(mat) || !nrow(mat)) {
    stop("Invalid or empty cached matrix: ", active_cache_paths[k])
  }
  if (!identical(colnames(mat), rownames(library_cell_metadata[[active_libraries[k]]]))) {
    stop("Cached cell names do not match metadata for ", sample_manifest$library_id[active_libraries[k]])
  }
  library_nnz[k] <- Matrix::nnzero(mat)
  all_features <- union(all_features, rownames(mat))
  rm(mat)
  if (k %% 5L == 0L) gc(verbose = FALSE)
}
if (!length(all_features) || anyDuplicated(all_features)) {
  stop("Feature identifiers are missing or duplicated before sparse matrix assembly.")
}

total_nnz <- sum(as.double(library_nnz))
if (!is.finite(total_nnz) || total_nnz < 1 || total_nnz > .Machine$integer.max) {
  stop("Unexpected total nonzero count; sparse matrix assembly cannot continue safely.")
}
i_values <- integer(total_nnz)
j_values <- integer(total_nnz)
x_values <- numeric(total_nnz)
position <- 1L
column_offset <- 0L
for (k in seq_along(active_libraries)) {
  mat <- readRDS(active_cache_paths[k])
  entries <- Matrix::summary(mat)
  end_position <- position + nrow(entries) - 1L
  if (nrow(entries) != library_nnz[k] || end_position > total_nnz) {
    stop("Nonzero count changed while combining ", sample_manifest$library_id[active_libraries[k]])
  }
  feature_map <- match(rownames(mat), all_features)
  if (anyNA(feature_map)) stop("Feature mapping failed while combining ", sample_manifest$library_id[active_libraries[k]])
  i_values[position:end_position] <- feature_map[entries$i]
  j_values[position:end_position] <- as.integer(entries$j + column_offset)
  x_values[position:end_position] <- entries$x
  position <- end_position + 1L
  column_offset <- column_offset + ncol(mat)
  rm(mat, entries, feature_map)
  if (k %% 5L == 0L) gc(verbose = FALSE)
}
if (position != total_nnz + 1L || column_offset != length(all_cell_names)) {
  stop("Sparse matrix assembly did not consume all cached counts and cells.")
}

combined_counts <- Matrix::sparseMatrix(
  i = i_values,
  j = j_values,
  x = x_values,
  dims = c(length(all_features), length(all_cell_names)),
  dimnames = list(all_features, all_cell_names)
)
rm(i_values, j_values, x_values, library_nnz, active_cache_paths)
gc(verbose = FALSE)
merged_pmc <- CreateSeuratObject(
  counts = combined_counts, project = "Barrington_PMC", min.cells = 0,
  meta.data = cell_metadata
)
rm(combined_counts, library_cell_metadata, cell_metadata, all_cell_names)
gc(verbose = FALSE)
mt_features <- grep("^mt-", rownames(merged_pmc), ignore.case = TRUE, value = TRUE)
hb_features <- grep("^Hb[ab]", rownames(merged_pmc), ignore.case = TRUE, value = TRUE)
if (!length(mt_features)) stop("No mitochondrial features detected in the merged count matrix.")
if (!length(hb_features)) stop("No hemoglobin features detected in the merged count matrix.")
merged_pmc[["percent.mt"]] <- PercentageFeatureSet(merged_pmc, features = mt_features)
merged_pmc[["percent.hb"]] <- PercentageFeatureSet(merged_pmc, features = hb_features)

qc_plot <- VlnPlot(merged_pmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb"), group.by = "sex", ncol = 4, pt.size = 0)
ggsave("results/figures/qc_vlnplot_raw.png", qc_plot, width = 14, height = 4.8, dpi = 300, bg = "white")

qc_keep <- with(merged_pmc[[]], nFeature_RNA >= 400 & nFeature_RNA <= 10000 &
                  percent.mt <= 10 & percent.hb <= 5)
qc_summary <- data.frame(
  metric = c("cells_after_200_feature_prefilter", "cells_retained_at_publication_qc"),
  cells = c(ncol(merged_pmc), sum(qc_keep))
)
write.table(qc_summary, "results/tables/qc_cell_counts.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)
pmc_filtered <- subset(merged_pmc, cells = colnames(merged_pmc)[qc_keep])
pmc_filtered <- NormalizeData(pmc_filtered, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
pmc_filtered <- FindVariableFeatures(pmc_filtered, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
features <- VariableFeatures(pmc_filtered)
pmc_filtered <- ScaleData(pmc_filtered, features = features, verbose = FALSE)
pmc_filtered <- RunPCA(pmc_filtered, features = features, npcs = 30, seed.use = 42, verbose = FALSE)
pmc_filtered <- FindNeighbors(pmc_filtered, dims = 1:20, verbose = FALSE)
pmc_filtered <- FindClusters(pmc_filtered, resolution = 0.4, random.seed = 42, verbose = FALSE)
pmc_filtered <- RunUMAP(pmc_filtered, dims = 1:20, seed.use = 42, verbose = FALSE)

saveRDS(pmc_filtered, "data/processed/pmc_barrington_seurat.rds")
p1 <- DimPlot(pmc_filtered, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle("Barrington's nucleus clusters")
p2 <- DimPlot(pmc_filtered, reduction = "umap", group.by = "sex") + ggtitle("Sample sex")
ggsave("results/figures/umap_clusters_sex.png", p1 + p2, width = 14, height = 6, dpi = 300, bg = "white")
writeLines(capture.output(sessionInfo()), "run_logs/sessionInfo.txt")
message("QC, clustering, and sample-level metadata saved.")
