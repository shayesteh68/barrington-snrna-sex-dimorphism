#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(GEOquery))
raw_dir <- "data/raw"
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

GEOquery::getGEOSuppFiles("GSE226809", baseDir = raw_dir, makeDirectory = FALSE)
archives <- list.files(raw_dir, pattern = "\\.tar(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
for (archive in archives) {
  message("Extracting ", basename(archive))
  utils::untar(archive, exdir = raw_dir)
}

sample_files <- list.files(raw_dir, pattern = "^GSM[0-9]+_Mm_PMC_.*_dge\\.tsv\\.gz$", full.names = TRUE)
if (!length(sample_files)) {
  stop("GEO files were fetched, but expected PMC DGE files were not found in data/raw/. Inspect the GEO supplementary archive.")
}
message("Found ", length(sample_files), " PMC DGE files in ", raw_dir, ".")
