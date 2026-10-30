#!/usr/bin/env Rscript

raw_dir <- "data/raw"
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

archive <- file.path(raw_dir, "GSE226809_RAW.tar")
archive_entries <- function(path) {
  if (!file.exists(path)) return(character())
  tryCatch(utils::untar(path, list = TRUE), error = function(e) character())
}

valid_archive <- function(path) {
  entries <- archive_entries(path)
  dge_entries <- grep("^GSM[0-9]+_.*_dge\\.tsv\\.gz$", entries, value = TRUE)
  pmc_entries <- grep("^GSM[0-9]+_Mm_PMC_.*_dge\\.tsv\\.gz$", entries, value = TRUE)
  length(entries) == 72L && length(dge_entries) == 72L && length(pmc_entries) == 30L
}

if (!valid_archive(archive)) {
  if (file.exists(archive)) {
    stop("Existing GEO archive is incomplete or invalid: ", archive,
         ". Preserve it, move it aside, and rerun the downloader.")
  }
  archive_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE226nnn/GSE226809/suppl/GSE226809_RAW.tar"
  curl_bin <- Sys.which("curl")
  if (nzchar(curl_bin)) {
    partial_archive <- paste0(archive, ".part")
    status <- system2(curl_bin, c(
      "--location", "--fail", "--retry", "10", "--retry-all-errors",
      "--continue-at", "-", "--output", partial_archive, archive_url
    ))
    if (status != 0L || !valid_archive(partial_archive)) {
      stop("GEO archive download did not pass the 72-file/30-PMC archive checks. ",
           "A partial download is kept at ", partial_archive, " for a resumable retry.")
    }
    if (!file.rename(partial_archive, archive)) {
      stop("Downloaded archive is valid but could not be moved into place: ", archive)
    }
  } else {
    suppressPackageStartupMessages(library(GEOquery))
    GEOquery::getGEOSuppFiles("GSE226809", baseDir = raw_dir, makeDirectory = FALSE)
    if (!valid_archive(archive)) {
      stop("GEOquery download did not pass the 72-file/30-PMC archive checks.")
    }
  }
}

message("Extracting validated GEO archive: ", basename(archive))
utils::untar(archive, exdir = raw_dir)
sample_files <- list.files(raw_dir, pattern = "^GSM[0-9]+_Mm_PMC_.*_dge\\.tsv\\.gz$", full.names = TRUE)
if (length(sample_files) != 30L) {
  stop("Expected 30 PMC DGE files after extraction; found ", length(sample_files), ".")
}
message("Found all 30 PMC DGE files in ", raw_dir, ".")
