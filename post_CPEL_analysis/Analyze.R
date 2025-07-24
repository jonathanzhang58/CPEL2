#!/usr/bin/env Rscript

# Load required libraries and source analysis.R
suppressPackageStartupMessages({
  library(Gmisc)
  library(rtracklayer)
  library(data.table)
})

source("analysis.R")

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript Analyze.R <input_file> <output_dir>")
}
input_file <- args[1]
outdir <- args[2]

# Ensure output directory exists
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}
# Check if input file contains "cpeltdm_" and read accordingly
if (grepl("cpeltdm_", input_file)) {
   cat("Running read.differential.stat...\n")
  data <- read.differential.stat(input_file)
  cat("Saving results...\n")
  diff_stat_file <- file.path(outdir, "differential_stat.rds")
  saveRDS(data, file.path(outdir, paste0(tools::file_path_sans_ext(basename(input_file)), "_diff_stat.rds")))
  cat(paste0("Differential stat saved to: ", diff_stat_file, "\n"))
} else {
   cat("Running read.sample.stat...\n")
  data <- read.sample.stat(input_file)
  cat("Saving results...\n")
  sample_stat_file <- file.path(outdir, "sample_stat.rds")
  saveRDS(data, file.path(outdir, paste0(tools::file_path_sans_ext(basename(input_file)), "_sample_stat.rds")))
  cat(paste0("Sample stat saved to: ", file.path(outdir, paste0(tools::file_path_sans_ext(basename(input_file)), "_sample_stat.rds")), "\n"))
}
