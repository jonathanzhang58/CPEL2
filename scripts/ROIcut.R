#!/usr/bin/env Rscript








# Ensure BiocManager is installed (for Bioconductor packages)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Mus.musculus")

# List of Bioconductor packages you need
bioc_pkgs <- c(
  "bsseq",
  "GenomicRanges",
  "AnnotationDbi",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "rtracklayer",
  "Mus.musculus"
)

# Install any missing Bioconductor packages
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE)
}

# CRAN packages you need
cran_pkgs <- c("devtools")

# Install any missing CRAN packages
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
}

# Now load everything
library(bsseq)
library(GenomicRanges)
library(AnnotationDbi)
library(rtracklayer)
library(devtools)
# Try to load the package
library(Mus.musculus)

########################################################################################

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript ROIcut.R <directory>")
}

# Get the directory argument
input_dir <- args[1]

# Check if directory exists
if (!dir.exists(input_dir)) {
  stop("Directory does not exist: ", input_dir)
}

# Find .bed file in the directory
bed_files <- list.files(input_dir, pattern = "\\.bed$", full.names = TRUE)

if (length(bed_files) == 0) {
  stop("No .bed files found in directory: ", input_dir)
}

if (length(bed_files) > 1) {
  warning("Multiple .bed files found, using the first one: ", bed_files[1])
}



my_roi_file <- bed_files[1]
cat("Using ROI file:", my_roi_file, "\n")

# Import the ROI file
my_roi <- rtracklayer::import.bed(my_roi_file)
seqlevelsStyle(my_roi) <- "ncbi"

######################################################################

# Create chr subdirectory
chr_dir <- file.path(input_dir, "chr")
if (!dir.exists(chr_dir)) {
  dir.create(chr_dir, recursive = TRUE)
  cat("Created chr directory:", chr_dir, "\n")
}

base_dir <- chr_dir

for (i in 1:19){
  chr_x_rois <- my_roi[seqnames(my_roi) == i,]
  file_name1 = paste0(base_dir, "/myROIS_chr", i, ".bed")
  rtracklayer::export.bed(chr_x_rois, file_name1)
  cat("Exported chromosome", i, "ROIs to:", file_name1, "\n")
}

cat("ROI processing completed successfully!\n")


# Append "chr" prefix to chromosome numbers in all bed files
for (i in 1:19) {
  bed_file <- file.path(base_dir, paste0("myROIS_chr", i, ".bed"))
  
  # Read the bed file
  bed_data <- readLines(bed_file)
  
  # Add "chr" prefix to first column
  bed_data <- gsub(paste0("^", i), paste0("chr", i), bed_data)
  
  # Write back to file
  writeLines(bed_data, bed_file)
}
# Process each bed file to keep only first 4 columns
for (i in 1:19) {
  bed_file <- file.path(base_dir, paste0("myROIS_chr", i, ".bed"))
  
  # Read the bed file
  bed_data <- read.table(bed_file, sep="\t", stringsAsFactors=FALSE)
  
  # Keep only first 4 columns
  bed_data <- bed_data[, 1:4]
  
  # Write back to file with tab separation
  write.table(bed_data, bed_file, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}


cat("ROI transformation completed successfully!\n")