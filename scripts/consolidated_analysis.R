#!/usr/bin/env Rscript

# Consolidated analysis script that combines Analyze.R, merge.R, and filter_sum.R
# to eliminate unnecessary RDS file I/O operations

# Load required libraries
suppressPackageStartupMessages({
  library(Gmisc)
  library(rtracklayer)
  library(data.table)
})

# Get script directory and source analysis functions
script_dir <- dirname(sys.frame(1)$ofile)
if (is.null(script_dir) || script_dir == "") {
  script_dir <- getwd()
}
source(file.path(script_dir, "analysis.R"))

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript consolidated_analysis.R <input_directory> <p_value> [output_base_dir]")
}

in_dir <- args[1]
p_val <- as.numeric(args[2])
output_base_dir <- if (length(args) >= 3) args[3] else in_dir

# Get current date for output directories
DATE <- format(Sys.time(), "%y_%m_%d_%H_%M_%S")

# Create output directories
rds_outdir <- file.path(output_base_dir, paste0("RDS_files", DATE))
concat_outdir <- file.path(output_base_dir, paste0("concatenated_output", DATE))
filter_outdir <- file.path(concat_outdir, paste0("filtered_output", DATE))

dir.create(rds_outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(concat_outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(filter_outdir, recursive = TRUE, showWarnings = FALSE)

cat("Processing input directory:", in_dir, "\n")
cat("P-value threshold:", p_val, "\n")

# =============================================================================
# STEP 1: ANALYZE (equivalent to Analyze.R) - Process bedGraph files
# =============================================================================

# Initialize storage for concatenated data (equivalent to what merge.R does)
concatenated_data <- list()

# Process each subdirectory
subdirs <- list.dirs(in_dir, recursive = FALSE)
for (subdir in subdirs) {
  dir_name <- basename(subdir)
  sub_outdir <- file.path(rds_outdir, paste0(dir_name, "_out"))
  dir.create(sub_outdir, recursive = TRUE, showWarnings = FALSE)
  
  # Process each bedGraph file in the subdirectory
  bedgraph_files <- list.files(subdir, pattern = "\\.bedGraph$", full.names = TRUE)
  
  for (bedgraph_file in bedgraph_files) {
    if (file.exists(bedgraph_file)) {
      cat("Processing", bedgraph_file, "-> in-memory processing\n")
      
      # Process the bedGraph file (equivalent to Analyze.R logic)
      if (grepl("cpeltdm_", bedgraph_file)) {
        cat("Running read.differential.stat...\n")
        data <- read.differential.stat(bedgraph_file)
        var_name <- paste0(tools::file_path_sans_ext(basename(bedgraph_file)), "_diff_stat")
        
        # Save individual RDS file (to maintain same file output as original)
        rds_filename <- file.path(sub_outdir, paste0(var_name, ".rds"))
        saveRDS(data, rds_filename)
        cat("Saved individual RDS:", rds_filename, "\n")
        
      } else {
        cat("Running read.sample.stat...\n")
        data <- read.sample.stat(bedgraph_file)
        var_name <- paste0(tools::file_path_sans_ext(basename(bedgraph_file)), "_sample_stat")
        
        # Save individual RDS file (to maintain same file output as original)
        rds_filename <- file.path(sub_outdir, paste0(var_name, ".rds"))
        saveRDS(data, rds_filename)
        cat("Saved individual RDS:", rds_filename, "\n")
      }
      
      # =============================================================================
      # STEP 2: MERGE (equivalent to merge.R) - Concatenate in memory
      # =============================================================================
      
      # Instead of writing and reading RDS files, concatenate in memory
      if (var_name %in% names(concatenated_data)) {
        if (inherits(concatenated_data[[var_name]], "GRanges")) {
          concatenated_data[[var_name]] <- c(concatenated_data[[var_name]], data)
        } else {
          concatenated_data[[var_name]] <- rbind(concatenated_data[[var_name]], data)
        }
        cat("Concatenated:", var_name, "with", bedgraph_file, "\n")
      } else {
        concatenated_data[[var_name]] <- data
        cat("Loaded:", var_name, "from", bedgraph_file, "\n")
      }
    }
  }
}

# =============================================================================
# STEP 3: Save concatenated data (equivalent to merge.R output)
# =============================================================================

cat("Saving concatenated data to RDS files...\n")
for (obj_name in names(concatenated_data)) {
  obj <- concatenated_data[[obj_name]]
  out_file <- file.path(concat_outdir, paste0(obj_name, ".rds"))
  write.csv(obj, file = sub("\\.rds$", ".bedGraph", out_file), row.names = FALSE)
  cat("Saved concatenated:", obj_name, "to", out_file, "\n")
}

# =============================================================================
# STEP 4: FILTER (equivalent to filter_sum.R) - Filter and export to CSV
# =============================================================================

cat("Filtering data with p-value threshold:", p_val, "\n")

# Process only differential stat data (cpeltdm files)
for (obj_name in names(concatenated_data)) {
  if (grepl("cpeltdm", obj_name)) {
    cat("Processing differential data:", obj_name, "\n")
    
    gr <- concatenated_data[[obj_name]]
    
    # Filter: keep only rows with p.adj <= p_val and not NA
    keep <- !is.na(gr$p.adj) & gr$p.adj <= p_val
    gr_filtered <- gr[keep]
    
    cat("Filtered", sum(keep), "regions out of", length(gr), "total regions\n")
    
    # Convert to data.frame for export
    df <- as.data.frame(gr_filtered)
    
    # Write to CSV (human-readable)
    out_file <- file.path(filter_outdir, paste0(obj_name, "_filtered.csv"))
    write.csv(df, out_file, row.names = FALSE)
    cat("Filtered and exported:", out_file, "\n")
  }
}

cat("Analysis complete!\n")
cat("Individual RDS files saved to:", rds_outdir, "\n")
cat("Concatenated RDS files saved to:", concat_outdir, "\n")
cat("Filtered CSV files saved to:", filter_outdir, "\n")