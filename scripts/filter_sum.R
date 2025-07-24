# Get command line arguments for input and output directories
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
# Get command line arguments for input and output directories
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
   print(args)
  stop("Usage: Rscript filter_sum.R <input_directory> <output_directory>")
}

# Set input and output directories from command line arguments
concat_dir <- args[1]
out_dir <- args[2]
p_val <- args[3]

# Create output directory if it doesn't exist
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# List all RDS files starting with 'cpeltdm' in the directory
rds_files <- list.files(concat_dir, pattern = "^cpeltdm.*\\.rds$", full.names = TRUE)

for (rds_file in rds_files) {
  # Load the GRanges object
  gr <- readRDS(rds_file)
  
  # Filter: keep only rows with p.adj >= 0.1 and not NA
  keep <- !is.na(gr$p.adj) & gr$p.adj <= p_val
  gr_filtered <- gr[keep]
  
  # Convert to data.frame for export
  df <- as.data.frame(gr_filtered)
  
  # Write to CSV (human-readable)
  out_file <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(rds_file)), "_filtered.csv"))
  write.csv(df, out_file, row.names = FALSE)
  cat("Filtered and exported:", out_file, "\n")
}

}

# Set input and output directories from command line arguments
concat_dir <- args[1]
out_dir <- args[2]
p_val <- args[3]

# Create output directory if it doesn't exist
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# List all RDS files starting with 'cpeltdm' in the directory
rds_files <- list.files(concat_dir, pattern = "^cpeltdm.*\\.rds$", full.names = TRUE)

for (rds_file in rds_files) {
  # Load the GRanges object
  gr <- readRDS(rds_file)
  
  # Filter: keep only rows with p.adj >= 0.1 and not NA
  keep <- !is.na(gr$p.adj) & gr$p.adj <= p_val
  gr_filtered <- gr[keep]
  
  # Convert to data.frame for export
  df <- as.data.frame(gr_filtered)
  
  # Write to CSV (human-readable)
  out_file <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(rds_file)), "_filtered.csv"))
  write.csv(df, out_file, row.names = FALSE)
  cat("Filtered and exported:", out_file, "\n")
}
