# Get command line arguments for input and output directories
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript statistics.R <input_directory> <output_directory>")
}

# Create output directory if it doesn't exist
outdir <- args[2]
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# Set the main directory from command line argument
main_dir <- args[1]

# List all subdirectories
subdirs <- list.dirs(main_dir, recursive = FALSE)
# pull just the folder names
dirnames      <- basename(subdirs)               # → "chr_1_out", "chr_10_out", …

# build the exact ordering you want
desired_levels <- c(
  paste0("chr_", 1:19, "_out"),  # chr_1_out, chr_2_out, …, chr_19_out
  "chr_X_out", "chr_Y_out", "chr_M_out"
)

# make a factor with those levels (absent ones are ignored)
dir_f <- factor(dirnames, levels = desired_levels)

# reorder the full paths by that factor
subdirs <- subdirs[order(dir_f)]
# Print sorted order for verification
cat("Sorted chromosome order:\n")
print(basename(subdirs))


for (subdir in subdirs) {
  # List all .rds files in the subdirectory
  rds_files <- list.files(subdir, pattern = "\\.rds$", full.names = TRUE)
  
  for (rds_file in rds_files) {
    # Get the variable name: filename without directory or extension
    var_name <- tools::file_path_sans_ext(basename(rds_file))
    cat("Accessing:", rds_file, "from", subdir, "\n")
    # Read the RDS file
    new_data <- readRDS(rds_file)
    
if (exists(var_name, envir = .GlobalEnv)) {
  old_data <- get(var_name, envir = .GlobalEnv)

  if (inherits(old_data, "GRanges")) {
    combined_data <- c(old_data, new_data)
  } else {
    combined_data <- rbind(old_data, new_data)
  }

  assign(var_name, combined_data, envir = .GlobalEnv)
  cat("Concatenated:", var_name, "with", rds_file, "\n")
} else {
  assign(var_name, new_data, envir = .GlobalEnv)
  cat("Loaded:", var_name, "from", rds_file, "\n")
}

    
  }
}


# Get all objects in global environment, excluding specified variables
excluded_vars <- c("args", "combined_data", "main_dir", "new_data", "old_data", 
                  "outdir", "rds_file", "rds_files", "subdir", "subdirs", "var_name", "excluded_vars", "dirnames", "desired_levels", "dir_factor")
all_objects <- ls(envir = .GlobalEnv)[!ls(envir = .GlobalEnv) %in% excluded_vars]

# Save each concatenated data frame to RDS file
for (obj_name in all_objects) {
  # Get the object
  obj <- get(obj_name, envir = .GlobalEnv)
  
  # Create output filename
  out_file <- file.path(outdir, paste0(obj_name, ".rds"))
  
  # Save to RDS
  saveRDS(obj, file = out_file)
  cat("Saved:", obj_name, "to", out_file, "\n")
}
