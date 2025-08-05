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

for (subdir in subdirs) {
  # List all .rds files in the subdirectory
  rds_files <- list.files(subdir, pattern = "\\.rds$", full.names = TRUE)
  
  for (rds_file in rds_files) {
    # Get the variable name: filename without directory or extension
    var_name <- tools::file_path_sans_ext(basename(rds_file))
    cat("Accessing:", rds_file, "from", subdir, "\n")
    # Read the RDS file
    new_data <- readRDS(rds_file)
    
    # If the variable already exists, concatenate (row-bind) the new data
    if (exists(var_name, envir = .GlobalEnv)) {
      old_data <- get(var_name, envir = .GlobalEnv)
      # Try to row-bind, if possible
      combined_data <- tryCatch({
        rbind(old_data, new_data)
      }, error = function(e) {
        warning(paste("Could not rbind", var_name, ":", e$message))
        old_data
      })
      assign(var_name, combined_data, envir = .GlobalEnv)
    } else {
      assign(var_name, new_data, envir = .GlobalEnv)
    }
    cat("Loaded:", var_name, "from", rds_file, "\n")
  }
}


# Get all objects in global environment, excluding specified variables
excluded_vars <- c("args", "combined_data", "main_dir", "new_data", "old_data", 
                  "outdir", "rds_file", "rds_files", "subdir", "subdirs", "var_name", "excluded_vars")
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
