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
library(Mus.musculus)

# Get the path of the currently executing script
script_path <- commandArgs(trailingOnly = FALSE)
script_path <- script_path[grep("--file=", script_path)]
if (length(script_path) > 0) {
  script_dir <- dirname(sub("--file=", "", script_path))
} else {
  script_dir <- getwd()  # Fallback
}

# Source the filtering function (relative to this script's directory)
source(file.path(script_dir, "filtering.R"))

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript partitioning.R <tile_width> <landmark> <bs_obj>")
}

tile_width <- as.numeric(args[1])
landmark <- args[2]
bs_obj <- readRDS(args[3])

# Validate landmark parameter
valid_landmarks <- c("promoters", "shores", "islands")
if (!landmark %in% valid_landmarks) {
  stop("Invalid landmark parameter. Must be one of: ", paste(valid_landmarks, collapse = ", "))
}


# TODO: 
# Define sample groups by extracting from bs_obj
sampList <- list(
  Old   = c("O-1","O-2","O-3","O-4","O-5"),
  Young = c("Y-1","Y-2","Y-3","Y-4","Y-5")
)

# Extract high-coverage CpGs for each group
cgs_gr_list <- lapply(sampList, function(samps) {
  tmp <- bs_obj[, sampleNames(bs_obj) %in% samps]
  gr  <- granges(tmp)
  cov <- getCoverage(tmp)
  mcols(gr) <- cov
  gr[rowSums(as.matrix(mcols(gr)) >= 10) >= 5, ]
})

# Get CpGs that pass coverage requirements in all groups
cgs_in_all_groups <- Reduce(intersect, cgs_gr_list)
seqlevelsStyle(cgs_in_all_groups) <- "UCSC"

# Extract the appropriate landmark regions
# Get project root (two levels up from this script)
project_root <- dirname(dirname(script_dir))
cpg_islands_file <- file.path(project_root, "CpG_islands_mm10.txt")

landmark_data <- filtering(landmark, cpg_islands_file)
regions <- landmark_data$regions
genes <- landmark_data$genes

# Tile the regions
regions_tiled <- unlist(tile(regions, width = tile_width))

# Count overlaps of CpGs in each tile
cgs_in_tiles <- countOverlaps(regions_tiled, cgs_in_all_groups)

# Keep tiles with at least 3 well-covered CpGs
regions_tiled_3cgs <- regions_tiled[cgs_in_tiles >= 3]

# For each tile, find the nearest gene and get its symbol
nearest_gene_idx <- nearest(regions_tiled_3cgs, genes, select = "arbitrary") 
nearest_gene_symbol <- genes$gene_symbol[nearest_gene_idx]

# Calculate distance from each tile to nearest gene
dist_to_gene <- distanceToNearest(regions_tiled_3cgs, genes)

# Keep only tiles that are within 1kb of a gene
keep_idx <- queryHits(dist_to_gene)[subjectHits(dist_to_gene) & mcols(dist_to_gene)$distance <= 1000]
regions_tiled_3cgs <- regions_tiled_3cgs[keep_idx]


# Print summary statistics
cat("Total kb covered:", sum(width(regions_tiled_3cgs)/1000), "kb\n")
cat("Number of regions:", length(regions_tiled_3cgs), "\n")

# Generate output filename based on landmark type
if (landmark == "promoters") {
  output_filename <- paste0("CPEL_regions_3CGs_Passing10xCovatLeast5Samps_", tile_width, "genefilteredbpTiledProms_nearGenes.bed")
} else {
  output_filename <- paste0("CPEL_regions_3CGs_Passing10xCovatLeast5Samps_", tile_width, "bpTiled", toupper(substring(landmark, 1, 1)), substring(landmark, 2), ".bed")
}

# Export to BED file
export.bed(regions_tiled_3cgs, output_filename)

cat("Results exported to:", output_filename, "\n") 