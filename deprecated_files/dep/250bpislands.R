library(bsseq)
library(GenomicRanges)
library(rtracklayer)

# 1. Load your BSseq object
bs_obj <- readRDS("/dcs07/afeinber/data/personal/jzhan/CPEL2/BS.sorted_by_condition.RDS")

# 2. Split by sample group and extract high-coverage CpGs
sampList <- list(
  Old   = c("O-1","O-2","O-3","O-4","O-5"),
  Young = c("Y-1","Y-2","Y-3","Y-4","Y-5")
)
cgs_gr_list <- lapply(sampList, function(samps) {
  tmp <- bs_obj[, sampleNames(bs_obj) %in% samps]
  gr  <- granges(tmp)
  cov <- getCoverage(tmp)
  mcols(gr) <- cov
  gr[rowSums(as.matrix(mcols(gr)) >= 10) >= 5, ]
})
cgs_in_all_groups <- Reduce(intersect, cgs_gr_list)
seqlevelsStyle(cgs_in_all_groups) <- "UCSC"




# Load CpG islands
cpg_islands <- import.bed("/dcs07/afeinber/data/personal/jzhan/CpG_islands_mm10.txt", trackLine = FALSE)
cpg_islands$name <- "island"
   
# Get 2kb shores upstream and downstream
# shore1 <- flank(cpg_islands, 2000)
# shore2 <- flank(cpg_islands, 2000, start = FALSE)

# Combine overlapping shore regions
cpg_shores <- cpg_islands

# Optional: restrict to UCSC-style seqlevels and sort
seqlevelsStyle(cpg_shores) <- "UCSC"
cpg_shores <- sortSeqlevels(cpg_shores)


args <- commandArgs(trailingOnly = TRUE)
tile_width <- as.numeric(args[1])
cpg_shores_tiled <- unlist(tile(cpg_shores, width = tile_width))


# # Count overlaps of CpGs in each tile
# cgs_in_tiles <- countOverlaps(cpg_shores_tiled, cgs_in_all_groups)

# # Keep tiles with at least 3 well-covered CpGs
# cpg_shores_tiled_3cgs <- cpg_shores_tiled[cgs_in_tiles >= 3]

# # You should check how many promoters are covered, how many kb are covered, 
# # and how many regions you have
# # Then you can adjust coverage/CpG number filters to be more or less stringent

# # Print total kb covered
# cat("Total kb covered:", sum(width(cpg_shores_tiled_3cgs)/1000), "kb\n")
# # Print number of regions
# cat("Number of regions:", length(cpg_shores_tiled_3cgs), "\n") 


export.bed(cpg_shores_tiled, "CPEL_regions_3CGs_Passing10xCovatLeast5Samps_250bpTiledCpGShores.bed")


# library(GenomicRanges)
# library(rtracklayer)
# library(bsseq)

# # 1. Import CpG islands and build ±2kb shores
# cpg_islands <- import.bed(
#   "/Users/oscarcamacho/.../CpG_islands_mm10.txt", 
#   trackLine = FALSE
# )
# seqlevelsStyle(cpg_islands) <- "UCSC"

# shore_up   <- flank(cpg_islands, 2000)                 # upstream 2kb
# shore_down <- flank(cpg_islands, 2000, start=FALSE)   # downstream 2kb
# shores     <- reduce(c(shore_up, shore_down))         # merge overlaps

# # 2. Tile each shore into 250 bp windows
# shore_tiles <- unlist(tile(shores, width=250))

# # 3. Standardize to your BSseq object’s style, and load your filtered CpGs
# bs_obj <- readRDS("…/BS.sorted_by_condition.RDS")
# # (re-use your previous code to produce `cgs_in_all_groups`)

# # 4. Count overlaps and filter to windows with ≥3 high-coverage CpGs
# hits_per_tile <- countOverlaps(shore_tiles, cgs_in_all_groups)
# good_tiles    <- shore_tiles[hits_per_tile >= 3]

# # 5. (Optional) collapse adjacent good tiles back into larger blocks
# good_regions  <- reduce(good_tiles)

# # 6. Export
# export.bed(good_tiles,   "shore_tiles_250bp.bed")
# export.bed(good_regions, "shore_regions_merged.bed")
