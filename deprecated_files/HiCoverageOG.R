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



###########################
library(bsseq)
library(GenomicRanges)
library(AnnotationDbi)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(rtracklayer)
library(devtools)
library(org.Mm.eg.db)
#######################################

# Load gene annotations for mouse
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes <- genes(txdb)
genes_gr_symbol <- sort(genes, by=~gene_id)

# Get gene symbols
out <- AnnotationDbi::select(org.Mm.eg.db, keys=genes_gr_symbol$gene_id, 
                            columns=c("SYMBOL"), keytype="ENTREZID")
idx <- match(genes_gr_symbol$gene_id, out$ENTREZID)
genes_gr_symbol$gene_symbol <- out$SYMBOL[idx]
genes_gr_symbol <- sort(genes_gr_symbol)

#######################################


subdivide_granges <- function(gr, width = 100) {
  gr_list <- lapply(seq_along(gr), function(i) {
    range <- gr[i]
    starts <- seq(start(range), end(range), by = width)
    ends <- pmin(starts + width - 1, end(range))
    GRanges(seqnames = seqnames(range),
            ranges = IRanges(start = starts, end = ends),
            strand = strand(range))
  })
  do.call(c, gr_list)
}

#######################################




# Load in bsseq object with methylation calls
bs_obj <- readRDS("/dcs07/afeinber/data/personal/jzhan/CPEL2/BS.sorted_by_condition.RDS")

cgs_gr_withCov <- granges(bs_obj)
cov <- getCoverage(bs_obj)
mcols(cgs_gr_withCov) <- cov

# Set CpG coverage filter here
# Using 10x coverage as in publication
cov_filter <- 1

#----

# If you have multiple genotypes in the bsseq object you should filter them separately
## So make a list with the sample names in each genotype

sampList <- list(

  "Old" = c("O-1", "O-2", "O-3", "O-4", "O-5"),

  "Young" = c("Y-1", "Y-2", "Y-3", "Y-4", "Y-5")

)
desired_samps <- do.call("c", sampList)
#----

# Only analyze autosomal CpGs
cgs_gr_autosomes <- cgs_gr_withCov[ seqnames(cgs_gr_withCov) %in% 1:19 ]
mcols(cgs_gr_autosomes) <- mcols(cgs_gr_autosomes)[ , colnames(mcols(cgs_gr_autosomes)) %in% desired_samps ]

# Check coverage for each CpG in each sample
cgs_pass_ind <- data.matrix(mcols( cgs_gr_autosomes )) >= cov_filter
#-----
# Required function
get_rle_indices <- function(rle_vec_curr){
  # Print parameters
  cat("Parameters:\n")
  cat("Coverage filter:", cov_filter, "\n")
  cat("Sample groups:\n")
  for(group in names(sampList)) {
    cat("  ", group, ":", paste(sampList[[group]], collapse=", "), "\n") 
  }
  cat("Minimum CpGs per region: 3\n")
  cat("Maximum gap between CpGs: 50bp\n")
  cat("Chromosomes analyzed: 1-19\n")
  
  # Check if RLE result is empty
  if (length(rle_vec_curr$values) == 0) {
    cat("Warning: RLE result is empty, returning empty data frame\n")
    return(data.frame(start = integer(0), end = integer(0), val = logical(0)))
  }
  
  # Turn rle indices into a convenient dataframe containing start/end/value
  # used to construct CPEL analysis regions 
  tmp_vals <- rle_vec_curr$values
  tmp_ends <- cumsum(rle_vec_curr$lengths)
  tmp_starts <- c(1, head(tmp_ends + 1, -1))
  inds <- data.frame("start" = tmp_starts, "end" = tmp_ends, "val" = tmp_vals)
  
  return(inds)
  
}
#----- 
# GENERATE REGIONS FOR REQUIRING CpGs TO PASS COVERAGE IN EVERY SAMPLE
# Get CpGs that pass coverage in every sample ( = every column of matrix)
cg_pass_ind_AllSamps <- rowSums(cgs_pass_ind)
cg_pass_ind_AllSamps <- cg_pass_ind_AllSamps == ncol(cgs_pass_ind)
cat("Number of CpGs in autosomes: ", length(cgs_gr_autosomes), "\n")
cat("Matrix dimensions: ", dim(cgs_pass_ind), "\n")

rle_vec_AllSamps <- rle(cg_pass_ind_AllSamps)
rle_inds_AllSamps <- get_rle_indices(rle_vec_AllSamps) 
# returns dataframe containing start, end, and true/false from rles

# Check if we have any regions that pass
if (nrow(rle_inds_AllSamps) == 0 || !any(rle_inds_AllSamps$val)) {
  cat("No regions found that pass coverage filter in all samples\n")
} else {
  # Get indices of CpGs in each candidate region
  rle_inds_pass_AllSamps <- rle_inds_AllSamps[ rle_inds_AllSamps$val ]
  rle_inds_pass_AllSamps <- apply(rle_inds_pass_AllSamps, 1, function(x){ seq(x[1], x[2]) })

  # And filter only to include candidate regions with at least 3 CpGs
  rle_inds_pass_AllSamps <- Filter(function(x) length(x) >= 3, rle_inds_pass_AllSamps)
  
  if (length(rle_inds_pass_AllSamps) > 0) {
    names(rle_inds_pass_AllSamps) <- paste0("region", 1:length(rle_inds_pass_AllSamps))


    # Now assign a region indicator to each of the CpGs in the CpG grange
    ind_vec <- rep("None", nrow(mcols(cgs_gr_autosomes)))
    all_indices <- unlist(rle_inds_pass_AllSamps)
    all_labels <- rep(names(rle_inds_pass_AllSamps), lengths(rle_inds_pass_AllSamps))
    ind_vec[ all_indices ] <- all_labels
    cgs_gr_autosomes$CG_region_status <- ind_vec

    # And filter to CpGs that have passed
    candidate_cgs <- cgs_gr_autosomes[ cgs_gr_autosomes$CG_region_status != "None" ]


    # And generate candidate regions
    # We must split before reducing, otherwise it may reduce across areas that contain CpGs that were filtered out
    # That is also why we assign a region name to each CG
    candidate_cgs <- split(candidate_cgs, f=candidate_cgs$CG_region_status)
    candidate_regions <-  GenomicRanges::reduce(candidate_cgs, min.gap=50)
    candidate_regions <- unlist(candidate_regions)
    candidate_regions$region <- names(candidate_regions)
    candidate_regions <- sort(candidate_regions)


    # Split large regions into multiple regions (see subdivideGRanges)
    # Gives better granularity for CPEL analysis
    candidate_regions_split <- subdivide_granges(candidate_regions, 100)

    # And finally filter by number of CpGs in region - another check to make sure regions have at least 3 CpGs
    cgs_in_roi <- countOverlaps(candidate_regions_split, cgs_gr_autosomes)
    candidate_regions_final <- candidate_regions_split[ cgs_in_roi >= 3 ]

    # Can also subset to candidate regions to promoters
    seqlevelsStyle(genes_gr_symbol) <- "ncbi"
    proms <- promoters(genes_gr_symbol, upstream=2000, downstream=2000)
    candidate_regions_final_promoters <- subsetByOverlaps(candidate_regions_final, proms)

    export_file <- paste0("CPEL_analysis_regions_", cov_filter, "xcoverage_in_all_samples.bed")
    export.bed(candidate_regions_final, export_file)

    export_file <- paste0("CPEL_analysis_regions_Promoters_", cov_filter, "xcoverage_in_all_samples.bed")
    export.bed(candidate_regions_final_promoters, export_file)
    
    cat("Successfully exported regions for all samples analysis\n")
  } else {
    cat("No regions with at least 3 CpGs found for all samples analysis\n")
  }
}

#------ 
# GENERATE REGIONS FOR REQUIRING CpGs TO PASS COVERAGE IN AT LEAST 4 SAMPLES PER GROUP

# Split matrix indicating coverage filter pass status
# into multiple matrices, one per AML subtype/normal
cg_pass_ind_perGroup <- lapply(sampList, function(x){ cgs_pass_ind[, colnames(cgs_pass_ind) %in% x] })

# Get CpGs that pass coverage in at least 4 samples in each group
cg_pass_ind_atLeast4 <- lapply(cg_pass_ind_perGroup, rowSums)
cg_pass_ind_atLeast4 <- lapply(cg_pass_ind_atLeast4, function(x){ x >= 4 })


rle_vecs_atLeast4 <- lapply(cg_pass_ind_atLeast4, rle)
rle_inds_atLeast4 <- lapply(rle_vecs_atLeast4, get_rle_indices)

# Check if we have any regions that pass
has_valid_regions <- sapply(rle_inds_atLeast4, function(x) nrow(x) > 0 && any(x$val))

if (all(has_valid_regions)) {
  # Get indices of CpGs in each candidate region
  rle_inds_pass_atLeast4 <- lapply(rle_inds_atLeast4, function(x){ x[x$val, ]})
  rle_inds_pass_atLeast4 <- lapply(rle_inds_pass_atLeast4, function(x){
    apply(x, 1, function(y){ seq(y[1], y[2]) })
  })
  # And filter only to include candidate regions with at least 3 CpGs
  rle_inds_pass_atLeast4 <- lapply(rle_inds_pass_atLeast4, function(x){
    tmp <- Filter(function(y) length(y) >= 3, x)
    names(tmp) <- paste0("region", 1:length(tmp))
    return(tmp)
  })


  # Now assign a region indicator to each of the CpGs in the CpG grange
  # With one instance per subtype
  needed_vecs <- lapply(rle_inds_pass_atLeast4, function(x){
    ind_vec <- rep("None", nrow(mcols(cgs_gr_autosomes)))
    
    all_indices <- unlist(x)
    all_labels <- rep(names(x), lengths(x))
    
    ind_vec[ all_indices ] <- all_labels
    return(ind_vec)
  })

  cgs_gr_byGroup <- lapply(sampList, function(x){
    tmp <- cgs_gr_autosomes
    mcols(tmp) <- mcols(tmp)[ , colnames(tmp) %in% x ]
    return(tmp)
  })

  ## and filter those to CpGs that have passed in the given group
  candidate_cgs_gr_byGroup <- mapply(function(gr, inds){
    gr$CG_region_status <- inds
    tmp <- gr[ gr$CG_region_status != "None" ]
    return(tmp)
  }, cgs_gr_byGroup, needed_vecs, SIMPLIFY=F)


  # Now get CpGs that have passed in all groups
  candidate_cgs_gr_allGroups <- Reduce(subsetByOverlaps, candidate_cgs_gr_byGroup)

  if (length(candidate_cgs_gr_allGroups) > 0) {
    candidate_regions_CovAtLeast4 <- split(candidate_cgs_gr_allGroups, f=candidate_cgs_gr_allGroups$CG_region_status)
    candidate_regions_CovAtLeast4 <-  GenomicRanges::reduce(candidate_regions_CovAtLeast4, min.gap=50)
    candidate_regions_CovAtLeast4 <- unlist(candidate_regions_CovAtLeast4)
    candidate_regions_CovAtLeast4 <- sort(candidate_regions_CovAtLeast4)
    candidate_regions_CovAtLeast4$region <- paste0("region", 1:length(candidate_regions_CovAtLeast4))


    # Split large regions into multiple regions (see subdivideGRanges)
    # Gives better granularity for CPEL analysis
    candidate_regions_CovAtLeast4_split <- subdivide_granges(candidate_regions_CovAtLeast4, 100)

    # And finally filter by number of CpGs in region - another check to make sure regions have at least 3 CpGs
    cgs_in_roi <- countOverlaps(candidate_regions_CovAtLeast4_split, cgs_gr_autosomes)
    candidate_regions_CovAtLeast4_final <- candidate_regions_CovAtLeast4_split[ cgs_in_roi >= 3 ]

    # Now subset to candidate regions to promoters
    seqlevelsStyle(genes_gr_symbol) <- "ncbi"
    proms <- promoters(genes_gr_symbol, upstream=2000, downstream=2000)
    candidate_regions_CovAtLeast4_final_promoters <- subsetByOverlaps(candidate_regions_CovAtLeast4_final, proms)


    export_file <- paste0("CPEL_analysis_regions_", cov_filter, "xcoverage_in_AtLeast4SamplesPerGroup.bed")
    export.bed(candidate_regions_CovAtLeast4_final, export_file)

    export_file <- paste0("CPEL_analysis_regions_Promoters_", cov_filter, "xcoverage_in_AtLeast4SamplesPerGroup.bed")
    export.bed(candidate_regions_CovAtLeast4_final_promoters, export_file)
    
    cat("Successfully exported regions for at least 4 samples per group analysis\n")
  } else {
    cat("No overlapping regions found between groups for at least 4 samples per group analysis\n")
  }
} else {
  cat("Not all groups have valid regions for at least 4 samples per group analysis\n")
  for (i in seq_along(has_valid_regions)) {
    if (!has_valid_regions[i]) {
      cat("Group", names(has_valid_regions)[i], "has no valid regions\n")
    }
  }
}