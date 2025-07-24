# Load libraries (make sure conda_R/4.4 is loaded)
library(bsseq)
library(GenomicRanges)
library(AnnotationDbi)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(rtracklayer)
library(devtools)
library(exomeCopy)
library(org.Mm.eg.db)

# Load BSseq object
bs_obj <- readRDS("BS.sorted_by_condition.RDS")

# Extract GRanges with coverage
cgs_gr_withCov <- granges(bs_obj)
cov <- getCoverage(bs_obj)
mcols(cgs_gr_withCov) <- cov

# Define sample groups
sampList <- list(
  "Old" = c("O-1", "O-2", "O-3", "O-4", "O-5"),
  "Young" = c("Y-1", "Y-2", "Y-3", "Y-4", "Y-5")
)
desired_samps <- unlist(sampList)

# Filter autosomal CpGs only
cgs_gr_autosomes <- cgs_gr_withCov[seqnames(cgs_gr_withCov) %in% paste0("chr", 1:19)]
mcols(cgs_gr_autosomes) <- mcols(cgs_gr_autosomes)[, colnames(mcols(cgs_gr_autosomes)) %in% desired_samps]

# Check coverage filter (10x)
cov_filter <- 10
cgs_pass_ind <- data.matrix(mcols(cgs_gr_autosomes)) >= cov_filter

# Helper function
get_rle_indices <- function(rle_vec_curr) {
  tmp_vals <- rle_vec_curr$values
  tmp_ends <- cumsum(rle_vec_curr$lengths)
  tmp_starts <- c(1, head(tmp_ends + 1, -1))
  inds <- data.frame("start" = tmp_starts, "end" = tmp_ends, "val" = tmp_vals)
  return(inds)
}

# --- Regions where ALL samples pass ---
cg_pass_ind_AllSamps <- rowSums(cgs_pass_ind) == ncol(cgs_pass_ind)
rle_vec <- rle(cg_pass_ind_AllSamps)
rle_inds <- get_rle_indices(rle_vec)
pass_inds <- rle_inds[rle_inds$val, ]
pass_inds <- apply(pass_inds, 1, function(x) seq(x[1], x[2]))
pass_inds <- Filter(function(x) length(x) >= 3, pass_inds)
names(pass_inds) <- paste0("region", seq_along(pass_inds))

ind_vec <- rep("None", length(cgs_gr_autosomes))
ind_vec[unlist(pass_inds)] <- rep(names(pass_inds), lengths(pass_inds))
cgs_gr_autosomes$CG_region_status <- ind_vec

candidate_cgs <- cgs_gr_autosomes[cgs_gr_autosomes$CG_region_status != "None"]
candidate_cgs <- split(candidate_cgs, f=candidate_cgs$CG_region_status)
candidate_regions <- reduce(candidate_cgs, min.gap=50)
candidate_regions <- unlist(candidate_regions)
candidate_regions$region <- names(candidate_regions)
candidate_regions <- sort(candidate_regions)
candidate_regions_split <- subdivideGRanges(candidate_regions, 100)
cgs_in_roi <- countOverlaps(candidate_regions_split, cgs_gr_autosomes)
candidate_regions_final <- candidate_regions_split[cgs_in_roi >= 3]

# Promoter annotation
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes <- genes(txdb)
chrs <- paste0("chr", 1:19)
genes <- keepSeqlevels(genes, c(chrs, "chrX", "chrY"), pruning.mode="coarse")
genes <- sortSeqlevels(genes)
out <- AnnotationDbi::select(org.Mm.eg.db, keys=as.character(gene_id(geneRanges=genes)), keytype="ENTREZID", columns="SYMBOL")
idx <- match(genes$gene_id, out$ENTREZID)
genes$symbol <- out$SYMBOL[idx]
seqlevelsStyle(genes) <- "UCSC"
proms <- promoters(genes, upstream=2000, downstream=2000)
candidate_regions_final_promoters <- subsetByOverlaps(candidate_regions_final, proms)

# Export
export.bed(candidate_regions_final, paste0("CPEL_analysis_regions_", cov_filter, "x_all_samples.bed"))
export.bed(candidate_regions_final_promoters, paste0("CPEL_analysis_regions_promoters_", cov_filter, "x_all_samples.bed"))

# --- Regions where >=4 samples per group pass ---
cg_pass_perGroup <- lapply(sampList, function(samps) cgs_pass_ind[, colnames(cgs_pass_ind) %in% samps])
cg_pass_groupPass <- lapply(cg_pass_perGroup, function(x) rowSums(x) >= 4)
rle_perGroup <- lapply(cg_pass_groupPass, rle)
rle_inds_group <- lapply(rle_perGroup, get_rle_indices)
pass_inds_group <- lapply(rle_inds_group, function(df) apply(df[df$val, ], 1, function(y) seq(y[1], y[2])))
pass_inds_group <- lapply(pass_inds_group, function(x) {
  y <- Filter(function(z) length(z) >= 3, x)
  names(y) <- paste0("region", seq_along(y))
  return(y)
})

needed_vecs <- lapply(pass_inds_group, function(idx_list) {
  v <- rep("None", length(cgs_gr_autosomes))
  v[unlist(idx_list)] <- rep(names(idx_list), lengths(idx_list))
  return(v)
})

cgs_gr_byGroup <- lapply(sampList, function(samps) {
  tmp <- cgs_gr_autosomes
  mcols(tmp) <- mcols(tmp)[, colnames(tmp) %in% samps]
  return(tmp)
})

candidate_cgs_byGroup <- mapply(function(gr, v) {
  gr$CG_region_status <- v
  gr[gr$CG_region_status != "None"]
}, cgs_gr_byGroup, needed_vecs, SIMPLIFY=FALSE)

cgs_allGroups <- Reduce(subsetByOverlaps, candidate_cgs_byGroup)
candidate_regions <- split(cgs_allGroups, f=cgs_allGroups$CG_region_status)
candidate_regions <- reduce(candidate_regions, min.gap=50)
candidate_regions <- unlist(candidate_regions)
candidate_regions$region <- paste0("region", seq_along(candidate_regions))
candidate_regions <- sort(candidate_regions)
candidate_regions_split <- subdivideGRanges(candidate_regions, 100)
cgs_in_roi <- countOverlaps(candidate_regions_split, cgs_gr_autosomes)
candidate_regions_final <- candidate_regions_split[cgs_in_roi >= 3]
candidate_regions_final_promoters <- subsetByOverlaps(candidate_regions_final, proms)

export.bed(candidate_regions_final, paste0("CPEL_analysis_regions_", cov_filter, "x_4samp_min.bed"))
export.bed(candidate_regions_final_promoters, paste0("CPEL_analysis_regions_promoters_", cov_filter, "x_4samp_min.bed"))
