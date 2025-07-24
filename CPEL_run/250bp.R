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


#################################################################

bs_obj <- readRDS("/dcs07/afeinber/data/personal/jzhan/CPEL2/BS.sorted_by_condition.RDS")

# If you have multiple genotypes in the bsseq object you should filter them separately
## So make a list with the sample names in each genotype
sampList = list("Old" = c("O-1", "O-2", "O-3", "O-4", "O-5"), 
                "Young" = c("Y-1", "Y-2", "Y-3", "Y-4", "Y-5")) # etc

# Get one CpG grange for each genotype 
## where the mcols are the coverage values of each sample
cgs_gr_list_withCov <- lapply(sampList, 
  function(samps) {
    tmp <- bs_obj[, sampleNames(bs_obj) %in% samps ]
    grObj <- granges(tmp)
    covObj <- getCoverage( tmp)
    mcols(grObj) <- covObj
  return(grObj)
})


# Find CpGs with at least 10x coverage in at least 5 samples (in each group = element of list)
## You can play around with the coverage filter and sample number requirement
## I don't recommend dropping coverage below 5x though
cgs_w_10xcov_in_5_samps <- lapply( cgs_gr_list_withCov, function(x){
  x[ rowSums( as.matrix( mcols(x) ) >= 10 ) >= 5, ] 
})

# Get CpGs that pass this requirement in all groups
cgs_in_all_groups <- Reduce(intersect, cgs_w_10xcov_in_5_samps)
seqlevelsStyle(cgs_in_all_groups) <- "UCSC"
# Generate tiled promoters
txdb <-TxDb.Mmusculus.UCSC.mm10.knownGene # You may change promoter annotation as desired. This is the one I use.
genesAll<-genes(txdb)
chrsOfInterest=paste("chr",1:19,sep="")

bpBeforeTSS=2000
bpPastTSS=2000

genes<-keepSeqlevels(genesAll,c(chrsOfInterest, "chrM"),pruning.mode="coarse")
genes <- sortSeqlevels(genes)

out <- AnnotationDbi::select(Mus.musculus, key=as.character(genes$gene_id), keytype="ENTREZID", columns=c("SYMBOL"))
genes_gr_symbol <-sort(genes, by=~gene_id)
idx <- match(genes_gr_symbol$gene_id, out$ENTREZID)
genes_gr_symbol$gene_symbol <- out$SYMBOL[  idx ]
genes_gr_symbol <- sort(genes_gr_symbol)

seqlevelsStyle(genes_gr_symbol) <- "UCSC"
proms <- promoters(genes_gr_symbol, upstream=2000, downstream=2000)
proms_tiled <- tile(proms, width=250) # you can adjust the width of tiles too
proms_tiled <- unlist( proms_tiled) 


# Now check each tile to see if it passes the requirements:
# Contains least 3 CpGs 
# With CpGs passing the coverage filter in all groups (coverage filter already applied)

## first get number of cgs passing the coverage requirements in each tile
cgs_in_tiles <- countOverlaps(proms_tiled, cgs_in_all_groups) 
## then filter to have at least 3 cgs in each (play around - try 3, 4, 5?)
proms_tiled_3cgs_10x <- proms_tiled[cgs_in_tiles>=3]

# You should check how many promoters are covered, how many kb are covered, 
# and how many regions you have
# Then you can adjust coverage/CpG number filters to be more or less stringent

# How many promoters are covered
length(subsetByOverlaps(proms, proms_tiled_3cgs_10x))
# How many kb
sum(width(proms_tiled_3cgs_10x)/1000)
# How many regions 
length(proms_tiled_3cgs_10x)

seqlevels(proms_tiled)
seqlevels(cgs_in_all_groups)


export.bed(proms_tiled_3cgs_10x, "CPEL_regions_3CGs_Passing10xCovatLeast5Samps_250bpTiledProms.bed")