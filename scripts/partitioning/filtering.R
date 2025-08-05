# Function to extract genomic landmarks based on the landmark parameter
filtering <- function(landmark, cpg_islands_file = "data/CpG_islands_mm10.txt") {
  # Load required libraries
  library(GenomicRanges)
  library(AnnotationDbi)
  library(rtracklayer)
  library(Mus.musculus)
  
  # Load CpG islands (needed for all landmark types)
  if (!file.exists(cpg_islands_file)) {
    stop("CpG islands file not found: ", cpg_islands_file)
  }
  cpg_islands <- import.bed(cpg_islands_file, trackLine = FALSE)
  cpg_islands$name <- "island"
  seqlevelsStyle(cpg_islands) <- "UCSC"
  
  if (landmark == "promoters") {
    # Generate promoters
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    genesAll <- genes(txdb)
    chrsOfInterest <- c(paste("chr", 1:19, sep=""), "chrX", "chrY")
    
    genes <- keepSeqlevels(genesAll, chrsOfInterest, pruning.mode="coarse")
    genes <- sortSeqlevels(genes)
    
    out <- AnnotationDbi::select(Mus.musculus, key=as.character(genes$gene_id), keytype="ENTREZID", columns=c("SYMBOL"))
    genes_gr_symbol <- sort(genes, by=~gene_id)
    idx <- match(genes_gr_symbol$gene_id, out$ENTREZID)
    genes_gr_symbol$gene_symbol <- out$SYMBOL[idx]
    genes_gr_symbol <- sort(genes_gr_symbol)
    
    seqlevelsStyle(genes_gr_symbol) <- "UCSC"
    
    # Extract CpG island shores for promoters
    shore1 <- flank(cpg_islands, 2000)
    shore2 <- flank(cpg_islands, 2000, start = FALSE)
    shore1_2 <- GenomicRanges::reduce(c(shore1, shore2))
    
    return(list(regions = shore1_2, genes = genes_gr_symbol))
    
  } else if (landmark == "shores") {
    # Extract CpG island shores
    shore1 <- flank(cpg_islands, 2000)
    shore2 <- flank(cpg_islands, 2000, start = FALSE)
    cpg_shores <- GenomicRanges::reduce(c(shore1, shore2))
    
    # Optional: restrict to UCSC-style seqlevels and sort
    seqlevelsStyle(cpg_shores) <- "UCSC"
    cpg_shores <- sortSeqlevels(cpg_shores)
    
    return(list(regions = cpg_shores, genes = NULL))
    
  } else if (landmark == "islands") {
    # Use CpG islands directly
    seqlevelsStyle(cpg_islands) <- "UCSC"
    cpg_islands <- sortSeqlevels(cpg_islands)
    
    return(list(regions = cpg_islands, genes = NULL))
    
  } else {
    stop("Invalid landmark parameter. Must be one of: 'promoters', 'shores', 'islands'")
  }
} 