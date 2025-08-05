library(rtracklayer)
library(GenomicRanges)

# Import CpG island BED (or .txt in BED format)
cpg_islands <- import.bed("/dcs07/afeinber/data/personal/jzhan/CpG_islands_mm10.txt", trackLine = FALSE)

# Name all islands (optional)
cpg_islands$name <- "island"

# Create 2kb upstream (shore1) and downstream (shore2) flanks
shore1 <- flank(cpg_islands, 2000)
shore2 <- flank(cpg_islands, 2000, start = FALSE)

# Combine and reduce overlaps
shore1_2 <- GenomicRanges::reduce(c(shore1, shore2))

# Ensure consistent seqlevels
seqlevelsStyle(shore1_2) <- "UCSC"


export.bed(shore1_2, "CpG_shores_2kb_merged.bed")
