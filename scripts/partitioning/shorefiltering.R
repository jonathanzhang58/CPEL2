cpg_islands<-import.bed("/Users/oscarcamacho/Documents/Feinberg_lab/skin_aging/CpG_islands/CpG_islands_mm10.txt",trackLine = FALSE)
cpg_islands$name="island"
cpg_islands
###############################################################
#             Extract CpG island shores
###############################################################
# extract the shore defined by 2000 bp upstream of cpg islands
shore1=flank(cpg_islands, 2000)
# extract the shore defined by 2000 bp downstream of cpg islands
shore2=flank(cpg_islands,2000,FALSE)
# perform intersection and combine the shores where they overlap
shore1_2=GenomicRanges::reduce(c(shore1,shore2))