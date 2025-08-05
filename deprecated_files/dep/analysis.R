library(Gmisc)
library(rtracklayer)
library(data.table)

read.differential.stat <- function(file_in, comparison=NA){
  
  #cat('processing:',file_in,'\n')
  CPEL_in=import.bedGraph(file_in)
  if(length(CPEL_in)>0){
    
    colnames(elementMetadata(CPEL_in))= c('score','p.val','p.adj')
    if(all(seqlevels(CPEL_in)==gsub('chr','',seqlevels(CPEL_in)))){
      seqlevels(CPEL_in)=paste('chr',seqlevels(CPEL_in),sep='')
    }
    seqlevelsStyle(CPEL_in) <- "ucsc"
    # 
    # #fit  bedGraph reads, import.bedGraph will remove 1 from start
    start(CPEL_in)=start(CPEL_in)-1
    
    if (!is.na(comparison)){
      CPEL_in$Comparison = comparison
    }
    
    return(CPEL_in)
  }
}

read.sample.stat <-function(fn, comparison=NA, replicate="all"){
  fn_sub=gsub('_merged|.bedGraph|_sorted','', basename(fn))
  
  
  # print(strsplit(fn_sub, "_"))
  sample = strsplit(fn_sub, "_")[[1]][1]
  stat_type = strsplit(fn_sub, "_")[[1]][2]
  
  file_in=fn
  
  #cat('processing:',file_in,'\n')
  CPEL_in=import.bedGraph(file_in)
  if(length(CPEL_in)>0){
    
    if(all(seqlevels(CPEL_in)==gsub('chr','',seqlevels(CPEL_in)))){seqlevels(CPEL_in)=paste('chr',seqlevels(CPEL_in),sep='')}
    
    seqlevelsStyle(CPEL_in) <- "ucsc"
    # #fit  bedGraph reads, import.bedGraph will remove 1 from start
    start(CPEL_in)=start(CPEL_in)-1
    
    CPEL_in$Sample = sample
    CPEL_in$stat_type=stat_type
    
    if (!is.na(comparison)){
      CPEL_in$Comparison = comparison
    }
    
    return(CPEL_in)
  }
}