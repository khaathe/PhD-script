rank.table <- function(deseq.result, alpha = 0.05){
  rank <- deseq.result
  rank <- rank[!is.na(rank$padj),] #Delete NAs
  rank <- rank[rank$padj<alpha,] #Filter significant result
  rank <- rank[order(rank$log2FoldChange, decreasing = TRUE), ] #Decreasing sort on log2FoldChange
  rank <- rank["log2FoldChange"] #Selecting log2FoldChange column
  return(rank)
}