setwd("/home/spinicck/PhD/")

######### Remove every variable in memory
rm(list = ls())

######## Load Library needed
library(fgsea)
library(dplyr)

#IMPORTANT!! We must choose a seed to ensure the analysis is reproducible. This is because GSEA uses a 
#permutation test to span the null distribution of the statistic. Hence, each time we run it, results are 
#expected to change slightly.
set.seed(1)

reactome.gmt <- gmtPathways("Data/gene-set/reactome_pathways_symbol.gmt")
bioplanet.gmt <- gmtPathways("Data/gene-set/bioplanet_pathways_symbol.gmt")

# run gsea analysis with reactome and bioplanet GMT files for a list of DE analysis result
run.gsea <- function(f, reactome.gmt, bioplanet.gmt, create.ranked.list){
  message("###### GSEA Analysis for : ", f)
  ranked.list <- create.ranked.list(f)
  res <- list()
  message("\tRunning Reactome Enrichment ...")
  res[["reactome"]] <- fgsea(reactome.gmt, ranked.list, minSize=15, maxSize=500, nperm = 1000)
  message("\tDone !")
  message("\tRunning Bioplanet Enrichment ...")
  res[["bioplanet"]] <- fgsea(bioplanet.gmt, ranked.list, minSize=15, maxSize=500, nperm = 1000)
  message("\tDone !")
  return(res)
}

# Save a list of Generic Enrichment Map to tab delimited txt files
# Files name conttain the tool used for differential expression analysis, which patient
# has been compared to the control (analysis), the source of the gmt file (db.source)
save.gem.list.to.txt <- function(gem.list, save.dir, de.analysis){
  for (analysis in names(gem.list) ) {
    for (db.source in names(gem.list[[analysis]]) ){
      f <- paste0(save.dir, "gsea_", de.analysis, "_",analysis, "_", db.source, ".csv")
      gem.dt <- as.data.frame(gem.list[[analysis]][[db.source]]) %>% mutate(leadingEdge=sapply(leadingEdge, function(e){ paste0(e, collapse = ", ")}))
      write.table(gem.dt, file = f, row.names = F, col.names = T, sep = "\t")
    }
  }
}

# Remove path and file extension from the file name and return only the analysis name
# which is in the form : patientID_vs_control
get.analysis.name <- function(file.name){
  analysis.name <- sub("\\.\\w+", "", file.name, perl = T, ignore.case = T)
  analysis.name <- sub(".+/(deseq2|limma)_", "", analysis.name, perl = T, ignore.case = T)
  return(analysis.name)
}

######### Run GSEA for DESeq2 result
deseq2.res.dir <- "Result/PDCL/deseq2/"
deseq2.res.files <- list.files(deseq2.res.dir)
deseq2.res.files <- paste0(deseq2.res.dir, deseq2.res.files)
deseq2.gsea <- lapply(deseq2.res.files[1], run.gsea, reactome.gmt, bioplanet.gmt, function(f) { 
  x <- read.table(f, header = T, sep = ",")
  x <- as.data.frame(x) %>% filter(padj<0.1 & !is.na(padj)) %>% transmute(rank=log2FoldChange, gene=X) %>% arrange(desc(rank))
  rnk <- x$rank
  names(rnk) <- x$gene
  rnk
})
names(deseq2.gsea) <- sapply(deseq2.res.files, get.analysis.name)
saveRDS(deseq2.gsea, file = "Result/PDCL/gsea_deseq2_genes.rds")
save.gem.list.to.txt(deseq2.gsea, "Result/PDCL/gsea/deseq2/", "deseq2")

######### Run GSEA for Limma result
limma.res.dir <- "Result/PDCL/limma/"
limma.res.files <- list.files(limma.res.dir)
limma.res.files <- paste0(limma.res.dir, limma.res.files)
limma.gsea <- lapply(limma.res.files, run.gsea, reactome.gmt, bioplanet.gmt, function(f) { 
  x <- read.table(f, header = T, sep = ",")
  x <- as.data.frame(x) %>% filter(adj.P.Val<0.1 & !is.na(adj.P.Val)) %>% transmute(rank=logFC, gene=X) %>% arrange(desc(rank))
  rnk <- x$rank
  names(rnk) <- x$gene
  rnk
})
names(limma.gsea) <- sapply(limma.res.files, get.analysis.name)
saveRDS(limma.gsea, file = "Result/PDCL/gsea_limma_genes.rds")
save.gem.list.to.txt(limma.gsea, "Result/PDCL/gsea/limma/", "limma")
