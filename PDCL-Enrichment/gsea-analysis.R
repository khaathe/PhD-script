######### Remove every variable in memory
rm(list = ls())

######## Load Library needed
library(fgsea)
library(stringr)

######### Source all functions
rfun_dir <- "~/Bureau/EnrichmentAnalysis/Script/R/function/"
for(f in list.files(rfun_dir)) {
  source(paste0(rfun_dir, f, ""))
}
rm(f, rfun_dir)

#IMPORTANT!! We must choose a seed to ensure the analysis is reproducible. This is because GSEA uses a 
#permutation test to span the null distribution of the statistic. Hence, each time we run it, results are 
#expected to change slightly.
set.seed(1)

reactome.gmt <- gmtPathways("/home/spinicck/PhD/Data/gene-set/reactome_pathways_symbol.gmt")
bioplanet.gmt <- gmtPathways("/home/spinicck/PhD/Data/gene-set/bioplanet_pathways_symbol.gmt")

######### Load all DE Analysis Results
deseq2.res.dir <- "/home/spinicck/PhD/Data/PDCL/deseq2-analysis/"
deseq2.res.files <- list.files(deseq2.res.dir)
deseq2.res.files <- paste0(deseq2.res.dir, deseq2.res.files)
limma.res.dir <- "/home/spinicck/PhD/Data/PDCL/limma/"
limma.res.files <- list.files(limma.res.dir)
limma.res.files <- paste0(limma.res.dir, limma.res.files)
res.de.files <- c(deseq2.res.files, limma.res.files)

####### Function for creating a ranked list of genes
create.ranked.gene.list <- function(genes.dt){
  genes.dt <- na.omit(genes.dt)
  genes.dt <- genes.dt[genes.dt$padj<1,]
  ranked <- genes.dt$padj
  names(ranked) <- row.names(genes.dt)
  return(ranked)
}

######### Run GSEA Analysis
gseares <- lapply(res.de.files, function(f){
  message("###### GSEA Analysis for : ", f)
  genes.dt <- read.table(f, header = T, row.names = 1, sep = ",")
  ranked.list <- create.ranked.gene.list(genes.dt)
  res <- list()
  message("\tRunning Reactome Enrichment ...")
  res[["reactome"]] <- fgsea(reactome.gmt, ranked.list, minSize=15, maxSize=500, nperm = 1000)
  message("\tDone !")
  message("\tRunning Bioplanet Enrichment ...")
  res[["bioplanet"]] <- fgsea(bioplanet.gmt, ranked.list, minSize=15, maxSize=500, nperm = 1000)
  message("\tDone !")
  return(res)
})

analysis.name <- sub("\\.\\w+", "", res.de.files, perl = T)
analysis.name <- sub("/.+/(deseq2|limma)_", "", analysis.name, perl = T, ignore.case = T)
names(gseares) <- analysis.name

######### Save Result to Enrichment Files
save.dir <- "/home/spinicck/PhD/Data/PDCL/gsea/"
for (analysis in names(gseares) ) {
  for (db.source in names(gseares[[analysis]]) ){
    f <- paste0(save.dir, "gsea_", analysis, "_", db.source, ".csv")
    enrichment.res <- gseares[[analysis]][[db.source]]
    enrichment.res$GeneList <- sapply(enrichment.res$leadingEdge, function(x){ paste0(x, collapse = ", ")})
    enrichment.res <- enrichment.res[, !"leadingEdge" ]
    write.table(enrichment.res, file = f, row.names = F, col.names = T, sep = "\t")
  }
}
