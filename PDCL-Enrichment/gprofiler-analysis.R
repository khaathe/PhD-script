setwd("/home/spinicck/PhD/")

######### Remove every variable in memory
rm(list = ls())

######## Load Library needed
library(gprofiler2)
library(dplyr)

######### Upload GMT File before running G:Profiler Enrichment Analysis
# This portion is commented because the file are already uploaded, use the gmt token instead
# upload_GMT_file(gmtfile = "/home/spinicck/PhD/Data/gene-set/reactome_pathways_symbol.gmt")
# upload_GMT_file(gmtfile = "/home/spinicck/PhD/Data/gene-set/bioplanet_pathways_symbol.gmt")

######### Define GMT token to use for Enrichment Analysis
reactome.gmt.token <- "gp__HCf5_1G7e_SCI"
bioplanet.gmt.token <- "gp__7tZx_iEOw_ros"

# run gost analysis with reactome and bioplanet GMT files for a list of DE analysis result
run.gost <- function(f, reactome.gmt, bioplanet.gmt, filtering.function){
  message("###### G:Profiler Analysis for : ", f)
  genes.symbol <- read.table(f, header = T, row.names = 1, sep = ",")
  genes.symbol <- filtering.function(genes.symbol)
  genes.symbol <- row.names(genes.symbol)
  res <- list()
  res[["reactome"]] <- gost(query = genes.symbol, organism = reactome.gmt, user_threshold = 0.05, 
                            significant = F,  correction_method = "fdr")
  res[["bioplanet"]] <- gost(query = genes.symbol, organism = bioplanet.gmt, user_threshold = 0.05,
                             significant = F, correction_method = "fdr")
  return(res)
}

# Convert a list of Gost Analysis result to a list of Data Frame which is similar to Generic Enrichment Map (GEM)
convert.gost.to.gem <- function(res){
  lapply(res, function(x){
    x$result %>% dplyr::rename(pathway_name = term_id,pathway_id = term_name, pathway_size = term_size, p_value_adjusted = p_value) %>%
      select(pathway_id,pathway_name, p_value_adjusted,source, pathway_size, query_size,intersection_size)
  })
}

# Save a list of Generic Enrichment Map to tab delimited txt files
# Files name conttain the tool used for differential expression analysis, which patient
# has been compared to the control (analysis), the source of the gmt file (db.source)
save.gem.list.to.txt <- function(gem.list, save.dir, de.analysis){
  for (analysis in names(gem.list) ) {
    for (db.source in names(gem.list[[analysis]]) ){
      f <- paste0(save.dir, "gprofiler_", de.analysis, "_",analysis, "_", db.source, ".csv")
      write.table(gem.list[[analysis]][[db.source]], file = f, row.names = F, col.names = T, sep = "\t")
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

######### RUn G:Profiler for DESeq2 result
deseq2.res.dir <- "Result/PDCL/deseq2/"
deseq2.res.files <- list.files(deseq2.res.dir)
deseq2.res.files <- paste0(deseq2.res.dir, deseq2.res.files)
deseq2.gostres <- lapply(deseq2.res.files[1], run.gost, reactome.gmt.token, bioplanet.gmt.token, function(x) { 
    as.data.frame(x) %>% filter(padj<0.1 & !is.na(padj))
})
names(deseq2.gostres) <- sapply(deseq2.res.files, get.analysis.name)
saveRDS(deseq2.gostres, file = "Result/gost_deseq2_genes.rds")
deseq2.enrichment <- lapply(deseq2.gostres, convert.gost.to.gem)
save.gem.list.to.txt(deseq2.enrichment, "Result/PDCL/gprofiler/deseq2/", "deseq2")


######### RUn G:Profiler for Limma result
limma.res.dir <- "Result/PDCL/limma/"
limma.res.files <- list.files(limma.res.dir)
limma.res.files <- paste0(limma.res.dir, limma.res.files)
limma.gostres <- lapply(limma.res.files, run.gost, reactome.gmt.token, bioplanet.gmt.token, function(x) { 
  as.data.frame(x) %>% filter(adj.P.Val<0.1 & !is.na(adj.P.Val))
})
names(limma.gostres) <- sapply(limma.res.files, get.analysis.name)
saveRDS(limma.gostres, file = "Result/gost_limma_genes.rds")
limma.enrichment <- lapply(limma.gostres, convert.gost.to.gem)
save.gem.list.to.txt(limma.enrichment, "Result/PDCL/gprofiler/limma/", "limma")
