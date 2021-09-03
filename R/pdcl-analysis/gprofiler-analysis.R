######### Remove every variable in memory
rm(list = ls())

######## Load Library needed
library(gprofiler2)

######### Source all functions
rfun_dir <- "~/Bureau/EnrichmentAnalysis/Script/R/function/"
for(f in list.files(rfun_dir)) {
  source(paste0(rfun_dir, f, ""))
}
rm(f, rfun_dir)

######### Upload GMT File before running G:Profiler Enrichment Analysis
# This portion is commented because the file are already uploaded, use the gmt token instead
# upload_GMT_file(gmtfile = "/home/spinicck/PhD/Data/gene-set/reactome_pathways_symbol.gmt")
# upload_GMT_file(gmtfile = "/home/spinicck/PhD/Data/gene-set/bioplanet_pathways_symbol.gmt")

######### Define GMT token to use for Enrichment Analysis
reactome.gmt.token <- "gp__HCf5_1G7e_SCI"
bioplanet.gmt.token <- "gp__7tZx_iEOw_ros"

######### Load all DE Analysis Results
deseq2.res.dir <- "/home/spinicck/PhD/Data/PDCL/deseq2-analysis/"
deseq2.res.files <- list.files(deseq2.res.dir)
deseq2.res.files <- paste0(deseq2.res.dir, deseq2.res.files)
limma.res.dir <- "/home/spinicck/PhD/Data/PDCL/limma/"
limma.res.files <- list.files(limma.res.dir)
limma.res.files <- paste0(limma.res.dir, limma.res.files)
res.de.files <- c(deseq2.res.files, limma.res.files)

######### Run G:Profiler Enrichment Analysis
gostres <- lapply(res.de.files, function(f){
  message("###### G:Profiler Analysis for : ", f)
  genes.symbol <- read.table(f, header = T, row.names = 1, sep = ",")
  genes.symbol <- genes.symbol[genes.symbol$padj<1,]
  genes.symbol <- row.names(genes.symbol)
  res <- list()
  res[["reactome"]] <- gost(query = genes.symbol, organism = reactome.gmt.token, user_threshold = 0.05, 
                                significant = F,  correction_method = "fdr")
  res[["bioplanet"]] <- gost(query = genes.symbol, organism = bioplanet.gmt.token, user_threshold = 0.05,
                                 significant = F, correction_method = "fdr")
  return(res)
})

analysis.name <- sub("\\.\\w+", "", res.de.files, perl = T)
analysis.name <- sub("/.+/(deseq2|limma)_", "", analysis.name, perl = T, ignore.case = T)
names(gostres) <- analysis.name

######### Convert G:Profiler result into more friendly enrichment result
# We keep the Gost result as is in case for later use
enrichment.res <- lapply(gostres,function(res){
  lapply(res, function(x){
    x$result %>% dplyr::rename(pathway_name = term_id,pathway_id = term_name, pathway_size = term_size, p_value_adjusted = p_value) %>%
      select(pathway_id,pathway_name, p_value_adjusted,source, pathway_size, query_size,intersection_size)
  })
})
names(enrichment.res) <- analysis.name

######### Save the more friendly into more useful text file
save.dir <- "/home/spinicck/PhD/Data/PDCL/gprofiler/"
for (analysis in names(enrichment.res) ) {
  for (db.source in names(enrichment.res[[analysis]]) ){
    f <- paste0(save.dir, "gprofiler_", analysis, "_", db.source, ".csv")
    write.table(enrichment.res[[analysis]][[db.source]], file = f, row.names = F, col.names = T, sep = "\t")
  }
}
