setwd("/home/spinicck/PhD/")

######### Remove every variable in memory
rm(list = ls())

######## Load Library needed
library(gprofiler2)
library(dplyr)

######### Upload GMT File before running G:Profiler Enrichment Analysis
# This portion is commented because the file are already uploaded, use the gmt token instead
# reactome.gmt.token <- upload_GMT_file(gmtfile = "Data/gene-set/jorge-gmt/reactome_gmt_symbol_no_unwanted_categ_no_less_10.gmt")
# bioplanet.gmt.token <- upload_GMT_file(gmtfile = "Data/gene-set/jorge-gmt/bioplanet_gmt_symbol_no_unwanted_categ_no_less_10.gmt")

######### Define GMT token to use for Enrichment Analysis
reactome.gmt.token <- "gp__Nqmx_yxm7_37U"
bioplanet.gmt.token <- "gp__vB7V_DT57_yNY"

run.gost <- function(genes.symbol, reactome.gmt, bioplanet.gmt){
  message("Running Gost ...")
  res <- list()
  res[["reactome"]] <- gost(query = genes.symbol, organism = reactome.gmt, user_threshold = 0.05, 
                            significant = F,  correction_method = "fdr", evcodes = T)
  res[["bioplanet"]] <- gost(query = genes.symbol, organism = bioplanet.gmt, user_threshold = 0.05,
                             significant = F, correction_method = "fdr", evcodes = T)
  return(res)
}

# run gost analysis with reactome and bioplanet GMT files for a list of DE analysis result
# the function will read a file and filter the gene according to a filtering function
run.gost.on.file <- function(f, reactome.gmt, bioplanet.gmt, filtering.function){
  message("###### G:Profiler Analysis for : ", f)
  genes.symbol <- read.table(f, header = T, row.names = 1, sep = ",")
  genes.symbol <- filtering.function(genes.symbol)
  genes.symbol <- row.names(genes.symbol)
  res <- run.gost(genes.symbol, reactome.gmt, bioplanet.gmt)
  return(res)
}

# Convert a list of Gost Analysis result to a list of Data Frame containing the column in the Generic Enrichment Map (GEM)
# as described in G:Profiler documentation (https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html#enrichment-analysis)
convert.gost.to.gem <- function(res){
  lapply(res, function(x){
    x$result %>% dplyr::rename(pathway_id = term_id, description = term_name, genes = intersection) %>%
      mutate(fdr = p_value, phenotype = 1) %>%
      select(pathway_id,description, p_value, fdr, phenotype, genes)
  })
}

# Convert a list of Gost Analysis result to a list of Data Frame with the same column defined in Jorge Scripts
convert.gost.to.enrichment <- function(res){
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
      f <- paste0(save.dir, "gprofiler_", de.analysis, "_",analysis, "_", db.source, ".txt")
      write.table(gem.list[[analysis]][[db.source]], file = f, row.names = F, col.names = T, sep = "\t", quote = F)
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
deseq2.gostres <- lapply(deseq2.res.files, run.gost.on.file, reactome.gmt.token, bioplanet.gmt.token, function(x) { 
    as.data.frame(x) %>% filter(padj<0.1 & !is.na(padj))
})
names(deseq2.gostres) <- sapply(deseq2.res.files, get.analysis.name)
saveRDS(deseq2.gostres, file = "Result/PDCL/gost_deseq2_genes.rds")
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
saveRDS(limma.gostres, file = "Result/PDCL/gost_limma_genes.rds")
limma.enrichment <- lapply(limma.gostres, convert.gost.to.gem)
save.gem.list.to.txt(limma.enrichment, "Result/PDCL/gprofiler/limma/", "limma")

######### RUn G:Profiler for PENDA result
penda.res.file <- "Data/PDCL/results_combine.csv"
penda.res <- read.table(penda.res.file, sep = ";", header = T, row.names = 1, check.names = F)
pdcl.samples.names <- names(penda.res)
pdcl.samples.names <- sub(".genes.results", "", pdcl.samples.names, ignore.case = T)
names(penda.res) <- pdcl.samples.names
all.pdcl.gene <- row.names(penda.res)
penda.gost.res <- lapply(penda.res, function(patient){
  gene.list <- all.pdcl.gene[patient != 0]
  res <- run.gost(gene.list, reactome.gmt.token, bioplanet.gmt.token)
  return(res)
})
saveRDS(penda.gost.res, file = "Result/PDCL/gost_penda_genes.rds")
penda.enrichment <- lapply(penda.gost.res, convert.gost.to.gem)
save.gem.list.to.txt(penda.enrichment, "Result/PDCL/gprofiler/penda/", "penda")
