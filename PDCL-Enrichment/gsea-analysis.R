######## Load Library needed
library(fgsea)
library(dplyr)
library(tidyr)

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
      f <- paste0(save.dir, "gsea_", de.analysis, "_",analysis, "_", db.source, ".txt")
      write.table(gem.list[[analysis]][[db.source]], file = f, row.names = F, col.names = T, sep = "\t", quote=F)
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

# Convert a list of GSEA Analysis result to a list of Data Frame containing the column in the Generic Enrichment Map (GEM)
# as described in GEM documentation (https://enrichmentmap.readthedocs.io/en/latest/FileFormats.html#generic-results-files)
convert.gsea.to.gem <- function(res){
  lapply(res, function(x){
    x %>% dplyr::rename(pathway_id = pathway, p_value= pval, fdr=padj) %>%
      mutate(phenotype = sign(ES), description = "", genes = sapply(leadingEdge, function(e){ paste0(e, collapse = ",")} ) ) %>%
      select(pathway_id,description, p_value, fdr, phenotype, genes) %>%
      tidyr::drop_na()
  })
}
