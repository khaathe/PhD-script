######## Load Library needed
library(dplyr)
library(tidyr)

# Filter G:Profiler result to keep only pathways deregulated, and turn the result in GEM format
filter.gprofiler.pathways <- function(gost){
  gost.gem <- gost$result %>% dplyr::rename(pathway_id = term_id, description = term_name, genes = intersection) %>%
    mutate(fdr = p_value, phenotype = 1) %>%
    select(pathway_id,description, p_value, fdr, phenotype, genes) %>%
    filter(fdr < 0.1)
  message("Number of deregulated pathways for G:Profiler: ", nrow(gost.gem))
  return(gost.gem)
}

# Filter GSEA result to keep only pathways deregulated, and turn the result in GEM format
filter.gsea.pathways <- function(gsea){
  gsea.gem <- gsea %>% dplyr::rename(pathway_id = pathway, p_value= pval, fdr=padj) %>%
    mutate(phenotype = sign(ES), description = "", genes = sapply(leadingEdge, function(e){ paste0(e, collapse = ",")} ) ) %>%
    select(pathway_id,description, p_value, fdr, phenotype, genes) %>%
    tidyr::drop_na() %>%
    filter(fdr < 0.1)
  message("Number of deregulated pathways for GSEA: ", nrow(gsea.gem))
  return(gsea.gem)
}

# Return the deregulated pathways found in both G:Profiler+PENDA and GSEA+DESeq2 analysis
# in one enrichment method (bioplanet, reactome, ...)
get.common.deregulated.pathways <- function(gost.penda, gsea.deseq2){
  common.pathways <- inner_join(gost.penda, gsea.deseq2, by = "pathway_id") %>%
    dplyr::rename(description = description.x, genes_gprofiler = genes.x, genes_gsea = genes.y, p_value_gprofiler = p_value.x, p_value_gsea = p_value.y,
           fdr_gsea = fdr.y, phenotype_gsea = phenotype.y) %>%
    dplyr::select(pathway_id, description, p_value_gprofiler, p_value_gsea, fdr_gsea, 
           phenotype_gsea, genes_gprofiler, genes_gsea)
  message("Numer of common deregulated pathways: ", nrow(common.pathways))
  return(common.pathways)
}

# Return all the deregulated pathways found in both G:Profiler+PENDA and GSEA+DESeq2 analysis
get.all.common.deregulated.pathways <- function(p, gost.penda, gsea.deseq2, database){
  message("###### Getting Common Pathways for: ", p)
  element <- list()
  for (db in database){
    message("### Treating: ", db, " database")
    gost.gem <- filter.gprofiler.pathways(gost.penda[[p]][[db]])
    gsea.gem <- filter.gsea.pathways(gsea.deseq2[[p]][[db]])
    element[[db]] <- get.common.deregulated.pathways(gost.gem, gsea.gem)
  }
  return(element)
}

# Save the common deregulated pathways to tab delimited txt files
save.result <- function(pathways.list, save.dir){
  for (patient in names(pathways.list) ) {
    for (db.source in names(pathways.list[[patient]]) ){
      if (nrow(pathways.list[[patient]][[db.source]]) > 0) {
        f <- paste0(save.dir, "common_deregulated_pathways_", patient, "_", db.source, ".txt")
        write.table(pathways.list[[patient]][[db.source]], file = f, row.names = F, col.names = T, sep = "\t", quote = F) 
      }
    }
  }
}
