######## Load Library needed
library(gprofiler2)
library(dplyr)

run.gost <- function(genes.symbol, domain, reactome.gmt, bioplanet.gmt){
  message("Running Gost ...")
  res <- list()
  res[["reactome"]] <- gost(query = genes.symbol, organism = reactome.gmt, user_threshold = 0.05, 
                            significant = F,  correction_method = "fdr", evcodes = T, custom_bg = domain)
  res[["bioplanet"]] <- gost(query = genes.symbol, organism = bioplanet.gmt, user_threshold = 0.05,
                             significant = F, correction_method = "fdr", evcodes = T, custom_bg = domain)
  return(res)
}

# run gost analysis with reactome and bioplanet GMT files for a list of DE analysis result
# the function will read a file and filter the gene according to a filtering function
run.gost.on.file <- function(f, reactome.gmt, bioplanet.gmt, filtering.function){
  message("###### G:Profiler Analysis for : ", f)
  genes.symbol <- read.table(f, header = T, row.names = 1, sep = ",")
  domain <- row.names(genes.symbol)
  genes.symbol <- filtering.function(genes.symbol)
  genes.symbol <- row.names(genes.symbol)
  res <- run.gost(genes.symbol, domain, reactome.gmt, bioplanet.gmt)
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
