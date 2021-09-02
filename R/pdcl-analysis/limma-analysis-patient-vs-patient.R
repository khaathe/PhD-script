######### Remove every variable in memory
rm(list = ls())

######## Load Library needed
library(limma)
library(edgeR)
library(readxl)

######### Source all functions
rfun_dir <- "~/Bureau/EnrichmentAnalysis/Script/R/function/"
for(f in list.files(rfun_dir)) {
  source(paste0(rfun_dir, f, ""))
}
rm(f, rfun_dir)

######### Load PDCL DataBase
# data type in the excel sheet, we do not want R to guess
col.type <- c(rep("text", 2), rep("numeric", 25)) 
pdcl.tpm <- read_xlsx(path = "/home/spinicck/PhD/Data/PDCL/lignees_TPM_genes_ALL_PDCL.xlsx", sheet = 1, col_names = T, col_types = col.type)
pdcl.tpm <- as.data.frame(pdcl.tpm)
row.names(pdcl.tpm) <- pdcl.tpm$gene_id
pdcl.tpm <- pdcl.tpm[, -c(1:4)] # remove column that do not contain tpm expression values
nb.samples <- ncol(pdcl.tpm)
p.to.p.de.result <- list() # used to store result for later processing

######### Differential Expression Analysis for each patient against each other
save.dir <- "/home/spinicck/PhD/Data/PDCL/limma/"
for ( i in 1:(nb.samples-1) ){
  for (j in (i+1):nb.samples){
    p.to.p.tpm <- pdcl.tpm[,c(i,j)]
    fit <- lmFit(p.to.p.tpm)
    fit <- eBayes(fit)
    res.name <- paste0(colnames(p.to.p.tpm)[1], "_vs_", colnames(p.to.p.tpm)[2])
    p.to.p.de.result[[res.name]] <- fit
    write.csv(topTable(fit, number = nrow(pdcl.tpm)), file = paste0(save.dir,res.name, ".csv"))
  }
}
