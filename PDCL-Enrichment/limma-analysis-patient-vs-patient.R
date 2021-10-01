setwd("/home/spinicck/PhD/")

######### Remove every variable in memory
rm(list = ls())

######## Load Library needed
library(limma)
library(readxl)

######### Load PDCL DataBase
# data type in the excel sheet, we do not want R to guess
col.type <- c(rep("text", 2), rep("numeric", 25)) 
pdcl.tpm <- read_xlsx(path = "Data/PDCL/lignees_TPM_genes_ALL_PDCL.xlsx", sheet = 1, col_names = T, col_types = col.type)
pdcl.tpm <- as.data.frame(pdcl.tpm)
row.names(pdcl.tpm) <- pdcl.tpm$gene_id
pdcl.tpm <- pdcl.tpm[, -c(1:4)] # remove column that do not contain tpm expression values
pdcl.tpm <- pdcl.tpm[rowSums(pdcl.tpm) != 0, ] # remove genes with no expression across all patient
nb.samples <- ncol(pdcl.tpm)
p.to.p.de.result <- list() # used to store result for later processing

######### Differential Expression Analysis for each patient against each other
# Normally we should normalize expression values by using either cpm or voom method of limma
# But here epxression values are already normalize using TPM
# Some forums suggest that differential expression analysis should be executed with raw count values
# as fpkm and tpm normalize away some interesting informations (ex: sequencing depth)
# Thus result here should be taken carefully
save.dir <- "Result/PDCL/limma/"
for ( i in 1:(nb.samples-1) ){
  for (j in (i+1):nb.samples){
    p.to.p.tpm <- pdcl.tpm[,c(i,j)]
    res.name <- paste0("limma_" ,colnames(p.to.p.tpm)[1], "_vs_", colnames(p.to.p.tpm)[2])
    message("Limma Analysis : ", res.name)
    fit <- lmFit(p.to.p.tpm)
    fit <- eBayes(fit)
    p.to.p.de.result[[res.name]] <- fit
    de.gene.table <- topTable(fit, number = nrow(pdcl.tpm), adjust.method = "BH")
    # Rename the column so that the logFoldChange, p-value and p-value adjusted column names are the same
    # between DESeq2 and Limma for further analysis
    colnames(de.gene.table) <- c("log2FoldChange","AveExpr","t","pvalue","padj","B")
    write.csv(de.gene.table, file = paste0(save.dir,res.name, ".csv"))
  }
}
