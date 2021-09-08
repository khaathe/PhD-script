setwd("/home/spinicck/PhD/")

######### Remove every variable in memory
rm(list = ls())

######## Load Library needed
library(limma)
library(edgeR)

######### Load controls
astrocyte.file <- "Data/PDCL/astrocyte_GSE109001_counts_unique_genes_samples_filtered.txt"
astrocyte.dt <- read.table(astrocyte.file, header = T, sep = "\t", as.is = c("id", "symbol"))
control.samples.names <- c("AF22_NES_Astro_Br1_d29_37_S46","AF22_NES_Astro_Br2_d29_38_S56","AF22_NES_Astro_Br3_d29_39_S66",
             "CCF.STTG1_p24_Br1_S16","CCF.STTG1_p24_Br2_S17","CCF.STTG1_p24_Br3_S18",
             "CDIAstrocytes_p2_Br1_S19","CDIAstrocytes_p2_Br2_S20","CDIAstrocytes_p2_Br3_S21",
             "phaAstrocyte_p2_Br1_S1","phaAstrocyte_p2_Br2_S2","phaAstrocyte_p2_Br3_S3")

######### Load PDCL DataBase
## Load the PDCL count table
# Disabled check.names otherwise it will an X as column names aren't valid according to R
pdcl.count.table <- read.table("/home/spinicck/PhD/Data/PDCL/lignees_count_genes_PDCL.txt", sep = "\t", header = T, as.is = "gene_id", check.names = F)
pdcl.samples.name <- c("4339-p21","4371-p37","5706-p14","6190-p43","6240-p12","7015-p17","7060-p18",
                           "7142-p14","N13-1300","N13-1520-p9","N14-0072","N14-0870","N14-1208","N14-1525",
                           "N15_0460","N15_0516","N15-0385","N15-0661","N16_0535","N16-0240")

######### Prepare Object needed by DESeq2 for Analysis
count.data <- merge(astrocyte.dt, pdcl.count.table, by.x = "symbol", by.y = "gene_id")
row.names(count.data) <- count.data$symbol
count.data <- count.data[,-2:-1]
group <- as.factor(c(rep("control", length(control.samples.names)), rep("pdcl", length(pdcl.samples.name))))
limma.res.list <- list()

######### Differential Analysis

# Perform a differential analysis for each patient in the database against the others
for ( patient in pdcl.samples.name ){
  message("#### Limma Analysis for : ", patient)
  samples.index <- which(colnames(count.data) %in% c(control.samples.names, patient))
  dge <- DGEList(counts = count.data[,samples.index], group = group[samples.index])
  keep <- filterByExpr(dge)
  dge <- dge[keep, keep.lib.sizes=F]
  dge <- calcNormFactors(dge)
  logCpm <- cpm(dge, log = T, priorCount=3)
  fit <- lmFit(logCpm)
  fit <- eBayes(fit, robust = T)
  limma.res.list[[patient]] <- fit # store the result in a list for later data processing
  res.file <- paste0("Result/PDCL/limma/limma_", patient, "_vs_control.csv")
  message("Writing result to file : ", res.file)
  write.csv(topTable(fit, adjust.method = "BH", number = nrow(logCpm)), file = res.file)
  message("Done !")
}

saveRDS(limma.res.list, file = "Result/PDCL/limma_pdcl_vs_control.rds")
