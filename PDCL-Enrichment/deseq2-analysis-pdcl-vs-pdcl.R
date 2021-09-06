######### Remove every variable in memory
rm(list = ls())

######## Load Library needed
library(DESeq2)

######### Source all functions
rfun_dir <- "~/Bureau/EnrichmentAnalysis/Script/R/function/"
for(f in list.files(rfun_dir)) {
  source(paste0(rfun_dir, f, ""))
}
rm(f, rfun_dir)

######### Load PDCL DataBase

# Prepare count table with all PDCL samples and write a txt file
# Uncoment this code only to create a count table or update it 
# pdcl.db.dir <- "/home/spinicck/PhD/Data/PDCL/data_Annabelle/"
# pdcl.files.name <- list.files(path = pdcl.db.dir)
# pdcl.samples.name <- sub(".genes.results", "", pdcl.files.name)
# pdcl.count.table <- read.table(paste0(pdcl.db.dir, pdcl.files.name[1]), header = T, sep = "\t", as.is = "gene_id")
# pdcl.count.table <- pdcl.count.table[, c("gene_id", "expected_count")]
# colnames(pdcl.count.table)[2] <- pdcl.samples.name[1]
# for ( i in  2:length(pdcl.files.name) ){
#   count.table <- read.table(paste0(pdcl.db.dir, pdcl.files.name[i]), header = T, sep = "\t",as.is = "gene_id")
#   count.table <- count.table[, c("gene_id", "expected_count")]
#   colnames(count.table)[2] <- pdcl.samples.name[i]
#   pdcl.count.table <- merge(pdcl.count.table, count.table)
# }
# write.table(pdcl.count.table, file = "/home/spinicck/PhD/Data/PDCL/lignees_count_genes_PDCL.txt",
#             sep = "\t", row.names = F, col.names = T)

# Load the PDCL count table
# Disabled check.names otherwise it will an X as column names aren't valid according to R
pdcl.count.table <- read.table("/home/spinicck/PhD/Data/PDCL/lignees_count_genes_PDCL.txt", 
                               sep = "\t", header = T, as.is = "gene_id", check.names = F, row.names = 1)
# R use the type numeric instead of integer when reading the count table which throw an error later with DESeq2
# Here we convert all numeric column into integers
for (c in colnames(pdcl.count.table)) { pdcl.count.table[, c] <- as.integer(pdcl.count.table[,c]) }
rm(c)

#Prepare metadata for DESeqDataSet Object
col.data <- data.frame(row.names = colnames(pdcl.count.table), condition=factor(rep("pdcl", ncol(pdcl.count.table)), levels = c("selected", "pdcl")))
deseq2.res.list <- list()

######### Differential Analysis

# Perform a differential analysis for each patient in the database against the others
for ( patient in colnames(pdcl.count.table) ){
  col.data.copy <- col.data
  col.data.copy[patient, "condition"] <- "selected"
  dds <- DESeqDataSetFromMatrix(countData = pdcl.count.table, colData = col.data.copy, design = ~condition)
  dds <- DESeq(dds)
  res <- results(dds, alpha = 0.05)
  deseq2.res.list[[patient]] <- res # store the result in a list for later data propcessing
  write.csv(as.data.frame(res), file = paste0("/home/spinicck/PhD/Data/PDCL/deseq2-analysis/deseq2_", patient, "_vs_pdcl.csv"))
}
rm(patient, dds, col.data.copy, res)

