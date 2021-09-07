setwd("/home/spinicck/PhD/")

######### Remove every variable in memory
rm(list = ls())

######## Load Library needed
library(DESeq2)

######### Load controls
astrocyte.file <- "Data/PDCL/Control/astrocyte/GSE109001_counts.txt"
astrocyte.dt <- read.table(astrocyte.file, header = T, sep = "\t", as.is = c("id", "symbol"))
# Only keep those samples
control.samples.names <- c("AF22_NES_Astro_Br1_d29_37_S46","AF22_NES_Astro_Br2_d29_38_S56","AF22_NES_Astro_Br3_d29_39_S66",
             "CCF.STTG1_p24_Br1_S16","CCF.STTG1_p24_Br2_S17","CCF.STTG1_p24_Br3_S18",
             "CDIAstrocytes_p2_Br1_S19","CDIAstrocytes_p2_Br2_S20","CDIAstrocytes_p2_Br3_S21",
             "phaAstrocyte_p2_Br1_S1","phaAstrocyte_p2_Br2_S2","phaAstrocyte_p2_Br3_S3")
astrocyte.dt <- astrocyte.dt[, c("symbol",control.samples.names)]
colnames(astrocyte.dt) <- c("gene_id", control.samples.names)
# Remove non unique gene names from the dataset
non.uniq.gene.id <- unique(astrocyte.dt$gene_id[!isUnique(astrocyte.dt$gene_id)])
astrocyte.dt <- astrocyte.dt[ !(astrocyte.dt$gene_id %in% non.uniq.gene.id), ]

######### Load PDCL DataBase

# Prepare count table with all PDCL samples and write a txt file
# Uncoment this code only to create a count table or update it 
# pdcl.db.dir <- "/home/spinicck/PhD/Data/PDCL/data_Annabelle/"
# pdcl.files.name <- list.files(path = pdcl.db.dir)
# pdcl.samples.name <- sub(".genes.results", "", pdcl.files.name)
# pdcl.count.table <- read.table(paste0(pdcl.db.dir, pdcl.files.name[1]), header = T, sep = "\t", as.is = "gene_id")
# pdcl.count.table <- pdcl.count.table[, c("gene_id", "expected_count")]
# pdcl.count.table$expected_count <- as.integer(pdcl.count.table$expected_count)
# colnames(pdcl.count.table)[2] <- pdcl.samples.name[1]
# for ( i in  2:length(pdcl.files.name) ){
#   count.table <- read.table(paste0(pdcl.db.dir, pdcl.files.name[i]), header = T, sep = "\t",as.is = "gene_id")
#   count.table$expected_count <- as.integer(count.table$expected_count)
#   count.table <- count.table[, c("gene_id", "expected_count")]
#   colnames(count.table)[2] <- pdcl.samples.name[i]
#   pdcl.count.table <- merge(pdcl.count.table, count.table)
# }
# write.table(pdcl.count.table, file = "/home/spinicck/PhD/Data/PDCL/lignees_count_genes_PDCL.txt",
#             sep = "\t", row.names = F, col.names = T)


## Load the PDCL count table
# Disabled check.names otherwise it will an X as column names aren't valid according to R
pdcl.count.table <- read.table("/home/spinicck/PhD/Data/PDCL/lignees_count_genes_PDCL.txt", sep = "\t", header = T, as.is = "gene_id", check.names = F)
pdcl.samples.name <- colnames(pdcl.count.table)[-1]

######### Prepare Object needed by DESeq2 for Analysis
count.data <- merge(astrocyte.dt, pdcl.count.table)
row.names(count.data) <- count.data$gene_id
count.data <- as.matrix(count.data[,-1])
condition.vect <- c(rep("control", length(control.samples.names)), rep("pdcl", length(pdcl.samples.name)))
col.data <- data.frame(sample = c(control.samples.names, pdcl.samples.name), condition=condition.vect)
col.data$condition <- as.factor(col.data$condition)
deseq2.res.list <- list()

######### Differential Analysis

# Perform a differential analysis for each patient in the database against the others
for ( patient in pdcl.samples.name ){
  message("#### DESeq2 Analysis for : ", patient)
  dataset.samples <- c(control.samples.names, patient)
  dataset.matrix <- count.data[,dataset.samples]
  dataset.coldata <- col.data[ col.data$sample %in% dataset.samples, ]
  dds <- DESeqDataSetFromMatrix(countData = dataset.matrix, colData = dataset.coldata, design = ~condition)
  dds <- DESeq(dds)
  res <- results(dds, alpha = 0.05)
  deseq2.res.list[[patient]] <- res # store the result in a list for later data propcessing
  res.file <- paste0("Result/PDCL/deseq2/deseq2_", patient, "_vs_control.csv")
  message("Writing result to file : ", res.file)
  write.csv(as.data.frame(res), file = res.file)
  message("Done !")
}
