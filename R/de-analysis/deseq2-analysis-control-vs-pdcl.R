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

######### Load controls
astrocyte.file <- "/home/spinicck/PhD/Data/PDCL/Control/astrocyte/GSE109001_counts.txt"
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
# Use this function if the file isn't already created or need update
prepare.pdcl.count.table <- function (){
  pdcl.db.dir <- "/home/spinicck/PhD/Data/PDCL/data_Annabelle/"
  pdcl.files.name <- list.files(path = pdcl.db.dir)
  pdcl.samples.name <- sub(".genes.results", "", pdcl.files.name)
  pdcl.count.table <- read.table(paste0(pdcl.db.dir, pdcl.files.name[1]), header = T, sep = "\t", as.is = "gene_id")
  pdcl.count.table <- pdcl.count.table[, c("gene_id", "expected_count")]
  colnames(pdcl.count.table)[2] <- pdcl.samples.name[1]
  for ( i in  2:length(pdcl.files.name) ){
    count.table <- read.table(paste0(pdcl.db.dir, pdcl.files.name[i]), header = T, sep = "\t",as.is = "gene_id")
    count.table <- count.table[, c("gene_id", "expected_count")]
    colnames(count.table)[2] <- pdcl.samples.name[i]
    pdcl.count.table <- merge(pdcl.count.table, count.table)
  }
  write.table(pdcl.count.table, file = "/home/spinicck/PhD/Data/PDCL/lignees_count_genes_PDCL.txt",
              sep = "\t", row.names = F, col.names = T)
}
## Load the PDCL count table
# Disabled check.names otherwise it will an X as column names aren't valid according to R
pdcl.count.table <- read.table("/home/spinicck/PhD/Data/PDCL/lignees_count_genes_PDCL.txt", sep = "\t", header = T, as.is = "gene_id", check.names = F)
pdcl.count.table <- lapply(pdcl.count.table, function(x){ if(is.numeric(x)){ return(as.integer(x)) } else {return(x)} } ) 
pdcl.count.table <- as.data.frame(pdcl.count.table)

######### Prepare DESeqDataSet object for Analysis
count.data <- merge(astrocyte.dt, pdcl.count.table)
row.names(count.data) <- count.data$gene_id
count.data <- as.matrix(count.data[,-1])
control.samples <- colnames(astrocyte.dt[,-1])
pdcl.samples <- colnames(pdcl.count.table[,-1])
col.data <- data.frame(row.names = c(control.samples, pdcl.samples), condition=c(rep("control", length(control.samples)), rep("pdcl", length(pdcl.samples))))
col.data$condition <- as.factor(col.data$condition)
dds <- DESeqDataSetFromMatrix(countData = count.data, colData = col.data, design = ~condition)

####### Run differential expression analysis
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)
write.csv(as.data.frame(res), file = "/home/spinicck/PhD/Data/PDCL/pdcl_vs_control.csv")
