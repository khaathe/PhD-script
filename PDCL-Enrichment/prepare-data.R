setwd("/home/spinicck/PhD/")

######### Remove every variable in memory
rm(list = ls())

######### Load Library (if needed)

######### Create a Control Raw Count Table
# Data can be retrieved from GEO using the accession number : GSE109001
# Gene expressions raw count comes from in vitro cultived astrocytes
astrocyte.file <- "Data/PDCL/Control/astrocyte/GSE109001_counts.txt"
astrocyte.dt <- read.table(astrocyte.file, header = T, sep = "\t", as.is = c("id", "symbol"))
# Samples to keep
control.samples.names <- c("AF22_NES_Astro_Br1_d29_37_S46","AF22_NES_Astro_Br2_d29_38_S56","AF22_NES_Astro_Br3_d29_39_S66",
             "CCF.STTG1_p24_Br1_S16","CCF.STTG1_p24_Br2_S17","CCF.STTG1_p24_Br3_S18",
             "CDIAstrocytes_p2_Br1_S19","CDIAstrocytes_p2_Br2_S20","CDIAstrocytes_p2_Br3_S21",
             "phaAstrocyte_p2_Br1_S1","phaAstrocyte_p2_Br2_S2","phaAstrocyte_p2_Br3_S3")
# Filter to keep only gene ids, symbols and the samples of interest
astrocyte.dt <- astrocyte.dt[, c("id", "symbol",control.samples.names)]
# All ids are unique bu not the genes symol. Here we remove genes with duplicated symbols
# as they can caused trouble later (for example during DESeq2 analysis).
non.uniq.gene.symbol <- unique(astrocyte.dt$symbol[!isUnique(astrocyte.dt$symbol)])
astrocyte.dt <- astrocyte.dt[ !(astrocyte.dt$symbol %in% non.uniq.gene.symbol), ]
write.table(astrocyte.dt, file = "Data/PDCL/astrocyte_GSE109001_counts_unique_genes_samples_filtered.txt",
            sep = "\t", row.names = F, col.names = T)

######### Create PDCL Raw Count Table
# Loop through all files in a data directory and merge them (by gene symbol)
# to create a data table containing all the genes raw counts for each patient
# The resulting DataTable is saved in a tab delimited txt file.
pdcl.db.dir <- "Data/PDCL/data_Annabelle/"
pdcl.files.name <- list.files(path = pdcl.db.dir)
pdcl.samples.name <- sub(".genes.results", "", pdcl.files.name)
pdcl.count.table <- read.table(paste0(pdcl.db.dir, pdcl.files.name[1]), header = T, sep = "\t", as.is = "gene_id")
pdcl.count.table <- pdcl.count.table[, c("gene_id", "expected_count")]
pdcl.count.table$expected_count <- as.integer(pdcl.count.table$expected_count)
colnames(pdcl.count.table)[2] <- pdcl.samples.name[1]
for ( i in  2:length(pdcl.files.name) ){
  count.table <- read.table(paste0(pdcl.db.dir, pdcl.files.name[i]), header = T, sep = "\t",as.is = "gene_id")
  count.table$expected_count <- as.integer(count.table$expected_count)
  count.table <- count.table[, c("gene_id", "expected_count")]
  colnames(count.table)[2] <- pdcl.samples.name[i]
  pdcl.count.table <- merge(pdcl.count.table, count.table)
}
write.table(pdcl.count.table, file = "Data/PDCL/lignees_count_genes_PDCL.txt",
            sep = "\t", row.names = F, col.names = T)
