# Remove precedent variables
rm(list=ls())

# Set working Dir
setwd("~/PhD/")

# Function computing the null distribution
# Perform nb.times random sampling without replacement of nb.genes from the list 
# and calculate the number of deregulated genes in the sample.
# return a vector containing the number of deregulated genes
# for each sampling
compute.null.distribution <- function(gene.vect, nb.times, nb.genes){
  
}

# Read a gmt file
# Return a list object indexing using the id of the pathway (1st column e). 
# Each element in the list contains the description of the pathway (2nd
# column) and the list of genes.
read.gmt <- function(gmt.file, delim="\t"){
  gmt.list <- list()
  for (l in strsplit(readLines(con = gmt.file), delim, perl=T) ) {
    gmt.list[[l[1]]] <- list("description"=l[2], "genes"= as.vector(l[-2:-1]) )
  }
  return(gmt.list)
}

# Perform the home made test.
# Compute a null distribution and calculate the quantile of 1-alpha -> C statistic.
# Calculate number of deregulated in the pathway -> T statistic.
# If C>T, the null hypothesis is rejected with a risk alpha.
# Return a data.frame containing the following column:
# - pathway_id : ID of the pathway (1st column of gmt file)
# - pathway_description : Description of the pathway (2nd column of gmt file)
# - T statistic : Quantile of order 1-alpha of the null distribution
# - C statistic : Number of deregulated genes
# - result : 0 if H0 accepted, 1 if H0 rejected
# - nb_deregulated : number of deregulated genes in the pathway
homemade.test <- function(gene.vect, nb.times, nb.genes, gmt, alpha) {
  
}

# Read penda result and delete the ".genes.results" string from the column names
penda.res <- read.csv2("Data/PDCL/results_combine.csv", check.names = F)
colnames(penda.res) <- gsub("\\.genes\\.results", "", colnames(penda.res), perl=T)
gmt <- read.gmt("Data/gene-set/jorge-gmt/reactome_gmt_symbol_no_unwanted_categ_no_less_10.gmt")

test.result <- list()
for (patient in colnames(penda.res)[1] ) {
  
}