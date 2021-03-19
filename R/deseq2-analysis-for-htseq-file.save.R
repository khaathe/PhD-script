rm(list = ls())

########## Init Packages
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(httr)
library(jsonlite)
library(data.table)
library(fgsea)
library(data.table)
library(ggplot2)

######### Source all functions
rfun_dir <- "~/Bureau/EnrichmentAnalysis/Script/R/function/"
for(f in list.files(rfun_dir)) {
  source(paste0(rfun_dir, f, ""))
}
rm(f, rfun_dir)

########## Init Working Directory
wDir = "/home/spinicck/Bureau/EnrichmentAnalysis/R"
setwd(wDir)

######## Init Constants
htseqFile <- "/home/spinicck/Bureau/EnrichmentAnalysis/Data/TCGA-GBM/tcga-gbm-sample-table.txt"
htseqDir <- "/home/spinicck/Bureau/EnrichmentAnalysis"
resultDir <- "/home/spinicck/Bureau/EnrichmentAnalysis/Result/deseq2/tcga-gbm/"
dataDir <- "/home/spinicck/Bureau/EnrichmentAnalysis/Data/"

####### Reading gene informations
col.class <- c("character", "character", 	"character",	"character",	"character",	"logical",	"logical",	"logical",	"logical")
ginf <- read.csv(file = paste0(dataDir, "TCGA-GBM/gene-info.csv", ""), header = T, colClasses = col.class, row.names = 1)

######## Creating DataSet
data.htseq <- read.table(file = htseqFile, header = TRUE, sep="\t")
sampleTable <- data.frame(sampleName = data.htseq$sampleName, 
                          fileName = data.htseq$fileName, 
                          condition = data.htseq$condition)
sampleTable$condition <- factor(sampleTable$condition)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = htseqDir,
                                       design=~condition) #design=~condition si j'ai plusieurs conditions

####### Vector ensembl id to gene name
id2name <- gsub(row.names(ddsHTSeq), pattern = "\\.\\d+", replacement = "")
names(id2name) <- id2name
id2name[as.character(ginf$Gene.stable.ID)] <- as.character(ginf$Gene.name)

######## Differential gene expression Analysis
dds <- DESeq(ddsHTSeq)
plotCounts(dds, gene = "ENSG00000100644.15", intgroup = "condition")

####### Extracting result
alpha = 0.05
resultsNames(dds)
res <- results(dds, name = "condition_tumor_vs_solid_tissue_normal", alpha = alpha)

######## Shrinkage
#utiliser apeglm ou ashr => biais plus faible
resShrinked <- lfcShrink(dds, coef="condition_tumor_vs_solid_tissue_normal", type="apeglm") 
plotMA(resShrinked)

######## Ordering and filtering Data
rank <- res
# Suppression du num de version à la fin des id
row.names(rank) <- sub(x = row.names(rank), pattern = "\\.\\d+", replacement = "")
row.names(rank) <- id2name[row.names(rank)]
rank <- rank[!is.na(rank$padj),] # Suppression des NA
rank <- rank[rank$padj<alpha,]
rank <- rank[order(rank$log2FoldChange, decreasing = TRUE), ] # Trie décroissant sur le "log2FoldChange"
rank <- rank["log2FoldChange"] #On conserve uniquement la colonne "log2FoldChange"

######## Values transformation
vsd <- vst(dds, blind = TRUE) # Variance Stabilizing Transformation

####### Creating expression table
exprTable <- assay(vsd)
row.names(exprTable) <- gsub(row.names(exprTable), pattern = "\\.\\d+", replacement = "")
row.names(exprTable) <- id2name[row.names(exprTable)]

######## Savind result in files
write.table(exprTable, file = paste0(resultDir, "tcga-gbm_expr.txt", sep = ""), sep="\t")
write.table(rank, file = paste0(resultDir, "tumor-vs-solid-tissue-normal_expr-change.rnk", sep = ""), sep = "\t", col.names = FALSE)
de.result <- res
row.names(de.result) <- sub(x = row.names(de.result), pattern = "\\.\\d+", replacement = "")
row.names(de.result) <- ginf[row.names(de.result), "Gene.name"]
de.result <- de.result[!is.na(de.result$padj),]
de.result <- de.result[de.result$padj<alpha,]
write.csv(de.result, file = paste0(resultDir, "tcga_gbm_vs_normal.de_analysis.csv"))

######## Performing GSEA Analysis
gmt.file <- paste0(dataDir, "gene-set/metabolism_pathways.gmt", "")
pathways <- gmtPathways(gmt.file)
rnk.vect <- setNames(rank$log2FoldChange, row.names(rank))
fgseaRes <- fgsea(pathways = pathways, stats =  rnk.vect, minSize = 0, maxSize = 500, nperm = 1000)

####### PCA
gToKeep <- row.names( res[!is.na(res$pvalue),] )
t_vsd <- as.data.frame(t(assay(vsd)))
t_vsd <- t_vsd[, gToKeep]
pca <- prcomp(t_vsd, center = TRUE, scale. = TRUE)
dt <- data.frame( PC1 = pca$x[,1], PC2= pca$x[,2], condition=vsd$condition)
ggplot(data = dt, aes(PC1, PC2, colour = condition) ) + geom_point() + ggtitle("PCA TCGA-GBM Genes")
genes.ids <- colnames(t_vsd)
genes.ids <- gsub(pattern = "\\.\\d+", replacement = "", x = genes.ids)
col.vect.metabolism <- ginf[genes.ids, -1:-4]
pc.rotation <- data.frame(PC1 = pca$rotation[,1],
                          PC2 = pca$rotation[,2], 
                          row.names = genes.ids)
pc.rotation <- cbind(pc.rotation, col.vect.metabolism)
index <- which(ginf[row.names(pc.rotation), "Cellular.response.hypoxia"])
lab <- ginf[row.names(pc.rotation), "Gene.name"]
lab[-index] <- ""
hypoxia <- as.factor(ginf[row.names(pc.rotation), "Cellular.response.hypoxia"])
ggplot(data = pc.rotation, aes(PC1, PC2, label = lab)) + 
  ggtitle("Correlation Circle") + 
  geom_point(aes(colour = hypoxia, shape = hypoxia)) + 
  geom_text(colour = "black")
