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
library(hexbin)
library(RColorBrewer)

######## Load Session
load("/home/spinicck/Bureau/EnrichmentAnalysis/Result/tcga-gbm/tcga-gbm-session.RData")

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
gene.info <- read.csv(file = paste0(dataDir, "TCGA-GBM/gene-info.csv", ""), header = T, colClasses = col.class, row.names = 1)
rm(col.class)

###### Vector for converting gene id to gene name
id.to.name <- gene.info$Gene.name
names(id.to.name) <- gene.info$Gene.stable.ID
id.to.name[is.na(id.to.name)] <- names(id.to.name[is.na(id.to.name)])

######## Creating DataSet
ddsHTSeq <- create.deseq2.dataset.from.htseq.count(htseqFile, htseqDir, design = ~condition)

######## Differential gene expression Analysis
dds <- DESeq(ddsHTSeq)
row.names(dds) <- sub(x = row.names(dds), pattern = "\\.\\d+", replacement = "")

####### Extracting result
alpha = 0.05
resultsNames(dds)
res <- results(dds, name = "condition_tumor_vs_solid_tissue_normal", alpha = alpha)

######## Shrinkage
#utiliser apeglm ou ashr => biais plus faible
resShrinked <- lfcShrink(dds, coef="condition_tumor_vs_solid_tissue_normal", type="apeglm") 
plotMA(resShrinked)

######## Creating Differential Expression table
de.genes <- res
de.genes$Gene.Name <- id.to.name[row.names(de.genes)]
de.genes <- de.genes[!is.na(de.genes$padj),]
de.genes <- de.genes[de.genes$padj<alpha,]

######## Creating rank table
rank <- rank.table(de.result, alpha = alpha)
row.names(rank) <- id.to.name[row.names(rank)]

######## Values transformation
d.norm <- vst(dds, blind = TRUE) # Variance Stabilizing Transformation

####### Creating expression table
exprTable <- assay(d.norm)
row.names(exprTable) <- id.to.name[row.names(exprTable)]

####### PCA
keep.genes <- row.names(de.genes)
d.norm.t <- as.data.frame(t(assay(d.norm)))
d.norm.t <- d.norm.t[, keep.genes]
#Perfom PCA
pca <- prcomp(d.norm.t, center = T, scale. = T)

#Plot PCA
plot.pca( pca.x = pca$x, col.factor = d.norm$condition, title = "PCA TCGA-GBM Differentially Expressed Genes")
plot.pc.contribution(sdev = pca$sdev)
plot.sum.pc.contribution(sdev = pca$sdev)

#Plot Coefficient Density
plot.pc.density(pca.rotation = pca$rotation, num.bin = 40)

#Plot principal component coefs for hypoxia
hypoxia.factor <- as.factor(gene.info[keep.genes, "Cellular.response.hypoxia"])
hypoxia.label <- gene.info[keep.genes, "Gene.name"]
hypoxia.label[!gene.info[keep.genes, "Cellular.response.hypoxia"]] <- ""
plot.correlation.circle(pca.rotation = pca$rotation,
                        col.factor = hypoxia.factor,
                        lab = hypoxia.label,
                        title = "Coef Principal Component PCA DE Genes (hypoxia)",
                        legend.title = "Involved in Cellular Response to Hypoxia")

#Plot principal component coefs for The Citric Acid Cycle
tca.factor <- as.factor(gene.info[keep.genes, "Tca.cycle.and.respiratory.electron.transport"])
tca.label <- gene.info[keep.genes, "Gene.name"]
tca.label[!gene.info[keep.genes, "Tca.cycle.and.respiratory.electron.transport"]] <- ""
plot.correlation.circle(pca.rotation = pca$rotation, 
                        col.factor = tca.factor, 
                        lab = tca.label, 
                        title = "Coef Principal Component PCA DE Genes (TCA Cycle)",
                        legend.title = "Involved in The Citric Acid Cycle")

#Plot principal component coefs for Glycolysis
glycolysis.factor <- as.factor(gene.info[keep.genes, "Glycolysis"])
glycolysis.label <- gene.info[keep.genes, "Gene.name"]
glycolysis.label[!gene.info[keep.genes, "Glycolysis"]] <- ""
plot.correlation.circle(pca.rotation = pca$rotation, 
                        col.factor = glycolysis.factor, 
                        lab = glycolysis.label, 
                        title = "Coef Principal Component PCA DE Genes (Glycolysis)",
                        legend.title = "Involved in Glycolysis")

######## Savind result in files
write.table(exprTable, file = paste0(resultDir, "tcga-gbm_expr.txt", sep = ""), sep="\t")
write.table(rank, file = paste0(resultDir, "tumor-vs-solid-tissue-normal_expr-change.rnk", sep = ""), sep = "\t", col.names = FALSE)
write.csv(de.genes, file = paste0(resultDir, "tcga_gbm_vs_normal.de_analysis.csv"))

######## Performing GSEA Analysis
gmt.file <- paste0(dataDir, "gene-set/metabolism_pathways.gmt", "")
pathways <- gmtPathways(gmt.file)
rnk.vect <- setNames(rank$log2FoldChange, row.names(rank))
fgseaRes <- fgsea(pathways = pathways, stats =  rnk.vect, minSize = 0, maxSize = 500, nperm = 1000)