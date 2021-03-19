########## Init working dir
wDir = "/home/spinicck/Bureau/EnrichmentAnalysis/R"
setwd(wDir)

########## Construction table metadonnées
normal_dir_path <- "/home/spinicck/Bureau/EnrichmentAnalysis/Data/TCGA-GBM/solid-tissue-normal/htseq-counts/"
tumor_dir_path <- "/home/spinicck/Bureau/EnrichmentAnalysis/Data/TCGA-GBM/tumor/htseq-count/"
#normal_dir_path <- "/home/spinicck/Bureau/EnrichmentAnalysis/Data/TCGA-BRCA/solid-tissue-normal/"
#tumor_dir_path <- "/home/spinicck/Bureau/EnrichmentAnalysis/Data/TCGA-BRCA/primary-tumor/"
normal_sample <- list.files(normal_dir_path)
tumor_sample <- list.files(tumor_dir_path)
metadata <- data.frame( file = c(normal_sample, tumor_sample))
metadata$dir <- c( rep(normal_dir_path, times= length(normal_sample)), rep(tumor_dir_path, times = length(tumor_sample)) )
metadata$type <- c( rep("normal", times= length(normal_sample)), rep("tumor", times = length(tumor_sample)) )

########## Writting expression values in csv
source('~/Bureau/EnrichmentAnalysis/Script/R/create-expression-data-frame.R')
expr_file_list <- paste0(metadata$dir, metadata$file, sep = "")
expr_dt <- create.expression.data.frame(expr_file_list) 
# /!\ les fichiers htseq-count possède des lignes supplémentaires à la fin du fichier
write.csv(expr_data_frame, file = "tcga-gbm_expression_values.csv")

######## PCA
t_expr_dt <- as.data.frame(t(expr_dt))
expr_pca <- prcomp(t_expr_dt)
plot(expr_pca$x[,1], expr_pca$x[,2], xlab = "PC1", ylab = "PC2", main = "PCA Expression", col=as.factor(metadata[metadata$file[row.names(expr_pca$x)], "type"]), pch=1)
legend("topright", legend = unique(metadata$type), col = unique(as.factor(metadata$type)), pch=1)
