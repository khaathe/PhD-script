plot.pca <- function(pca.x, col.factor, title = "PCA", legend.title = "condition", pc.num = c(1,2)){
  pca.dt <- data.frame( PC1 = pca.x[,pc.num[1]], PC2= pca.x[,pc.num[2]], color=col.factor)
  axes.label <- paste0("PC", pc.num, "")
  ggplot(data = pca.dt, aes(PC1, PC2, colour = color) ) + 
    geom_point() + 
    ggtitle(title) +
    labs(colour = legend.title)+
    xlab(axes.label[1]) +
    ylab(axes.label[2])
}

plot.pc.contribution <- function(sdev, title = "Contribution of each Principal Components"){
  var.percentage.dt <- data.frame(pc.num = c(1:length(sdev)), percentage = sdev/sum(sdev))
  ggplot(data = var.percentage.dt, aes(pc.num, percentage)) +
    geom_col() + 
    labs(title = title, x = "Principal Component Number", y = "Percentage of Variance Explained")
}

plot.sum.pc.contribution <- function(sdev, title = "Sum of the Contribution Principal Components"){
  var.percent <- sdev/sum(sdev)
  percent.sum <- 0
  var.sum.percent <-c()
  for (percent in var.percent) {
    percent.sum <- percent.sum + percent
    var.sum.percent <- c(var.sum.percent, percent.sum)
  }
  var.explained.dt <- data.frame(pc.num = c(1:length(var.percent)), sum.percentages = var.sum.percent)
  ggplot(data = var.explained.dt, aes(pc.num, sum.percentages)) +
    geom_col() + 
    labs(title = title, x = "Principal Component Number", y = "Percentage of Variance Explained") +
    geom_abline(slope = 0.0, intercept = 0.80, colour = "red")
}

plot.correlation.circle <- function(pca.rotation, col.factor=NULL, lab="", 
                                    title = "Correlation Circle", legend.title = "condition", 
                                    check.overlap = F, show.legend = NA, v.adjust = 0, 
                                    pc.num = c(1,2)){
  pc.rotation.dt <- data.frame(PC1 = pca.rotation[,pc.num[1]], PC2 = pca.rotation[,pc.num[2]] )
  axes.label <- paste0("PC", pc.num, "")
  ggplot(data = pc.rotation.dt, aes(PC1, PC2, label = lab)) + 
    ggtitle(title) + 
    geom_point(aes(colour = col.factor), show.legend = show.legend) + 
    geom_text(colour = "black", check_overlap = check.overlap, vjust=0, nudge_y = v.adjust) +
    labs(colour = legend.title) +
    xlab(axes.label[1]) +
    ylab(axes.label[2])
}

ggplot.pc.density <- function(pca.rotation, num.bin = 40, pc.num = c(1,2)){
  hbin <- data.frame(PC1 = pca.rotation[,pc.num[1]], PC2 = pca.rotation[,pc.num[2]])
  axes.label <- paste0("PC", pc.num, "")
  ggplot(hbin, aes(PC1, PC2)) +
    ggtitle(paste0("Principal Components Coefficient Density (xbin=", num.bin, ")", collapse = "")) +
    geom_hex(bins = num.bin) + 
    scale_fill_distiller(palette = "Spectral")
    xlab(axes.label[1]) +
    ylab(axes.label[2])
}

plot.pc.density <- function(pca.rotation, num.bin = 40, pc.num = c(1,2)){
  hbin <- hexbin(pca.rotation[,pc.num[1]], pca.rotation[,pc.num[2]], xbins = num.bin)
  spectral.col <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
  plot(hbin, colramp = spectral.col, 
       xlab = paste0("PC", pc.num[1], ""), ylab = paste0("PC", pc.num[2], ""), 
       main=paste0("Principal Components Coefficient Density (xbin=", num.bin, ")", collapse = ""))
}