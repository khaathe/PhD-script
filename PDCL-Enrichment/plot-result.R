######## Load Library needed
library(ggplot2)

collapse.data <- function(gost.penda.gem, gsea.deseq2.gem, pdcl.samples.name, db.source){
  collapsed <- data.frame()
  for (pdcl in pdcl.samples.name){
    for (db in db.source) {
      gost <- gost.penda.gem[[pdcl]][[db]]
      gsea <- gsea.deseq2.gem[[pdcl]][[db]]
      join <- gost %>% 
        full_join(gsea, by = "pathway_id") %>%
        rename(p_value_gprofiler = p_value.x, p_value_gsea = p_value.y, fdr_gprofiler = fdr.x, fdr_gsea = fdr.y) %>%
        select(pathway_id, p_value_gprofiler, fdr_gprofiler, p_value_gsea, fdr_gsea) %>%
        mutate(pdcl = pdcl, db = db) 
      collapsed <- rbind(collapsed, join)
    }
  }
  collapsed
}

collapse.data.by.method <- function(gost.penda.gem, gsea.deseq2.gem, pdcl.samples.name, db.source){
  collapsed <- data.frame()
  for (pdcl in pdcl.samples.name){
    for (db in db.source) {
      gost <- gost.penda.gem[[pdcl]][[db]]
      gost$method <- "gprofiler"
      gsea <- gsea.deseq2.gem[[pdcl]][[db]]
      gsea$method <- "gsea"
      join <- bind_rows(gost, gsea) %>%
        mutate(pdcl = pdcl, db = db) 
      collapsed <- rbind(collapsed, join)
    }
  }
  collapsed
}

plot.gprofiler.vs.gsea <- function(x, threshold = 0.1){
  data <- na.omit(x)
  data$is.common <- ( data$fdr_gprofiler<threshold & data$fdr_gsea<threshold )
  ggplot(
    data = data, 
    aes(
      x = -log10(fdr_gprofiler), 
      y = -log10(fdr_gsea),
      colour = is.common
    )
  ) +
    geom_point(
      size = 0.5
    ) +
    labs(
      title = "pvalue G:Profiler against GSEA",
      subtitle = paste0("Threshold alpha : ", threshold),
      x = "-log10(pvalue G:Profiler)",
      y = "-log10(pvalue GSEA)",
      colour = "Common"
    )
}

plot.count.pathways <- function(x, threshold = 0.1){
  data <- x
  data <- na.omit(data)
  data$is.common <- ( data$fdr_gprofiler<threshold & data$fdr_gsea<threshold )
  gprofiler <- data[(data$fdr_gprofiler < threshold),]
  gprofiler$method <- "gprofiler"
  gsea <- data[(data$fdr_gsea < threshold),]
  gsea$method <- "gsea"
  data <- rbind(gprofiler, gsea)
  ggplot(
    data = data, 
    aes(x = pdcl, fill = is.common )
  ) +
    geom_bar() +
    labs(
      title = "Number of enriched pathways for each enrichment method",
      subtitle = paste0("Threshold alpha : ", threshold),
      x = "PDCL",
      y = "Count",
      fill = "Common"
    ) +
    facet_wrap(vars(method)) +
    theme(axis.text.x = element_text(angle = 90) ) 
}

