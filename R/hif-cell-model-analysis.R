library(plotly)
library(stringr)
library(data.table)
library(deSolve)

rm(list=ls())
source('~/PhD/Code/R/function/model-gene-level.R')
param.dir <- '/home/spinicck/PhD/Data/param-model/parameter/'
save.dir <- '/home/spinicck/PhD/Data/param-model/genes_level/'
li.param.dir <- '/home/spinicck/PhD/Data/param-model/li_and_wang/param_li_and_wang/'
regulator.csv <- '~/PhD/Data/param-model/parameter/regulator.csv'

genes.names <- c("Akt","AMPK","cMyc","HIF.1","mTOR","NOX","p53","PDK",
                "PI3K","PTEN","RAS","SOD","VEGF","GluT1","HK","G6PD.6PGD",
                "GPI","PFKFB2.3","PFK.1","ALD","TPI","GAPDH","PGK","PHGDH",
                "PGAM", "ENO","PKM2","PDH","ACC","LDH", "Glucose","G6P","F6P",
                "FBP", "G3P","DHAP","1,3BPG","3PG","2PG","PEP", "Pyruvate",
                "Lactate", "R5P","F2,6BP","Serine","Citrate","AMP","ADP", "ATP",
                "NAD","NADH","complex2","ROS")

all.param <- load.all.param(param.dir)

# Read Li Param CSV whithout Gene Name as Row and Column Label
# all.li.param.old <- load.all.param("PhD/Data/param-model/li_and_wang/param_li_and_wang_save/", is.matlab.tab = T)

# Add Gene Names as Row and Column Label in CSV file
# add.parameters.name(all.li.param.old, li.param.dir, genes.names)

all.li.param <- load.all.param(li.param.dir)

gene.regulator <- load.regulator(regulator.csv)

ode.dir = '~/PhD/Data/param-model/li_and_wang/result_ode_liand_wang/result_ode3_without_noise/paramset_1/'
ode.files <- list.files(ode.dir)
ode.files <- grep(".*\\.csv", ode.files, perl = TRUE, value=TRUE)
ode.files <- str_sort(ode.files, numeric = T)
all.li.result <- list()
for (i in 1:length(ode.files) ){
  all.li.result[[i]] <- fread(paste0(ode.dir, ode.files[i]), header = T, data.table = F)
  # colnames(all.li.result[[i]]) <- c("T", genes.names)
  # write.csv(all.li.result[[i]], file = paste0(ode.dir, ode.files[i]), row.names = F)
}

y0.idx <- 3
li.res.idx <- 3

index <- c(4,8,28,30)
# index <- c(1:30)

y0 <- all.li.result[[y0.idx]] [1, index]
y0 <- as.numeric(y0)
names(y0) <- genes.names[index]

gene.to.compute <- genes.names[index]
constant.gene <- genes.names[-1*index]

all.genes.level <- compute.all.genes.levels.with.solver(
  all.param, 
  gene.to.compute, 
  gene.regulator, 
  constant.gene = constant.gene,
  initial.level = y0,
  time = 7200,
  time.step = 6
)

box.li <- plot_ly(all.li.result[[li.res.idx]], type="box") %>% 
  add_trace(y = ~HIF.1, name = "HIF.1 Level") %>%
  add_trace(y = ~LDH, name = "LDH Level") %>%
  add_trace(y = ~PDK, name = "PDK Level") %>%
  add_trace(y = ~PDH, name = "PDH Level") %>% 
  add_trace(y = ~ATP, name = "ATP Level") %>% 
  add_trace(y = ~Lactate, name = "Lactate Level") %>%
  layout(title = paste0("Box Plot Li and Wang Results ", li.res.idx), yaxis=list(title="Level"))
box.li

plot.li <- plot_ly(all.li.result[[li.res.idx]], x= ~T, y = ~HIF.1, name = "HIF.1 Level", type="scatter", mode = "lines") %>% 
  add_trace(y = ~LDH, name = "LDH Level") %>%
  add_trace(y = ~PDK, name = "PDK Level") %>%
  add_trace(y = ~PDH, name = "PDH Level") %>% 
  add_trace(y = ~ATP, name = "ATP Level") %>% 
  add_trace(y = ~Lactate, name = "Lactate Level") %>% 
  layout(title = paste0("Line Plot Li and Wang Results ", li.res.idx), yaxis=list(title="Level"))
plot.li


box.my.model <- plot_ly(all.genes.level[[1]], type="box") %>%
  add_trace(y = ~HIF.1, name = "HIF.1 Level") %>%
  add_trace(y = ~LDH, name = "LDH Level") %>%
  add_trace(y = ~PDK, name = "PDK Level") %>%
  add_trace(y = ~PDH, name = "PDH Level") %>%
  layout(title = paste0("Box Plot Model Results ", y0.idx))
box.my.model

plot.my.model <- plot_ly(all.genes.level[[1]], x= ~T, y = ~HIF.1, name = "HIF.1 Level", type="scatter", mode = "lines") %>% 
  add_trace(y = ~LDH, name = "LDH Level") %>%
  add_trace(y = ~PDK, name = "PDK Level") %>%
  add_trace(y = ~PDH, name = "PDH Level") %>% 
  layout(title = paste0("Line Plot My Model Results ", y0.idx), yaxis=list(title="Level"))
plot.my.model


one.g.level <- data.frame(row.names=c(1:1201))
for ( i in 1:10){
  one.g.level[1:1201, i] <- all.li.result[[i]]$HIF.1  
}
boc.one.g <- plot_ly(one.g.level, type="box") 
for ( i in 1:10 ){
  boc.one.g <-  boc.one.g %>% add_trace(y = one.g.level[1:1201, i], name = i) 
}
boc.one.g <-  boc.one.g %>% layout(title = paste0("One Gene Level ", i), yaxis=list(title="Level"))
boc.one.g
