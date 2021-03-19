create.deseq2.dataset.from.htseq.count = function (htseq.file, htseq.dir, ...){
  data.htseq <- read.table(file = htseq.file, header = TRUE, sep="\t")
  sample.table <- data.frame(sampleName = data.htseq$sampleName, 
                            fileName = data.htseq$fileName, 
                            condition = data.htseq$condition)
  sample.table$condition <- factor(sample.table$condition)
  
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sample.table,
                                         directory = htseq.dir,
                                         design=...)
  return(ddsHTSeq)
}

create.deseq2.dataset.from.files <- function(file.list, condition, ...){
  cts <- create.expression.data.frame(file.list)
  coldata <- data.frame(condition = condition)
  coldata$condition <- factor(coldata$condition)
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ...)
  return(dds)
}

create.expression.data.frame <- function(file.list){
  expression.data <- read.expression.file(file_path = file.list[1])
  for (f in file.list[-1]){
    to.merge <- read.expression.file(f)
    if ( sum(expression.data$gene %in% to.merge$gene) != length(expression.data$gene)){
      msg <- paste0("the file ", f, " has different rows")
      warning(msg)
    }
    newCol <- colnames(to.merge)[1]
    expression.data <- cbind(expression.data, to.merge[,2])
  }
  expression.data <- as.data.frame(expression.data)
  row.names(expression.data) <- expression.data$gene
  expression.data <- expression.data[,-1]
  return(expression.data)
}

read.expression.file <- function(file_path){
  col.name <- gsub( x = file_path, pattern = ".*/", replacement = "")
  expr.dt <- fread(file = file_path, header = FALSE)
  colnames(expr.dt) <- c("gene", col.name)
  return(expr.dt)
}