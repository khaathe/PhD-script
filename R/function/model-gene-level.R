compute.hif.level <- function(pO2min = 0.0, pO2max = 0.21, by = 0.001 ){
  hif.levels.vect <- c()
  for (o in seq(pO2min, pO2max, by) ) {
    h <- 0.0085*exp( 2.5*(1-o))
    hif.levels.vect <- c(hif.levels.vect, h)
  }
  hif.levels.vect <- sort(hif.levels.vect)
  return(hif.levels.vect)
}

compute.all.genes.levels <- function(param, gene.names, gene.regulators, 
                                     time = 10, time.step = 0.1,
                                     paramset.num=1, save.dir=NULL, 
                                     initial.level = 0.0, constant.gene=NULL){
  all.genes.levels <- list()
  t.vect <- seq(from = 0.0, to = time, by = time.step)
  N = length(t.vect)
  for ( i in paramset.num ){
    levels.computed <- data.frame(T = t.vect)
    levels.computed[1, gene.names] = initial.level
    for ( j in 2:N ){
      for (gene in gene.names){
        regulator.names <-gene.regulators[[gene]]
        variables.regulators <- regulator.names[!(regulator.names %in% constant.gene)]
        precedent.gene.levels <- as.numeric(levels.computed[j-1, c(gene, variables.regulators)])
        names(precedent.gene.levels) <- c(gene, variables.regulators)
        h =  compute.h.value( precedent.gene.levels, gene, regulator.names, param, i, constant.gene)
        A = param[['A']][[i]] [gene,1]
        D = param[['A']][[i]] [gene,1]
        X = levels.computed[j-1, gene]
        levels.computed[j, gene] = X + time.step * (A * h - D * X)
      }
    }
    all.genes.levels[[i]] <- levels.computed
    if ( !is.null(save.dir) ){
      f <- paste0(save.dir,'result_paramset_', i ,'.csv')
      write.csv(genes.level, f, row.names = F)
    }
  }
  return(all.genes.levels)
}

compute.all.genes.levels.with.solver <- function(param, gene.names, gene.regulators, 
                                     initial.level,
                                     time = 10, time.step = 0.1,
                                     paramset.num=1, save.dir=NULL, 
                                     constant.gene=NULL){
  all.genes.levels <- list()
  t.vect <- seq(from = 0.0, to = time, by = time.step)
  N = length(t.vect)
  for ( i in paramset.num ){
    parameters = list(param = param, 
                      paramset = i,
                      gene.names = gene.names,
                      gene.regulator = gene.regulator,
                      constant.gene = constant.gene
    )
    desolve.out <- ode(y = initial.level, times = t.vect, func = model.odes, parms = parameters, method = 'euler')
    all.genes.levels[[i]] <- as.data.frame(desolve.out)
    colnames(all.genes.levels[[i]])[1] <- "T"
    if ( !is.null(save.dir) ){
      f <- paste0(save.dir,'result_paramset_', i ,'.csv')
      write.csv(genes.level, f, row.names = F)
    }
  }
  return(all.genes.levels)
}

model.odes <- function(t, y0, parameters) {
  dX <- list()
  p <- parameters[['param']]
  paramset <- parameters[['paramset']]
  gene.names <- parameters[['gene.names']]
  gene.regulators <- parameters[['gene.regulator']]
  constant.gene <- parameters[['constant.gene']]
  for (i in 1:length(y0)) {
    gene <- names(y0)[i]
    A = p [['A']] [[paramset]] [gene,1]
    D = p [['D']] [[paramset]] [gene,1]
    regulator.names <- gene.regulators[[gene]]
    variables.regulators <- regulator.names[!(regulator.names %in% constant.gene)]
    precedent.gene.levels <- y0[c(gene, variables.regulators)]
    h =  compute.h.value( precedent.gene.levels, gene, regulator.names, p, paramset, constant.gene)
    dX[i] <- A * h - D* y0[i]
  }
  list(as.numeric(dX))
}



compute.h.value = function(gene.levels, gene, regulator.names, param, paramset.num, 
                           constant.gene.names=NULL) {
  h = 1.0
  variable.genes <- regulator.names[!(regulator.names %in% constant.gene.names)]
  constant.genes <- regulator.names[regulator.names %in% constant.gene.names]
  for ( r in variable.genes ){
    y <- as.numeric(gene.levels[r])
    n <- param [["n"]] [[paramset.num]] [r, gene]
    s <- param [["S"]] [[paramset.num]] [r, gene]
    gamma <- param [["gamma"]] [[paramset.num]] [r, gene]
    h = h * (s^n/(s^n+y^n))+gamma*(y^n/(s^n+y^n))
  }
  for ( r in constant.genes ){
    h <- h * param [["gamma"]] [[paramset.num]] [r, gene]
  }
  return(h)
}

load.all.param <- function(param.dir, is.matlab.tab =F){
  param.name <- list.dirs(param.dir, full.names = F)
  param.name <- param.name[param.name != ""]
  all.param <- list()
  for (p in param.name) {
    dir <- paste0(param.dir, p)
    files <- paste0(dir, "/", list.files(dir))
    for ( i in 1:length(files) ){
      if ( is.matlab.tab ){
        param <- read.table(files[i], sep = " ")
      } else {
        param <- read.csv(files[i], row.names = 1) 
      }
      all.param[[p]] [[i]] <- as.data.frame(param)
    }
  }
  return(all.param)
}

load.regulator <- function(regular.file){
  regulator.tab <- read.table(regular.file, sep = ",", fill = T, colClasses = 'character')
  regulator.list <- list()
  for (i in 1:nrow(regulator.tab) ){
    # Loop over all the genes in the csv from the second column and return a vector without names
    reg <- sapply(regulator.tab[i,2:ncol(regulator.tab)], function(x){ return(x) }, USE.NAMES = F) 
    reg <- reg[ reg != ""]
    regulator.list[[ regulator.tab[i,1] ]] <- reg
  }
  return(regulator.list)
}

add.parameters.name <- function (all.li.param,li.param.dir, genes.names) {
  #Write HILL's parameter
  for (var in c("A", "D", "S", "n", "gamma", "lx") ) {
    for (i in 1:length(all.li.param[[var]]) ){
      file <- paste0(li.param.dir, var, "/", var ,"_", i, "_param.csv")
      df <- all.li.param [[var]] [[i]]
      n <- nrow(df)
      colname.vect <- c()
      if (nrow(df) == 53){
        df$V1 <- genes.names
      }
      if ( ncol(df ) == 2 ){
        colnames(df) <- c("X", var)
      } else if ( ncol(df) == 54 ){
        colnames(df) <- c("X", genes.names)
      }
      write.csv(x = df, file = file, row.names = FALSE)
    }
  }
  
  #Write kinetic parameter
  file <- paste0(li.param.dir, "/kinetic_param.csv")
  all.var.idx <- names(all.li.param) %in% c("A", "D", "S", "n", "gamma", "lx")
  all.var <- names(all.li.param)[!all.var.idx]
  param.df <- data.frame(row.names = c(1:18) )
  param.df[1:18, all.var] <- NA
  for (var in all.var ) {
    for (i in 1:length(all.li.param[[var]]) ){
      df <- all.li.param [[var]] [[i]]
      param.df[i, var] <- df[1,2]
    }
  }
  write.csv(x = param.df, file = file)
}