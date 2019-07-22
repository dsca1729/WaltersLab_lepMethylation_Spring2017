#dependencies
library(mclust)
#library(xlsx)
#library(stringr)

paths <- c() #add file paths here

#data input and preprocessing
read_input_data <- function()
{
  dataList <- list()
  fileNames <- list.files()
  dataList <- lapply(fileNames, read.delim)
  names(dataList) <- fileNames
  
  return(dataList)
}

filter_genelen <- function(x)
{
  x <- subset(x, x$gene_length >= 200)
}

filter_oe <- function(x)
{
  x <- subset(x, x$CpGOE <= 3.5)  
}

filter_input <- function(dataList)
{
  dataList <- lapply(dataList, na.omit)
  dataList <- lapply(dataList, filter_genelen)
  dataList <- lapply(dataList, filter_oe)
  
  return(dataList)
}

#density functions
gen_density_data <- function (dataList) {
  density_data <- list()
  for (x in dataList) {
    dens_x <- densityMclust(x$CpGOE, G=1:2) 
    density_data <- append(density_data, list(dens_x))
  }
  names(density_data) <- gen_dens_names()
  return(density_data)
}

mean_table <- function (dataList, ...) {
 names <- names(dataList)
 dim <- length(dataList)
 df.mean <- data.frame(mean1 = 1:dim, mean2 = 1:dim, row.names = names)
 for(i in 1:dim) {
   com.mean1 <- paste0("df.mean[i, 1] <- dataList$'", names[i], "'$parameters$mean[1]")
   com.mean2 <- paste0("df.mean[i, 2] <- dataList$'", names[i], "'$parameters$mean[2]")
   eval(parse(text=com.mean1))
   eval(parse(text=com.mean2))
 }
 return(df.mean)
}

variance_table <- function (dataList, ...) {
  names <- names(dataList)
  dim <- length(dataList)
  df.variance <- data.frame(var1 = 1:dim, var2 = 1:dim, row.names = names)
  for(i in 1:dim) {
    com.var1 <- paste0("df.variance[i, 1] <-dataList$'", names[i], "'$parameters$variance$sigmasq[1]")
    com.var2 <- paste0("df.variance[i, 2] <-dataList$'", names[i], "'$parameters$variance$sigmasq[2]")
    eval(parse(text=com.var1))
    eval(parse(text=com.var2))
  }
  return(df.variance)
}


gen_plots <- function(dataList, file.title, ...) {
  names <- gen_dens_names()
  command.pdf <- gen_pdf_command(file.title)
  eval(parse(text=command.pdf))
  for(i in 1:length(dataList)) {
    command.plots <- gen_plot_command(dataList, names, i)
    eval(parse(text=command.plots))
  }
  dev.off()
}
plot.dens <- function(clus, title, ...) {
  oe.range <- range(clus$data)
  oe.x <- seq(from = oe.range[1], to = oe.range[2], length.out = 500)
  
  sd1 <- sqrt(clus$parameters$variance$sigmasq[1])
  sd2 <- ifelse(clus$parameters$variance$modelName == "E", sd1, sqrt(clus$parameters$variance$sigmasq[2]) )
  
  norm1 <- dnorm(x = oe.x, mean = clus$parameters$mean[1], sd = sd1 ) * clus$parameters$pro[1]
  norm2 <- dnorm(x = oe.x, mean = clus$parameters$mean[2], sd = sd2 ) * clus$parameters$pro[2]
  
  hist(clus$data, breaks = 50, freq=F, xlab = "CpG O/E", col = "grey80", border = "grey80", main = title, ...)
  lines(x = oe.x, y = norm1, col = "blue", lwd = 2)
  lines(x = oe.x, y = norm2, col = "orange", lwd = 2)
}

#helper functions
dirTraverse <- function (pathvec, ...) {
  for (i in 1:length(pathvec)) {
    setwd(pathvec[i])
    print(pathvec[i])
  }
}

gen_dens_names <- function()
{
  fnames <- list.files()
  nameVec <- list()
  for(i in 1:length(fnames)) {
    nameVec <- append(nameVec, paste0("dclus_", fnames[i])) 
  }
  names(nameVec) <- nameVec
  return(nameVec)
}

gen_plot_command <- function(dataList,nameVec,index, ...) {
  com.name <- paste0("$'", nameVec[index], "'")
  com.call <- paste0("plot.dens(dataList", com.name,", '", nameVec[index], "'", ")")
  return(com.call)
}

gen_pdf_command <- function(file.title,...) {
  command <- paste0("pdf('",file.title,"-histograms.pdf')")
  return (command)
}

export.tables <- function(file.name, means, vars, ...) {
  file.means <- paste0(file.name, "mean-data.txt")
  file.vars <- paste0(file.name, "var-data.txt")
  write.table(x=means, file = file.means, sep = "\t")
  write.table(x= vars, file = file.vars, sep = "\t")
}



main <- function() {
  #for (i in 1:length(paths)){
  #  setwd(paths[i])
  #  l.cpg <- read_input_data()
  #  l.cpg <- filter_input(l.cpg)
  #  l.clus <- gen_density_data(l.cpg)
  #  df.mean <- mean_table(l.clus)
  #  df.var <- variance_table(l.clus)
    
  #  if(i == 1) {
  #    dataset <- c("cpgoe")
      
  #  }
  #  if(i == 2) {
  #    dataset <- c("transcriptome")
  #  }
  #  if(i == 3) {
  #    dataset <- c("kawahara")
  #  }
  #  if(i == 4) {
  #    dataset <- c("bazinet")
  #  }
  #  gen_plots(l.clus, dataset)
  #  export.tables(dataset, means = df.mean, vars = df.var)
    
  #}
  setwd(paths[1])
  l.cpg <- read_input_data()
  l.cpg <- filter_input(l.cpg)
  clus.cpg <- gen_density_data(l.cpg)
  mean.cpg <- mean_table(clus.cpg)
  var.cpg <- variance_table(clus.cpg)
  gen_plots(clus.cpg, "cpgoe-data")
  export.tables("cpgoe-data", means=mean.cpg, vars = var.cpg)
  
  setwd(paths[3])
  l.kawa <- read_input_data()
  l.kawa <- filter_input(l.kawa)
  clus.kawa <- gen_density_data(l.kawa)
  mean.kawa <- mean_table(clus.kawa)
  var.kawa <- variance_table(clus.kawa)
  gen_plots(clus.kawa, "kawahara-data")
  export.tables("kawahara-data", means = mean.kawa, vars = var.kawa)
  
  setwd(paths[4])
  l.baz <- read_input_data()
  l.baz <- filter_input(l.baz)
  clus.baz <- gen_density_data(l.baz)
  mean.baz <- mean_table(clus.baz)
  var.baz <- variance_table(clus.baz)
  gen_plots(clus.baz, "bazinet-data")
  export.tables("bazinet-data", means = mean.baz, vars = var.baz)
}