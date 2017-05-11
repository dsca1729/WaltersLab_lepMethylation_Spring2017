#Script assumes working directory (wd) is contains the entirety of data-sets of interest
fileNames <- list.files()

#orgNames <- c("A_Transitella", "B_Anynana1", "B_Anynana2", "Bombyx_Mori1", "Bombyx_Mori ASM", "C._Cecrops", "C._Supressalis", "D_Plexippus1", "D. Plexippus3", "H. Erato", "H. Melpomene1", "H. Melpomene2", "L. Accius", "M. Sexta", "M. Cinxia", "O. Brumata", "P. Glaucus", "P. Machaon", "P. Polytes", "P. Xuthus Pap", "P. Xuthus Pxut", "P. Sennae", "P. Interpunctella", "P. Xylostella") 

library(mclust)

#generate and plot individual density clusters
fit_and_plot <- function(x)
{
  #z <- densityMclust(x$CpGOE, G=1:2)
  plot.dens(clus = z)
}

filter_geneLength <- function(x)
{
  x <- subset(x, x$gene_length >= 200)
}





#Density Plotting Function
plot.dens <- function(clus, ...) {
  oe.range <- range(clus$data)
  oe.x <- seq(from = oe.range[1], to = oe.range[2], length.out = 500)
  
  sd1 <- sqrt(clus$parameters$variance$sigmasq[1])
  sd2 <- ifelse(clus$parameters$variance$modelName == "E", sd1, sqrt(clus$parameters$variance$sigmasq[2]) )
  
  norm1 <- dnorm(x = oe.x, mean = clus$parameters$mean[1], sd = sd1 ) * clus$parameters$pro[1]
  norm2 <- dnorm(x = oe.x, mean = clus$parameters$mean[2], sd = sd2 ) * clus$parameters$pro[2]
  
  hist(clus$data, breaks = 50, freq=F, xlab = "CpG O/E", col = "grey80", border = "grey80", ...)
  lines(x = oe.x, y = norm1, col = "blue", lwd = 2)
  lines(x = oe.x, y = norm2, col = "orange", lwd = 2)
  
}

#read in all files for usage.
dataList <- lapply(fileNames, read.delim)
names(dataList) <- fileNames
dataList <- lapply(dataList, na.omit)
dataList <- lapply(dataList, filter_geneLength)


#generating density clusters
#pdf("histograms_kawahara.pdf")
clustData <- lapply(dataList, fit_and_plot)
#dev.off()


