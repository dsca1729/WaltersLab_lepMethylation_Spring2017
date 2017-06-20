
#setwd("C:/Documents/GitHub/WaltersLab_lepMethylation_Spring2017/input_files/CpGOE_data")
#setwd("C:/Documents/GitHub/WaltersLab_lepMethylation_Spring2017/input_files/transcriptome_data")
#setwd("C:/Documents/GitHub/WaltersLab_lepMethylation_Spring2017/input_files/kawahara_data")
#setwd("C:/Documents/GitHub/WaltersLab_lepMethylation_Spring2017/input_files/bazinet_data")

fileNames <- list.files()

kawahara_Names <- c("Pterodecta felderi", "Notoplusia minuta", "Therinia lactucina", "Morpheis mathani", "Megalopyge tharops", "Aepytus sp.", "Dalcera abrasa", "Lacosoma ludolpha", "Adhemarius daphne", "Nothus lunus", "Artace sp.", "Semomesia campanea", "Zeuzerodes maculata", "Eacles sp.", "Macrosoma sp.", "Hylephila phyleus", "Hemiargus ceraunus", "Lantanophaga pusillidactyla", "Macaria distribuaria", "Nemoria lixaria", "Urodus parvula", "Megathymus yuccae", "Pseudothyris sepulchralis", "Phyllocnistis citrella", "Cnaphalocrocis medinalis", "Grapholita dimorpha", "Papilio glaucus", "Papilio polytes", "Thubana sp.", "Anigraea sp.", "Lyssa zampa", "Dichomeris sp.", "Nemophora sp.", "Manoba major", "Eudocima salaminia", "Xylophanes tersa")
refGenome_Names <- c("Amyelois transitella", "Bicyclus anynana nBA", "Bicyclus anynana v1.2", "Bombyx mori", "Bombyx mori ASM", "Calycopis cecrops", "Chilo suppressalis", "Danaus plexippus", "Danaus plexippus v3", "Heliconius erato", "Heliconius melpomene", "Heliconius melpomene 2", "Lerema accius", "Manduca sexta", "Melitaea cinxia", "Operophtera brumata", "Papilio glaucus", "Papilio machaon", "Papilio polytes", "Papilio xuthus 1", "Papilio xuthus 2", "Phoebis sennae", "Plodia interpunctella", "Plutella xylostella")
bazinet_Names <- c("Dryadaula visaliella", "Palaephatus luteolus", "Phymatopus californicus", "Thyridopteryx ephemeraeformis", "Tineola bisselliella", "Tischeria quercitella", "Acanthopteroctetes unifascia", "Eudarcia simulatricella", "Ptyssoptera sp. AB-2015", "Acanthopteroctetes unifascia", "Acanthopteroctetes unifascia", "Agathiphaga queenslandensis", "Agathiphaga queenslandensis", "Agathiphaga queenslandensis", "Andesiana lamellata", "Azaleodes micronipha", "Dyseriocrania griseocapitella", "Dryaduala visaliella", "Enteucha acetosae", "Eudarcia simulatricella", "Heterobathmia pseuderiocrania", "Lophocorona astiptica", "Metaphatus ochraceus","Metaphatus ochraceus","Metaphatus ochraceus", "Palaephatus luteolus", "Phymatopus californicus", "Palaephatus nielseni", "Thyridopterix ephemeraeformis", "Tineola bisselliella", "Neopseustis meyricki", "Pseudopostega quadristigella", "Tischeria quercitella", "Tegeticula yuccasella", "Ptyssoptera sp. AB-2015")
trans_Names <- c("Chilo suppressalis", "Heliconius melpomene", "Melitaea cinxia", "Manduca sexta", "Plutella Xylostella")

library(mclust)

#generate and plot individual density clusters
fit_and_plot <- function(x)
{
  z <- densityMclust(x$CpGOE, G=1:2)
  plot.dens(z)
}

filter_geneLength <- function(x)
{
  x <- subset(x, x$gene_length >= 100)
}

filter_OE <- function(x)
{
  x <- subset(x, x$CpGOE <= 3.5)  
}

tab_mean <- function(x)
{
  p <- clustData2[[x]]$parameters$mean
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
dataList2 <- lapply(fileNames, read.delim)
names(dataList2) <- fileNames
dataList2 <- lapply(dataList2, na.omit)
dataList2 <- lapply(dataList2, filter_geneLength)
dataList2 <- lapply(dataList2, filter_OE)


#generating density clusters
pdf("histograms_transcriptomes.pdf")
lapply(dataList2, fit_and_plot)
dev.off()

#generating mean tables
trans_mean <- lapply(1:5, tab_mean)
kawahara_mean <- lapply(1:, tab_mean)
bazinet_mean <- lapply(1:35, tab_mean)
