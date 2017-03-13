library(mclust)

#Functions
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
Atrans <- read.delim(file="Amyelois_transitella_v1_-_cds.fa.out.txt", as.is=TRUE)
Banynana1 <- read.delim(file="Bicyclus_anynana_nBa.0.1_-_cds.fa.out.txt", as.is=TRUE)
Banynana2 <- read.delim(file="Bicyclus_anynana_v1.2_-_cds.fa.out.txt" , as.is=TRUE)
Bmori1 <- read.delim(file="Bombyx_mori_-_cds.fa.out.txt", as.is=TRUE)
BmoriAsm <- read.delim(file="Bombyx_mori_ASM15162v1_-_cds.fa.out.txt", as.is=TRUE) 
Ccecro <- read.delim(file="Calycopis_cecrops_v1.1_-_cds.fa.out.txt", as.is=TRUE)
Csupp <- read.delim(file="Chilo_suppressalis_CsuOGS1.0_-_cds.fa.out.txt", as.is=TRUE)
Dplex1 <- read.delim(file="Danaus_plexippus_-_cds.fa.out.txt", as.is=TRUE)
Dplex3 <- read.delim(file="Danaus_plexippus_v3_-_cds.fa.out.txt", as.is=TRUE)
Herato <- read.delim(file="Heliconius_erato_v1_-_cds.fa.out.txt", as.is=TRUE)
Hmelpo1 <-read.delim(file="Heliconius_melpomene_-_cds.fa.out.txt", as.is=TRUE)
Hmelpo2 <-read.delim(file="Heliconius_melpomene_Hmel2_-_cds.fa.out.txt", as.is=TRUE)
Laccius<- read.delim(file="Lerema_accius_v1.1_-_cds.fa.out.txt", as.is=TRUE) 
Msex <- read.delim(file="Manduca_sexta_Msex_1.0_-_cds.fa.out.txt", as.is=TRUE) 
Mcinx <-read.delim(file="Melitaea_cinxia_-_cds.fa.out.txt", as.is=TRUE)
Obrum <- read.delim(file="Operophtera_brumata_v1_-_cds.fa.out.txt", as.is=TRUE)
Pglauc <- read.delim(file="Papilio_glaucus_v1.1_-_cds.fa.out.txt", as.is=TRUE)
Pmach <- read.delim(file="Papilio_machaon_Pap_ma_1.0_-_cds.fa.out.txt", as.is=TRUE)
Ppoly <- read.delim(file="Papilio_polytes_Ppol_1.0_-_cds.fa.out.txt", as.is=TRUE)
Pxuth1<- read.delim(file="Papilio_xuthus_Pap_xu_1.0_-_cds.fa.out.txt", as.is=TRUE) 
Pxuth2<- read.delim(file="Papilio_xuthus_Pxut_1.0_-_cds.fa.out.txt", as.is=TRUE)
Psennae<- read.delim(file="Phoebis_sennae_v1.1_-_cds.fa.out.txt", as.is=TRUE) 
Pinter<- read.delim(file="Plodia_interpunctella_v1_-_cds.fa.out.txt", as.is=TRUE)
Pxylo<- read.delim(file="Plutella_xylostella_DBM_FJ_v1.1_-_cds.fa.out.txt", as.is=TRUE)

#Removing NA's
Banynana1 <- na.omit(Banynana1)
Ccecro <- na.omit(Ccecro)
Dplex3 <- na.omit(Dplex3)
Hmelpo1 <- na.omit(Hmelpo1)
Mcinx <- na.omit(Mcinx)
Msex <-na.omit(Msex)
Pinter <- na.omit(Pinter)
Psennae <- na.omit(Psennae)
Pxuth1 <- na.omit(Pxuth1)

#Generating Density Clusters
components <- 1:2

Atrans.clust <- densityMclust(Atrans$CpGOE, G=components)
Banynana1.clust <- densityMclust(Banynana1$CpGOE, G=components)
Banynana2.clust <- densityMclust(Banynana2$CpGOE, G=components)
Bmori1.clust <- densityMclust(Bmori1$CpGOE, G=components)
BmoriAsm.clust <- densityMclust(BmoriAsm$CpGOE, G=components)
Ccecro.clust <- densityMclust(Ccecro$CpGOE, G= components)
Csupp.clust <- densityMclust(Csupp$CpGOE, G= components)
Dplex1.clust <- densityMclust(Dplex1$CpGOE, G= components)
Dplex3.clust <- densityMclust(Dplex3$CpGOE, G= components)
Herato.clust <- densityMclust(Herato$CpGOE, G= components)
Hmelpo1.clust <- densityMclust(Hmelpo1$CpGOE, G= components)
Hmelpo2.clust <- densityMclust(Hmelpo2$CpGOE, G= components)
Laccius.clust <- densityMclust(Laccius$CpGOE, G= components)
Msex.clust <- densityMclust(Msex$CpGOE, G= components)
Mcinx.clust <- densityMclust(Mcinx$CpGOE, G= components)
Obrum.clust <- densityMclust(Obrum$CpGOE, G= components)
Pglauc.clust <- densityMclust(Pglauc$CpGOE, G= components)
Pinter.clust <- densityMclust(Pinter$CpGOE, G= components)
Pmach.clust <- densityMclust(Pmach$CpGOE, G= components)
Ppoly.clust <- densityMclust(Ppoly$CpGOE, G= components)
Psennae.clust <- densityMclust(Psennae$CpGOE, G= components)
Pxuth1.clust <- densityMclust(Pxuth1$CpGOE, G= components)
Pxuth2.clust <- densityMclust(Pxuth2$CpGOE, G= components)
Pinter.clust <- densityMclust(Pinter$CpGOE, G= components)
Pxylo.clust <- densityMclust(Pxylo$CpGOE, G= components)


#Generating Plots

plot.dens(clus = Atrans.clust, main= "A. Transitella")
plot.dens(clus = Banynana1.clust, main= "B. Anynana_nBa")
plot.dens(clus = Banynana2.clust, main= "B. Anynana")
plot.dens(clus = Bmori1.clust, main= "B. Mori")
plot.dens(clus = BmoriAsm.clust, main="B. Mori ASM15162")
plot.dens(clus = Ccecro.clust, main = "C. Cecrops")
plot.dens(clus = Csupp.clust, main="C. Suppressalis")
plot.dens(clus = Dplex1.clust, main="D. Plexippus")
plot.dens(clus = Dplex3.clust, main="D. Plexippus v3")
plot.dens(clus = Herato.clust, main="H. Erato")
plot.dens(clus = Hmelpo1.clust, main="H. Melpomene")
plot.dens(clus = Hmelpo2.clust, main="H. Melpomene Hmel2")
plot.dens(clus = Laccius.clust, main= "L. Accius")
plot.dens(clus = Msex.clust, main= "M. Sexta")
plot.dens(clus = Mcinx.clust, main= "M. Cinxia")
plot.dens(clus = Obrum.clust, main= "O. Brumata")
plot.dens(clus = Pglauc.clust, main= "P. Glaucus")
plot.dens(clus = Pmach.clust, main= "P. Machaon")
plot.dens(clus = Ppoly.clust, main= "P. Polytes")
plot.dens(clus = Pxuth1.clust, main= "P. Xuthus Pap_xu")
plot.dens(clus = Pxuth2.clust, main= "P. Xuthus Pxut")
plot.dens(clus = Psennae.clust, main="P. Sennae")
plot.dens(clus = Pinter.clust, main="P. Interpunctella")
plot.dens(clus = Pxylo.clust, main="P. Xylostella")

