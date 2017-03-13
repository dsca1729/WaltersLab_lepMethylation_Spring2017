Atrans <- read.delim(file = "Amyelois_transitella_v1_-_cds.fa.out.txt", as.is = TRUE)
Bmori <- read.delim(file = "Bombyx_mori_-_cds.fa.out.txt", as.is = TRUE)

hist(Atrans$CpGOE, breaks = 100)
hist(Bmori$CpGOE, breaks = 100)


components <- 1:2
Atrans.clust <- densityMclust(Atrans$CpGOE, G=components) 
Bmori.clust <- densityMclust(Bmori$CpGOE, G=components) 

summary(Atrans.clust)
summary(Bmori.clust)

plot(Atrans.clust, what = "density", data = Atrans$CpGOE, breaks = 50)
#plot(Atrans.clust, what = "BIC")
plot(Bmori.clust, what = "density", data = Atrans$CpGOE, breaks = 50)
par(mfrow = c(1,2))
plot(Bmori.clust, what = "diagnostic")

#plot(Bmori.clust, what = "BIC")




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

par(mfrow = c(1,2))
plot.dens(clus = Bmori.clust, main = "B. mori")

plot.dens(clus = Atrans.clust, main = "A. transitella")




# Below here were preliminary data simulations to understand how Mclust worked...

data(acidity)
data(faithful)
mod4 <- densityMclust(acidity)
mod4

x <- rnorm(1000)
y <- rnorm(300, mean = 3, sd = 3)

xy <- c(x,y)

hist(xy, breaks = 50)

modx <- densityMclust(x)
summary(modx)
plot(modx, what = "density", data = x)

mod <-densityMclust(xy)
summary(mod)
plot(mod, what = "density", data = x)
plot(mod)

