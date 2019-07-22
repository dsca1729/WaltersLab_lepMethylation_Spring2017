# import dependencies----------------------------------------------------------
library(mclust)

# data input/helper variables--------------------------------------------------
list_filenames <- list.files()
names(list_filenames) <- list_filenames
index <- length(list_filenames)
data_density <- list()
# produce----------------------------------------------------------------------
data_cpg <- lapply(list_filenames, read.delim)


gen_density_clus <- function(data, index_vector, ...) {
  for (i in index_vector) {
    print (i)
    print(data[12])
  }
}
nameVec <- names(dataList2)
prop_list <- list()
mean_list <- vector("list", 24)
var_list <- list()
mean_frame <- data.frame(x = 1:24, y = 1:2, row.names = nameVec)
var_frame <- data.frame(x = 1:24, y = 1:2, row.names = nameVec )
#obj_list <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

for (i in 1:24) {
 dens_obj <- densityMclust(dataList2[[i]]$CpGOE, G=1:2)
 #append(mean_list[i], values = c(dens_obj$parameters$mean[1], dens_obj$parameters$mean[2])) #
 #prop_list[i] <- paste(dens_obj$parameters$mean[1], dens_obj$parameters$mean[2], sep="\t")
 mean_frame[i,1] <- dens_obj$parameters$mean[1]
 mean_frame[i,2] <- dens_obj$parameters$mean[2]
 var_frame[i,1] <- dens_obj$parameters$variance$sigmasq[1]
 var_frame[i,2] <- dens_obj$parameters$variance$sigmasq[2]
 #var_list[i] <- paste("variance1_", names(dataList2[i]), " ", dens_obj$parameters$variance$sigmasq[1], '\t', "var2_", names(dataList2[i]), dens_obj$parameters$variance$sigmasq[2], '\n', sep="")
 #obj_list[i] <- as.list.data.frame(dens_obj)
}