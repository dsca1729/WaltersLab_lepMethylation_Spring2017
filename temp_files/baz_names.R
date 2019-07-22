#==========================Dependencies========================================
library(stringr)
source("~/GitHub/WaltersLab_lepMethylation_Spring2017/analysis_scripts/oe_clust_v0.R")

#=============Function Definitions=============================================
get.seqinfo <- function(dataList, ...) {
  df.seq <- read.delim(dataList, sep = "\t", nrows = 1, as.is = TRUE)
  seqinfo <- df.seq$seq_info1seq_info2
  return(seqinfo)
}

get.taxa <- function(string) {
  i.start <- 21 #uniform file name length ensures all names begin at 21st index position
  loc <- str_locate(string, "transcribed")
  i.end <- loc[1]-1  #sets cutoff point to start
  taxa <- str_sub(string, start = i.start, end=i.end)
  return(taxa)
} 
gen.taxalist <- function(list, fun, exp = FALSE, fname = "", ...) {
  list.out <- lapply(list, fun)
  names(list.out) <- list.out
  if((exp == TRUE) && (fname != "")) {
    write.table(list.out, file = fname, row.names = FALSE, col.names = FALSE)
  }
  return(list.out)
}
#====================Script Body===============================================
files <- list.files()
l.seqs <- lapply(files, get.seqinfo)
taxa.list <- lapply(l.seqs, get.taxa)
#write.table(taxa.list, file = "bazinet_taxa.txt", row.names = FALSE, col.names = FALSE)

#setwd("~/GitHub/WaltersLab_lepMethylation_Spring2017/analysis_scripts/input_files/bazinet_data")
l.baz <- read_input_data()
l.baz <- filter_input(l.baz)
clus.baz <- gen_density_data(l.baz)
mean.baz <- mean_table(clus.baz)
var.baz <- variance_table(clus.baz)
gen_plots(clus.baz, "bazinet-data", taxa.list)
export.tables("bazinet-data", means = mean.baz, vars = var.baz)