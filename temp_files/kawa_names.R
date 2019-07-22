library(stringr)
source("~/GitHub/WaltersLab_lepMethylation_Spring2017/analysis_scripts/oe_clust_v0.R")

df.taxa <- read.delim("~/GitHub/WaltersLab_lepMethylation_Spring2017/temp_files/final_kawa_clean.txt", header = TRUE, as.is = TRUE)
taxa.list <- df.taxa$Taxa


l.kawa <- read_input_data()
l.kawa <- filter_input(l.kawa)
clus.kawa <- gen_density_data(l.kawa)
names(clus.kawa) <- taxa.list
mean.kawa <- mean_table(clus.kawa)
var.kawa <- variance_table(clus.kawa)
gen_plots(clus.kawa, "kawahara-named", taxa.list)
export.tables("kawahara-named", means = mean.kawa, vars = var.kawa)