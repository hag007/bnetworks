# source("https://bioconductor.org/biocLite.R")
# biocLite("prot2D")
.libPaths("/specific/netapp5/gaga/hagailevi/evaluation/Renv")
library(fdrtool, lib.loc  = "/specific/netapp5/gaga/hagailevi/evaluation/Renv")
library(st, lib.loc  = "/specific/netapp5/gaga/hagailevi/evaluation/Renv")
# library(prot2D, lib.loc  = "/specific/netapp5/gaga/hagailevi/evaluation/Renv")
# library(prot2D, lib.loc  = "/specific/netapp5/gaga/hagailevi/evaluation/Renv")
# library(prot2D, lib.loc  = "/specific/netapp5/gaga/hagailevi/evaluation/Renv")
library(prot2D, lib.loc  = "/specific/netapp5/gaga/hagailevi/evaluation/Renv")

# emp_file_name="/specific/netapp5/gaga/hagailevi/evaluation/output/emp_fdr/MAX/emp_diff_ERS_1_bionet_passed_oob.tsv"
# discrete_interval=0.001
# data<-read.delim(emp_file_name, row.names = 1)
# pvals[pvals==0]=discrete_interval
# pvals<-data["emp_pval"]

pvals<-pvals[!is.na(pvals)]

res<-robust.fdr(pvals, sides = 1, discrete = T, use8 = T)

result<-res$q
# print(pvals)
