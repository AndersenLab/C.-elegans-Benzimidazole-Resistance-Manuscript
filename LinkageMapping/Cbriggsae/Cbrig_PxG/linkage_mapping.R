#devtools::install_github("AndersenLab/linkagemapping")
library(linkagemapping)
library(dplyr)
library(broom)
library(tidyr)

setwd("~/Dropbox/Andersenlab/LabFolders/Mostafa/Cbrig_mappings/")
#pheno <- read.csv("~/Dropbox/Andersen\ Lab/LabFolders/Mostafa/Cbrig_rawsorter/CbRIAILphenotypes_pr.csv")
#pheno <- readRDS("~/Dropbox/AndersenLab/LabFolders/PastMembers/Tyler/ForTrip/RIAILs2_processed.rds")

load("~/Dropbox/Andersenlab/LabFolders/Mostafa/Cbrig_rawsorter/cb_regressed.Rda")
pheno <- regressed

#get cross object
data("AF16xHK104cross")
cross <- AF16xHK104cross

# Merge the cross object and the phenotype data
cross <- mergepheno(cross, pheno)

#FDR
map <- fsearch(cross, permutations = 1000, doGPU = FALSE,threshold = "FDR")
#GWER
map <- fsearch.manmarker(cross, permutations = 1000, doGPU = FALSE, threshold = "GWER")

# Annotate the LOD scores
annotatedlods <- annotate_lods(map, cross)

#save gwer mapping
cb_mapping_gwer <- annotatedlods
save(cb_mapping_gwer, file="~/Dropbox/Andersenlab/LabFolders/Mostafa/Cbrig_mappings/cb_mapping_gwer.Rda")

#save fdr mapping
cb_mapping_fdr <- annotatedlods
save(cb_mapping_fdr, file="~/Dropbox/Andersenlab/LabFolders/Mostafa/Cbrig_mappings/cb_mapping_fdr.Rda")



