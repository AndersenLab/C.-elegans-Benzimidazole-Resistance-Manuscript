#Script developed by Tyler, all functions are within linkagemapping
#load the linkagemapping package

devtools::install_github("AndersenLab/linkagemapping")
library("linkagemapping")
library(dplyr)
library(ggplot2)
library(Hmisc)
#insert cross data
data("N2xCB4856cross")
cross <- N2xCB4856cross

bzeds <- c("albendazole","fenbendazole-15","fenbendazole-30","mebendazole","thiabendazole-125","thiabendazole-625")

#insert pheno data
pheno <- readRDS("~/Dropbox/AndersenLab/LabFolders/Mostafa/Cel_Linkage_Mapping/RIAILs2_processed.rds") %>%
#filter(grepl("monepantel", condition))
pheno <- readRDS("~/Dropbox/AndersenLab/LabFolders/Mostafa/Cel_Linkage_Mapping/RIAILs1_processed.rds") %>%
filter(condition %in% bzeds) 
#filter(grepl("abamectin", condition))
#filter(grepl("albendazole", condition))
#filter(grepl("fenbendazole-15", condition))
#filter(grepl("fenbendazole-30", condition))
#filter(grepl("mebendazole", condition))
#filter(grepl("thiabendazole-125", condition))
#filter(grepl("thiabendazole-625", condition))

  
#create completed cross object with pheno data set
mapcross <- mergepheno(cross, pheno, set = 2)

# get counts
temp <- mapcross$pheno %>%
  filter(set == 2) %>%
  select(-strain,-set) %>%
length(which(rowSums(is.na(temp))==ncol(temp))) #67
nrow(temp) #359

#map from subset of phenotype or load map object if mapping has been saved
map <- linkagemapping::fsearch(mapcross, phenotype="monepantel", permutations = 100, threshold = "GWER")
#map all df phenotypes (benzids)
map <- linkagemapping::fsearch(mapcross, phenotype=NULL, permutations = 100, threshold = "GWER")

an_map <- annotate_lods(lods = map,cross = mapcross)
#save(an_map, file="Monepantel_an_map.Rda")
#save(an_map, file="Abamectin_an_map.Rda")
#save(an_map, file="Albendazole_an_map.Rda")
#save(an_map, file="Fenbendazole-15_an_map.Rda")
#save(an_map, file="Fenbendazole-30_an_map.Rda")
#save(an_map, file="Mebendazole_an_map.Rda")
#save(an_map, file="Thiabendazole-125_an_map.Rda")
#save(an_map, file="Thiabendazole-625_an_map.Rda")
#save(an_map, file="Benzimidazoles_an_map.Rda")

setwd("~/Dropbox/AndersenLab/LabFolders/Mostafa/Linkage_Mapping/")
load("Albendazole_an_map.Rda")
q75.EXT <- filter(an_map, trait == "albendazole.q75.EXT")
q90.TOF <- filter(an_map, trait == "albendazole.q90.TOF")
mean.EXT <- filter(an_map, trait == "albendazole.mean.EXT")
n <- filter(an_map, trait == "albendazole.n")
norm.n <- filter(an_map, trait == "albendazole.norm.n")
lodplot(q75.EXT)
lodplot(q90.TOF)
lodplot(mean.EXT)
lodplot(n)
lodplot(norm.n)

load("Abamectin_an_map.Rda")
q75.EXT <- filter(an_map, trait == "abamectin.q75.EXT")
q90.TOF <- filter(an_map, trait == "abamectin.q90.TOF")
q75.TOF <- filter(an_map, trait == "abamectin.q75.TOF")
mean.TOF <- filter(an_map, trait == "abamectin.mean.TOF")
mean.EXT <- filter(an_map, trait == "abamectin.mean.EXT")
n <- filter(an_map, trait == "abamectin.n")
norm.n <- filter(an_map, trait == "abamectin.norm.n")


lodplot(mean.TOF)
pxgplot(cross = mapcross, map = mean.TOF)

lodplot(mean.EXT)
pxgplot(cross = mapcross, map = mean.EXT)

#### look at variants in interval (cegwas package)
library(cegwas)
#V2
V2regiondb <- snpeff("V:5894814-7482174", severity = c("HIGH", "MODERATE"), elements = c("exon"),
                   long = TRUE, remote = FALSE, impute = FALSE)
V2filterdb <- V2regiondb %>% 
  filter(strain == "CB4856", transcript_biotype == "Coding")

length(unique(V2filterdb$gene_id)) #355
V2_list <- unique(V2filterdb$gene_id)
write(V2_list, file = "V2list.txt",
      ncolumns = 1,
      append = FALSE, sep = "\n")


V1regiondb <- snpeff("V:2675410-3081026", severity = c("HIGH", "MODERATE"), elements = c("exon"),
                   long = TRUE, remote = FALSE, impute = FALSE)
V1filterdb <- V1regiondb %>% 
  filter(strain == "CB4856", transcript_biotype == "Coding")

length(unique(V1filterdb$gene_id)) #109
V1_list <- unique(V1filterdb$gene_id)
write(V1_list, file = "V1list.txt",
      ncolumns = 1,
      append = FALSE, sep = "\n")

#send off V1/V2filterdb and unique gene_id list
