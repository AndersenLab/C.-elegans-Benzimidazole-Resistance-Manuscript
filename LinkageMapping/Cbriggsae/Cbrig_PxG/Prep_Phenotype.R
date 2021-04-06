library(stringr)
library(dplyr)
library(reshape2)

setwd("~/Dropbox/Andersen\ Lab/LabFolders/Mostafa/Cbrig_PxG/")
pheno <- read.csv("~/Dropbox/Andersen\ Lab/LabFolders/Mostafa/Cbrig/PrunedProcessedCbRIAILphenotypes_r.csv") 
#pheno <- read.csv("~/Dropbox/Andersen\ Lab/LabFolders/Mostafa/Cbrig/PrunedProcessedCbRIAILphenotypes_p1.csv") 
#pheno <- read.csv("~/Dropbox/Andersen\ Lab/LabFolders/Mostafa/Cbrig/PrunedProcessedCbRIAILphenotypes_p2.csv") 

#original version (missing parental strains)
phenoa <- pheno %>% 
  filter(!(strain %in% c("AF16", "HK104"))) %>% 
  gather(trait, value, resid.cv.EXT:resid.var.yellow) %>%
  mutate(n.trait = paste(drug, trait, sep="_")) %>% 
  select(strain, n.trait, value) %>% 
  rename(trait = n.trait) %>%
  filter(! is.na(value))
Lphen_cb <- phenoa
save(Lphen_cb, file="Lphen_cb.Rdata")

#include parental strains (AF and HK)
phenob <- pheno %>% 
  gather(trait, value, resid.cv.EXT:resid.var.yellow) %>%
  mutate(n.trait = paste(drug, trait, sep="_")) %>% 
  select(strain, n.trait, value) %>% 
  rename(trait = n.trait) %>%
  filter(! is.na(value))
Lphen_cb_wp <- phenob
save(Lphen_cb_wp, file="Lphen_cb_wp.Rdata")
load("~/Dropbox/Andersen\ Lab/LabFolders/Mostafa/Cbrig_PxG/Lphen_cb_wp.Rdata")



  
