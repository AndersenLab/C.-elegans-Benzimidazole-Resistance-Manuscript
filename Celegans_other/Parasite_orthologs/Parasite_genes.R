library(dplyr)
library(ggplot2)
library(cegwas)
db <- get_db()

setwd("~/Dropbox/Andersenlab/LabFolders/Mostafa/Manuscripts/Benzimidazoles/Github-BZ/Celegans_other/Parasite_orthologs")

#### C. elegans (takes all QTL - finds overlapping set and all genes with variants in intervals)
#outputs ceQTLvars.csv

#Load TableS1 (C. elegans QTL)
ceQTL = read.csv("~/Dropbox/Andersenlab/LabFolders/Mostafa/Manuscripts/Benzimidazoles/Github-BZ/LinkageMapping/Celegans/Ce_Linkage.csv", header = TRUE) %>%
  select(-X)

#Reduced set: check for overlapping QTL for each drug 
ceQTLR <- ceQTL %>%
  group_by(Drug,Chr) %>%
  arrange(Drug,Chr,Left,Right)
j <- 1
ceQTLR$unique <- NA
for(i in 1:nrow(ceQTLR)){
  ceQTLR$unique[i] <- j
  ifelse((ceQTLR$Left[i+1] - ceQTLR$Right[i]) <= 1000 & ceQTLR$Chr[i] == ceQTLR$Chr[i+1] & ceQTLR$Drug[i] == ceQTLR$Drug[i+1], j <- j, j <- j+1)
}
ceQTLR <- ceQTLR %>%
  mutate(Left_min = min(Left)) %>%
  mutate(Right_max= max(Right)) %>%
  select(Drug,Chr,Left_min,Right_max) %>% 
  ungroup() %>%
  group_by(Drug,Chr,Left_min,Right_max) %>% 
  distinct() %>%
  ungroup()

ceQTLvars <- data.frame(Drug=character(),
                 QTLregion=character(), 
                 Gene_ID=character(), 
                 Feature_ID=character(), 
                 Gene_Name=character(), 
                 stringsAsFactors=FALSE) 

for(i in 1:length(ceQTLR)){
  Drug_temp <- ceQTLR$Drug[i]
  Chr_temp <- ceQTLR$Chr[i]
  L_temp <- ceQTLR$Left_min[i]
  R_temp <- ceQTLR$Right_max[i]
  region_temp <- paste0(Chr_temp,":",L_temp,"-",R_temp)
  regiondb <- snpeff(region_temp, severity = c("HIGH", "MODERATE"), elements = c("exon"),
                     long = TRUE, remote = FALSE, impute = FALSE)
  filterdb <- regiondb %>%
    filter(strain == "CB4856", transcript_biotype == "Coding") %>%
    select(gene_id,feature_id,gene_name) %>%
    distinct(gene_id,.keep_all=TRUE) %>%
    mutate(Drug = Drug_temp, QTLregion = region_temp) %>%
    select(Drug,QTLregion,Gene_ID = gene_id,Feature_ID = feature_id, Gene_Name = gene_name)
  ceQTLvars <- rbind(ceQTLvars,filterdb)
}
ceQTLvars <- ceQTLvars %>%
  mutate(Gene_Name = ifelse(is.na(Gene_Name),"NA",Gene_Name))
  
write.csv(ceQTLvars, file = "ceQTLvars.csv")

