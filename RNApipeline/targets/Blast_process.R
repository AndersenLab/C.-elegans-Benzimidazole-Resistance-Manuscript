library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

setwd("RNApipeline/targets")

N2unique_targets <- read.csv("blast/N2_blast.txt", header=F, sep="\t") 
colnames(N2unique_targets) <- c("Query","Chr","PID","AlnLen","Mismatch","Gaps","Qstart","Qend","Start","End","Evalue","Bitscore")
N2unique_targets <- N2unique_targets %>% mutate(Source = "N2unique")

CBunique_targets <- read.csv("blast/CB_blast.txt", header=F, sep="\t") 
colnames(CBunique_targets) <- c("Query","Chr","PID","AlnLen","Mismatch","Gaps","Qstart","Qend","Start","End","Evalue","Bitscore")
CBunique_targets <- CBunique_targets %>% mutate(Source = "CBunique")

All_targets <- rbind(N2unique_targets,CBunique_targets) %>%
  dplyr::select(-Gaps,-Evalue,-Bitscore) %>%
  dplyr::mutate(Query2 = Query) %>%
  tidyr::separate(Query, c("Del", "Query"), sep = ":", remove = TRUE) %>%
  dplyr::select(-Del) %>%
  tidyr::separate(Query, c("QChr","QLeft", "QRight"), sep = "_", remove = TRUE) %>%
  dplyr::mutate(Query = Query2) %>% 
  dplyr::select(-Query2) %>%
  #Filter for piRNAs originating from N2 QTL region (N2unique and US) N2 coords
  #Filter for piRNAs originating from CB QTL region (CBunique) CB coords
  dplyr::mutate(Delete = NA) %>%
  dplyr::mutate(Delete = ifelse(QChr == "IV" & Source == "N2unique" & QLeft >= 15570000 & QRight <=16010000,"NO",Delete)) %>%
  dplyr::mutate(Delete = ifelse(QChr == "IV" & Source == "CBunique" & QLeft >= 15338314 & QRight <=15772091,"NO",Delete)) %>%
  dplyr::filter(Delete == "NO") %>% 
  dplyr::select(-Delete) %>%
  dplyr::mutate(Delete = NA) %>%
  dplyr::mutate(Delete = ifelse(Chr == "IV" & Start > 15570000 & End <=16010000,"YES","NO")) %>%   #Filter out targets mapping to QTL region
  dplyr::filter(Delete == "NO") %>% 
  dplyr::select(-Delete)
All_targets <- All_targets %>% #rearrange columns
  dplyr::select(Query,Source,QChr,QLeft,QRight,Qstart,Qend,PID,AlnLen,Mismatch,Chr,Start,End) %>%
  dplyr::select(-QChr,-QLeft,-QRight,-Qstart,-Qend) %>%
  dplyr::mutate(TChr = Chr, TStart = Start, TEnd = End) %>% 
  dplyr::select(-Chr,-Start,-End)
save(All_targets, file = "Blast_targets.Rda")

load("Blast_targets.Rda")

#cross-reference putative piRNAs with known piRNAs
