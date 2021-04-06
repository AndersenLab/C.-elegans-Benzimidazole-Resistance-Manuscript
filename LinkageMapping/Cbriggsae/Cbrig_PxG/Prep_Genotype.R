library(dplyr)
library(reshape2)

setwd("~/Dropbox/Andersen\ Lab/LabFolders/Mostafa/Cbrig_PxG/")
load("~/Dropbox/Andersen\ lab/LabFolders/Mostafa/Cbrig/AF16xHK104_RIAILs.RData")

df1 <- data.frame(AF16xHK104cross$pheno$id) %>% rename(strain = AF16xHK104cross.pheno.id) #
#df2 <- data.frame(cbind(AF16xHK104cross$geno$I$data)) # single chromosome
df2 <- data.frame(cbind(AF16xHK104cross$geno$I$data,AF16xHK104cross$geno$II$data,AF16xHK104cross$geno$III$data,AF16xHK104cross$geno$IV$data,AF16xHK104cross$geno$V$data,AF16xHK104cross$geno$X$data)) #
df <- cbind(df1,df2)

Lgeno_cb <- melt(df) %>% rename(geno = value, chr_pos = variable)
Lgeno_cb$chr_pos <- substring(Lgeno_cb$chr_pos, 2)
save(Lgeno_cb, file="Lgeno_cb.Rdata")
#load("~/Dropbox/Andersen\ Lab/LabFolders/Mostafa/Cbrig_PxG/Lphen_cb.Rdata")