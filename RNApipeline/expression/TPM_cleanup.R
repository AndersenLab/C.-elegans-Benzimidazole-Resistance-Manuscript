library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

########################################################
####### TPM Counts for mRNA and 22G expression #########
########################################################
setwd("~/GitHub/BZlocal/expression/")

#Load in stringtie normalized abundances for R
RN2A <- read.csv("~/GitHub/BZlocal/expression/Rexpression/N2A/N2A_abund.tab", header=T, sep="\t") %>%
  mutate(N2_A_R_TPM = TPM) %>% select(-FPKM,-TPM)
RN2B <- read.csv("~/GitHub/BZlocal/expression/Rexpression/N2B/N2B_abund.tab", header=T, sep="\t") %>%
  mutate(N2_B_R_TPM = TPM) %>% select(Gene.ID,N2_B_R_TPM) 
RN2C <- read.csv("~/GitHub/BZlocal/expression/Rexpression/N2C/N2C_abund.tab", header=T, sep="\t") %>%
  mutate(N2_C_R_TPM = TPM) %>% select(Gene.ID,N2_C_R_TPM) 
RN2D <- read.csv("~/GitHub/BZlocal/expression/Rexpression/N2D/N2D_abund.tab", header=T, sep="\t") %>%
  mutate(N2_D_R_TPM = TPM) %>% select(Gene.ID,N2_D_R_TPM) 
RCBA <- read.csv("~/GitHub/BZlocal/expression/Rexpression/CBA/CBA_abund.tab", header=T, sep="\t") %>%
  mutate(CB_A_R_TPM = TPM) %>% select(Gene.ID,CB_A_R_TPM) 
RCBB <- read.csv("~/GitHub/BZlocal/expression/Rexpression/CBB/CBB_abund.tab", header=T, sep="\t") %>%
  mutate(CB_B_R_TPM = TPM) %>% select(Gene.ID,CB_B_R_TPM) 
RCBC <- read.csv("~/GitHub/BZlocal/expression/Rexpression/CBC/CBC_abund.tab", header=T, sep="\t") %>%
  mutate(CB_C_R_TPM = TPM) %>% select(Gene.ID,CB_C_R_TPM) 
RCBD <- read.csv("~/GitHub/BZlocal/expression/Rexpression/CBD/CBD_abund.tab", header=T, sep="\t") %>%
  mutate(CB_D_R_TPM = TPM) %>% select(Gene.ID,CB_D_R_TPM) 
R_TPM <- merge(RN2A,RN2B,by="Gene.ID")
R_TPM <- merge(R_TPM,RN2C,by="Gene.ID")
R_TPM <- merge(R_TPM,RN2D,by="Gene.ID")
R_TPM <- merge(R_TPM,RCBA,by="Gene.ID")
R_TPM <- merge(R_TPM,RCBB,by="Gene.ID")
R_TPM <- merge(R_TPM,RCBC,by="Gene.ID")
R_TPM <- merge(R_TPM,RCBD,by="Gene.ID")
colnames(R_TPM)[1] <- c("WB_ID")
colnames(R_TPM)[2] <- c("Gene_ID")
R_TPM <- R_TPM %>% select(-Gene_ID,-Coverage)

#Load in stringtie normalized abundances for G
GN2A <- read.csv("~/GitHub/BZlocal/expression/Gexpression/N2AP/N2AP_abund.tab", header=T, sep="\t") %>%
  mutate(N2_A_G_TPM = TPM) %>% select(-FPKM,-TPM)
GN2B <- read.csv("~/GitHub/BZlocal/expression/Gexpression/N2BP/N2BP_abund.tab", header=T, sep="\t") %>%
  mutate(N2_B_G_TPM = TPM) %>% select(Gene.ID,N2_B_G_TPM) 
GN2C <- read.csv("~/GitHub/BZlocal/expression/Gexpression/N2CP/N2CP_abund.tab", header=T, sep="\t") %>%
  mutate(N2_C_G_TPM = TPM) %>% select(Gene.ID,N2_C_G_TPM) 
GN2D <- read.csv("~/GitHub/BZlocal/expression/Gexpression/N2DP/N2DP_abund.tab", header=T, sep="\t") %>%
  mutate(N2_D_G_TPM = TPM) %>% select(Gene.ID,N2_D_G_TPM) 
GCBA <- read.csv("~/GitHub/BZlocal/expression/Gexpression/CBAP/CBAP_abund.tab", header=T, sep="\t") %>%
  mutate(CB_A_G_TPM = TPM) %>% select(Gene.ID,CB_A_G_TPM) 
GCBB <- read.csv("~/GitHub/BZlocal/expression/Gexpression/CBBP/CBBP_abund.tab", header=T, sep="\t") %>%
  mutate(CB_B_G_TPM = TPM) %>% select(Gene.ID,CB_B_G_TPM) 
GCBC <- read.csv("~/GitHub/BZlocal/expression/Gexpression/CBCP/CBCP_abund.tab", header=T, sep="\t") %>%
  mutate(CB_C_G_TPM = TPM) %>% select(Gene.ID,CB_C_G_TPM) 
GCBD <- read.csv("~/GitHub/BZlocal/expression/Gexpression/CBDP/CBDP_abund.tab", header=T, sep="\t") %>%
  mutate(CB_D_G_TPM = TPM) %>% select(Gene.ID,CB_D_G_TPM) 
G_TPM <- merge(GN2A,GN2B,by="Gene.ID")
G_TPM <- merge(G_TPM,GN2C,by="Gene.ID")
G_TPM <- merge(G_TPM,GN2D,by="Gene.ID")
G_TPM <- merge(G_TPM,GCBA,by="Gene.ID")
G_TPM <- merge(G_TPM,GCBB,by="Gene.ID")
G_TPM <- merge(G_TPM,GCBC,by="Gene.ID")
G_TPM <- merge(G_TPM,GCBD,by="Gene.ID")
colnames(G_TPM)[1] <- c("WB_ID")
colnames(G_TPM)[2] <- c("Gene_ID")
G_TPM <- G_TPM %>% select(-Gene_ID,-Coverage)

#merge mRNA (R) and 22G (G) dfs
TPMf <- merge(R_TPM,G_TPM,by="WB_ID") %>% select(1:13,18:25)
colnames(TPMf)[2] <- c("Chr")
colnames(TPMf)[3] <- c("Strand")
colnames(TPMf)[4] <- c("Start")
colnames(TPMf)[5] <- c("End")

save(TPMf,file="TPM.Rda")