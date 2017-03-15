#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
#biocLite("ReportingTools")

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(DESeq2)
library(ReportingTools)

setwd("~/GitHub/BZlocal/diffexp")

##########################################
# mRNA and 22G DESEQ2
##########################################

RcountData = read.delim("~/GitHub/BZlocal/diffexp/Rdiffexp/gene_count_matrix.csv",sep=',', row.names = 1)
colnames(RcountData)[1] <- c("WB_ID")
GcountData = read.delim("~/GitHub/BZlocal/diffexp/Gdiffexp/gene_count_matrix.csv",sep=',', row.names = 1)
colnames(GcountData)[1] <- c("WB_ID")

samples <- data.frame(row.names=c("CBA","CBB","CBC","CBD","N2A","N2B","N2C","N2D"), condition=as.factor(c(rep("CB",4),rep("N2",4))))

RcountMatrix <- DESeqDataSetFromMatrix(RcountData, colData=samples, design =~condition, tidy = FALSE, ignoreRank = FALSE)
GcountMatrix <- DESeqDataSetFromMatrix(GcountData, colData=samples, design =~condition, tidy = FALSE, ignoreRank = FALSE)
#Remove all 0s: nrow(countMatrix) 46766  -> 17399
#countMatrix <- countMatrix[ rowSums(counts(countMatrix)) > 1, ]


## Pre-DE analysis: log-transform and look at sample distances
library("pheatmap")
library("RColorBrewer")
Rrld <- rlog(RcountMatrix, blind=FALSE) #head(assay(Rrld), 3)
Grld <- rlog(GcountMatrix, blind=FALSE) #head(assay(Grld), 3)
RsampleDists <- dist(t(assay(Rrld)))
RsampleDists
GsampleDists <- dist(t(assay(Grld)))
GsampleDists
#euclidean distances of log-transformed data (to estimate contribution across larger geneset)
RsampleDistMatrix <- as.matrix( RsampleDists )
GsampleDistMatrix <- as.matrix( GsampleDists )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
Rheatmap <- pheatmap(RsampleDistMatrix,
        clustering_distance_rows=RsampleDists,
        clustering_distance_cols=RsampleDists,
        col=colors)
Gheatmap <- pheatmap(GsampleDistMatrix,
         clustering_distance_rows=GsampleDists,
         clustering_distance_cols=GsampleDists,
         col=colors)

## RUN DESEQ
RcountMatrixDES <- DESeq(RcountMatrix)
Rres <- results(RcountMatrixDES)
GcountMatrixDES <- DESeq(GcountMatrix)
Gres <- results(GcountMatrixDES)
#mcols(Rres, use.names=TRUE)
#summary(Rres, alpha=0.05)
#summary(Gres, alpha=0.05)
#plotMA(Rres, alpha = 0.05, main="MA-plot")
#plotMA(Gres, alpha = 0.05, main="MA-plot")

Rres.ordered <- Rres[order(Rres$padj),]
Rres.df <- as.data.frame(Rres.ordered)
Gres.ordered <- Gres[order(Gres$padj),]
Gres.df <- as.data.frame(Gres.ordered)

head(Rres.df)
head(Gres.df)

#Combine the dataframes
Rres.df_comb <- Rres.df %>%
  mutate(WB_ID = rownames(Rres.df))
colnames(Rres.df_comb) <- c("RbaseMean","Rlog2FC","RlfcSE","Rstat","Rpvalue","Rpadj","WB_ID")

Gres.df_comb <- Gres.df %>%
  mutate(WB_ID = rownames(Gres.df))
colnames(Gres.df_comb) <- c("GbaseMean","Glog2FC","GlfcSE","Gstat","Gpvalue","Gpadj","WB_ID")

exp_merged <- merge(Rres.df_comb,Gres.df_comb,by="WB_ID")
exp_merged <- exp_merged %>%
  select(WB_ID, Rlog2FC, Rpadj, Glog2FC, Gpadj)

save(exp_merged,file="Diffexp.Rda")

#ignore below

plot <- ggplot(exp_merged, aes(x = Glog2FC, y=Rlog2FC))+
  geom_point(aes(), size = 2, alpha = 0.75) +
  ylab("22G log2FC") + xlab("mRNA log2FC") + theme(axis.ticks = element_blank()) +
  theme_bw() +
  theme(axis.ticks.y = element_blank()) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
  theme(legend.position="none")
plot


write.csv(exp_merged, file="DES2results.csv")


########################################################
####### TPM Counts for mRNA and 22G expression #########
########################################################

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

