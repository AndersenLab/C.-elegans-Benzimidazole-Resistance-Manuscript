library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

setwd("RNApipeline/analyses")

#load DEseq differential expression output
load("../diffexp/Diffexp.Rda")

#load strintie TPM expression values
load("../expression/TPM.Rda")

#merge DE and TPM data and reorder columns 
DE_TPM <- merge(exp_merged,TPMf,by="WB_ID") %>%
  select(1,6:9,2:5,10:25)
save(DE_TPM, file = "DE_TPM.Rda")

#annotate DE_TPM
load("DE_TPM.Rda")
DE_TPM_ann <- DE_TPM %>%
  mutate(Annot = NA)
for (i in 1:nrow(DE_TPM_ann)){
  gene <- DE_TPM_ann$WB_ID[i]
  print(gene)
  req <- httr::GET(paste("http://www.wormbase.org/rest/widget/gene/",gene,
                         "/overview", sep = ""), config = httr::content_type_json())
  raw_text <- httr::content(req, "parsed")
  Sum_text <- raw_text$fields$concise_description$data$text
  
  # get gene location=======================================================
  req <- httr::GET(paste("http://www.wormbase.org/rest/widget/gene/", "WBGene00138721", "/location", sep = ""),
                   config = httr::content_type_json())
  loc = gsub("\\.\\.", "-", (httr::content(req, "parsed"))$fields$genomic_position$data[[1]]$label)
  #tempdf <- dplyr::data_frame(GeneID = gene, Gene_location=loc, Summary=Sum_text)
  if (length(Sum_text >0)){
    DE_TPM_ann$Annot[i] <- Sum_text
  }
}
save(DE_TPM_ann, file = "DE_TPM_ann.Rda")

##### START HERE
setwd("RNApipeline/analyses")

#load differential expression call and TPM values for mRNA and 22G expression + annotations
load("DE_TPM_ann.Rda")

#### filter based on p-value, fold-change...etc
DE_TPMf <- DE_TPM_ann %>%
  #Filter p-value (<= 0.05)
  dplyr::filter(Rpadj <= 0.05) %>%
  dplyr::filter(Gpadj <= 0.05) %>%
  #Filter FC (>=1.5-fold)
  #dplyr::filter(abs(Rlog2FC) > 0.5849) %>%
  dplyr::filter(abs(Glog2FC) > 0.5849) %>%
  #Filter expression in both CB and N2 (min ~ 1 TPM per sample)
  dplyr::filter(N2_A_R_TPM >= 1 & N2_B_R_TPM >= 1 & N2_C_R_TPM >= 1 & N2_D_R_TPM >= 1) %>%
  dplyr::filter(CB_A_R_TPM >= 1 & CB_B_R_TPM >= 1 & CB_C_R_TPM >= 1 & CB_D_R_TPM >= 1) %>%
  #filter anti-correlation
  #dplyr::filter((sign(Rlog2FC)+sign(Glog2FC)) == 0) %>%
  #Populate anti-correlation column
  dplyr::mutate(Anticorr = ifelse(sign(Rlog2FC)+sign(Glog2FC) == 0,"YES","NO")) %>%
  #Populate direction of effect column (N2 mRNA low (Rlog2FC < 0))
  dplyr::mutate(Direction = ifelse(Rlog2FC < 0,"CB4856-high","N2-high"))
save(DE_TPMf, file = "DE_TPMf.Rda")

############ 
#### Look for overlapping targets (blast-predicted) for unique 21mers
#overlap between DE_TMPf (Chr,Start,End) and All_targets (TChr,TStart,TEnd)

load("DE_TPMf.Rda") 

#load predicted blast targets for unique piRNAs (for unique piRNAs within QTL interval)
load("../targets/Blast_targets.Rda")

DE_TPMf_blast <- DE_TPMf
DE_TPMf_blast$Query <- 'NA'
DE_TPMf_blast$Source <- 'NA'
DE_TPMf_blast$PID <- 'NA'
DE_TPMf_blast$AlnLen <- 'NA'
DE_TPMf_blast$Mismatch <- 'NA'
DE_TPMf_blast$TChr <- 'NA'
DE_TPMf_blast$TStart <- 'NA'
DE_TPMf_blast$TEnd <- 'NA'

for (i in 1:nrow(DE_TPMf_blast)){
  print(i)
  for (j in 1:nrow(All_targets)){
    if(DE_TPMf_blast$Chr[i] == All_targets$TChr[j] && DE_TPMf_blast$Start[i] < All_targets$TStart[j] && DE_TPMf_blast$End[i] > All_targets$TEnd[j]){
      DE_TPMf_blast$Query[i] <- All_targets$Query[j]
      DE_TPMf_blast$Source[i] <- All_targets$Source[j]
      DE_TPMf_blast$PID[i] <- All_targets$PID[j]
      DE_TPMf_blast$AlnLen[i] <- All_targets$AlnLen[j]
      DE_TPMf_blast$Mismatch[i] <- All_targets$Mismatch[j]
      DE_TPMf_blast$TChr[i] <- All_targets$TChr[j]
      DE_TPMf_blast$TStart[i] <- All_targets$TStart[j]
      DE_TPMf_blast$TEnd[i] <- All_targets$TEnd[j]
    }
  }
}
save(DE_TPMf_blast, file = "DE_TPMf_blast.Rda")


#### OUTPUTS 
setwd("RNApipeline/analyses")

### Table: Supplementary Table 5

load("DE_TPMf_blast.Rda")
write.csv(DE_TPMf_blast, file = "DE_TPMf_blast.csv")

### PLOT: Figure 3C

load("DE_TPMf_blast.Rda")
volcano_data <- DE_TPMf_blast

volcano_plot <- ggplot(volcano_data)+
  aes(x=Glog2FC, y=Rlog2FC)+
  geom_point(colour = "gray", size = 1.5, alpha = 0.75)+
  geom_point(data= subset(volcano_data, Anticorr == "YES" & Direction == "CB4856-high"), colour = "red", size = 2, alpha = 0.6)  +
  geom_hline(aes(yintercept=0),linetype="dashed") +
  geom_vline(aes(xintercept=0),linetype="dashed") +
  theme_bw() +
  #ylim(-2.5,2.5) + 
  #xlim(-4,4) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(size=14, face="bold", color="black"),
        axis.text.y = element_text(size=12, face="bold", color="black"),
        axis.title.x = element_text(size=14, face="bold", color="black"),
        axis.title.y = element_text(size=14, face="bold", color="black"),
        strip.text.x = element_text(size=12, face="bold", color="black"),
        strip.text.y = element_text(size=12, face="bold", color="black"),
        plot.title = element_text(size=0, face="bold")) +
  xlab("log2(22G RNA Ratio)") + ylab("log2(mRNA Ratio)")
volcano_plot

  