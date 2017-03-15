#library(devtools)
#install_github("AndersenLab/easysorter")
library(easysorter)
library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(cowplot)
library(data.table)
library(ggrepel)
library(viridis)

setwd("~/Dropbox/Andersenlab/LabFolders/Mostafa/Manuscripts/Benzimidazoles/Github-BZ/Celegans_other/QTLnarrow")

########################
#### PLOT: NIL NARROWING
########################
### load and plot
load("D_NILs.Rda")
regressed$strain <- factor(regressed$strain, 
                           levels=c("SnuIR3","SnuIR6","SnuIR7","SnuIR9","SnuIR10","SnuIR12"))

strainsplot <- c("SnuIR3","SnuIR6","SnuIR7","SnuIR9","SnuIR10","SnuIR12")

alb <- subset(regressed, condition=="Albendazole") %>%
  filter(trait == "q75.TOF")%>%
  filter(strain %in% strainsplot)

NILphenoplot <- ggplot(alb)+
  aes(x = strain, y = phenotype, fill = strain)+
  geom_boxplot(outlier.size= NA, alpha = 0.5, lwd=1)+
  geom_jitter(size =2.5, alpha=0.5, width =0.35)+
  theme_bw()+
  scale_colour_manual(values = c("N2" = "orange","CB" = "blue","SnuIR3" = "black", "SnuIR4" = "black", "SnuIR5" = "black", "SnuIR6" = "black", "SnuIR7" = "black","SnuIR9" = "black","SnuIR10" = "black","SnuIR11" = "black","SnuIR12" = "black","SnuIR13" = "black")) +
  scale_fill_manual(values = c("N2" = "orange","CB" = "blue","SnuIR3" = "black", "SnuIR4" = "black", "SnuIR5" = "black", "SnuIR6" = "black", "SnuIR7" = "black","SnuIR9" = "black","SnuIR10" = "black","SnuIR11" = "black","SnuIR12" = "black","SnuIR13" = "black")) +
  theme(axis.text.x = element_text(size=12, face="bold", color="black"),
        axis.text.y = element_text(size=14, face="bold", color="black"),
        axis.title.x = element_text(size=14, face="bold", color="black"),
        axis.title.y = element_text(size=16, face="bold", color="black"),
        strip.text.x = element_text(size=16, face="bold", color="black"),
        strip.text.y = element_text(size=16, face="bold", color="black"),
        plot.title = element_text(size=14, face="bold")) +
  labs(x = "", y = "Regressed Length (µm)", title = " ") + 
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  theme(legend.position = "none") +
  coord_flip() 
NILphenoplot

strains <- c("snuIR3","snuIR6","snuIR7","snuIR9","snuIR10","snuIR12")
plotRIAILgenotypes <- function(NILs){
  genos <- data.frame(fread("20141214_availableNILs.csv",
                            header = T))%>%
    select(strain, nil_genotype,chr,start,stop)%>%
    filter(strain %in% strains)
  genos$strain <- factor(genos$strain, levels = c("snuIR3","snuIR6","snuIR7","snuIR9","snuIR10","snuIR12"))
  genoplot <-  ggplot(genos) +
    aes(y=strain,color = ifelse(nil_genotype=="N2","N2",
                                ifelse(nil_genotype=="CB4856","CB",NA)))+
    geom_segment(mapping= aes(x = start/1e6, xend = stop/1e6, y=strain, yend=strain), alpha=1, size = 2.5)+
    geom_segment(mapping= aes(x = 15.2, xend = 16.1, y = "snuIR6", yend = "snuIR6"), colour = "blue", alpha = 0.2, size =2.5) +
    geom_segment(mapping= aes(x = 15.2, xend = 16.1, y = "snuIR7", yend = "snuIR7"), colour = "blue", alpha = 0.2, size =2.5) +
    geom_segment(mapping= aes(x = 15.2, xend = 16.1, y = "snuIR9", yend = "snuIR9"), colour = "blue", alpha = 0.2, size =2.5) +
    geom_segment(mapping= aes(x = 15.2, xend = 16.1, y = "snuIR10", yend = "snuIR10"), colour = "blue", alpha = 0.2, size =2.5) +
    geom_segment(mapping= aes(x = 15.2, xend = 16.1, y = "snuIR12", yend = "snuIR12"), colour = "blue", alpha = 0.2, size =2.5) +
    geom_segment(mapping= aes(x = 15.2, xend = 16.1, y = "snuIR3", yend = "snuIR3"), colour = "orange", alpha = 0.2, size =2.5) +
    geom_segment(mapping= aes(x = start/1e6, xend = stop/1e6, y=strain, yend=strain), alpha=1, size = 2.5)+
    scale_color_manual(values=c("N2"="orange","CB"="blue"),name="Genotype")+
    facet_grid(.~chr, scales = "free")+
    theme_bw() +
    theme(axis.text.x = element_text(size=12, face="bold", color="black"),
          axis.text.y = element_text(size=12, face="bold", color="black"),
          axis.title.x = element_text(size=14, face="bold", color="black"),
          axis.title.y = element_text(size=16, face="bold", color="black"),
          strip.text.x = element_text(size=16, face="bold", color="black"),
          strip.text.y = element_text(size=16, face="bold", color="black"),
          plot.title = element_text(size=0, face="bold")) +
    geom_vline(aes(xintercept = 15466666/1e6), size =0.75)+  
    geom_vline(aes(xintercept = 15907521/1e6), size =0.75)+  
    geom_vline(aes(xintercept = 15570359/1e6), color = "gray", linetype =2, size =1)+
    #geom_vline(aes(xintercept = 15731721/1e6), color = "gray", linetype =2, size =1)+  
    geom_vline(aes(xintercept = 15645536/1e6), color = "gray", linetype =2, size =1)+ 
    #geom_text(mapping=aes(x=start/1e6, y=strain, label=start), size=4, vjust=-1.5, hjust=0.5) +
    #geom_text(mapping=aes(x=stop/1e6, y=strain, label=stop), size=4, vjust=-1.5, hjust=0.5) +
    theme(legend.position = "none") + #added for lab meeting, mz
    labs(x = "Chr IV Genomic Position (Mb)", y = "") + 
    theme(strip.background = element_blank(),
          strip.text.x = element_blank()) +
    scale_y_discrete(labels = c("LJ1203","LJ1206","LJ1207","LJ1209","LJ1210","LJ1212")) +
    xlim(15.2,16.1)
}
NILgenoplot <- plotRIAILgenotypes(c(strains))
NILgenoplot

top_row <- plot_grid(NILgenoplot + theme(plot.margin=unit(c(0,0,0,0),"cm")), NILphenoplot + theme(plot.margin=unit(c(-0.6,0,0,0),"cm")), labels = c("A","B"), rel_widths = c(1.6, 0.9), vjust = -0.2, scale = 0.95)
top_row


########################
#### PLOT: NARROWED REGION VARIANT ANNOTATION
########################

#####################
#### LOAD Annotations
#####################

gff3 <- read.csv("gff3_all2.txt", header=F)
colnames(gff3) <- c("chrom","start","end","strand","WBid","seqname","biotype")

other_list <- c("ncRNA","tRNA","snoRNA","pseudogene")
interval <- gff3 %>%
  filter(chrom == "IV") %>%
  filter() %>%
  filter(biotype != "exon") %>%
  mutate(biotype2 = biotype) %>% #switch to 3 categories (protein,piRNA,other)
  mutate(biotype2 = ifelse(biotype %in% other_list, "Other", paste0(biotype))) %>%
  filter(start > 15570359 & end < 15645536)
interval$biotype2 <- factor(interval$biotype2, levels = c("Other","piRNA","protein_coding"))

#get protein exons (15570359 - 15645536)
protein_list <- gff3 %>%
  filter(chrom == "IV") %>%
  filter(start > 15570359 & end < 15645536) %>%
  filter(biotype == "protein_coding")
protein_list <- protein_list$seqname
protein_list_grepl <- paste(protein_list, collapse = '[a-z]|')
protein_list_grepl <- paste0(protein_list_grepl,'[a-z]')
protein_exons1 <- gff3 %>%
  filter(chrom == "IV") %>%
  filter(start > 15570359 & end < 15645536) %>%
  filter(seqname %in% protein_list)
protein_exons2 <- gff3 %>%
  filter(chrom == "IV") %>%
  filter(start > 15570359 & end < 15645536) %>%
  filter(grepl(protein_list_grepl,seqname))
protein_exons <- rbind(protein_exons1,protein_exons2) %>%
  filter(biotype != "protein_coding")

#get piRNAs (15570359 - 15645536)
piRNAs <- gff3 %>%
  filter(chrom == "IV") %>%
  filter(start > 15570359 & end < 15645536) %>%
  filter(biotype == "piRNA")


#####################
#### LOAD variants
#####################

#For figure generation start here
load("piCBvars_small.Rda")
CBvars_small <- CBvars_small %>%
  filter(FILTER == "PASS")

# #2. look at all variants (regardless of impact) that overlap with exons
# #unique positions of all variants (361 total)
# varpos <- unique(CBvars_small$POS) 
# #overlap with exons in updated gff (protein_exons)
# for (i in varpos){
#   for (j in 1:nrow(protein_exons)){
#     if(i > protein_exons$start[j] && i < protein_exons$end[j]){
#       print(paste0("Pos: ",i,", start: ",protein_exons$start[j],", end: ",protein_exons$end[j],", strand: ",protein_exons$strand[j],", seqname: ",protein_exons$seqname[j]))
#     }
#   }
# }
# #[1] "Pos: 15579097, start: 15579086, end: 15579139, strand: +, seqname: Y105C5A.3"
# #[1] "Pos: 15585361, start: 15585225, end: 15585500, strand: +, seqname: Y105C5A.508"
# #[1] "Pos: 15594290, start: 15593960, end: 15594913, strand: -, seqname: Y105C5A.5"
# #[1] "Pos: 15596217, start: 15595692, end: 15596645, strand: +, seqname: Y105C5A.6"

#coding exon vars
CBvars_c <- CBvars_small %>% 
  filter(transcript_biotype == "Coding") %>%
  filter(POS > 15570359 & POS < 15645536) %>%
  filter(impact == "HIGH" | impact == "MODERATE") %>%
  #group_by(gene_id) %>%
  distinct(POS, .keep_all = TRUE)

# #noncoding vars within ncRNA ('Other') exons (within 100, ignore)
# CBvars_nc <- CBvars_small %>% 
#   filter(transcript_biotype == "Noncoding") %>%
#   filter(POS > 15570359 & POS < 15645536) %>%
#   #filter(impact == "LOW") %>%
#   #filter(impact == "MODIFIER") %>%
#   #group_by(gene_id) %>%
#   distinct(POS, .keep_all = TRUE)

#Cross reference all vcf variants with piRNA-filtered gff for overlap
# varpos <- unique(CBvars_small$POS)  #positions of all variants (361)
# #overlap with piRNAs in updated gff (piRNAs)
# piRNAs$exonvar <- 'NA'
# for (i in varpos){
#   for (j in 1:nrow(piRNAs)){
#     if(i >piRNAs$start[j] && i < piRNAs$end[j]){
#       piRNAs$exonvar[j] <- 1
#       #print(paste0("Pos: ",i,", start: ",piRNAs$start[j],", end: ",piRNAs$end[j],", strand: ",piRNAs$strand[j],", seqname: ",piRNAs$seqname[j]))
#     }
#   }
# }
# #overlap with US (50 nts)
# piRNAs <- piRNAs %>%
#   mutate(USl = ifelse(strand == "+", start-51, end+1)) %>%
#   mutate(USr = ifelse(strand == "+", start-1, end+50))
# piRNAs$USvar <- 'NA'
# for (i in varpos){
#   for (j in 1:nrow(piRNAs)){
#     if(i >piRNAs$USl[j] && i < piRNAs$USr[j]){
#       piRNAs$USvar[j] <- 1
#       #print(paste0("Pos: ",i,", start: ",piRNAs$start[j],", end: ",piRNAs$end[j],", strand: ",piRNAs$strand[j],", seqname: ",piRNAs$seqname[j]))
#     }
#   }
# }

#change to Mb
interval <- interval %>%
  mutate(start = start/1e6, end = end/1e6)
CBvars_c <- CBvars_c %>%
  mutate(POS = POS/1e6)
piRNAs <- piRNAs %>%
  mutate(POSexon = (start+end)/(2*1e6)) #%>%
#mutate(POSus = (USl+USr)/(2*1e6))

plot_annotation <- ggplot(interval)+
  geom_segment(aes(x=start, y=biotype2, xend=end, yend=biotype2, colour=biotype2), size=12)+
  #plot protein variants with predicted consequences
  geom_point(data=CBvars_c, aes(x=POS, y=3), size = 4, shape = 21, colour = "black", fill = "black", alpha = 1) +
  geom_point(data=CBvars_c, aes(x=POS, y=3), size = 2.5, shape = 21, colour = "royalblue2", fill = "royalblue2", alpha = 1) +
  #plot US piRNAs variants
  #geom_point(data=subset(piRNAs,USvar == 1), aes(x=POSus, y=2.1), size = 4, shape = 21, colour = "black", fill = "black", alpha = 0.8) +
  #geom_point(data=subset(piRNAs,USvar == 1), aes(x=POSus, y=2.1), size = 2.5, shape = 21, colour = "yellow", fill = "yellow", alpha = 1) +
  #plot unique piRNAs (piRNA variants, includes 2 in C)
  #geom_point(data=subset(piRNAs,exonvar == 1), aes(x=POSexon, y=1.9), size = 4, shape = 21, colour = "black", fill = "black", alpha = 0.8) +
  #geom_point(data=subset(piRNAs,exonvar == 1), aes(x=POSexon, y=1.9), size = 2.5, shape = 21, colour = "royalblue2", fill = "royalblue2", alpha = 1) +
  theme_bw() +  
  labs(x = "", y = "") + 
  theme_bw(base_size = 12) +
  #scale_fill_hue(c=45, l=80) +
  #scale_color_viridis(discrete=TRUE) +
  scale_colour_manual(values = c("protein_coding" = "grey70","piRNA" = "grey70","Other" = "grey70")) +
  scale_fill_manual(values = c("protein_coding" = "black","piRNA" = "grey46","Other" = "grey46")) +
  geom_label_repel(data=subset(interval, biotype2 == "protein_coding"),aes(x=(start+end)/2, y=biotype2, label = seqname), 
                   fill = "grey55", 
                   segment.color = "grey",
                   size = 2.5, fontface = 'bold', color = "white",
                   segment.size = 0.35,
                   box.padding = unit(0.4, "lines"),point.padding = unit(0.15, "lines"),
                   force = 3.5, max.iter = 2e5,
                   nudge_x = ifelse(interval$strand == "+",0,0),
                   nudge_y = ifelse(interval$strand == "+",0.5,-0.5))+
  theme(axis.text.x = element_text(size=10, face="bold", color="black"),
        axis.text.y = element_text(size=12, face="bold", color="black"),
        axis.title.x = element_text(size=14, face="bold", color="black"),
        axis.title.y = element_text(size=12, face="bold", color="black"),
        strip.text.x = element_text(size=14, face="bold", color="black"),
        strip.text.y = element_text(size=14, color="black", face="bold", angle = 0),
        plot.title = element_text(size=14, face="bold")) +
  theme(legend.position = c(0.5,0.075),
        legend.direction="horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.key = element_rect(colour = "NA"),
        legend.key.width = unit(0.5, "in"),
        legend.key.height = unit(0.2, "in"),
        legend.background = element_rect(
          fill = "NA",
          #color = "black",
          size = 0.15)) +
  labs(x = "Chr IV Genomic Position (Mb)", y = "") + 
  theme(legend.position="none") +
  xlim(15.57,15.65) +
  scale_y_discrete(labels=c("piRNA" = "piRNA", "protein_coding" = "Protein", "Other" = "Other"))
plot_annotation
#NOTE: geom_repel needs fixed

###############
#### FINAL PLOT (NIL narrowing + protein-coding annotation)
###############

final_plot <- plot_grid(top_row, plot_annotation + theme(plot.margin=unit(c(1,0.4,0.2,0.2),"cm")), nrow=2, labels = c('','C'), rel_heights = c(1.25,1), vjust = 2)
final_plot
ggsave(plot = final_plot, filename = "QTLinterval.tiff", width = 15, height = 8)


########################
#### SUPPLEMENTAL FILES
########################

#####################
####### (snpeff/WS255 variant annotation across QTL region )
#####################

#IV:13500000-17200000 (piRNA cluster)
#IV:15466666-15907521 (QTL) 
#IV:15570359-15769804 (Narrow1)
#IV:15570359-15645536 (Narrow2)

#TableS2A 15466666-15907521 (QTL region) 
#QTL interval 'coding'
#QTL other 'coding' #ignore
#QTL piRNA 'coding'
#QTL piRNA 'upstream'

#ALL PROTEIN-CODING VARIANTS IN QTL AND NARROWED INTERVAL
load("piCBvars_big.Rda")
CBvars_c_big <- CBvars_big %>%
  filter(FILTER == "PASS") %>%
  filter(transcript_biotype == "Coding") %>%
  filter(impact == "HIGH" | impact == "MODERATE") %>%
  distinct(POS, .keep_all = TRUE) %>%
  mutate(NIL_narrowed_interval = ifelse(POS >= 15570359 & POS <= 15645536,"YES","NO")) %>%
  select(CHROM,POS,REF,ALT,effect,impact,gene_name,gene_id,feature_id,exon_intron_rank,nt_change,aa_change,protein_position,NIL_narrowed_interval)
colnames(CBvars_c_big) <- c("Chr","Pos", "Ref", "Alt","Effect","Impact","Gene_Name","Gene_ID", "Feature_ID", "Exon_Intron_Rank", "NT_change","AA_change","Prot_position","NIL_narrowed_interval")
write.csv(CBvars_c_big, file = "Ce_QTL_protvars.csv")

#Includes all piRNAS in interval
gff3 <- read.csv("gff3_all2.txt", header=F)
colnames(gff3) <- c("chrom","start","end","strand","WBid","seqname","biotype")
piRNAs <- gff3 %>%
  filter(chrom == "IV") %>%
  filter(start > 15466666 & end < 15907521) %>% #QTL interval
  filter(biotype == "piRNA") %>%
  mutate(NIL_narrowed_interval = ifelse(start >= 15570359 & end <= 15645536,"YES","NO"))
colnames(piRNAs) <- c("Chr","Start", "End", "Strand","Gene_ID", "Feature_ID", "Biotype", "NIL_narrowed_interval")
write.csv(piRNAs, file = "Ce_QTL_piRNAvars.csv")
