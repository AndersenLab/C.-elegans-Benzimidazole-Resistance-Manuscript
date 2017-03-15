library("linkagemapping")
library("easysorter")
library(dplyr)
library(ggplot2)
library(stringr)
#library(ggrepel)
library(gridExtra)
library(cowplot)
library(viridis)
library(tidyr)

setwd("~/Dropbox/Andersenlab/LabFolders/Mostafa/Manuscripts/Benzimidazoles/Github-BZ/LinkageMapping/Celegans")
load("Ce_BZ_GWER.Rda")

benzpeaks <- Ce_BZ_GWER %>%
  filter (var_exp > 0.05) %>%
  filter(!grepl("cv|iqr|q10|q25|var|norm.yellow|norm.red|norm.EXT|f.L1|f.L2L3|f.L4|f.ad", trait)) %>%
  mutate(trait = sub("fenbendazole.15", "Fenbendazole_15",trait)) %>%
  mutate(trait = sub("fenbendazole.30", "Fenbendazole_30",trait)) %>%
  mutate(trait = sub("thiabendazole.125", "Thiabendazole_125",trait)) %>%
  mutate(trait = sub("thiabendazole.625", "Thiabendazole_625",trait)) %>%
  mutate(trait = sub("albendazole", "Albendazole",trait)) %>%
  mutate(sp.cond = factor(str_split_fixed(trait, "\\.", 2)[,1], levels=c("Albendazole","Fenbendazole_30","Fenbendazole_15","Mebendazole", "Thiabendazole_625","Thiabendazole_125"))) %>%
  group_by(chr,trait,lod) %>% distinct(trait,.keep_all=TRUE) %>% ungroup() %>%
  mutate(trait = sub("Fenbendazole_15.","",trait)) %>%
  mutate(trait = sub("Fenbendazole_30.","",trait)) %>%
  mutate(trait = sub("Thiabendazole_125.","",trait)) %>%
  mutate(trait = sub("Thiabendazole_625.","",trait)) %>%
  mutate(trait = sub("Albendazole.","",trait)) %>%
  mutate(drug = sp.cond) %>% #Add drug column (group concentrations for each drug)
  mutate(drug = sub("Fenbendazole_15", "Fenbendazole",drug)) %>%
  mutate(drug = sub("Fenbendazole_30", "Fenbendazole",drug)) %>%
  mutate(drug = sub("Thiabendazole_125", "Thiabendazole",drug)) %>%
  mutate(drug = sub("Thiabendazole_625", "Thiabendazole",drug)) 
#Group traits
sizetraits <- c("median.TOF","median.EXT","mean.TOF","mean.EXT","q75.TOF","q75.EXT","q90.TOF","q90.EXT")
broodtraits <- c("n","norm.n")
fltraits <- c("mean.yellow","median.yellow","q75.yellow","q90.yellow","mean.red","median.red","q75.red","q90.red")
benzpeaks$traitgroup <- "NA"
benzpeaks <- benzpeaks %>%
  mutate(traitgroup = ifelse(trait %in% sizetraits,"Length",traitgroup)) %>%
  mutate(traitgroup = ifelse(trait %in% broodtraits,"Brood Size",traitgroup)) %>%
  mutate(traitgroup = ifelse(trait %in% fltraits,"Pumping",traitgroup))

####################################
### PLOT A (Reduced lineplot of QTL, grouped by drug)
####################################

#Reduced set: check for overlapping QTL for each drug 
benzpeaksR <- benzpeaks %>%
  group_by(sp.cond,chr) %>%
  arrange(sp.cond,chr,ci_l_pos,ci_r_pos)
j <- 1
benzpeaksR$unique <- NA
for(i in 1:nrow(benzpeaksR)){
   benzpeaksR$unique[i] <- j
#  ifelse((benzpeaksR$ci_r_pos[i] - benzpeaksR$ci_l_pos[i+1]) >= 0.01*(benzpeaksR$ci_r_pos[i+1] - benzpeaksR$ci_l_pos[i+1]) &  benzpeaksR$chr[i] == benzpeaksR$chr[i+1] &  benzpeaksR$sp.cond[i] == benzpeaksR$sp.cond[i+1], j <- j, j <- j+1)
   ifelse((benzpeaksR$ci_l_pos[i+1] - benzpeaksR$ci_r_pos[i]) <= 1000000 &  benzpeaksR$chr[i] == benzpeaksR$chr[i+1] &  benzpeaksR$sp.cond[i] == benzpeaksR$sp.cond[i+1], j <- j, j <- j+1)
}
benzpeaksR <- benzpeaksR %>%
  group_by(unique) %>%
  arrange(desc(lod)) %>%
  distinct(unique,.keep_all=TRUE)
# for faceting issues
benzpeaksR <- benzpeaksR %>%
  mutate(height = "") %>%
  mutate(height = ifelse(sp.cond == "Albendazole","12.5 μM",height)) %>%
  mutate(height = ifelse(sp.cond == "Fenbendazole_15","15 μM",height)) %>%
  mutate(height = ifelse(sp.cond == "Fenbendazole_30","30 μM",height)) %>%
  mutate(height = ifelse(sp.cond == "Thiabendazole_125","125 μM",height)) %>%
  mutate(height = ifelse(sp.cond == "Thiabendazole_625","62.5 μM",height)) 

benzpeaksR$drug <- factor(benzpeaksR$drug, 
                            levels=c("Thiabendazole","Fenbendazole","Albendazole"))
benzpeaksR$height <- factor(benzpeaksR$height,
                            levels=c("12.5 μM","15 μM","30 μM","62.5 μM","125 μM"))

#set chromosome boundaries
newrows <- benzpeaksR[1,] 
newrows[1,] = c(NA,"I",5000000,NA,0,NA,NA,NA,NA,NA,1,NA,14972282,"Albendazole","Albendazole","Length",NA,"12.5 μM")
newrows[2,] = c(NA,"II",5000000,NA,0,NA,NA,NA,NA,NA,1,NA,15173999,"Albendazole","Albendazole","Length",NA,"12.5 μM")
newrows[3,] = c(NA,"III",5000000,NA,0,NA,NA,NA,NA,NA,1,NA,13829314,"Albendazole","Albendazole","Length",NA,"12.5 μM")
newrows[4,] = c(NA,"IV",5000000,NA,0,NA,NA,NA,NA,NA,1,NA,17450860,"Albendazole","Albendazole","Length",NA,"12.5 μM")
newrows[5,] = c(NA,"V",5000000,NA,0,NA,NA,NA,NA,NA,1,NA,20914693,"Albendazole","Albendazole","Length",NA,"12.5 μM")
newrows[6,] = c(NA,"X",5000000,NA,0,NA,NA,NA,NA,NA,1,NA,17748731,"Albendazole","Albendazole","Length",NA,"12.5 μM")
newrows$ci_l_pos <- as.numeric(newrows$ci_l_pos)
newrows$ci_r_pos <- as.numeric(newrows$ci_r_pos)
newrows$pos <- as.numeric(newrows$pos)
newrows$lod <- as.numeric(newrows$lod)
#CE_I 14,972,282, #CE_II 15,173,999, #CE_III 13,829,314, #CE_IV 17,450,863, #CE_V 20,914,693, #CE_VI 17,748,731 

plot_benz <- ggplot(benzpeaksR)+
  aes(x=pos/1E6, y=height)+
  theme_bw() +
  scale_fill_viridis(name = "LOD") + scale_color_viridis(name = "LOD") +
  geom_segment(aes(x = ci_l_pos/1e6, y = height, xend = ci_r_pos/1e6, yend = height, color = lod), size = 2, alpha = 1) +
  geom_segment(data=newrows,aes(x = 0, y = height, xend = ci_r_pos/1e6, yend = height), size = 2, alpha = 0) +
  geom_point(aes(fill=lod),colour = "black",size = 5, alpha = 1, shape = 21)+
  theme(legend.title = element_text(size = 12, face = "bold")) +
  xlab("Genomic Position (Mb)") + ylab("") +
  theme(strip.background = element_rect(colour = "black", fill = "white",
                                        size = 0.75, linetype = "solid")) +
  theme(axis.text.x = element_text(size=10, face="bold", color="black"),
        axis.text.y = element_text(size=12, face="bold", color="black"),
        axis.title.x = element_text(size=14, face="bold", color="black"),
        axis.title.y = element_text(size=12, face="bold", color="black"),
        strip.text.x = element_text(size=14, face="bold", color="black"),
        strip.text.y = element_text(size=14, color="black", face="bold", angle = 0),
        plot.title = element_text(size=14, face="bold")) +
  facet_grid(drug ~ chr, scales = "free", space = "free")
plot_benz

####################################
### PLOT B / C (Linkage plot albendazole.q75.EXT + split at peak marker)
####################################

trait_input <- "albendazole.q75.TOF"
linkplot <- filter(Ce_BZ_GWER, trait == trait_input) %>%
  filter(iteration == 1)
peaks <- linkplot %>%
    filter(var_exp > 0)
linkplot <- mutate(linkplot,LOD2 = 0)
for(i in 1:nrow(peaks)){
    linkplot$LOD2 <- ifelse(linkplot$pos >= peaks$ci_l_pos[i] & linkplot$pos <= peaks$ci_r_pos[i] & linkplot$chr == peaks$chr[i], linkplot$lod, linkplot$LOD2)}
plot <- ggplot(data=linkplot,aes(x=pos/1e6, y = lod, group = factor(iteration))) +
    geom_line(aes(colour = factor(iteration)), size = 1) +
    theme_bw() +
    geom_hline(aes(yintercept=2.993084), color ="#F8766D", linetype="dashed", size=0.70) +
    geom_ribbon(aes(ymin = 0, ymax = LOD2),  alpha = 0.5) +
    ggtitle(paste0(trait_input)) +
    theme(legend.position = "none") +
    scale_color_viridis(discrete=TRUE) +
    labs(x = "Genomic Position (Mb)", y = "LOD", title = " ") + 
    theme(strip.background = element_rect(colour = "black", fill = "white",
                                        size = 0.75, linetype = "solid")) +
    theme(axis.text.x = element_text(size=10, face="bold", color="black"),
          axis.text.y = element_text(size=12, face="bold", color="black"),
          axis.title.x = element_text(size=14, face="bold", color="black"),
          axis.title.y = element_text(size=14, face="bold", color="black"),
          strip.text.x = element_text(size=14, face="bold", color="black"),
          strip.text.y = element_text(size=14, face="bold", color="black")) +
    facet_grid(. ~ chr, scales = "free_x", space= "free_x") +
    theme(legend.position = "none")
plot

### PxG Split
data("N2xCB4856cross")
cross <- N2xCB4856cross
#pheno <- readRDS("~/Dropbox/AndersenLab/LabFolders/Mostafa/Cel_Linkage_Mapping/RIAILs1_processed.rds")
pheno <- readRDS("~/Dropbox/AndersenLab/LabFolders/Mostafa/Manuscripts/Benzimidazoles/GitHub-BZ-prep/LinkageMapping/Celegans/RIAILs1_processed.rds")
mapcross <- mergepheno(cross, pheno, set = 2)
pgtrait <- filter(Ce_BZ_GWER, trait == "albendazole.q75.EXT")
pxgplot <- pxgplot2(cross = mapcross, map = pgtrait)
pxgplot

#sum(mapcross$pheno$set == 2) #359

### for manuscript (calculating parental differences)
# pheno_parents <- pheno %>%
#   filter(trait == "q75.EXT", condition == "albendazole") %>%
#   filter(strain == "N2" | strain == "CB4856") %>%
#   group_by(strain) %>% 
#   mutate(strain_mean = mean(phenotype))

########################
#### FINAL Celegans LINKAGE PLOT FIGURE (A,B,C)
########################

bottom_row <- plot_grid(plot, pxgplot, labels = c('B', 'C'), rel_widths = c(2.7, 1), scale = 0.95)
bottom_row
final_plot <- plot_grid(plot_benz, bottom_row, nrow=2, labels = c('A',''), rel_heights = c(0.9,1))
final_plot
ggsave(plot = final_plot, filename = "Ce_Linkage.tiff", width = 15, height = 8)

####################################
#### PxG plot functions
####################################

extract_genotype=function(cross){
  
  # Pull out the genotypes into snp x strain matrix
  genomat <- qtl::pull.geno(cross)
  class(genomat) <- "numeric"
  
  # Handle genotype encodings with heterozygous individuals
  # (encoded as 1, 2, 3) and without (encoded 1, 2)
  if(max(genomat, na.rm = TRUE) == 3) {
    genomat <- genomat - 2
  } else {
    genomat <- (genomat * 2) - 3
  }
  return(genomat)
}

pxgplot2 <- function(cross, map, parent="N2xCB4856") {
  peaks <- map %>% 
    dplyr::group_by(iteration) %>%
    dplyr::filter(!is.na(var_exp)) %>%
    dplyr::do(head(., n=1))
  
  if(nrow(peaks) == 0) {
    stop("No QTL identified")
  }
  
  uniquemarkers <- gsub("-", "\\.", unique(peaks$marker))
  colnames(cross$pheno) <- gsub("-", "\\.", colnames(cross$pheno))
  
  pheno <- cross$pheno %>%
    dplyr::select_(map$trait[1])
  geno <- data.frame(extract_genotype(cross)) %>%
    dplyr::select(which(colnames(.) %in% uniquemarkers)) %>%
    data.frame(., pheno)
  
  colnames(geno)[1:(ncol(geno)-1)] <- sapply(colnames(geno)[1:(ncol(geno)-1)],
                                             function(marker) {
                                               paste(unlist(peaks[peaks$marker == gsub("\\.","-",marker),c("chr", "pos")]),collapse = ":")
                                             })
  colnames(geno)[ncol(geno)] <- "pheno"
  split <- tidyr::gather(geno, marker, genotype, -pheno)
  
  split$genotype <- sapply(split$genotype, function(x){
    if(is.na(x)) {
      return(NA)
    }
    if(parent=="N2xCB4856") {
      if(x == -1) {
        "N2"
      } else {
        "CB4856"
      }
    } else if(parent=="LSJ2xN2") {
      if(x == -1) {
        "LSJ2"
      } else {
        "N2"
      }
    } else if(parent=="AF16xHK104") {
      if(x==-1) {
        "AF16"
      } else {
        "HK104"
      }
    }
  })
  
  split$genotype <- factor(split$genotype, 
                             levels=c("N2","CB4856"))
  ggplot2::ggplot(split) +
    ggplot2::geom_boxplot(ggplot2::aes(x = genotype, y = pheno, fill = genotype), outlier.size= NA, alpha = 0.5, lwd=1) +
    ggplot2::scale_fill_manual(values = c("N2" = "orange", "CB4856" = "blue")) +
    ggplot2::geom_jitter(ggplot2::aes(x = genotype, y = pheno), size =2.5, alpha=.5, width = 0.35) +
    ggplot2::facet_wrap(~ marker, ncol = 5) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=14, face="bold", color="black"),
                   axis.text.y = ggplot2::element_text(size=12, face="bold", color="black"),
                   axis.title.x = ggplot2::element_text(size=14, face="bold", color="black", vjust=-.3),
                   axis.title.y = ggplot2::element_text(size=16, face="bold", color="black"),
                   strip.text.x = ggplot2::element_blank(),
                   strip.text.y = ggplot2::element_text(size=14, face="bold", color="black"),
                   plot.title = ggplot2::element_text(size=14, face="bold", vjust = 1),
                   legend.position="none",
                   panel.background = ggplot2::element_rect(color="black",size=1.2)) +
    ggplot2::ggtitle(peaks$trait[1]) +
    ggplot2::labs(x = "", y = "Regressed Length (µm)", title = " ")
}


########################
#### SUPPLEMENTAL FILES
########################

####################
####### TABLE S1
#####################
####### TABLE of QTLs with varexp > 0.05
#####################

TableCe <- benzpeaks %>%
  filter(var_exp >= 0.05) %>%
  select(sp.cond,chr, pos, trait, lod, threshold, var_exp, eff_size, ci_l_pos, ci_r_pos)
colnames(TableCe) <- c("Drug","Chr","Peak Pos", "Trait", "Peak LOD","LOD Threshold", "Variance Explained", "Effect Size", "Left", "Right")
write.csv(TableCe, file = "Ce_Linkage.csv")

#####################
####### S FIGURE: correlations for all drugs
#####################

phenoS2 <- readRDS("~/Dropbox/AndersenLab/LabFolders/Mostafa/Cel_Linkage_Mapping/RIAILs1_processed.rds")
bzs <- c("albendazole","fenbendazole-15","fenbendazole-30","thiabendazole-625","thiabendazole-125")
sizetraits <- c("median.TOF","median.EXT","mean.TOF","mean.EXT","q75.TOF","q75.EXT","q90.TOF","q90.EXT")
broodtraits <- c("n","norm.n")
fltraits <- c("mean.yellow","median.yellow","q75.yellow","q90.yellow","mean.red","median.red","q75.red","q90.red")
traitfilter <- c(sizetraits,broodtraits,fltraits)
phenoS2 <- phenoS2 %>%
  filter(condition %in% bzs, trait %in% traitfilter) %>%
  mutate(condition = sub("albendazole", "Albendazole (12.5 uM)",condition)) %>%
  mutate(condition = sub("fenbendazole-15", "Fenbendazole (15 uM)",condition)) %>%
  mutate(condition = sub("fenbendazole-30", "Fenbendazole (30 uM)",condition)) %>%
  mutate(condition = sub("thiabendazole-625", "Thiabendazole (62.5 uM)",condition)) %>%
  mutate(condition = sub("thiabendazole-125", "Thiabendazole (125 uM)",condition)) %>%
  select(condition,strain,trait,phenotype)

alb <- correlationplots(Drug = "Albendazole (12.5 uM)")
fen15 <- correlationplots(Drug = "Fenbendazole (15 uM)")
fen30 <- correlationplots(Drug = "Fenbendazole (30 uM)")
thia625 <- correlationplots(Drug = "Thiabendazole (62.5 uM)")
thia125 <- correlationplots(Drug = "Thiabendazole (125 uM)")

Ce_corr <- plot_grid(alb, ggdraw(), fen15, fen30, thia625, thia125, nrow=3, ncol=2)
Ce_corr
ggsave(plot = Ce_corr, filename = "Ce_Correlation.pdf", width = 12, height = 14)

correlationplots <- function(Drug) {
  phenoS2s <- phenoS2 %>%
    filter(condition == Drug) %>%
    select(-condition) %>%
    group_by(strain,trait)%>%
    mutate(phenotype2 = mean(phenotype)) %>%
    select(-phenotype) %>%
    distinct(strain,trait,.keep_all=TRUE) %>%
    spread(trait,phenotype2)

  correlation.df <- reshape2::melt(cor(phenoS2s[,3:ncol(phenoS2s)], use = "pairwise.complete.obs", method = "spearman"))

    corr_plot <- ggplot(data=correlation.df, aes(Var1,Var2, fill=value^2))+
    geom_tile() +
    theme_bw() +
    scale_fill_gradient(low="red", high="white", space="Lab",name=expression(bold(italic("r^2")))) +
    theme(axis.text.x = element_text(size=10, color="black", angle=60, hjust=1),
          axis.text.y = element_text(size=10, color="black"),
          panel.background = element_rect(fill = '#E5E8E8')) +
  labs(title=paste(Drug), x=NULL, y=NULL)
  corr_plot
}


#####################
####### S FIGURE: reduced lineplot of QTL, grouped by drug and divided by trait groups
#####################

#run Figure A first to generate benzpeaks
benzpeaksS3 <- benzpeaks %>%
  group_by(sp.cond,chr,traitgroup) %>%
  arrange(sp.cond,chr,traitgroup,ci_l_pos,ci_r_pos)
j <- 1
benzkpeaksS3$unique <- NA
for(i in 1:nrow(benzpeaksS3)){
  benzpeaksS3$unique[i] <- j
  ifelse((benzpeaksS3$ci_l_pos[i+1] - benzpeaksS3$ci_r_pos[i]) <= 1000000 & benzpeaksS3$chr[i] == benzpeaksS3$chr[i+1] & benzpeaksS3$traitgroup[i] == benzpeaksS3$traitgroup[i+1] & benzpeaksS3$sp.cond[i] == benzpeaksS3$sp.cond[i+1], j <- j, j <- j+1)
}
benzpeaksS3 <- benzpeaksS3 %>%
  group_by(unique) %>%
  arrange(desc(lod)) %>%
  distinct(unique,.keep_all=TRUE) %>%
  ungroup()

#set chromosome boundaries
newrows <- benzpeaksS3[1,] 
newrows[1,] = c(NA,"I",5000000,NA,0,NA,NA,NA,NA,NA,1,NA,14972282,"Albendazole","Albendazole","Length",NA)
newrows[2,] = c(NA,"II",5000000,NA,0,NA,NA,NA,NA,NA,1,NA,15173999,"Albendazole","Albendazole","Length",NA)
newrows[3,] = c(NA,"III",5000000,NA,0,NA,NA,NA,NA,NA,1,NA,13829314,"Albendazole","Albendazole","Length",NA)
newrows[4,] = c(NA,"IV",5000000,NA,0,NA,NA,NA,NA,NA,1,NA,17450860,"Albendazole","Albendazole","Length",NA)
newrows[5,] = c(NA,"V",5000000,NA,0,NA,NA,NA,NA,NA,1,NA,20914693,"Albendazole","Albendazole","Length",NA)
newrows[6,] = c(NA,"X",5000000,NA,0,NA,NA,NA,NA,NA,1,NA,17748731,"Albendazole","Albendazole","Length",NA)
newrows$ci_l_pos <- as.numeric(newrows$ci_l_pos)
newrows$ci_r_pos <- as.numeric(newrows$ci_r_pos)
newrows$pos <- as.numeric(newrows$pos)
newrows$lod <- as.numeric(newrows$lod)
#CE_I 14,972,282, #CE_II 15,173,999, #CE_III 13,829,314, #CE_IV 17,450,863, #CE_V 20,914,693, #CE_VI 17,748,731 


benzpeaksS3 <- benzpeaksS3 %>%
  mutate(sp.cond = sub("Fenbendazole_15", "Fenbendazole (15 μM)",sp.cond)) %>%
  mutate(sp.cond = sub("Fenbendazole_30", "Fenbendazole (30 μM)",sp.cond)) %>%
  mutate(sp.cond = sub("Thiabendazole_125", "Thiabendazole (125 μM)",sp.cond)) %>%
  mutate(sp.cond = sub("Thiabendazole_625", "Thiabendazole (62.5 μM)",sp.cond)) %>%
  mutate(sp.cond = sub("Albendazole", "Albendazole (12.5 μM)",sp.cond))
newrows <- newrows %>%
  mutate(sp.cond = sub("Albendazole", "Albendazole (12.5 μM)",sp.cond))

benzpeaksS3$sp.cond <- factor(benzpeaksS3$sp.cond, 
                              levels=c("Thiabendazole (125 μM)","Thiabendazole (62.5 μM)","Fenbendazole (30 μM)","Fenbendazole (15 μM)","Albendazole (12.5 μM)"), order = TRUE)
benzpeaksS3$traitgroup <- factor(benzpeaksS3$traitgroup, 
                                 levels=c("Pumping","Brood Size","Length"), order = TRUE)
newrows$sp.cond <- factor(newrows$sp.cond, 
                              levels=c("Thiabendazole (125 μM)","Thiabendazole (62.5 μM)","Fenbendazole (30 μM)","Fenbendazole (15 μM)","Albendazole (12.5 μM)"), order = TRUE)
newrows$traitgroup <- factor(newrows$traitgroup, 
                                 levels=c("Pumping","Brood Size","Length"), order = TRUE)

plot_S3 <- ggplot(benzpeaksS3)+
  aes(x=pos/1E6, y=traitgroup)+
  theme_bw() +
  scale_fill_viridis(name = "LOD") + scale_color_viridis(name = "LOD") +
  geom_segment(aes(x = ci_l_pos/1e6, y = traitgroup, xend = ci_r_pos/1e6, yend = traitgroup, color = lod), size = 2, alpha = 1) +
  geom_segment(data=newrows,aes(x = 0, y = traitgroup, xend = ci_r_pos/1e6, yend = traitgroup), size = 2, alpha = 0) +
  geom_point(aes(fill=lod),colour = "black",size = 5, alpha = 1, shape = 21)+
  theme(legend.title = element_text(size = 12, face = "bold")) +
  xlab("Genomic Position (Mb)") + ylab("") +
  theme(strip.background = element_rect(colour = "black", fill = "white",
                                        size = 0.75, linetype = "solid")) +
  theme(axis.text.x = element_text(size=10, face="bold", color="black"),
        axis.text.y = element_text(size=12, face="bold", color="black"),
        axis.title.x = element_text(size=14, face="bold", color="black"),
        axis.title.y = element_text(size=12, face="bold", color="black"),
        strip.text.x = element_text(size=14, face="bold", color="black"),
        strip.text.y = element_text(size=14, color="black", face="bold", angle = 0),
        plot.title = element_text(size=14, face="bold")) +
  facet_grid(sp.cond ~ chr, scales = "free", space = "free")
plot_S3
ggsave(plot = plot_S3, filename = "Ce_Linkage_Grouped.tiff", width = 15, height = 6)

