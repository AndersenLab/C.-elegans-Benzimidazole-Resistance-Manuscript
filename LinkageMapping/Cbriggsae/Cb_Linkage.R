library("linkagemapping")
library("easysorter")
library(dplyr)
library(ggplot2)
library(stringr)
#library(ggrepel)
library(gridExtra)
library(cowplot)
library(viridis)
options(stringsAsFactors = FALSE)

setwd("LinkageMapping/Cbriggsae")
load("Cb_BZ_GWER.Rda")

cb_mapping_gwer <- Cb_BZ_GWER %>%
  mutate(pos2 = as.numeric(str_split_fixed(marker, "\\_", 2)[,2])) %>%
  mutate(ci_l_pos2 = as.numeric(str_split_fixed(ci_l_marker, "\\_", 2)[,2])) %>%
  mutate(ci_r_pos2 = as.numeric(str_split_fixed(ci_r_marker, "\\_", 2)[,2]))
  
benzpeaks <- cb_mapping_gwer %>%
  filter (var_exp > 0.05) %>%
  filter(!grepl("cv|iqr|q10|q25|var|norm.yellow|norm.red|norm.EXT|green|f.L1|f.L2L3|f.L4|f.ad", trait)) %>%
  mutate(sp.cond = factor(str_split_fixed(trait, "\\.", 2)[,1], levels=c("abamectin","albendazole", "fenbendazole","thiabendazole"))) %>%
  filter(sp.cond != "abamectin") %>%
  mutate(sp.cond = sub("albendazole", "Albendazole",sp.cond)) %>%
  mutate(sp.cond = sub("fenbendazole", "Fenbendazole",sp.cond)) %>%
  mutate(sp.cond = sub("thiabendazole", "Thiabendazole",sp.cond)) %>%
  group_by(chr,trait,lod) %>% distinct(trait,.keep_all=TRUE) %>% ungroup() %>%
  mutate(trait = sub("albendazole.","",trait)) %>%
  mutate(trait = sub("thiabendazole.","",trait)) %>%
  mutate(trait = sub("fenbendazole.","",trait))
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
  arrange(sp.cond,chr,ci_l_pos2,ci_r_pos2)
j <- 1
benzpeaksR$unique <- NA
for(i in 1:nrow(benzpeaksR)){
  benzpeaksR$unique[i] <- j
  ifelse((benzpeaksR$ci_l_pos2[i+1] - benzpeaksR$ci_r_pos2[i]) <= 1000000 &  benzpeaksR$chr[i] == benzpeaksR$chr[i+1] &  benzpeaksR$sp.cond[i] == benzpeaksR$sp.cond[i+1], j <- j, j <- j+1)
}
benzpeaksR <- benzpeaksR %>%
  group_by(unique) %>%
  arrange(desc(lod)) %>%
  distinct(unique,.keep_all=TRUE)
# for faceting issues
benzpeaksR <- benzpeaksR %>%
  mutate(height = "") %>%
  mutate(height = ifelse(sp.cond == "Albendazole","25 μM",height)) %>%
  mutate(height = ifelse(sp.cond == "Fenbendazole","30 μM",height)) %>%
  mutate(height = ifelse(sp.cond == "Thiabendazole","40 μM",height)) %>%
  mutate(drug = sp.cond) 
benzpeaksR$drug <- factor(benzpeaksR$drug, 
                          levels=c("Thiabendazole","Fenbendazole","Albendazole"))
#change to roman numerals
benzpeaksR$chr <- factor(benzpeaksR$chr, levels = c("1","2","3","4","5","6"), ordered = T)
levels(benzpeaksR$chr) <- c("I","II","III","IV","V","X")

#set chromosome boundaries (use last marker as estimate)
newrows <- benzpeaksR[1,] 
newrows[1,] = c(NA,"I",NA,NA,0,NA,NA,NA,NA,NA,NA,NA,NA,5000000,1,15400454,"Albendazole","Length",NA,"25 μM","Albendazole")
newrows[2,] = c(NA,"II",NA,NA,0,NA,NA,NA,NA,NA,NA,NA,NA,5000000,1,16536341,"Albendazole","Length",NA,"25 μM","Albendazole")
newrows[3,] = c(NA,"III",NA,NA,0,NA,NA,NA,NA,NA,NA,NA,NA,5000000,1,14576522,"Albendazole","Length",NA,"25 μM","Albendazole")
newrows[4,] = c(NA,"IV",NA,NA,0,NA,NA,NA,NA,NA,NA,NA,NA,5000000,1,17269500,"Albendazole","Length",NA,"25 μM","Albendazole")
newrows[5,] = c(NA,"V",NA,NA,0,NA,NA,NA,NA,NA,NA,NA,NA,5000000,1,19209427,"Albendazole","Length",NA,"25 μM","Albendazole")
newrows[6,] = c(NA,"X",NA,NA,0,NA,NA,NA,NA,NA,NA,NA,NA,5000000,1,21485079,"Albendazole","Length",NA,"25 μM","Albendazole")
newrows$ci_l_pos2 <- as.numeric(newrows$ci_l_pos2)
newrows$ci_r_pos2 <- as.numeric(newrows$ci_r_pos2)
newrows$pos2 <- as.numeric(newrows$pos2)
newrows$lod <- as.numeric(newrows$lod)

####################################
### PLOT reduced set
####################################

plot_benz <- ggplot(benzpeaksR)+
  aes(x=pos2/1E6, y=height)+
  theme_bw() +
  scale_fill_viridis(name = "LOD") + scale_color_viridis(name = "LOD") +
  geom_segment(aes(x = ci_l_pos2/1E6, y = height, xend = ci_r_pos2/1E6, yend = height, color = lod), size = 2, alpha = 1) +
  geom_segment(data=newrows,aes(x = ci_l_pos2/1E6, y = height, xend = ci_r_pos2/1E6, yend = height), size = 2, alpha = 0) +
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
### PLOT B / C (Linkage plot fenbendazole.n + split at peak marker)
####################################

trait_input <- "fenbendazole.n"
linkplot <- filter(cb_mapping_gwer, trait == trait_input) %>%
  filter(iteration == 1)
peaks <- linkplot %>%
  filter(var_exp > 0)
linkplot <- mutate(linkplot,LOD2 = 0)
for(i in 1:nrow(peaks)){
  linkplot$LOD2 <- ifelse(linkplot$pos2 >= peaks$ci_l_pos2[i] & linkplot$pos2 <= peaks$ci_r_pos2[i] & linkplot$chr == peaks$chr[i], linkplot$lod, linkplot$LOD2)}

linkplot$chr <- factor(linkplot$chr, levels = c("1","2","3","4","5","6"), ordered = T)
levels(linkplot$chr) <- c("I","II","III","IV","V","X")

plot <- ggplot(data=linkplot,aes(x=pos2/1e6, y = lod, group = factor(iteration))) +
  geom_line(aes(colour = factor(iteration)), size = 1) +
  theme_bw() +
  geom_hline(aes(yintercept=2.7513), color ="#F8766D", linetype="dashed", size=0.70) +
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
data("AF16xHK104cross")
cross <- AF16xHK104cross
load("Github-BZ-prep/LinkageMapping/Cbriggsae/cb_regressed.Rda")
pheno <- regressed
mapcross <- mergepheno(cross, pheno)
pgtrait <- filter(cb_mapping_gwer, trait == "fenbendazole.n")
pxgplot <- pxgplot2(cross = mapcross, map = pgtrait)
pxgplot

#mapcross$geno #167 strains

## numbers for paper (average brood size difference)
# broodavg <- regressed %>%
#   filter(strain == "AF16" | strain == "HK104", trait == "norm.n") %>%
#   group_by(strain) %>%
#   mutate(mean_phenotype = mean(phenotype))

########################
#### FINAL Cbriggsae LINKAGE PLOT FIGURE (A,B,C)
########################

bottom_row <- plot_grid(plot, pxgplot, labels = c('B', 'C'), rel_widths = c(2.7, 1), scale = 0.95)
bottom_row
final_plot <- plot_grid(plot_benz, bottom_row, nrow=2, labels = c('A',''), rel_heights = c(0.9,1))
final_plot
ggsave(plot = final_plot, filename = "Cb_Linkage.tiff", width = 15, height = 8)

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

pxgplot2 <- function(cross, map, parent="AF16xHK104") {
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
  geno <- data.frame(extract_genotype(cross))
  colnames(geno) <- gsub("X", "", colnames(geno)) #Mostafa fix
  geno <- geno %>%
    dplyr::select(which(colnames(.) %in% uniquemarkers)) %>%
    data.frame(., pheno)
  
  colnames(geno)[1:(ncol(geno)-1)] <- sapply(colnames(geno)[1:(ncol(geno)-1)],
                                             function(marker) {paste(unlist(peaks[peaks$marker == gsub("\\.","-",marker),c("chr", "pos")]),collapse = ":")})
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
                           levels=c("AF16","HK104"))
  ggplot2::ggplot(split) +
    ggplot2::geom_boxplot(ggplot2::aes(x = genotype, y = pheno, fill = genotype), outlier.size= NA, alpha = 0.5, lwd=1) +
    ggplot2::scale_fill_manual(values = c("AF16" = "indianred", "HK104"= "gold")) +
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
    ggplot2::labs(x = "", y = "Brood size", title = " ")
}



########################
#### SUPPLEMENTAL FILES
########################

####################
####### TABLE of QTLs with varexp > 0.05
#####################

TableCb <- benzpeaks %>%
  filter(var_exp >= 0.05) %>%
  separate(marker, into = c("Chr", "Pos"), sep = "_") %>%
  select(sp.cond,chr,Pos,trait,lod,threshold,var_exp,eff_size,ci_l_pos2, ci_r_pos2)
colnames(TableCb) <- c("Drug","Chr","Peak Pos", "Trait", "Peak LOD","LOD Threshold", "Variance Explained", "Effect Size", "Left", "Right")
TableCb$Chr <- factor(TableCb$Chr, levels = c("1","2","3","4","5","6"), ordered = T)
levels(TableCb$Chr) <- c("I","II","III","IV","V","X")
write.csv(TableCb, file = "Cb_Linkage.csv")

#####################
####### S FIGURE: correlations for all drugs
#####################

load("Github-BZ-prep/LinkageMapping/Cbriggsae/cb_regressed.Rda")
phenoS5 <- regressed %>% ungroup()
bzs <- c("albendazole","fenbendazole","thiabendazole")
sizetraits <- c("median.TOF","median.EXT","mean.TOF","mean.EXT","q75.TOF","q75.EXT","q90.TOF","q90.EXT")
broodtraits <- c("n","norm.n")
fltraits <- c("mean.yellow","median.yellow","q75.yellow","q90.yellow","mean.red","median.red","q75.red","q90.red")
traitfilter <- c(sizetraits,broodtraits,fltraits)
phenoS5 <- phenoS5 %>%
  filter(condition %in% bzs, trait %in% traitfilter) %>%
  mutate(condition = sub("albendazole", "Albendazole (25 uM)",condition)) %>%
  mutate(condition = sub("fenbendazole", "Fenbendazole (30 uM)",condition)) %>%
  mutate(condition = sub("thiabendazole", "Thiabendazole (40 uM)",condition)) %>%
  select(condition,strain,trait,phenotype)

alb <- correlationplots(Drug = "Albendazole (25 uM)")
fen <- correlationplots(Drug = "Fenbendazole (30 uM)")
thia <- correlationplots(Drug = "Thiabendazole (40 uM)")

Cb_corr <- plot_grid(alb, ggdraw(), fen, ggdraw(), thia, ggdraw(), nrow=3, ncol=2)
Cb_corr
ggsave(plot = Cb_corr, filename = "Cb_Correlation.pdf", width = 12, height = 14)

correlationplots <- function(Drug) {
  phenoS5s <- phenoS5 %>%
    filter(condition == Drug) %>%
    select(-condition) %>%
    group_by(strain,trait)%>%
    mutate(phenotype2 = mean(phenotype)) %>%
    select(-phenotype) %>%
    distinct(strain,trait,.keep_all=TRUE) %>%
    spread(trait,phenotype2)
  
  correlation.df <- reshape2::melt(cor(phenoS5s[,3:ncol(phenoS5s)], use = "pairwise.complete.obs", method = "spearman"))
  
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
benzpeaksS6 <- benzpeaks %>%
  group_by(sp.cond,chr,traitgroup) %>%
  arrange(sp.cond,chr,traitgroup,ci_l_pos,ci_r_pos)
j <- 1
benzkpeaksS6$unique <- NA
for(i in 1:nrow(benzpeaksS6)){
  benzpeaksS6$unique[i] <- j
  ifelse((benzpeaksS6$ci_l_pos[i+1] - benzpeaksS6$ci_r_pos[i]) <= 1000000 & benzpeaksS6$chr[i] == benzpeaksS6$chr[i+1] & benzpeaksS6$traitgroup[i] == benzpeaksS6$traitgroup[i+1] & benzpeaksS6$sp.cond[i] == benzpeaksS6$sp.cond[i+1], j <- j, j <- j+1)
}
benzpeaksS6 <- benzpeaksS6 %>%
  group_by(unique) %>%
  arrange(desc(lod)) %>%
  distinct(unique,.keep_all=TRUE) %>%
  ungroup()

benzpeaksS6$chr <- factor(benzpeaksS6$chr, levels = c("1","2","3","4","5","6"), ordered = T)
levels(benzpeaksS6$chr) <- c("I","II","III","IV","V","X")
#set chromosome boundaries
newrows <- benzpeaksS6[1,] 
newrows[1,] = c(NA,"I",NA,NA,0,NA,NA,NA,NA,NA,NA,NA,NA,5000000,1,15400454,"Albendazole","Length",NA)
newrows[2,] = c(NA,"II",NA,NA,0,NA,NA,NA,NA,NA,NA,NA,NA,5000000,1,16536341,"Albendazole","Length",NA)
newrows[3,] = c(NA,"III",NA,NA,0,NA,NA,NA,NA,NA,NA,NA,NA,5000000,1,14576522,"Albendazole","Length",NA)
newrows[4,] = c(NA,"IV",NA,NA,0,NA,NA,NA,NA,NA,NA,NA,NA,5000000,1,17269500,"Albendazole","Length",NA)
newrows[5,] = c(NA,"V",NA,NA,0,NA,NA,NA,NA,NA,NA,NA,NA,5000000,1,19209427,"Albendazole","Length",NA)
newrows[6,] = c(NA,"X",NA,NA,0,NA,NA,NA,NA,NA,NA,NA,NA,5000000,1,21485079,"Albendazole","Length",NA)
newrows$ci_l_pos2 <- as.numeric(newrows$ci_l_pos2)
newrows$ci_r_pos2 <- as.numeric(newrows$ci_r_pos2)
newrows$pos2 <- as.numeric(newrows$pos2)
newrows$lod <- as.numeric(newrows$lod)

benzpeaksS6 <- benzpeaksS6 %>%
  mutate(sp.cond = sub("Fenbendazole", "Fenbendazole (30 μM)",sp.cond)) %>%
  mutate(sp.cond = sub("Thiabendazole", "Thiabendazole (40 μM)",sp.cond)) %>%
  mutate(sp.cond = sub("Albendazole", "Albendazole (25 μM)",sp.cond))

newrows <- newrows %>%
  mutate(sp.cond = sub("Albendazole", "Albendazole (25 μM)",sp.cond))

benzpeaksS6$sp.cond <- factor(benzpeaksS6$sp.cond, 
                              levels=c("Thiabendazole (40 μM)","Fenbendazole (30 μM)","Albendazole (25 μM)"), order = TRUE)
benzpeaksS6$traitgroup <- factor(benzpeaksS6$traitgroup, 
                                 levels=c("Pumping","Brood Size","Length"), order = TRUE)
newrows$sp.cond <- factor(newrows$sp.cond, 
                          levels=c("Thiabendazole (40 μM)","Fenbendazole (30 μM)","Albendazole (25 μM)"), order = TRUE)
newrows$traitgroup <- factor(newrows$traitgroup, 
                             levels=c("Pumping","Brood Size","Length"), order = TRUE)

plot_S6 <- ggplot(benzpeaksS6)+
  aes(x=pos2/1E6, y=traitgroup)+
  theme_bw() +
  scale_fill_viridis(name = "LOD") + scale_color_viridis(name = "LOD") +
  geom_segment(aes(x = ci_l_pos2/1e6, y = traitgroup, xend = ci_r_pos2/1e6, yend = traitgroup, color = lod), size = 2, alpha = 1) +
  geom_segment(data=newrows,aes(x = 0, y = traitgroup, xend = ci_r_pos2/1e6, yend = traitgroup), size = 2, alpha = 0) +
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
plot_S6
ggsave(plot = plot_S6, filename = "Cb_Linkage_Grouped.tiff", width = 15, height = 6)

### NEW BRIGGSAE FIGURE

bottom_row <- plot_grid(plot, pxgplot, labels = c('B', 'C'), rel_widths = c(2.7, 1), scale = 0.95)
bottom_row
final_plot <- plot_grid(plot_S6, bottom_row, nrow=2, labels = c('A',''), rel_heights = c(0.9,1))
final_plot
ggsave(plot = final_plot, filename = "Cb_Linkage2.tiff", width = 15, height = 8)


#####################
####### S Table: C. briggsae variants for fenbendazole interval
#####################

### PREPARE VARIANT FILE FOR BRIGGSAE FENBENDAZOLE QTL INTERVAL
setwd("LinkageMapping/Cbriggsae")

region <- snpeff2("IV:2563732-3447467", 
                  severity = c("LOW", "MODERATE", "HIGH", "MODIFIER"), 
                  elements = c("CDS", "five_prime_UTR", "exon","intron", "three_prime_UTR"),
                  remote = FALSE,
                  impute = FALSE,
                  vcf = "~/Dropbox/Andersenlab/LabFolders/Mostafa/Manuscripts/Benzimidazoles/Briggsae/Variants/briggsae.WS256.snpeff.vcf.gz")
Briggsae_vars <- region %>%
  dplyr::filter(strain == "HK104", GT!="REF")
save(Briggsae_vars, file="Briggsae_vars.Rda")

#load function below before running code above
snpeff2 <- function(...,
                    severity = c("HIGH","MODERATE"),
                    elements = c("exon"),
                    long = TRUE,
                    remote = FALSE,
                    impute = FALSE,
                    vcf = NA) {
  
  regions <- unlist(list(...))
  
  # Allow user to specify 'ALL'
  if ("ALL" %in% severity) {
    severity <-  c("LOW", "MODERATE", "HIGH", 'MODIFIER')
  }
  if ("ALL" %in% elements) {
    elements <- c("CDS", "five_prime_UTR", "exon", "intron", "three_prime_UTR")
  }
  
  # Ensure that bcftools is available:
  bcftools_version <- as.double(stringr::str_extract(readLines(pipe("bcftools --version"))[1], "[0-9]+\\.[0-9]+"))
  if(is.na(bcftools_version) | bcftools_version < 1.2) {
    stop("bcftools 1.2+ required for this function")
  }
  
  results <- suppressWarnings(lapply(regions, function(query) {
    # Save region as query
    
    # Fix region specifications
    query <- gsub("\\.\\.", "-", query)
    query <- gsub(",", "", query)
    
    # Resolve region names
    if (!grepl("(I|II|III|IV|V|X|MtDNA).*", query)) {
      elegans_gff <- get_db()
      # Pull out regions by element type.
      region <- paste((dplyr::bind_rows(lapply(elements, function(e) {
        dplyr::collect(dplyr::filter(elegans_gff, locus == query | gene_id == query | sequence_name == query, type_of == e) %>%
                         dplyr::select(chrom, start, end, gene_id, biotype, type_of, locus, sequence_name) %>%
                         dplyr::distinct( .keep_all = TRUE))
      })) %>%
        dplyr::summarize(chrom = chrom[1], start = min(start), end = max(end)) %>%
        dplyr::mutate(region_format = paste0(chrom, ":", start, "-", end)) %>%
        dplyr::select(region_format) %>%
        dplyr::distinct(.keep_all = TRUE))$region_format, collapse = ",")
      if (stringr::str_length(regions[[1]]) == 0) {
        message(paste0(query, " not found."))
        region <- NA
      }
    } else {
      region <- query
    }
    
    if(is.na(vcf)) {
      vcf_path <- get_vcf(remote = remote, impute = impute)
    } else {
      vcf_path <- vcf
      gene_ids <- NA
    }
    
    sample_names <- readr::read_lines(suppressWarnings(pipe(paste("bcftools","query","-l",vcf_path))))
    base_header <- c("CHROM", "POS", "REF","ALT","FILTER")
    ANN_header = c("allele", "effect", "impact",
                   "gene_name", "gene_id", "feature_type", 
                   "feature_id", "transcript_biotype","exon_intron_rank",
                   "nt_change", "aa_change", "cDNA_position/cDNA_len", 
                   "protein_position", "distance_to_feature", "error", "extra")
    # If using long format provide additional information.
    format <- "'%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t%ANN[\\t%TGT]\\n'"
    if (long == T) {
      format <- "'%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t%ANN[\\t%TGT!%FT!%DP!%DP4!%SP!%HP]\\n'"
    }
    command <- paste("bcftools","query","--regions", region, "-f", format ,vcf_path)
    if (!is.na(region)) {
      message(paste0("Query: ", query, "; region - ", region, "; "))
      result <- try(dplyr::tbl_df(data.table::fread(command, col.names = c(base_header, "ANN", sample_names ), sep = "\t")), silent = TRUE)
      if(!grepl("^Error.*", result[[1]][1])) {
        tsv <- result %>%
          dplyr::mutate(REF = ifelse(REF==TRUE, "T", REF), # T nucleotides are converted to 'true'
                        ALT = ifelse(ALT==TRUE, "T", ALT))
      } else {
        tsv <- as.data.frame(NULL)
      }
      # If no results are returned, stop.
      if (typeof(tsv) == "character" | nrow(tsv) == 0) {
        warning("No Variants")
        NA
      } else {
        tsv <-  dplyr::mutate(tsv, ANN=strsplit(ANN,",")) %>%
          tidyr::unnest(ANN) %>%
          tidyr::separate(ANN, into = ANN_header, sep = "\\|") %>%
          dplyr::select(one_of(c(base_header, ANN_header)), everything(), -extra) %>%
          dplyr::mutate(gene_name = as.character(gene_ids[gene_name])) %>%
          dplyr::mutate(query = query, region = region) %>%
          dplyr::select(CHROM, POS, query, region, everything())
        
        tsv <-  dplyr::filter(tsv, impact %in% severity) 
        if (nrow(tsv) == 0) {
          message(paste("No Results for", region, "after filtering"))
        }
        if (long == FALSE) {
          tsv
        } else {
          tsv <- tidyr::gather_(tsv, "strain", "GT", names(tsv)[23:length(tsv)])  %>%
            tidyr::separate(GT, into=c("a1","a2", "FT", "DP", "DP4", "SP", "HP"), sep="/|\\||\\!", remove=T) %>%
            dplyr::mutate(a1=ifelse(a1 == ".", NA, a1)) %>%
            dplyr::mutate(a2=ifelse(a2 == ".", NA, a2)) %>%
            dplyr::mutate(GT = NA) %>%
            dplyr::mutate(GT = ifelse(a1 == REF & a2 == REF & !is.na(a1), "REF",GT)) %>%
            dplyr::mutate(GT = ifelse(a1 != a2 & !is.na(a1), "HET",GT)) %>%
            dplyr::mutate(GT = ifelse(a1 == a2 & a1 != REF & !is.na(a1), "ALT",GT)) %>%
            dplyr::select(CHROM, POS, strain, REF, ALT, a1, a2, GT, FT, FILTER, DP, DP4, SP, HP, everything()) %>%
            dplyr::arrange(CHROM, POS) 
        }
        tsv
      }
    }
  }
  ))
  results <- do.call(rbind, results)
  if (!"CHROM" %in% names(results)){
    stop("No Results")
  }
  results <- results %>% dplyr::filter(!is.na(CHROM))
  results
}

### S TABLE Variants for QTL interval
load("LinkageMapping/Cbriggsae/QTLvariants/Briggsae_vars.Rda")
Briggsae_vars <- Briggsae_vars %>%
  filter(FILTER == "PASS") %>%
  filter(transcript_biotype == "protein_coding") %>%
  filter(impact == "HIGH" | impact == "MODERATE") %>%
  distinct(POS, .keep_all = TRUE) %>%
  select(CHROM,POS,REF,ALT,effect,impact,gene_id,feature_id,exon_intron_rank,nt_change,aa_change,protein_position)
colnames(Briggsae_vars) <- c("Chr","Pos", "Ref", "Alt","Effect","Impact","Gene_ID", "Feature_ID", "Exon_Intron_Rank", "NT_change","AA_change","Prot_position")
write.csv(Briggsae_vars, file = "Cb_QTL_vars.csv")

#unique(Briggsae_vars$Gene_ID) 275 variants / 84 unique genes