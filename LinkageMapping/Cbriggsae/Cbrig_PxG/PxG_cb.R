library(stringr)
library(dplyr)
library(ggplot2)
library(tidyr)

setwd("~/Dropbox/Andersen\ Lab/LabFolders/Mostafa/Cbrig_PxG/")

load("~/Dropbox/Andersen\ Lab/LabFolders/Mostafa/Cbrig_PxG/Lgeno_cb.Rdata")

#no parental
load("~/Dropbox/Andersen\ Lab/LabFolders/Mostafa/Cbrig_PxG/Lphen_cb.Rdata")

#code to plot RIAIL genotypes
Lgeno_cb2 <- Lgeno_cb %>%
  separate(chr_pos, into = c("chr", "pos"), sep = "_")
Lgeno_cb2$chr <- as.integer(Lgeno_cb2$chr)
Lgeno_cb2$pos <- as.numeric(Lgeno_cb2$pos)

plot <- Lgeno_cb2 %>%  #filter(chr == 4) %>%
  ggplot(.,aes(x=pos/1e6,y=strain, colour = as.factor(geno), group = strain)) +
  geom_line(size = 1.5) + # geom_point(shape = 21, size = 0.75) +
  facet_grid(.~chr, scales = "free_x")
plot

#plot frequency at each position (group by chr,pos,geno, then summarize n())
plot2 <- Lgeno_cb2 %>%
  group_by(chr,pos,geno) %>% 
  summarise(count = n()/167) %>%
  filter(geno == 2) %>%
  ggplot(.,aes(x=pos/1e6,y=count*100, colour = as.factor(geno), group = as.factor(geno))) + 
  geom_line(size = 1.5, alpha = 0.8) + geom_point(size = 1.5) +
  facet_grid(.~chr, scales = "free_x") + ylab("% geno 2") + ylim(0,100)
plot2

plot2b <- Lgeno_cb2 %>%
  group_by(chr,pos,geno) %>% 
  summarise(count = n()/167) %>%
  filter(geno == 2) %>%
  ggplot(.,aes(x=pos/1e6,y=count*100, colour = as.factor(chr), group = as.factor(chr))) + 
  geom_line(size = 1.75, alpha = 0.9) + #geom_point(size = 1.5) +
  ylab("% geno 2") + ylim(0,100)
plot2b


#### playing around ... recombination hotspots

Lgeno_temp <- Lgeno_cb2 %>%
  group_by(strain,chr) %>%
  arrange(pos) %>%
  mutate(chg = 0)

for(i in 2:nrow(Lgeno_temp)) {
  Lgeno_temp$chg[i] = ifelse(Lgeno_temp$strain[i] == Lgeno_temp$strain[i-1] & Lgeno_temp$chr[i] == Lgeno_temp$chr[i-1] & abs(Lgeno_temp$geno[i] - Lgeno_temp$geno[i-1]) > 0, "1", "0")
}

plot <- Lgeno_temp %>%  #filter(chr == 4) %>%
  ggplot(.,aes(x=pos/1e6,y=chg, colour = as.factor(chg), group = strain)) +
  geom_jitter(data=subset(Lgeno_temp, Lgeno_temp$chg > 0),size = 2.5, alpha = 0.5, position = position_jitter(width = 0)) +
  facet_grid(.~chr, scales = "free_x")  + ylab("recombination")
plot

plot2 <- Lgeno_temp %>%
  group_by(chr,pos,chg) %>% 
  summarise(count = n()/167) %>%
  filter(chg == 1) %>%
  ggplot(.,aes(x=pos/1e6,y=count*100, colour = as.factor(chg), group = as.factor(chg))) + 
  geom_line(size = 1, alpha = 0.8) + geom_point(size = 1) +
  facet_grid(.~chr, scales = "free_x") + ylab("# of recombination events")
plot2


#### plot splits

PxGlinkage_cb("4_2777733","fenbendazole_resid.n")
PxGlinkage_cb("6_1909385","resid.median.norm.yellow")

PxGlinkage_cb <- function(peakMarker, condtrt){
  temp <- str_split_fixed(condtrt,pattern="_",n=2)
  trt <- temp[,2]
  condition <- temp[,1]

  gen <- Lgeno_cb %>%
    filter(chr_pos == peakMarker)
  
  posit <- gen[1,5] #?? 
  
  Lphen_cb%>%
    filter(grepl(condition,trait))%>%
    filter(grepl(trt,trait))%>%
    left_join(.,gen,by="strain")%>%
    mutate(color = ifelse(strain=="AF16", "AF16",
                          ifelse(strain == "HK104", "HK104",
                                 ifelse(geno == 1, "RIAILs-AF16",
                                        ifelse(geno == 2, "RIAILs-HK104",0)))),
           label = ifelse(strain=="AF16", "AF16",
                          ifelse(strain == "HK104", "HK104",
                                 ifelse(geno == 1, "RIAILs-AF16",
                                        ifelse(geno == 2, "RIAILs-HK104",0)))))%>%
    filter(!is.na(color))%>%
    ggplot(.)+
    aes(x = color, y = value,fill = ifelse(strain=="AF16", "AF16",
                                           ifelse(strain == "HK104", "HK104",
                                                  ifelse(geno == 1, "RIAILs-AF16",
                                                         ifelse(geno == 2, "RIAILs-HK104",0)))),
        color = ifelse(strain=="AF16", "AF16",
                       ifelse(strain == "HK104", "HK104",
                              ifelse(geno == 1, "RIAILs-AF16",
                                     ifelse(geno == 2, "RIAILs-HK104",0)))))+
    scale_fill_manual(values = c("AF16"="orange","HK104"="blue","RIAILs-HK104"="blue","RIAILs-AF16"="orange"), name = "Genotype")+
    scale_color_manual(values = c("AF16"="black","HK104"="black","RIAILs-HK104"="#666666","RIAILs-AF16"="#666666"),name="Genotype")+
    geom_boxplot(outlier.shape=NA)+
    theme_bw()+
    theme(axis.text.x = element_text(size=16, face="bold", color="black"),
          axis.text.y = element_text(size=16, face="bold", color="black"),
          axis.title.x = element_text(size=20, face="bold", color="black"),
          axis.title.y = element_text(size=20, face="bold", color="black",vjust=1),
          strip.text.x = element_text(size=20, face="bold", color="black"),
          strip.text.y = element_text(size=20, face="bold", color="black"),
          plot.title = element_text(size=24, face="bold",vjust=1),
          legend.title = element_text(size=14),
          panel.border = element_rect(size=2))+
    #ylim(-1,1) +
    geom_jitter(alpha=.75)+
    labs(title=paste( trt,peakMarker,posit,sep="_"),y=trt, x="Genotype" )
}

