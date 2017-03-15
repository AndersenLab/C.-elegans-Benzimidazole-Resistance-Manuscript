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

setwd("~/Dropbox/Andersenlab/LabFolders/Mostafa/Manuscripts/Benzimidazoles/Github-BZ/Celegans_other/ben1")

########################
#### PLOT: Ben-1 mutant
########################

load("ben1.Rda")
regressed$strain <- factor(regressed$strain, 
                           levels=c("N2", "CB", "FX234"))
strainsplot <- c("N2", "CB", "FX234")

alb <- subset(regressed, condition=="Albendazole") %>%
  filter(trait == "q90.EXT")%>%
  filter(strain %in% strainsplot) 

plot_ben1 <-  ggplot(alb)+
  aes(x = strain, y = phenotype, fill = strain)+
  geom_boxplot(outlier.size= NA, alpha = 0.5, lwd=1)+
  geom_jitter(size =2.5, alpha=.5, width = 0.35)+
  theme_bw()+
  #facet_grid(.~assay) +
  theme(legend.position="none",
        panel.border = element_rect(size=1, colour = "black"))+
  labs(x = "Strain", y = "Regressed Length (µm)", title="Albendazole q75.EXT") +
  scale_fill_manual(values = c("N2" = "orange","CB" = "blue","ECA270" = "black","ECA271" = "black","ECA286" = "black","FX234" = "black")) +
  theme(axis.text.x = element_text(size=14, face="bold", color="black"),
        axis.text.y = element_text(size=12, face="bold", color="black"),
        axis.title.x = element_text(size=0, face="bold", color="black"),
        axis.title.y = element_text(size=16, face="bold", color="black"),
        strip.text.x = element_text(size=16, face="bold", color="black"),
        strip.text.y = element_text(size=16, face="bold", color="black"),
        plot.title = element_text(size=0, face="bold")) +
  scale_x_discrete(labels=c("N2" = expression(paste(bold("N2"))), 
                            "CB" = expression(paste(bold("CB4856"))),
                            "FX234" = expression(paste(bolditalic("ben-1")))))
plot_ben1

ggsave(plot = plot_ben1, filename = "ben1.tiff", width = 7, height = 5)


#ben-1 only
## get q75 stats for ben-1 vs N2-CB difference
#albN2 <- alb %>% filter(strain == "N2")
#albN2 <- albN2$phenotype
#albCB <- alb %>% filter(strain == "CB") 
#albCB <- albCB$phenotype
#t.test(albN2,albCB)

#33 - (-58) = 91
#139 - (-58) = 197