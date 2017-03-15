library(easysorter)
library(dplyr)
library(ggplot2)
library(tidyr)
library(R.utils)

#### Calculate and Plot Heritabilities
setwd("DoseResponses/Celegans")
load("Ce_BZ_DR.Rda")
source("DoseResponses/Celegans/Dose_Response_Fxns.R")
all_herits <- list()
dr_final <- data.frame(strain = character(), phenotype = numeric(), trait = character(), cond_dose = numeric(),stringsAsFactors=FALSE)

for(j in 1:length(unique(dr_df$trait))){ #for each trait
  herit_df <- dplyr::select(dr_df, strain, condition, dose_µM, phenotype, trait) %>%
    dplyr::filter(trait == unique(dr_df$trait)[j])
  herit_df$phenotype_r <- NA
  
  for (k in unique(herit_df$strain)){ #for each strain
    for (l in unique(herit_df$condition)){ #for each condition
      zero <- herit_df %>% filter(strain == k,condition == l,dose_µM == 0)
      zero <- mean(zero$phenotype) #0 dose data for trait,strain,condition
      herit_df <- herit_df %>%
        arrange(strain,condition,dose_µM) %>% 
        group_by(strain,condition,dose_µM) %>%
        mutate(phenotype_r = ifelse(strain == k & condition == l,phenotype - zero,phenotype_r)) %>% 
        ungroup()
    }
  }
  herit_df <- herit_df %>%
    dplyr::mutate(cond_dose = paste(condition, dose_µM, sep = "_")) %>%
    dplyr::select(-condition, -dose_µM) 
  
  #regressed vs unregressed phenotype (run for regressed)
  herit_df$phenotype <- herit_df$phenotype_r    
  herit_df <- herit_df %>% select(-phenotype_r)
  
  dr_final <- rbind(dr_final,herit_df)
  
  herit_list <- list()
  
  for(i in 1:length(unique(herit_df$cond_dose))){
    temp <- dplyr::filter(herit_df, cond_dose == unique(herit_df$cond_dose)[i]) %>%
      arrange(strain)
    h2 <- H2.calc(temp, boot = F)
    h2$phenotype <- unique(herit_df$cond_dose)[i]
    h2$trait <- as.character(unique(dr_df$trait)[j])
    herit_list[[i]] <- h2
  }
  
  all_herits[[j]] <- rbind_all(herit_list)
}

H2 <- dplyr::bind_rows(all_herits)
H2a <- tidyr::separate(H2, phenotype, into = c("condition", "dose_µM"), sep = "_")
H2a$dose_µM <- as.numeric(H2a$dose_µM)

dr_final <- dr_final %>% separate(cond_dose, into = c("condition", "dose_µM"), sep = "_")
dr_final$dose_µM <- as.numeric(dr_final$dose_µM)

dr_H2 <- left_join(dr_final,H2a, by = c("condition", "dose_µM", "trait")) %>%
  group_by(strain,trait,condition,dose_µM) %>%
  mutate(mean_phenotype = mean(phenotype)) %>%
  dplyr::filter(grepl("q75.EXT|norm.n|q75.yellow", trait))
  
plotdf <- dr_H2 %>%
  mutate(h2_labl = ifelse(condition == "Albendazole" & dose_µM == 12.5 & strain == "N2", "plot",
                          ifelse(condition == "Fenbendazole" & dose_µM == 30 | dose_µM == 15  & strain == "N2", "plot",
                                        ifelse(condition == "Thiabendazole" & dose_µM == 62.5 | dose_µM == 125 & strain == "N2", "plot",NA))))
plotdf$dose_µM <- factor(plotdf$dose_µM, ordered = T) 
plotdf$trait <- factor(plotdf$trait, levels = c("q75.EXT","norm.n","q75.yellow"), ordered = T)
levels(plotdf$trait) <- c("Length","Brood Size","Pharyngeal Pumping")

#DR Plot
plot <- ggplot(plotdf, aes(dose_µM,mean_phenotype, colour = strain)) +
  geom_point() + 
  geom_line(aes(group=strain), size = 1, alpha = 0.5) +
  theme_bw() +
  labs(y = "Δ Phenotype Value", x = "Dose (µM)", title = "") + 
  scale_color_manual(values = c("N2" = "orange", "CB4856" = "blue", "DL238" = "green", "JU258" = "cyan"))+
  facet_grid(trait~condition, scales = "free") +
  #facet_grid(~condition, scales = "free") +
  #geom_text(data=subset(plotdf,h2_labl =="plot"),aes(y=1, label = as.character(signif(H2,2))), color= "black", size = 3.5) +
  geom_text(data=subset(plotdf,h2_labl =="plot" & trait == "Length"),aes(y=-50, label = as.character(signif(H2,2))), color= "black", size = 3.5) +
  geom_text(data=subset(plotdf,h2_labl =="plot" & trait == "Brood Size"),aes(y=50, label = as.character(signif(H2,2))), color= "black", size = 3.5) +
  geom_text(data=subset(plotdf,h2_labl =="plot" & trait == "Pharyngeal Pumping"),aes(y=20, label = as.character(signif(H2,2))), color= "black", size = 3.5) +
  #geom_text(data=subset(plotdf,strain =="N2"),aes(y=0, label = as.character(signif(H2,2))), color= "red") +
  theme(legend.title = element_text( size=0, face="bold"),
        panel.background = element_rect(colour = "black", size=1),
        strip.background = element_rect(colour = "black", fill = "grey80",size = 1)) +
  theme(strip.background = element_rect(colour = "black", fill = "white",
                                          size = 0.5, linetype = "dotted")) +
  theme(strip.text.x = element_text(face="bold", colour = "black", size = 12,
                                      hjust = 0.5, vjust = 0.5),
        strip.text.y = element_text(face="bold", colour = "black", size = 10,
                                    hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.x = element_text(size=10, color="black"),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_text(size=12, face="bold", color="black"),
        axis.title.y = element_text(size=12, face="bold", color="black"))
plot

ggsave(plot = plot, filename = "Ce_DR.tiff", width = 12, height = 8)
