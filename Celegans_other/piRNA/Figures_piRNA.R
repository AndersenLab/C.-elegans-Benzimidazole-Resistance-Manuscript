#library(devtools)
#install_github("AndersenLab/easysorter")
library(easysorter)
library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(cowplot)

setwd("Celegans_other/piRNA")

########################
#### PLOT piRNA interval NILs
########################

load("Interval_NILs.Rda")
regressed$strain <- factor(regressed$strain, 
                           levels=c("N2","CB","ECA240","ECA241"))
strainsplot <- c("N2","CB","ECA240","ECA241")

alb <- subset(regressed, condition=="Albendazole") %>%
  filter(trait == "q75.EXT")%>%
  filter(strain %in% strainsplot) 
plot <- ggplot(alb)+
  aes(x = strain, y = phenotype, fill = strain)+
  geom_boxplot(outlier.size= NA, alpha = 0.5, lwd=1)+
  geom_jitter(size =2.5, alpha=.5, width = 0.35)+
  theme_bw()+
  #facet_grid(.~assay) +
  theme(legend.position="none",
        panel.border = element_rect(size=1, colour = "black"))+
  labs(x = "Strain", y = "Regressed Length (µm)", title="Albendazole q75.EXT") +
  scale_fill_manual(values = c("N2" = "orange","CB" = "blue","ECA240" = "black", "ECA241" = "black", "ECA242" = "black", "SnuIR3" = "black", "SnuIR4" = "black","SnuIR6" = "black","SnuIR7" = "black","SnuIR9" = "black","SnuIR10" = "black","SnuIR12" = "black")) +
  theme(axis.text.x = element_text(size=14, face="bold", color="black"),
        axis.text.y = element_text(size=12, face="bold", color="black"),
        axis.title.x = element_text(size=0, face="bold", color="black"),
        axis.title.y = element_text(size=16, face="bold", color="black"),
        plot.title = element_text(size=0, face="bold")) +
  scale_x_discrete(labels = c("N2","CB4856","ECA240","ECA241"))
plot

x1 <- .175
y1 <- .03
y2 <- 0.009
blck <- .12
spce <- .09
x2 <- x1+blck+spce
x3 <- x2+blck+spce
x4 <- x3+blck+spce
x5 <- x4+blck+spce
x6 <- x5+blck+spce
# chr IV length:	17493829

piRNA_NILs <- ggdraw() +
  draw_plot(plot, 0, .15, 1, .8) +
  draw_label("ChrIV", x = 0, y = .1, size = 12, hjust = -.16, fontface = 'bold')+
  geom_rect(aes(xmin = x1, xmax = x1+blck, ymin = y1*3, ymax = y1*3 +.02),
            colour = "black", fill = "orange")+
  geom_rect(aes(xmin = x2, xmax = x2+blck, ymin = y1*3, ymax = y1*3 +.02),
            colour = "black", fill = "blue")+
  geom_rect(aes(xmin = x3, xmax = x3+blck, ymin = y1*3, ymax = y1*3 +.02),
            colour = "black", fill = "orange")+
  geom_rect(aes(xmin = x3+blck*(13.37/17.49), xmax = x3+blck*(16.9/17.49), ymin = y1*3, ymax = y1*3 +.02),
            colour = "blue", fill = "blue")+
  geom_rect(aes(xmin = x4, xmax = x4+blck, ymin = y1*3, ymax = y1*3 +.02),
            colour = "black", fill = "blue")+
  geom_rect(aes(xmin = x4+blck*(13.37/17.49), xmax = x4+blck*(16.9/17.49), ymin = y1*3, ymax = y1*3 +.02),
            colour = "orange", fill = "orange")+
  geom_rect(aes(xmin = x1, xmax = x1+blck, ymin = y2*3, ymax = y2*3 +.02),
            colour = "black", fill = "orange")+
  geom_rect(aes(xmin = x2, xmax = x2+blck, ymin = y2*3, ymax = y2*3 +.02),
            colour = "black", fill = "blue")+
  geom_rect(aes(xmin = x3, xmax = x3+blck, ymin = y2*3, ymax = y2*3 +.02),
            colour = "black", fill = "orange")+
  geom_rect(aes(xmin = x4, xmax = x4+blck, ymin = y2*3, ymax = y2*3 +.02),
            colour = "black", fill = "blue")+
  # draw_label("13.37 - 17.24 Mb", x =  x3+blck*.5, y = y1*3 +.05, size = 12)+
  #  draw_label("13.37 - 17.24 Mb", x =  x4+blck*.5, y = y1*3 +.05, size = 12)+
  draw_label("Genome", x = 0, y = .04, hjust = -.1, size = 12, fontface = 'bold') 
piRNA_NILs



########################
#### PLOT Mutants / alternative with SD bars
########################

#bar-plot alternative (regression without linear model)
#se <- function(x) sqrt(var(x)/length(x))
se <- function(x) sd(x) #actually sd

load("GitHub-BZ-prep/Celegans_other/Mutants/UnRegressed_piRNAmut1ab.Rda")
biopruned <- unregressed %>% ungroup()

#pick strains 
strainslist <- c("N2", "CB", "ECA270","ECA271", "ECA286")

barplot <- biopruned %>%
  #filter(assay == "a") %>%
  filter(strain %in% strainslist) %>%
  select(condition,assay,strain,q75.EXT) %>%
  group_by(strain,condition) %>%
  mutate(mean_phenotype = ifelse(condition == "DMSO",mean(q75.EXT),0)) %>%
  ungroup() %>%
  group_by(strain) %>%
  mutate(mean_phenotype = ifelse(condition != "DMSO",max(mean_phenotype),mean_phenotype)) %>%
  ungroup() %>%
  filter(condition != "DMSO") %>%
  mutate(phenotype = q75.EXT - mean_phenotype) 
Assay_a_mean <- barplot %>% 
  filter(assay == "a")
Assay_a_mean <- mean(Assay_a_mean$phenotype) #-94.63392
Assay_b_mean <- barplot %>% 
  filter(assay == "b")
Assay_b_mean <- mean(Assay_b_mean$phenotype) #-233.0829
#-233.0829+94.63392 = 138.449 (subtract 138.449 from a)
barplot <- barplot %>%
  mutate(phenotype_R = ifelse(assay == "a", phenotype - 138.449, phenotype)) %>%
  group_by(strain) %>%
  mutate(Mean = mean(phenotype_R)) %>% ungroup() %>%
  mutate(SD = se(phenotype_R))%>% ungroup()

barplot$strain <- factor(barplot$strain, 
                         levels=c("N2", "CB", "ECA270","ECA271", "ECA286"))

strain_mod = lm(formula = phenotype_R ~ strain, data = barplot)
summary(strain_mod)
# Coefficients:
#                 Estimate    Std. Error t value Pr(>|t|)    
#   (Intercept)  -187.605      9.043 -20.747  < 2e-16 ***
#   strainCB     -129.641     12.577 -10.308  < 2e-16 ***
#   strainECA270   33.843     12.480   2.712  0.00702 ** 
#   strainECA271   28.241     12.480   2.263  0.02424 *  
#   strainECA286  -90.348     12.528  -7.212 3.32e-12 ***
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 68.27 on 357 degrees of freedom
# Multiple R-squared:  0.5113,	Adjusted R-squared:  0.5044 
# F-statistic: 74.69 on 5 and 357 DF,  p-value: < 2.2e-16
strain_anova = anova(strain_mod)
# Response: phenotype_R
# Df  Sum Sq Mean Sq F value    Pr(>F)    
# strain      5 1740666  348133  74.694 < 2.2e-16 ***
#   Residuals 357 1663902    4661                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

aov.out = aov(phenotype_R ~ strain, data=barplot)
TukeyHSD(aov.out)


N2test <- filter(barplot,strain == "N2")
N2test <- N2test$phenotype_R
CBtest <- filter(barplot,strain == "CB")
CBtest <- CBtest$phenotype_R
test270 <- filter(barplot,strain == "ECA270")
test270 <- test270$phenotype_R
test271 <- filter(barplot,strain == "ECA271")
test271 <- test271$phenotype_R
test286 <- filter(barplot,strain == "ECA286")
test286 <- test286$phenotype_R

t.test(N2test,CBtest) p-value < 2.2e-16
t.test(N2test,test270) p-value = 0.00845
t.test(N2test,test271) p-value = 0.01471
t.test(N2test,test286) p-value = 1.195e-11

library(scales)
plot_mutants_bar <-  ggplot(barplot)+
  aes(x = strain, y = phenotype_R, fill = strain)+
  geom_bar(stat = "summary", fun.y = "mean", color = "black", alpha = 0.75) +
  geom_errorbar(aes(ymax = Mean + SD, ymin= Mean - SD), width=0.2) +
  #geom_boxplot(outlier.size= NA, alpha = 0.5, lwd=1)+
  #geom_jitter(size =2.5, alpha=.5, width = 0.75)+
  theme_bw()+
  theme(legend.position="none",
        panel.border = element_rect(size=1, colour = "black"))+
  labs(x = "Strain", y = "Length change in albendazole (µm)") +
  scale_fill_manual(values = c("N2" = "orange","CB" = "blue","ECA270" = "orange","ECA271" = "orange","ECA286" = "orange","ECA352" = "blue","FX234" = "orange")) +
  scale_colour_manual(values = c("N2" = "orange","CB" = "blue","ECA270" = "orange","ECA271" = "orange","ECA286" = "orange","ECA352" = "blue","FX234" = "orange")) +
  theme(axis.text.x = element_text(size=12, face="bold", color="black"),
        axis.text.y = element_text(size=12, face="bold", color="black"),
        axis.title.x = element_text(size=0, face="bold", color="black"),
        axis.title.y = element_text(size=16, face="bold", color="black"),
        strip.text.x = element_text(size=16, face="bold", color="black"),
        strip.text.y = element_text(size=16, face="bold", color="black"),
        plot.title = element_text(size=0, face="bold")) +
  scale_x_discrete(labels=c("N2" = expression(paste(bold("N2"))), 
                            "CB" = expression(paste(bold("CB4856"))),
                            "ECA270" = expression(paste(bolditalic("alg-4(0); alg-3(0)"))), 
                            "ECA271" = expression(paste(bolditalic("ergo-1(0)"))),
                            "ECA286" = expression(paste(bolditalic("prg-1(0)")))))
plot_mutants_bar


#### FINAL PLOT
final_plot <- plot_grid(piRNA_NILs, plot_mutants_bar, volcano_plot, nrow=1, labels = c('A','B','C'), rel_widths = c(0.33,0.42,0.25),align = "h", scale =0.98)
final_plot
ggsave(plot = final_plot, filename = "Figure3.tiff", width = 18.5, height = 5)
