# Cross-species models, figures, tables for eBird-owls manuscript

# load packages and set working directory
library(tidyverse)
library(ggpubr)
setwd("~/mitchell/")

### Partial dependence figure (fig. 2)
tab_points = read_csv("pd_points.csv") # load model predictions at intervals
tab_plots = read_csv("pd_plots.csv") # load smoothed version of model predictions
# order, loop species
sp = sort(unique(tab_points$species))
for (i in 1:length(sp)){
  # calculation of GHOW metric - difference between species predicted
  # occurrence at GHOW occ. prob. of 0.1 to 0.5
  tabs = tab_plots[tab_plots$species==sp[i],] # data for single species
  fif = tabs[which.min(abs(tabs$x - .5)),]$y # y value at closest x to 0.5
  ten = tabs[which.min(abs(tabs$x - .1)),]$y # y value at closest x to 0.1
  dif = fif - ten
  # compiling from each loop
  if(i==1){difs = dif}else{difs = c(difs, dif)}
  tabs = tab_points[tab_points$species==sp[i],] # data for single species
  # plot - individual panels 
  plot = ggplot(tabs, aes(x,y)) +
    geom_point(color="gray20", size=.5) +
    geom_smooth(color="black", size=.5) +
    theme(text = element_text(size=7.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text.y = element_text(color = "black", size=8),
          axis.text.x = element_text(color = "black", size=8),
          legend.position="none") +
    ylim(0,1.6) +
    xlim(0,.8) +
    scale_x_continuous(breaks = c(0,.5)) +
    xlab("") +
    ylab("") +
    annotate("text", x = Inf, y = Inf, label = round(dif,2), 
             vjust = 2, hjust = 1.5, size=3.5) +
    ggtitle(sp[i])
  assign(paste0("sp.",i), plot)
}
# create table from species level outputs
mtable = data.frame(species=sp, metric=difs)
write_csv(mtable, "sp_metrics_0_v2.csv")
# create plot from panels
gga = ggarrange(sp.1, sp.2, sp.3, sp.4, sp.5, sp.6, sp.7,
                 sp.8, sp.9, sp.10, sp.11, sp.12, sp.13, sp.14,
                ncol=7, nrow=2)
ggsave("sp_pd_plots_v2.jpeg", gga, "jpeg", height=4, width=12, units="in")


### Traits analysis/figure (fig. 3)
sptable = read_csv("sp_metrics.csv") # note- habitat/geography/taxon fields assigned in excel
sptable = sptable %>% rename(metric = `GHOW metric`)
sptable = sptable[2:15,] # remove GHOW, as models will compare other species
# add log body mass data from Avonet
fulltable = read_csv("~/traits/avonet/TraitData/AVONET2_ebird.csv")
fulltable$`Scientific name` = gsub(" ","_",fulltable$Species2)
sptable = left_join(sptable, fulltable[,c("Scientific name","Mass"),
                                       by="Scientific name"])
# log-transform mass
sptable$logmass = log(sptable$Mass)


## Figure panels (fig 3)
# habitat
gg1 = ggplot(sptable, aes(Habitat, metric)) + 
  theme_bw() +
  theme(axis.text=element_text(size=13)) +
  geom_hline(yintercept=0, linetype="dashed", color="black") +
  geom_boxplot() + 
  ggtitle("Habitat") +
  ylab("Change in occurrence probability with GHOW") +
  xlab("")

# geography
gg2 = ggplot(sptable, aes(Geography, metric)) + 
  theme_bw() +
  theme(axis.text=element_text(size=13)) +
  geom_hline(yintercept=0, linetype="dashed", color="black") +
  geom_boxplot() + 
  ggtitle("Geography") +
  ylab(NULL) +
  xlab("")

# taxon
gg3 = ggplot(sptable, aes(Taxon, metric)) + 
  theme_bw() +
  theme(axis.text=element_text(size=13)) +
  geom_hline(yintercept=0, linetype="dashed", color="black") +
  geom_boxplot() + 
  ggtitle("Taxon") +
  ylab(NULL) +
  xlab("")

# body mass
gg4 = ggplot(sptable, aes(logmass, metric)) + 
  geom_point(col="gray60", pch=16, size=2, alpha=.3) +
  geom_smooth(color="black", method="lm") +
  theme_bw() +
  theme(axis.text=element_text(size=13)) +
  ggtitle("Body mass") +
  ylab(NULL) +
  xlab("")

# full model
# GLM
glm = glm(metric ~ Habitat + Geography + Taxon + logmass, data=sptable)
table = round(as.data.frame(summary(glm)$coefficients),3)
write.csv(table, "glm_full.csv")
# ANOVA for evaluation of categorical predictors
anovaglm = aov(glm)
table <- as.data.frame(summary(anovaglm)[[1]])
write.csv(table, "anova_full.csv")
# Tukey's HSD for comparison of categorical levels when 3+ levels
tukeyG = TukeyHSD(anovaglm, "Geography")$Geography
tukeyT = TukeyHSD(anovaglm, "Taxon")$Taxon
write.csv(tukeyG, "tukeys_geography.csv")
write.csv(tukeyT, "tukeys_taxon.csv")
# save figure
gga = ggarrange(gg1, gg2, gg3, gg4, ncol=4, nrow=1)
ggsave("sp_trait_boxplots.jpeg", gga, "jpeg", height=4, width=10.7, units="in")


### Model Metrics (fig. 4)
data = read_csv("sdm_birds_VM3.1.csv")
# Get change between models including or omitting GHOW occ. prob.
rm(out)
  for (i in 2:length(unique(data$sciname))){ # excluding ghow
    sci_name <- unique(data$sciname)[i]
    common <- unique(data$common_name)[i]
    sub_noghow <- data[data$sciname==sci_name & data$model=="without GHOW",]
    sub_ghow <- data[data$sciname==sci_name & data$model=="with GHOW",]
    difs <- sub_ghow[,9:65]-sub_noghow[,9:65]
        collect <- c(common, sci_name, as.numeric(difs))
        if(!exists("out")){out <- collect}else{out <- rbind(out, collect)}
  }
# process output data frame, save
colnames(out) <- c("common_name", "sciname", 
                   paste0("delta_",colnames(data)[9:65]))
rownames(out) <- NULL
out <- data.frame(out)
write_csv(out, "VM3.1_comp.csv")

# load data for  model and figure
out = read_csv("VM3.1_comp.csv")
# list metrics
metrics <- c("AUC","TSS","Sensitivity","Specificity")
# subset table - models without GHOW
datawo <- data[data$model=="without GHOW",]
# loop through four metrics
  for(metric in c("AUC","TSS","Sensitivity","Specificity")){
    # designate metric
    datawo$metric <- datawo[,metric][[1]]
    # violin plot panel for all models lacking GHOW
    gg <- ggplot(data=datawo, aes(x=model, y=metric)) +
      geom_hline(yintercept = median(datawo$metric), size = 1.5, color = "gray") +
      geom_violin(trim=F, alpha = .6, fill="transparent") +
      scale_fill_manual(values = "transparent")  +
      stat_summary(fun="median", geom="point", shape=21, size=4,
                   stroke=1.5, fill="white", col="black") +
      xlab(NULL) +
      ylab(metric) +
      theme_light() +
      theme(legend.position = "none",
            text = element_text(size = 25))
    
    # subset, designate metric for models with GHOW
    out = read_csv("VM3.1_comp.csv")
    out$model = "with GHOW"
    out$metric <- out[,paste0("delta_",metric)][[1]]
    
    # violin plot panel for all models with GHOW
    gg2 <- ggplot(data=out, aes(x=model, y=metric)) +
      geom_hline(yintercept = 0, size = 1.5, color = "gray") +
      geom_violin(trim=F, alpha = .6) +
      scale_fill_manual(values = "transparent") +
      stat_summary(fun="median", geom="point", shape=21, size=4,
                   stroke=1.5, fill="white", col="black") +
      ylab(expr(paste(Delta, " ", !!metric))) +
      xlab(NULL) +
      theme_light() +
      theme(legend.position = "none",
            text = element_text(size = 25))
    
    # save both
    assign(paste0("gg.",metric),gg)
    assign(paste0("gg2.",metric),gg2)
  }

# save
gg.ppm <- ggarrange(gg.AUC, gg2.AUC, gg.TSS, gg2.TSS,
                    gg.Sensitivity, gg2.Sensitivity, 
                    gg.Specificity, gg2.Specificity,
                    nrow=2, ncol=4, widths = c(2,2), labels=c(
                      "a)","","b)","","c)","","d)",""),
                    font.label=list(color="black",size=20))
ggsave("ppms_violin.jpeg", gg.ppm, "jpeg", height=8, width=16, units="in")

# means and errors for reported metrics in text
mean(out$delta_TSS)
sd(out$delta_TSS) / sqrt(14) # se
mean(out$delta_Sensitivity)
sd(out$delta_Sensitivity) / sqrt(14) # se

