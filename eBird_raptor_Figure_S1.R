#Species-specific rankings of model importance scores for eBird-owls manuscript (Figure S1)
#Plot shows covariates with top 10 highest importance scores for each species

#Load packages
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(tidytext)

#Set working directory
setwd("~/eBird_raptors/")

#Open dataframe with summary of importance scores
impdf <- read.csv("importance_adjusted_scores_GHOW_for_plotting.csv")

#Convert df to long format
impdf_long <- impdf %>%
  pivot_longer(
    cols = -COVARIATE,
    names_to = "SPECIES",
    values_to = "IMPORTANCE"
  )

#Get top 10 covariates per species
top10 <- impdf_long %>%
  group_by(SPECIES) %>%
  arrange(desc(IMPORTANCE), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

#Reorder for plotting highest to lowest importance score
top10_ordered <- top10 %>%
  group_by(SPECIES) %>%
  mutate(covariate_ordered = factor(COVARIATE, levels = COVARIATE[order(IMPORTANCE)]))

#Summarize per species
abow <- top10 %>%
  filter(SPECIES == "ABOW") #american barn owl

baow <- top10 %>%
  filter(SPECIES == "BAOW") #barred owl

bwha <- top10 %>%
  filter(SPECIES == "BWHA") #broad-winged hawk

coha <- top10 %>%
  filter(SPECIES == "COHA") #cooper's hawk

feha <- top10 %>%
  filter(SPECIES == "FEHA") #ferruginous hawk

leow <- top10 %>%
  filter(SPECIES == "LEOW") #long-eared owl

miki <- top10 %>%
  filter(SPECIES == "MIKI") #mississippi kite

noha <- top10 %>%
  filter(SPECIES == "NOHA") #northern harrier

pefa <- top10 %>%
  filter(SPECIES == "PEFA") #peregrine falcon

prfa <- top10 %>%
  filter(SPECIES == "PRFA") #prairie falcon

rsha <- top10 %>%
  filter(SPECIES == "RSHA") #red-shouldered hawk

rtha <- top10 %>%
  filter(SPECIES == "RTHA") #red-tailed hawk

seow <- top10 %>%
  filter(SPECIES == "SEOW") #short-eared owl

swha <- top10 %>%
  filter(SPECIES == "SWHA") #swainson's hawk

#Plot individual panels by species
##american barn owl
abow10 <- ggplot(abow, aes(x = IMPORTANCE, y = reorder(COVARIATE, IMPORTANCE))) +
  geom_col(fill = "cornflowerblue") +
  labs(x = "", y = "", title = "American Barn Owl") +
  theme(
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 22),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), #element_rect(color = "black", size = 0.5),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 16, angle = 0),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 20),
    text = element_text(size=10))

##barred owl
baow10 <- ggplot(baow, aes(x = IMPORTANCE, y = reorder(COVARIATE, IMPORTANCE))) +
  geom_col(fill = "cornflowerblue") +
  labs(x = "", y = "", title = "Barred Owl") +
  theme(
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 22),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), #element_rect(color = "black", size = 0.5),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 16, angle = 0),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 20),
    text = element_text(size=10))

#broad-winged hawk
bwha10 <- ggplot(bwha, aes(x = IMPORTANCE, y = reorder(COVARIATE, IMPORTANCE))) +
  geom_col(fill = "cornflowerblue") +
  labs(x = "", y = "", title = "Broad-winged Hawk") +
  theme(
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 22),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), #element_rect(color = "black", size = 0.5),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 16, angle = 0),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 20),
    text = element_text(size=10))

#cooper's hawk
coha10 <- ggplot(coha, aes(x = IMPORTANCE, y = reorder(COVARIATE, IMPORTANCE))) +
  geom_col(fill = "cornflowerblue") +
  labs(x = "", y = "", title = "Cooper's Hawk") +
  theme(
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 22),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), #element_rect(color = "black", size = 0.5),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 16, angle = 0),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 20),
    text = element_text(size=10))

##ferruginous hawk
feha10 <- ggplot(feha, aes(x = IMPORTANCE, y = reorder(COVARIATE, IMPORTANCE))) +
  geom_col(fill = "cornflowerblue") +
  labs(x = "", y = "", title = "Ferruginous Hawk") +
  theme(
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 22),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), #element_rect(color = "black", size = 0.5),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 16, angle = 0),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 20),
    text = element_text(size=10))

##long-eared owl
leow10 <- ggplot(leow, aes(x = IMPORTANCE, y = reorder(COVARIATE, IMPORTANCE))) +
  geom_col(fill = "cornflowerblue") +
  labs(x = "", y = "", title = "Long-eared Owl") +
  theme(
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 22),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), #element_rect(color = "black", size = 0.5),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 16, angle = 0),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 20),
    text = element_text(size=10))

##mississippi kite
miki10 <- ggplot(miki, aes(x = IMPORTANCE, y = reorder(COVARIATE, IMPORTANCE))) +
  geom_col(fill = "cornflowerblue") +
  labs(x = "", y = "", title = "Mississippi Kite") +
  theme(
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 22),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), #element_rect(color = "black", size = 0.5),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 16, angle = 0),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 20),
    text = element_text(size=10))

##northern harrier
noha10 <- ggplot(noha, aes(x = IMPORTANCE, y = reorder(COVARIATE, IMPORTANCE))) +
  geom_col(fill = "cornflowerblue") +
  labs(x = "", y = "", title = "Northern Harrier") +
  theme(
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 22),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), #element_rect(color = "black", size = 0.5),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 16, angle = 0),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 20),
    text = element_text(size=10))

##peregrine falcon
pefa10 <- ggplot(pefa, aes(x = IMPORTANCE, y = reorder(COVARIATE, IMPORTANCE))) +
  geom_col(fill = "cornflowerblue") +
  labs(x = "", y = "", title = "Peregrine Falcon") +
  theme(
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 22),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), #element_rect(color = "black", size = 0.5),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 16, angle = 0),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 20),
    text = element_text(size=10))

##prairie falcon
prfa10 <- ggplot(prfa, aes(x = IMPORTANCE, y = reorder(COVARIATE, IMPORTANCE))) +
  geom_col(fill = "cornflowerblue") +
  labs(x = "", y = "", title = "Prairie Falcon") +
  theme(
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 22),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), #element_rect(color = "black", size = 0.5),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 16, angle = 0),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 20),
    text = element_text(size=10))

##red-shouldered hawk
rsha10 <- ggplot(rsha, aes(x = IMPORTANCE, y = reorder(COVARIATE, IMPORTANCE))) +
  geom_col(fill = "cornflowerblue") +
  labs(x = "", y = "", title = "Red-shouldered Hawk") +
  theme(
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 22),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), #element_rect(color = "black", size = 0.5),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 16, angle = 0),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 20),
    text = element_text(size=10))

##red-tailed hawk
rtha10 <- ggplot(rtha, aes(x = IMPORTANCE, y = reorder(COVARIATE, IMPORTANCE))) +
  geom_col(fill = "cornflowerblue") +
  labs(x = "", y = "", title = "Red-tailed Hawk") +
  theme(
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 22),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), #element_rect(color = "black", size = 0.5),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 16, angle = 0),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 20),
    text = element_text(size=10))

##short-eared owl
seow10 <- ggplot(seow, aes(x = IMPORTANCE, y = reorder(COVARIATE, IMPORTANCE))) +
  geom_col(fill = "cornflowerblue") +
  labs(x = "", y = "", title = "Short-eared Owl") +
  theme(
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 22),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), #element_rect(color = "black", size = 0.5),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 16, angle = 0),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 20),
    text = element_text(size=10))

##swainson's hawk
swha10 <- ggplot(swha, aes(x = IMPORTANCE, y = reorder(COVARIATE, IMPORTANCE))) +
  geom_col(fill = "cornflowerblue") +
  labs(x = "", y = "", title = "Swainson's Hawk") +
  theme(
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 22),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), #element_rect(color = "black", size = 0.5),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 16, angle = 0),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 20),
    text = element_text(size=10))

#Create plot pane with all panels
grid <- plot_grid(abow10, baow10, bwha10, coha10, feha10,
                  leow10, miki10, noha10, pefa10, prfa10,
                  rsha10, rtha10, seow10, swha10, 
                  #labels = c('A','B','C','D','E','F','G','H','I','J','K','L','M','N'),
                  ncol = 4, nrow = 4)

#Add labels and adjust axis text
final_plot <- ggdraw() +
  draw_plot(grid) +
  draw_label("Adjusted Importance Score", x = 0.5, y = 0, vjust = -0.5, size = 20) +
  draw_label("Covariate", x = 0, y = 0.5, angle = 90, vjust = 1.5, size = 20)

#Save
ggsave("Figure_S1.png", final_plot, "png", height = 10, width = 14, units = "in")

