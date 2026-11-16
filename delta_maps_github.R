# Create maps that give delta probability of occurrence with/without GHOW (fig. 1f)

# load packages
library(raster)
library(tidyverse)
library(RColorBrewer)
library(plyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
# set working directory
setwd("/vast/palmer/home.mccleary/jc3893/mitchell/")

# load metrics, sort species
tab = read_csv("sp_metrics.csv")
sp = sort(tab$scientific)
# color spectrum
colors = brewer.pal(7, "Spectral")
cols = c("tomato1","lemonchiffon","transparent","steelblue3")

# world map base
world = ne_countries(scale='medium', returnclass = 'sf') %>%
  st_transform(crs=st_crs('+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs '))

# loop species
for (i in 1:length(sp)){
  id = tab[tab$scientific==sp[i],]$ID
  
  # relative occurrence rate maps
  # load spatial predictions without and with GHOW
  r_no =  raster(paste0("/vast/palmer/pi/jetz/from_gibbs/from_loomis/data/results_30x30/birds/VM3.1/",
                        sp[i],"/predictions/ROR/PO_predictions_",sp[i],"_",id,"_1.tif"))
  r_ghow = raster(paste0("/vast/palmer/pi/jetz/from_gibbs/from_loomis/data/results_30x30/birds/VM3.1/",
                  sp[i],"/predictions/ROR/PO_predictions_",sp[i],"_",id,"_2.tif"))
  # difference between models, trim to species range extent
  r_dif = r_ghow - r_no
  r_dif = trim(r_dif)
  # get difference values from raster, find max
  gv = getValues(r_dif)
  min = round_any(min(gv, na.rm=T), .1, floor)
  max = round_any(max(gv, na.rm=T), .1, ceiling)
  high = max(c(max, abs(min)))
  # plot map including specific value breaks
  jpeg(paste0(sp[i],"_delta.jpeg"), width = 800, height = 800, units = "px")
  plot(r_dif,  col=colors, breaks = c(high, high*.5, high*.25, high*.1, -high*.1,
                                      -high*.25, -high*.5, -high),
       axes=F)
  plot(world[1], add=T, col="transparent")
  dev.off()

}
