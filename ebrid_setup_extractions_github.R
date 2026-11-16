##### set up ebird data for SDMs based on existing datasets

# load packages
library(tidyverse)
library(raster)
library(vegan)

# combine existing seasonal ebird data, cut unneeded columns
# data originally compiled for Cohen 2025 Global Ecology and Biogeography 
# https://onlinelibrary.wiley.com/doi/10.1111/geb.13934
# refer to https://github.com/jeremy-cohen/grain_sdm for file assembly
setwd("/gpfs/gibbs/pi/jetz/data/temporary_datasets/scale/") # set working directory
# loop seasons
for(season in c("breeding","postbreeding_migration","nonbreeding","prebreeding_migration")){
  # load existing data
  ebird = read_csv(paste0("ebird_",season,"_CEA.csv"))
  # select needed columns only
  ebird = ebird[,c(1:10,238,241,688:731)]
  # compile iteratively
  if(season=="breeding"){ebirdfull = ebird}else{ebirdfull = rbind(ebirdfull, ebird)}
  # cleanup
  rm(ebird)
  gc()
}
# write master dataset
write_csv(ebirdfull, "ebird_mitchell.csv")

# stack crop rasters
setwd("/gpfs/gibbs/pi/jetz/data/environmental_datasets/CDL_2023/")
lf = list.files()
rasternames = c("Cash_Row_Crop","Cereals","Fallow_Idle",
                "Grazing_Grassland","Human_Veg_Fruit",
                "Orchard_Vineyard","Rice_Aquaculture")

# loop layers
for (i in 1:length(lf)){
  print(lf[i])
  # load raster
  r = raster(lf[i])
  # mean aggregate to ~1km
  r = raster::aggregate(r, 33, fun="mean")
  # transform to be consistent with ebird
  r = projectRaster(r, res=1000,
                    crs="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
  # write
  writeRaster(r, paste0(rasternames[i],".tif"), format="GTiff")
}

# load, stack layers
setwd("/gpfs/gibbs/pi/jetz/data/environmental_datasets/CDL_2023/")
lf = list.files()[8:14]
stack = stack()
# loop layers
for (i in 1:length(lf)){
  # load raster
  r = raster(lf[i])
  # stack
  stack = stack(stack, r)
}

# extract 7 crop rasters to ebird
setwd("/gpfs/gibbs/pi/jetz/data/temporary_datasets/scale/")
ebird = read_csv("ebird_mitchell.csv")
coordinates(ebird) = c("x","y")
ex = raster::extract(stack, ebird)
ebird = cbind(ebird, ex)
# NAs are outside the US and those obs can be removed
ebird = ebird[complete.cases(ebird$Rice_Aquaculture),]
ebird = data.frame(ebird)
ebird$optional = NULL
# renaming
names(ebird)[9:10] = c("Tyto_alba","Bubo_virginianus")
names(ebird)[11:54] <-  gsub('.{2}$', '', names(ebird)[11:54])
write_csv(ebird, "ebird_mitchell.csv")

# load and extract to prediction surface
pred_surface = read_csv("habitat_prediction-surface_1.csv")
coordinates(pred_surface) = c("x","y")
ex = raster::extract(stack, pred_surface)
pred_surface = cbind(pred_surface, ex)
# NAs are outside the US and those cells can be removed
pred_surface = pred_surface[complete.cases(pred_surface$Rice_Aquaculture),]
pred_surface = data.frame(pred_surface)
pred_surface$optional = NULL
# save
write_csv(pred_surface, "pred_surface_mitchell.csv")


### v3: add more raptors, include ldi and tighten crop rasters
setwd("/gpfs/gibbs/pi/jetz/data/temporary_datasets/scale/")
## ebird
# add species
ebird = read_csv("ebird_mitchell.csv")
species = c("Buteo platypterus", "Buteo lineatus", "Buteo jamaicensis",
            "Buteo regalis", "Buteo swainsoni", "Accipiter atricapillus",
            "Accipiter cooperii", "Ictinia mississippiensis",
            "Circus hudsonius", "Asio flammeus", "Asio otus", "Strix varia",
            "Falco peregrinus", "Falco mexicanus")
# loop seasons
for(season in c("breeding","postbreeding_migration","nonbreeding","prebreeding_migration")){
   old = read_csv(paste0("ebird_",season,"_CEA.csv"))
     names(old)
     wh = which(names(old) %in% species) 
     old = old[,c(1, wh)] # change to just an ID and species
     if(season=="breeding"){ebirdfull = old}else{ebirdfull = rbind(ebirdfull, old)}
     rm(old)
     gc()
}
# combine datasets
ebird = left_join(ebird, ebirdfull, by="checklist_id")
rm(ebirdfull)
gc()
# merge crop layers
ebird$row_crop = ebird$Rice_Aquaculture + ebird$Cash_Row_Crop +
                     ebird$Cereals + ebird$Human_Veg_Fruit
ebird$non_row_crop = ebird$Fallow_Idle + ebird$Grazing_Grassland +
                     ebird$Orchard_Vineyard
# Shannon habitat diversity
hab <- ebird[, 19:54]
ebird$shannon_diversity <- diversity(hab, index = "shannon")
# save
write_csv(ebird, "ebird_mitchell_v2.csv")


## prediction surface update for v3
pred_surface = read_csv("pred_surface_mitchell.csv")
# merge crop layers
pred_surface$row_crop = pred_surface$Rice_Aquaculture + pred_surface$Cash_Row_Crop +
                     pred_surface$Cereals + pred_surface$Human_Veg_Fruit
pred_surface$non_row_crop = pred_surface$Fallow_Idle + pred_surface$Grazing_Grassland +
                         pred_surface$Orchard_Vineyard
# Shannon diversity
hab <- pred_surface[, 9:44]
pred_surface$shannon_diversity <- diversity(hab, index = "shannon")
# save
write_csv(pred_surface, "pred_surface_mitchell_v2.csv")
