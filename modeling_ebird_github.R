# eBird modeling workflow
# v3- more raptor species, shannon diversity covariate, lumped crop layers

{library(maps)
  library(rgdal)
  library(raster)
  library(dismo)
  library(data.table)	
  library(scales)
  library(ranger)
  library(fields)  
  library(PresenceAbsence) 
  library(pdp)
  library(tidyverse)
  library(gbm)
  library(sf)
  library(rgeos)
  library(plyr)
  library(ff)
  library(ffbase)
  library(ggpubr)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(blockCV)
  
  sessionInfo()
  capabilities() 
}

###############################
#MODELS (RF)
###############################
{
  setwd("/vast/palmer/home.mccleary/jc3893/")
  version = "VM3.1"
  
  # Species names
  species = c("Bubo virginianus", "Bubo virginianus", # duplicate because we need annual and may-aug
              "Buteo platypterus", "Buteo lineatus", 
              "Buteo jamaicensis", "Buteo regalis", "Buteo swainsoni", 
              "Accipiter cooperii", "Ictinia mississippiensis", 
              "Circus hudsonius", "Tyto alba", 
              "Asio flammeus", "Asio otus", "Strix varia", 
              "Falco peregrinus", "Falco mexicanus")
  list = read_csv("30x30/species_list_2024.csv") %>%
    subset(Scientific %in% species)
  
  # shapefile
  countries = readOGR("30x30/shapefiles", "countries")
  countries_proj = spTransform(countries, CRS="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
  aba = countries[c(9,38),] %>%
    spTransform(countries, CRS="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ") %>%
    st_as_sf()
  
  {# grid for thinning
    bb = bbox(countries_proj)
    cs = c(5000,5000)
    cc = bb[, 1] + (cs/2)  # cell offset
    cd = ceiling(diff(t(bb))/cs)  # number of cells per direction
    grd = GridTopology(cellcentre.offset=cc, cellsize=cs, cells.dim=cd)
    grd
    sp_grd = SpatialGridDataFrame(grd,
                                  data=data.frame(id=1:prod(cd)),
                                  proj4string=CRS(proj4string(countries)))
    }
  
  # standard grid
  bbox = readOGR("30x30/shapefiles", "30x30BBox_NAmerica") %>%
    spTransform(countries, CRS="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ") %>%
    st_as_sf()
  grid = raster("30x30/grids/Global_reference_raster_1km_CEA.tif") %>%
    crop(bbox)
  
  # list migrant species for different treatment (may-aug data only)
  migrants = c("Buteo platypterus", "Buteo jamaicensis", 
               "Buteo regalis", "Buteo swainsoni", 
               "Ictinia mississippiensis", 
               "Circus hudsonius", "Asio flammeus", "Asio otus", 
               "Falco peregrinus", "Falco mexicanus")
  
  # directories
  dir = paste0("/vast/palmer/pi/jetz/from_gibbs/from_loomis/data/results_30x30/birds/",version,"/")
  dir.create(dir)
  dir.create(paste0(dir, "maps/"))
  dir.create(paste0(dir, "pointmaps/"))
  dir.create(paste0(dir, "pdp/"))
  setwd(dir)
}

# ----------------------------------------------------------------------
# RANGER
# ----------------------------------------------------------------------		
{# model GHOW, then each species with and without GHOW
  for(sp in c(1:16)){
    gc()
    # Establish names
    listsub = list[list$Scientific==species[sp],]
    common = as.character(listsub[,1])
    sci.name = as.character(listsub[,2])
    sci_name = gsub(" ","_",sci.name)
    code = as.character(listsub[,5])   
    print(c(sp, common, sci_name))
    # Expert range = modeling domain
    setwd("/vast/palmer/home.mccleary/jc3893/")
    try({
      total.range <- readOGR(unzip(
        paste0("30x30/range_polygons/",code,"-range-2020.gpkg.zip"),
        paste0(code,"-range-mr-2020.gpkg")), "range")
      # seasonal range for migrants
      if(sci.name %in% migrants){
        total.range = total.range[total.range$season_name=="breeding" |
                                    total.range$season_name=="resident" ,]
      }
    })
    # alternatives if first range doesn't load
    if(!exists("total.range") | length(total.range$season_name)==0){
      try({
        total.range <- readOGR(unzip(
          paste0("30x30/range_polygons/",code,"-range-2021.gpkg.zip"),
          paste0(code,"-range-mr-2021.gpkg")), "range")
        # Get seasonal range of species
        total.range <- total.range[total.range$season_name==ss |
                                 total.range$season_name=="resident" ,]
      })
    }
    # another alternative
    if(!exists("total.range") | length(total.range$season_name)==0){
      mol.range <- readOGR("/vast/palmer/home.mccleary/jc3893/30x30/range_polygons", 
                           "jetz_maps_v2020")
      total.range <- mol.range[mol.range$sciname==sci.name & 
                             (mol.range$season==2 |
                                mol.range$season==1),]
      }
    try(file.remove(paste0(code,"-range-mr-2020.gpkg")))
    
    # project, box range
    sprange = spTransform(total.range, CRS="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
    sprangebuff = buffer(sprange, 200000)
    box = bbox(sprangebuff)
    # plot(sprangebuff, col="bisque")
    # Back to SP dataframe
    range.df = data.frame(ID=1, row.names="buffer")
    sprangebuff = SpatialPolygonsDataFrame(sprangebuff, range.df, match.ID=F)
    # ST format for range
    rangest = st_as_sf(sprangebuff)
    ograngest = st_as_sf(sprange)
    sf::sf_use_s2(FALSE)
    
    # Predictors
    preds = c(
      "year","day_of_year","time_observations_started","duration_minutes",
      "effort_distance_km","number_observers","bio1","bio12",
      "tri","elev","pland_10_cropland_rainfed", 
      "pland_100_mosaic_tree_shrub","pland_11_cropland_rainfed", 
      "pland_110_mosaic_herbacious","pland_12_cropland_rainfed", 
      "pland_120_shrubland","pland_121_shrubland","pland_122_shrubland",
      "pland_130_grassland","pland_140_lichens_mosses","pland_150_sparse",
      "pland_152_sparse","pland_153_sparse","pland_160_flooded_freshwater",
      "pland_170_flooded_saltwater","pland_180_flooded_shrub", 
      "pland_190_urban","pland_20_cropland_irrigated","pland_200_barren",     
      "pland_201_barren","pland_202_barren","pland_210_water",
      "pland_220_ice","pland_30_mosaic_cropland","pland_40_mosaic_natural_veg",
      "pland_50_evergreen_broadleaf","pland_60_deciduous_broadleaf", 
      "pland_61_deciduous_broadleaf","pland_62_deciduous_broadleaf", 
      "pland_70_evergreen_needleleaf","pland_71_evergreen_needleleaf",
      "pland_72_evergreen_needleleaf","pland_80_deciduous_needleleaf",
      "pland_81_deciduous_needleleaf","pland_82_deciduous_needleleaf",
      "pland_90_mixed_forest", "evisum", "row_crop", "non_row_crop", "shannon_diversity")
    

    # directories
    setwd(dir)
    dir.create(paste0(sci_name,"/"))
    dir.create(paste0(sci_name,"/predictions/"))
    dir.create(paste0(sci_name,"/predictions/ROR/"))
    dir.create(paste0(sci_name,"/predictions/PA/"))
    dir.create(paste0(sci_name,"/points/"))
    
    # ebird - get cropped rows, columns needed
    # for first ghow run use data with no ghow information
    # for later runs use full or breeding season data depending on migrant status
    if(sp==1){
      ebird = read_csv("/vast/palmer/pi/jetz/from_gibbs/data/temporary_datasets/scale/ebird_mitchell_v2.csv")
      # adjust mean temp
      ebird$bio1 = ebird$bio1/10
    }else{
      if(sci.name %in% migrants){
      ebird = read_csv("/vast/palmer/pi/jetz/from_gibbs/data/temporary_datasets/scale/ebird_mitchell_ghow_breeding.csv")
    }else{
      ebird = read_csv("/vast/palmer/pi/jetz/from_gibbs/data/temporary_datasets/scale/ebird_mitchell_ghow_full.csv")
    }
    }
    # NA removal (only a few rows)
    ebird = ebird[complete.cases(ebird),]
    # readLines("/gpfs/gibbs/pi/jetz/data/temporary_datasets/scale/ebird_mitchell_ghow_v2.csv", 1)
    
    # fix species column names if first run
    if(sp==1){
    wh = which(colnames(ebird)=="Circus hudsonius"):which(colnames(ebird)=="Asio otus")
    colnames(ebird)[wh] = gsub(" ","_",colnames(ebird)[wh])
    }
    
    # for migrants OR second GHOW run, cut ebird to may-aug
    if (sci.name %in% migrants | sp==2){
      ebird = ebird[ebird$day_of_year>120 & ebird$day_of_year<244,]
    }
    
    # For non-GHOW, cut ebird to range
    if(sp>2){
      print("processing ebird data")
      ebird = st_as_sf(ebird, coords = c("x", "y"), crs = crs(rangest))
      ebird.int = st_intersects(ebird, rangest, sparse=F)
      ebird = ebird[as.vector(ebird.int),]
      #plot(ebird[sample.int(nrow(ebird), nrow(ebird)/20),], max.plot=1)
      ebird = data.frame(sf:::as_Spatial(ebird))
      ebird$optional = NULL
      names(ebird)[(ncol(ebird)-1):ncol(ebird)] = c("x", "y")
      rm(ebird.int)
    }
    
    {# Column for presence/absence
    so = ebird[,sci_name]
    ebird$species_observed = ifelse(so >= 1, 1, 0)
    # Presence and absence datasets for plot
    pres = ebird[ebird$species_observed==1,]
    abs = ebird[ebird$species_observed==0,]
    pres_st = st_as_sf(pres, coords = c("x", "y"), 
                       crs = crs(rangest))
    }
    
    { # Plot occurrence data (figs 1c-1d)
      world = ne_countries(scale='medium', returnclass = 'sf') %>%
        st_transform(crs="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
      lim = st_bbox(pres_st)            
      # plot
      sprange2 = st_as_sf(sprange)
      map = ggplot() +
        geom_sf(data=world, bg="gray95") + 
        geom_sf(data=sprange2, bg="bisque1") +
        coord_sf(xlim=lim[c(1,3)], ylim=c(3000000,lim[4]), expand=F) +
        theme_bw() +
        theme(text = element_text(size=13),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(color = "black"),
              axis.text.x = element_text(color = "black", size=13),
              axis.text.y = element_text(color = "black", size=13),
              plot.title = element_text(size = 22, 
                                        face = "bold", hjust=.5)) +
        geom_point(ebird, mapping=aes(x=x, y=y), col="gray50", size=.2) +
        geom_point(pres, mapping=aes(x=x, y=y), 
                   col="black", bg="white", size=.5) +
        xlab("") +
        ylab("")
      ggsave(paste0("pointmaps/",sci_name,"_point_map_with_abs_",sp,".jpeg"), 
             map, "jpeg", height=8, width=7, units="in")
      }
    
    # load prediction surface
    print("processing prediction surface")
    {# prediction surface
      if(sp==1){
        pred_surface = read_csv("/gpfs/gibbs/pi/jetz/data/temporary_datasets/scale/pred_surface_mitchell_v2.csv")
        # adjust mean temp
        pred_surface$bio1 = pred_surface$bio1/10
        }else{
        pred_surface = read_csv("/gpfs/gibbs/pi/jetz/data/temporary_datasets/scale/pred_surface_mitchell_ghow_v2.csv")
      }
      # Populate with values
      pred_surface$duration_minutes = 60
      pred_surface$effort_distance_km = 1
      pred_surface$number_observers = 1
      pred_surface$year = 2020
      pred_surface$optional = NULL
    }
    
    # for non-GHOW model, cut prediction surface to range
    if(sp>2){
      pred_surface = st_as_sf(pred_surface, coords = c("x", "y"),
                              crs = crs(sprange))
      surface.int = st_intersects(pred_surface, rangest, sparse=F)
      pred_surface = pred_surface[as.vector(surface.int),]
      pred_surface = data.frame(sf:::as_Spatial(pred_surface))
      pred_surface$optional = NULL
      names(pred_surface)[(ncol(pred_surface)-1):ncol(pred_surface)] =
        c("x", "y")
      pred_surface$twi = NULL
      pred_surface = pred_surface[complete.cases(pred_surface),]
      rm(surface.int)
      gc()
    }
    
    { # Find max hr of day for spatial predictions (MIGHT USE TIME TO DUSK)
      hrs = ebird[,c("time_observations_started","species_observed")] %>%
        group_by(round(time_observations_started)) %>%
        summarize_at(vars(species_observed), mean)
      pred_surface$time_observations_started = as.numeric(
        hrs[which.max(hrs$species_observed),1])
      # set day for spatial predictions 
      pred_surface$day_of_year = 75
    }
    
    
    # loop models with and without GHOW
    for (i in 1:2){
      # no with-GHOW model for GHOW
      try({if(sci.name=="Bubo virginianus" & i==2){STOP}
        # for models with GHOW, add GHOW covariate
        if(i==2 & sci.name %in% migrants){preds = c(preds, "gh")}
        if(i==2 & !(sci.name %in% migrants)){preds = c(preds, "ghow")}
        
        # Splits
        print("splitting and filtering data")
        ebirdpts = st_as_sf(ebird, coords = c("x", "y"), crs = crs(rangest))
        {#training and testing samples
          sb = cv_spatial(x = ebirdpts,
                          column = "species_observed",
                          rows_cols = c(15, 15),
                          k = 4,
                          hexagon = FALSE,
                          selection = "systematic")
        }
        # limit folds to those with sufficient train and test records
        whk = which(sb$records$test_1>10 & sb$records$train_1>5)
        # loop folds
        for (k in whk){
          print(k)
          train1 = ebird[sb$biomod_table[,k]==T,]
          test = ebird[sb$biomod_table[,k]==F,]
          test = test[complete.cases(test),]
          # random split OOB from train data
          index = runif(nrow(train1))
          train = train1[which(index>.25),]
          train = train[complete.cases(train),]
          train.oob = train1[which(index<=.25),]
          train.oob = train.oob[complete.cases(train.oob),]
          
          {# Spatiotemporal thinning
            train$id = NULL
            coordinates(train) = c("x", "y")
            crs(train) = crs(sp_grd)
            over = over(train, sp_grd)
            train = cbind(data.frame(train), over)
            # get cell id and week number for each checklist
            train$optional = NULL
            checklist_cell = train %>% 
              mutate(cell = id,
                     week = ceiling(day_of_year/7))
            # sample one checklist per grid cell per week
            # sample detection/non-detection independently 
            train = checklist_cell %>% 
              group_by(species_observed, year, week, cell) %>% 
              sample_n(size = 1) %>% 
              ungroup()
          }
          {# Assemble TrainOOB
            train.oob$id = NULL
            coordinates(train.oob) = c("x","y")
            crs(train.oob) = crs(sp_grd)
            over = over(train.oob, sp_grd)
            train.oob = cbind(data.frame(train.oob), over)
            train$optional = NULL
            checklist_cell_oob = train.oob %>% 
              mutate(cell = id,
                     week = ceiling(day_of_year/7))
            train.oob = checklist_cell_oob %>% 
              group_by(species_observed, year, week, cell) %>% 
              sample_n(size = 1) %>% 
              ungroup()
          }
          
          {# remove oversampling for occurrence model, limit size
            train = train[!duplicated(train$checklist_id), ]
            # Binary Occurrence Response, coded as numeric
            train$pres_abs = as.numeric(train[[which(names(train) == "species_observed")]] > 0)
            # balanced sampling
            pos_fraction = mean(as.numeric(train$pres_abs))
            # Check that positives ARE the minority class
            if (pos_fraction > 0.5) pos_fraction = 1 - pos_fraction
            # Ranger binary response model requires code response as factor
            train$pres_abs = as.factor(as.numeric(train$pres_abs))}
          
          # Save points alone
          if(k==max(whk)){
            write_csv(train[,c("y","x","pres_abs")], 
                      paste0(sci_name,"/points/pts_binary_occurrences_",
                             sci_name,".csv"))
          }
          
          {#MODEL
            # Model Formula, threads
            m_formula = paste("pres_abs ~", paste(preds, collapse = "+"))
            n_threads = 7
            # Balanced Ranger Occurrence Model
            rf_occ = ranger::ranger(
              formula =  m_formula, 
              num.trees = 100, 
              importance = "impurity",
              num.threads = n_threads,
              respect.unordered.factors = "order",
              always.split.variables = NULL,
              probability = TRUE,
              replace = TRUE, 
              sample.fraction = c(pos_fraction, pos_fraction),
              data = train)
          }
          
          # Save model
          save(rf_occ, file=paste0("rf_",sci_name,"_",sp,"_",i,".RData"))
          
          # Load model if needed
          load(file=paste0("rf_",sci_name,"_",sp,"_",i,".RData"))
          
          # Partial dependence plots (fig 1g)
          vars = c("bio1","bio12","row_crop","non_row_crop","shannon_diversity")
          xlabs = c("Mean Annual Temperature (°C)", "Total Annual Precipitation",
                    vars[3:5], "Great Horned Owl Occurrence Prob.")
          if(sci.name %in% migrants & i == 2){vars = c(vars, "gh")}
          if(!(sci.name %in% migrants) & i == 2){vars = c(vars, "ghow")}
          if((k == max(whk))){
            for (pr in 1:length(vars)){
              par = pdp::partial(rf_occ, pred.var=vars[pr],
                                 trim.outliers = TRUE, chull = TRUE, parallel = TRUE,
                                 grid.resolution = 30,  paropts = list(.packages = "ranger"))
              names(par)[1] = "var"
              ggpd = ggplot(data=par, aes(x=var, y=yhat)) +
                geom_point() +
                geom_smooth() +
                ylab("Predicted occurrence") +
                xlab(xlabs[pr]) +
                theme_light() +
                theme(legend.position = "none",
                      text = element_text(size = 12))
              assign(paste0("pd.",pr), ggpd)
            }
            if(i==1){pd.6 = ggplot()}
            # save
            gga = ggarrange(pd.1, pd.2, pd.3, pd.4, pd.5, pd.6, ncol=3, nrow=2)
            ggsave(paste0("pdp/pd_plots_",sci_name,"_",sp,"_",i,".jpeg"), 
                   gga, "jpeg", height=8, width=11, units="in")
            
          }
          
          {# Test set predictive performance metrics (PPMS)
            test$occ_pred = predict(
              rf_occ, data = test, 
              type = "response",
              num.threads = n_threads)$predictions[,2]
            # Select Optimal Threshold on train OOB data
            trainOOB.pred = predict(rf_occ, data=train.oob, type="response")
            pa.df = data.frame("space.holder.tag",
                               obs = train.oob$species_observed, 
                               ppp = trainOOB.pred$predictions[,2] )
            pa.metrics = presence.absence.accuracy(pa.df,
                                                   threshold = quantile(pa.df$ppp,
                                                                        probs = seq(from=0, to=1, length=1000), na.rm =T),
                                                   na.rm = T, st.dev = F)
            pa.metrics$PCC = pa.metrics$Kappa = NULL
            pa.metrics$TSS = (pa.metrics$sensitivity + pa.metrics$specificity) - 1
            optimal_thresh_position = which.max(pa.metrics$sensitivity +
                                                  pa.metrics$specificity)
          }
          {# Compute Test set PPMs
            test.pred = predict(rf_occ, data=test, type="response")
            pa.df = data.frame("space.holder.tag",
                               obs = test$species_observed, 
                               ppp = test.pred$predictions[,2] )
            pa.metrics = presence.absence.accuracy(pa.df,
                                                   threshold = pa.metrics$threshold[ optimal_thresh_position ] ,
                                                   na.rm = T, st.dev = F)
            thresh = pa.metrics$threshold
            pa.metrics["PCC"] = pa.metrics["Kappa"] = NULL
            pa.metrics["TSS"] = (pa.metrics$sensitivity + pa.metrics$specificity) - 1
            assign(paste0("pa.metrics.",k), pa.metrics)
          }
          
          # Spatial predictions
          pred_surface$occ_pred = predict(
            rf_occ, data = pred_surface, type = "response",
            num.threads = n_threads)$predictions[,2]
          pred_surface$occ_thresh = ifelse(
            pred_surface$occ_pred>pa.metrics$threshold,1,0)
          # Range size?
          new_surface = pred_surface[pred_surface$occ_thresh==1,]
          assign(paste0("range.",k), nrow(new_surface))
          
          if(k==max(whk)){
            {# Save rasterized predictions
              coordinates(pred_surface) = c("x","y")
              ras = rasterize(pred_surface, grid, pred_surface$occ_pred, fun='first')
              writeRaster(ras, paste0(sci_name,"/predictions/ROR/PO_predictions_",
                                      sci_name,"_",sp,"_",i,".tif"),
                          overwrite=T)
              pred_ones = pred_surface[pred_surface$occ_thresh==1,]
              ras = rasterize(pred_ones, grid, pred_ones$occ_thresh, fun='first')
              writeRaster(ras, paste0(sci_name,"/predictions/PA/Thresholded_predictions_",
                                      sci_name,"_",sp,"_",i,".tif"),
                          overwrite=T)
              pred_surface = data.frame(pred_surface)
              pred_ones = data.frame(pred_ones)
              occ.probs = pred_surface[,c("x","y","occ_pred")]
              write_csv(occ.probs, paste0(sci_name,"/predictions/ROR/PO_predictions_",
                                          sci_name,"_",sp,"_",i,".csv"))
              occ.probs = pred_ones[,c("x","y","occ_thresh")]
              write_csv(occ.probs, paste0(sci_name,"/predictions/PA/Thresholded_predictions_",
                                          sci_name,"_",sp,"_",i,".csv"))
              rm(pred_ones)}
            
            {# Plot maps (fig. 1e)
              jpeg(paste0("maps/PO_spatial_predictions_",
                          sci_name,"_",sp,"_",i,".jpeg"),
                   width=500, height=500, units="px")
              # Plot prob of occurrence
              plotres = 500
              quilt.plot(
                x = pred_surface$x,
                y = pred_surface$y,
                z = pred_surface$occ_pred,
                nrow = plotres,
                ncol = plotres,
                na.rm = T,
                main = paste(common,"-",sci.name))
              plot(countries_proj, col="transparent", border="black",
                   usePolypath=F, add=T)
              dev.off()
              # Plot thresholded occurrence
              jpeg(paste0("maps/Thresholded_spatial_predictions_",
                          sci_name,"_",sp,"_",i,".jpeg"),
                   width=500, height=500, units="px")
              quilt.plot(
                x = pred_surface$x,
                y = pred_surface$y,
                z = pred_surface$occ_thresh,
                nrow = plotres,
                ncol = plotres,
                na.rm = T,
                add.legend = F,
                col = c("transparent", "darkgreen"),
                main = paste(common,"-",sci.name))
              plot(countries_proj, col="transparent", border="black",
                   usePolypath=F, add=T)
              dev.off()
            } # end plot
            
            # Averaging across folds
            for (w in whk){
              pa.metrics.w = get(paste0("pa.metrics.",w))[2:6]
              if(w==min(whk)){pa.metrics.table=pa.metrics.w}else{
                pa.metrics.table=rbind(pa.metrics.table, pa.metrics.w)}
              range.w = get(paste0("range.",w))
              if(w==min(whk)){ras=range.w}else{ras=c(ras, range.w)}
            }
            pa.metrics.summary = colMeans(pa.metrics.table, na.rm=T)
            range.summary = mean(ras, na.rm=T)
            
            # If GHOW model, save occurrence probs for ebird and pred surface
            if(sp < 3){
            if(sp == 1){
              ebird$ghow = predict(
                rf_occ, data = ebird, type = "response",
                num.threads = n_threads)$predictions[,2]
              pred_surface$ghow = pred_surface$occ_pred}
              # save
              ebird$species_observed = NULL
              write_csv(ebird, "/gpfs/gibbs/pi/jetz/data/temporary_datasets/scale/ebird_mitchell_ghow_full.csv")
              write_csv(pred_surface, "/gpfs/gibbs/pi/jetz/data/temporary_datasets/scale/pred_surface_mitchell_ghow_v2.csv")
              if(sp == 2){
                ebird$gh = predict(
                  rf_occ, data = ebird, type = "response",
                  num.threads = n_threads)$predictions[,2]
                pred_surface$gh = pred_surface$occ_pred
              # save
              ebird$species_observed = NULL
              write_csv(ebird, "/gpfs/gibbs/pi/jetz/data/temporary_datasets/scale/ebird_mitchell_ghow_breeding.csv")
              write_csv(pred_surface, "/gpfs/gibbs/pi/jetz/data/temporary_datasets/scale/pred_surface_mitchell_ghow_v2.csv")
            }
          }
            
            {# Species values
              modelnames = c("without GHOW", "with GHOW")
              model.output = data.frame(t(unlist(c(common, sci_name, sp, version, 
                                                   modelnames[i], nrow(train), 
                                                   nrow(train[train$species_observed==1,]),
                                                   nrow(train[train$species_observed==0,]),
                                                   range.summary, 
                                                   pa.metrics.summary[4],
                                                   pa.metrics.summary[c(1,5,2,3)],
                                                   rf_occ$variable.importance))))
              colnames(model.output)= c("common_name","sciname","ID",
                                        "version","model","noPts","noPres",
                                        "noAbs","range_area",
                                        "AUC","SPSthreshold","TSS",
                                        "Sensitivity","Specificity",
                                        paste0(preds,"_importance"))
            }
          } # end k=max(whk) only 
        } # end k loop
        write_csv(data.frame(model.output), 
                  paste0(sci_name,"_complete_",sp,"_",i,".csv"))
      }) # end try
    } # end with/without GHOW loop
  } # end sp
} # end workflow



