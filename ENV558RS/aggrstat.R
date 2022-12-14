library(rLiDAR) 
library(lidR)
library(sp)
library(raster)
library(ForestTools)
library(rgdal)
library(stringr)
library(terra)
library(tidyr)

#Focal Stats
LASFiles <- list.files("Q:\\courseProject\\LASData")
plot_stats_50 <- data.frame(t(data.frame(sample = 1:10)))
colnames(plot_stats_50) <- c('elev_mean','slope_mean','aspect_mean','TreeCount','CrownArea_mean',
                           'CrownArea_SD', 'TreeH_mean', 'TreeH_min','TreeH_max','TreeH_SD')

for(i in 1:20){
  sitename <- substr(LASFiles[i],1,9)
  site_dtm <- raster(paste("Q:\\courseProject\\CourseProject\\scratch\\site", sitename, "_dtm.tif", sep = ""))
  #aggregate
  aggr_dtm <- aggregate(site_dtm, 50, fun=mean, expand=TRUE, na.rm=TRUE)
  crs(aggr_dtm) <- "+init=epsg:2056 +proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333 +k_0=1 +x_0=2600000 +y_0=1200000 +ellps=bessel +towgs84=674.374,15.056,405.346,0,0,0,0 +units=m +no_defs"
  #get polygon
  # dtm_poly <- rasterToPolygons(aggr_dtm, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE)
  # writeOGR(dtm_poly,paste("Q:\\courseProject\\CourseProject\\scratch\\site", sitename, "focdtm50.shp", sep = ""),'dtm_50',driver = "ESRI Shapefile", overwrite=T)
  # 
  #save value
  ag_dtm_value <- values(aggr_dtm)
  aggr_stats[1:4] <- c(mean(ag_dtm_value,na.rm=T), min(ag_dtm_value,na.rm=T), max(ag_dtm_value,na.rm=T), sd(ag_dtm_value,na.rm=T))
  #slope
  aggr_slope <- terra::terrain(aggr_dtm, "slope", unit='degrees')
  aggr_aspect <- terra::terrain(aggr_dtm, "aspect", unit='degrees')
  
  
  ###for vectors
  # Compute tree height statistics within a 50 m x 50 m cell grid
  crownA <- readOGR(paste("Q:\\courseProject\\CourseProject\\scratch\\site", sitename, "_crown.shp", sep = ""))
  gridStats <- sp_summarise(crownA, grid = 50, variables =  c("height", "crownArea"))
  # View layer names
  #names(gridStats)
  # [1] "TreeCount"       "heightMean"      "crownAreaMean"   "heightMedian"    "crownAreaMedian" "heightSD"        "crownAreaSD"    
  # [8] "heightMin"       "crownAreaMin"    "heightMax"       "crownAreaMax"   
  
  #loop through each pixel
  for (v in 1:400) {
    pix_stats <- 1:10
    pix_stats[1] <- values(aggr_dtm)[v]
    pix_stats[2] <- values(aggr_slope)[v]
    pix_stats[3] <- values(aggr_aspect)[v]
    pix_stats[4] <- values(gridStats[[1]])[v]
    pix_stats[5] <- values(gridStats[[3]])[v]
    pix_stats[6] <- values(gridStats[[7]])[v]
    pix_stats[7] <- values(gridStats[[2]])[v]
    pix_stats[8] <- values(gridStats[[8]])[v]
    pix_stats[9] <- values(gridStats[[10]])[v]
    pix_stats[10] <- values(gridStats[[6]])[v]
    
    pix_stats <- t(data.frame(pix_stats))
    colnames(pix_stats) <- c('elev_mean','slope_mean','aspect_mean','TreeCount','CrownArea_mean',
                                 'CrownArea_SD', 'TreeH_mean', 'TreeH_min','TreeH_max','TreeH_SD')
    plot_stats_50 <- rbind(plot_stats_50, pix_stats)
  }
  
  }
plot_stats_50 <- plot_stats_50[-1,]
write.table(plot_stats_50, file = "Q:\\courseProject\\R\\stats_50.csv", row.names = T, col.names = T)
aggr_nna<- drop_na(plot_stats_50)
write.table(aggr_nna, file = "Q:\\courseProject\\R\\stats_50_NNA.csv", row.names = T, col.names = T)
