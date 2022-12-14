#For course project 558
#Elsie Liu, Fall 2022
#install.packages("rLiDAR")
#install.packages("lidR")
library(rLiDAR) 
library(lidR)
library(sp)
library(raster)
library(ForestTools)
library(rgdal)
library(stringr)
library(terra)
#set working directory
setwd(dir = "Q:\\courseProject\\LASData")
#read in file list
LASFiles <- list.files()
#set empty tables and vectors
plot_stats <- data.frame(sample = 1:15)
row.names(plot_stats) <- c('elev_mean','elev_min','elev_max','elev_SD','slope_mean','slope_SD','aspect_mean','aspect_SD','TreeCount','CrownArea_mean',
                           'CrownArea_SD', 'TreeH_mean', 'TreeH_min','TreeH_max','TreeH_SD')
site_stats <- 1:15

#loop start
for(i in 1:length(LASFiles)){
  site_stats <- 1:15
  #read in LAS file
  sitename <- substr(LASFiles[i],1,9)
  site_pc <- lidR::readLAS(LASFiles[i])
  # Defining dtm function arguments
  cs_dtm <- 1 # output cellsize of the dtm
  # Creating dtm
  #site_dtm <- grid_terrain(site_ground_points, cs_dtm, knnidw()) 
  site_dtm <- rasterize_terrain(site_pc,algorithm = knnidw())
  #Height normalization
  site_hnorm <- normalize_height(site_pc, site_dtm)
  
  #export terrain stats
  dtm_value <- values(site_dtm)
  site_stats[1:4] <- c(mean(dtm_value), min(dtm_value), max(dtm_value), sd(dtm_value))
  slope <- terra::terrain(site_dtm, "slope", unit='degree')
  site_stats[5:6] <- c(mean(values(slope,na.rm=T)),sd(values(slope,na.rm=T)))
  aspect <- terra::terrain(site_dtm, "aspect", unit='degree')
  site_stats[7:8] <- c(mean(values(aspect,na.rm=T)),sd(values(aspect,na.rm=T)))
  
  
  #CHM
  # Defining chm function arguments
  cs_chm <- 1 # output cellsize of the chm # Creating chm
  site_chm <- grid_canopy(site_hnorm, cs_chm, p2r(na.fill = knnidw(k=3,p=2)))
  
  #Detecting tree tops
  #Defining tree detection function arguments
  ttop_lmf_ws <- 3 # treetop detection local maxima filter window size 
  ttop_lmf_hmin <- 2 # minimum height threshold
  site_ttops <- rLiDAR::FindTreesCHM(site_chm, ttop_lmf_ws, ttop_lmf_hmin)
  #transform into spatial points
  site_ttops <- SpatialPointsDataFrame(site_ttops[,1:2],site_ttops)
  #Export
  # Exporting dtm in TIFF format 
  dtm_dir <-paste("Q:\\courseProject\\CourseProject\\scratch\\site", sitename, "_dtm.tif", sep = "")
  writeRaster(site_dtm, dtm_dir,overwrite=TRUE) 
  # Exporting chm in TIFF format 
  chm_dir <- paste("Q:\\courseProject\\CourseProject\\scratch\\site", sitename, "_chm.tif", sep = "")
  writeRaster(site_chm,chm_dir, format = "GTiff",overwrite=TRUE)
  # Exporting tree tops in SHAPEFILE format
  ttop_dir <- paste("Q:\\courseProject\\CourseProject\\scratch\\site", sitename, "_ttops_ns.shp", sep = "")
  raster::shapefile(site_ttops, ttop_dir,overwrite=TRUE)
  
  #Crowns
  # Create crown map
  crowns <- mcws(treetops = site_ttops, CHM = site_chm, minHeight = 3, verbose = FALSE)
  # Create polygon crown map
  crownsPoly <- mcws(treetops = site_ttops, CHM = site_chm, 
                     format = "polygons", minHeight = 1.5, verbose = FALSE)
  #save crown layer
  crown_dir <- paste("Q:\\courseProject\\CourseProject\\scratch\\site", sitename, "_crown.shp", sep = "")
  writeOGR(crownsPoly, crown_dir, "treetops", driver = "ESRI Shapefile",overwrite=TRUE)
  
  #export tree stats
  site_stats[9:15] <- sp_summarise(crownsPoly, variables = c("crownArea", "height"))[c('TreeCount','crownAreaMean','crownAreaSD','heightMean','heightMin','heightMax','heightSD'),1]
  #append dataframe
  site_stats <- data.frame(site_stats)
  colnames(site_stats) <- sitename
  plot_stats <- cbind(plot_stats,site_stats)
}
plot_stats<- plot_stats[,-1]
write.table(plot_stats, file = "Q:\\courseProject\\R\\stats.csv", row.names = T, col.names = T)

###customized stats function
# quant98 <- function(x, ...) quantile(x, c(.98), na.rm = TRUE)
# #To have this function applied using sp_summarise, it must be put into a named list. 
# #Naming the functions in the list is needed for labeling the function’s outputs.
# # Create list of functions
# custFuns <- list(quant98, max)
# names(custFuns) <- c("98thQuantile", "Max")
# # Generate statistics for crown areas and tree heights
# sp_summarise(crownsPoly, variables = c("crownArea", "height"), statFuns = custFuns)

###Avalaible functions for dtm
# hist(site_dtm)
#plot(aspect,main="Aspect for CALI [degrees]", col=rainbow(25,alpha=0.7))
