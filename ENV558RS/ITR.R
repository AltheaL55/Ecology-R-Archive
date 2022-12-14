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
#set working directory
setwd(dir = "Q:\\courseProject")

#read in LAS file
site_2725_pc <- lidR::readLAS("swisssurface3d_2020_2725-1117_2056_5728\\2725_1117.las")
#boundary shapefile
#site_2725_2725_shp <- raster::shapefile()
#plot point cloud
#plot(site_2725_pc)

# Defining dtm function arguments
site_2725_dtm <- rasterize_terrain(site_2725_pc,algorithm = knnidw())
# Plot and visualize the dtm in 3d 
plot_dtm3d(site_2725_dtm)

#Height normalization
site_2725_hnorm <- normalize_height(site_2725_pc, site_2725_dtm)

#CHM
# Defining chm function arguments
cs_chm <- 1 # output cellsize of the chm # Creating chm
site_2725_chm <- grid_canopy(site_2725_hnorm, cs_chm, p2r(na.fill = knnidw(k=3,p=2)))
##Smoothing
#smooth_ws <- 5 #smooth filter window size
#site_2725_schm <- rLiDAR::CHMsmoothing(site_2725_chm, filter = "mean", ws=smooth_ws)
#plot(site_2725_schm)

#Detecting tree tops
#Defining tree detection function arguments
ttop_lmf_ws <- 3 # treetop detection local maxima filter window size 
ttop_lmf_hmin <- 2 # minimum height threshold
site_2725_ttops <- rLiDAR::FindTreesCHM(site_2725_chm, ttop_lmf_ws, ttop_lmf_hmin)
#transform into spatial points
site_2725_ttops <- SpatialPointsDataFrame(site_2725_ttops[,1:2],site_2725_ttops)
#plot together
#plot(site_2725_chm)
#plot(site_2725_ttops,add = T)

#Export
# Exporting dtm in TIFF format 
writeRaster(site_2725_dtm,"CourseProject/scratch/site_27251plot_dtm.tif") 
# Exporting chm in TIFF format 
writeRaster(site_2725_chm,"CourseProject/scratch/site_27251plot_chm.tif", format = "GTiff")
# Exporting tree tops in SHAPEFILE format 
raster::shapefile(site_2725_ttops,"CourseProject/scratch/site_2725plot_ttops_ns.shp",overwrite=TRUE)

#Crowns
# Create crown map
crowns <- mcws(treetops = site_2725_ttops, CHM = site_2725_chm, minHeight = 3, verbose = FALSE)
# Plot crowns
plot(crowns, col = sample(rainbow(50), length(unique(crowns[])), replace = TRUE), legend = FALSE, xlab = "", ylab = "", xaxt='n', yaxt = 'n')
# Create polygon crown map
crownsPoly <- mcws(treetops = site_2725_ttops, CHM = site_2725_chm, 
                   format = "polygons", minHeight = 1.5, verbose = FALSE)
#plot
#plot(crownsPoly, border = "blue", lwd = 0.5, add = TRUE)
#save
writeOGR(crownsPoly, "Q:\\courseProject\\CourseProject\\scratch\\site_2725_crown_FT.shp", "treetops", driver = "ESRI Shapefile")
