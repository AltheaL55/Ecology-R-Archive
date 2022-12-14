setwd("Q:\\courseProject\\CourseProject\\scratch\\")
library(raster)
library(rgdal)
library(ForestTools)
for (i in 1:length(LASFiles)){
  site_stats <- 1:15
  #read in LAS file
  sitename <- substr(LASFiles[i],1,9)
  site_dtm <- raster(paste("site", sitename, "_dtm.tif", sep = ""),crs=)
  crs(site_dtm) <- "+init=epsg:2056 +proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333 +k_0=1 +x_0=2600000 +y_0=1200000 +ellps=bessel +towgs84=674.374,15.056,405.346,0,0,0,0 +units=m +no_defs"
  #export terrain stats
  dtm_value <- values(site_dtm)
  site_stats[1:4] <- c(mean(dtm_value,na.rm=T), min(dtm_value,na.rm=T), max(dtm_value,na.rm=T), sd(dtm_value,na.rm=T))
  slope <- terra::terrain(site_dtm, "slope", unit='degrees')
  site_stats[5:6] <- c(mean(values(slope),na.rm=T),sd(values(slope),na.rm=T))
  aspect <- terra::terrain(site_dtm, "aspect", unit='degrees')
  site_stats[7:8] <- c(mean(values(aspect),na.rm=T),sd(values(aspect),na.rm=T))
  
  crownsPoly <- readOGR(paste("Q:\\courseProject\\CourseProject\\scratch\\site", sitename, "_crown.shp", sep = ""))
  #export tree stats
  site_stats[9:15] <- sp_summarise(crownsPoly, variables = c("crownArea", "height"))[c('TreeCount','crownAreaMean','crownAreaSD','heightMean','heightMin','heightMax','heightSD'),1]
  #append dataframe
  site_stats <- data.frame(site_stats)
  colnames(site_stats) <- sitename
  plot_stats <- cbind(plot_stats,site_stats)
  }
plot_stats<- plot_stats[,-1]
write.table(plot_stats, file = "Q:\\courseProject\\R\\stats.csv", row.names = T, col.names = T)




