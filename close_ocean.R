#funtion to find CODAP station IDs within a defined radius of a given USGS site code
library(geosphere)
library(tidyverse)
library(dataRetrieval)
#CODAP_unique <- distinct(CODAP_surface, Station_ID, .keep_all = TRUE) #sort CODAP surface data by unique Station IDs to more quickly sort Lat/Long
#max_dist <- 100*1000 #define max distance ocean station can be from USGS station, in meters

#Calculate the distances of all CODAP surface stations from USGS river station
#distm(c(lon1, lat1), c(lon2, lat2), fun = distHaversine)
#disttest <- distm(cbind(CODAP_unique$Longitude, CODAP_unique$Latitude), c(riv_site$dec_long_va, riv_site$dec_lat_va), fun = distHaversine) #distances in meters

close_ocean <- function(siteNo,max_dist) {
  riv_site <- whatWQPdata(siteNumber=siteNo)
  dists <- distm(cbind(CODAP_surface$Longitude, CODAP_surface$Latitude), c(riv_site$lon, riv_site$lat), fun = distHaversine) #distances in meters
  if(min(dists)<max_dist) {
    codap_close <- CODAP_surface$NEAA_ID[which(dists<max_dist)]
  } else {
    codap_close <- CODAP_surface$NEAA_ID[which(dists==min(dists))]
  }
    
  return(codap_close)
}
