# August 3 2022 code for NEAA analysis
# SRP
#setwd("C:/Users/spacella/OneDrive - Environmental Protection Agency (EPA)/R Code May 2022 SRP Review")
setwd("C:/Users/spacella/OneDrive - Environmental Protection Agency (EPA)/National Estuary Acidification Assessment/NEAA Analysis Files 240814/R Code May 2022 SRP Review")
#Download all river data and write it to folder
source("neaa_final_download_231028.R") #SRP checked 10/28/2023
source("neaa_corrections.R")#SRP checked 8/25/2022
source("ts_analysis_220816.R")
source("neaa_calcium.R")
source("export_neaa_csv_calcium.R") #writes calcium data to neaa_exports folder

#Download ocean data and organize
source("codap_load.R")
source("codap_organize.R")
source("export_ocean_csv.R")

# Make map of river basins
source("basin_maps.R")

# # Map of river and coastal ocean stations
# leaflet(data=siteINFO_final) %>% 
#   addProviderTiles("CartoDB.Positron") %>%
#   addCircleMarkers(~dec_long_va,~dec_lat_va,
#                    fillColor = "red", 
#                    radius=6, stroke=FALSE,
#                    fillOpacity = 0.8, opacity = 0.8,
#                    popup=~station_nm) %>%
#   addCircleMarkers(ocean_coords_all$Longitude,ocean_coords_all$Latitude,
#                    fillColor = "black", 
#                    radius=3, stroke=FALSE,
#                    fillOpacity = 0.8, opacity = 0.8,
#                    popup=ocean_coords_all$NEAA_ID) %>%
#   addCircleMarkers(ocean_coords$Longitude,ocean_coords$Latitude,
#                    fillColor = "blue", 
#                   radius=4, stroke=FALSE,
#                     fillOpacity = 0.8, opacity = 0.8,
#                    popup=ocean_coords$NEAA_ID)

 

# leaflet(data=siteINFO_final) %>% 
#   addProviderTiles("CartoDB.Positron") %>%
#   addCircleMarkers(~dec_long_va,~dec_lat_va,
#                    fillColor = "blue", 
#                    radius=6, stroke=FALSE,
#                    fillOpacity = 0.8, opacity = 0.8,
#                    popup=~station_nm) %>%
#   addCircleMarkers(siteINFO_recent$dec_long_va,siteINFO_recent$dec_lat_va,
#                    fillColor = "orange", 
#                    radius=6, stroke=FALSE,
#                    fillOpacity = 0.8, opacity = 0.8,
#                    popup=siteINFO_recent$station_nm) %>%
#   addCircleMarkers(siteINFO_2000$dec_long_va,siteINFO_2000$dec_lat_va,
#                    fillColor = "red", 
#                    radius=6, stroke=FALSE,
#                    fillOpacity = 0.8, opacity = 0.8,
#                    popup=siteINFO_2000$station_nm) 


