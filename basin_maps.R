## Create map of all site basins using info from: https://waterdata.usgs.gov/blog/nldi-intro/

library(dataRetrieval)
library(nhdplusTools)
library(sf)
library(knitr)

##  Functions to generate base maps
GetURL <- function(service, host = "basemap.nationalmap.gov") {
  sprintf("https://%s/ArcGIS/rest/services/%s/MapServer/tile/{z}/{y}/{x}",
          host, service)
}

get_base_map <- function(options = leaflet::leafletOptions()) {
  map <- leaflet::leaflet(options = options)
  grp <- c("USGS Topo", "USGS Imagery Only", "USGS Imagery Topo",
           "USGS Shaded Relief", "Hydrography")
  att <- paste0("<a href='https://www.usgs.gov/'>",
                "U.S. Geological Survey</a> | ",
                "<a href='https://www.usgs.gov/laws/policies_notices.html'>",
                "Policies</a>")
  map <- leaflet::addTiles(map, urlTemplate = GetURL("USGSTopo"),
                           group = grp[1], attribution = att)
  map <- leaflet::addTiles(map, urlTemplate = GetURL("USGSImageryOnly"),
                           group = grp[2], attribution = att)
  map <- leaflet::addTiles(map, urlTemplate = GetURL("USGSImageryTopo"),
                           group = grp[3], attribution = att)
  map <- leaflet::addTiles(map, urlTemplate = GetURL("USGSShadedReliefOnly"),
                           group = grp[4], attribution = att)
  map <- leaflet::addTiles(map, urlTemplate = GetURL("USGSHydroCached"),
                           group = grp[5], attribution = att)
  map <- leaflet::hideGroup(map, grp[5])
  opt <- leaflet::layersControlOptions(collapsed = FALSE)
  map <- leaflet::addLayersControl(map, baseGroups = grp[1:4],
                                   overlayGroups = grp[5], options = opt)
}




## test for USGS-01021050
nldiURLs <- list(site_data = "https://labs.waterdata.usgs.gov/api/nldi/linked-data/nwissite/USGS-01021050",
                 basin_boundary = "https://labs.waterdata.usgs.gov/api/nldi/linked-data/nwissite/USGS-01021050/basin",
                 UT = "https://labs.waterdata.usgs.gov/api/nldi/linked-data/nwissite/USGS-01021050/navigation/UT/flowlines?distance=999",
                 UM = "https://labs.waterdata.usgs.gov/api/nldi/linked-data/nwissite/USGS-01021050/navigation/UM/flowlines?distance=999",
                 DM = "https://labs.waterdata.usgs.gov/api/nldi/linked-data/nwissite/USGS-01021050/navigation/DM/flowlines?distance=999",
                 UTwqp = "https://labs.waterdata.usgs.gov/api/nldi/linked-data/nwissite/USGS-01021050/navigation/UT/wqp?distance=999",
                 DMwqp = "https://labs.waterdata.usgs.gov/api/nldi/linked-data/nwissite/USGS-01021050/navigation/DM/wqp?distance=999")

nldi_data <- list()

for(n in names(nldiURLs)) {
  nldi_data[n] <- list(sf::read_sf(nldiURLs[n][[1]]))
  print(paste(n, "is of class", class(nldi_data[[n]]), "and has", nrow(nldi_data[[n]]), "features"))
}

map <- get_base_map()

map <- leaflet::addPolygons(map,
                            data=nldi_data$basin_boundary,
                            color = "red",
                            fill = FALSE,
                            weight = 1,
                            opacity = 1)

print(map)


## Create loop to add all basins to map; plot stations as points if no basin boundary information available
map <- get_base_map()
counter <- 0
for(j in names_final){
  counter <- counter+1
  temp_site <- j
  siteINFO_temp <- readNWISsite(sub("USGS-", "", j))
  nldi_nwis <- list(featureSource = "nwissite", featureID = j)
  #site <- get_nldi_feature(nldi_nwis)
  basin <- get_nldi_basin(nldi_feature = nldi_nwis)
  if(is.null(basin) == "FALSE"){
    map <- leaflet::addPolygons(map,
                              data=basin,
                              color = "blue",
                              fill = "blue",
                              weight = 1,
                              opacity = 1)
    map <- leaflet::addCircleMarkers(map,siteINFO_temp$dec_long_va,siteINFO_temp$dec_lat_va,
                                     fillColor = "black", 
                                     radius=3, stroke=FALSE,
                                     fillOpacity = 0.8, opacity = 0.8,
                                     popup=siteINFO_temp$station_nm)
  }else {
    
    map <- leaflet::addCircleMarkers(map,siteINFO_temp$dec_long_va,siteINFO_temp$dec_lat_va,
                       fillColor = "red", 
                       radius=3, stroke=FALSE,
                       fillOpacity = 0.8, opacity = 0.8,
                       popup=siteINFO_temp$station_nm)
    
  }
}
print(map)

# ## Create loop to add all basins to map; plot stations as points if no basin boundary information available
# map <- get_base_map()
# counter <- 0
# for(j in names_recent){
#   counter <- counter+1
#   temp_site <- j
#   siteINFO_temp <- readNWISsite(sub("USGS-", "", j))
#   nldi_nwis <- list(featureSource = "nwissite", featureID = j)
#   #site <- get_nldi_feature(nldi_nwis)
#   basin <- get_nldi_basin(nldi_feature = nldi_nwis)
#   if(is.null(basin) == "FALSE"){
#     map <- leaflet::addPolygons(map,
#                                 data=basin,
#                                 color = "blue",
#                                 fill = "blue",
#                                 weight = 1,
#                                 opacity = 1)
#     map <- leaflet::addCircleMarkers(map,siteINFO_temp$dec_long_va,siteINFO_temp$dec_lat_va,
#                                      fillColor = "black", 
#                                      radius=3, stroke=FALSE,
#                                      fillOpacity = 0.8, opacity = 0.8,
#                                      popup=siteINFO_temp$station_nm)
#   }else {
#     
#     map <- leaflet::addCircleMarkers(map,siteINFO_temp$dec_long_va,siteINFO_temp$dec_lat_va,
#                                      fillColor = "red", 
#                                      radius=3, stroke=FALSE,
#                                      fillOpacity = 0.8, opacity = 0.8,
#                                      popup=siteINFO_temp$station_nm)
#     
#   }
# }
# print(map)
   