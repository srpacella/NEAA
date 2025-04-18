

# Clear workspace
rm(list = ls())

install.packages("gpclib")

write.delim <- function(df, file, quote=FALSE, row.names=FALSE,sep='\t',...){
  write.table(df,file,quote=quote,row.names=row.names,
              sep=sep,...)
}

library(dataRetrieval); library(dplyr); library(svDialogs); library(purrr); 
library(zoo); library(leaflet); library(mapview); library(sp); library(rgdal);
library(rgeos);library(maptools);library(ggplot2)
# Prompt parameter
#parameter <- as.character(dlgInput("Enter the parameter, ex: alk, pH, temp",
#                                Sys.info()["user"])$res)

setwd("/Users/spacella/Desktop/R Code May 2022 SRP Review")

# Timestamp is written to data file
time = Sys.time()

# Create parameter values
#pCode <- paste("pCode_",parameter,sep="")

# Specify output directory
drctry = paste("/Users/spacella/Desktop/R Code May 2022 SRP Review")

# Specify output txt file name
out_file <- paste("Preliminary_station_metadata_all",".txt",sep="")

# Specify output image file name
image_name <- paste("USGS_all", "_timeseries_stations.png", sep="")

# Set state parameters codes
coastal_states <- c("ME","NH","MA","RI","CT","NY","NJ","DE","PA","MD","VA","NC",
                    "SC","GA","FL","MS","AL","LA","TX","CA","OR","WA","AK")

# Search USGS NWIS database for all coastal states with these param codes
# Alk parameters include titration alkalinity and HCO3-
#if(parameter == "alk"){
#pCode <- c("00410", "00416", "00417", "00418", "00425", "00431", "00440",
#             "00449", "00450", "00453", "29801", "29802", "29803", "29804",
#             "39036", "39086", "39087", "90410", "95410", "95440", "99440",
#             "00409", "00411", "00413", "00415", "00419", "00421", "00430",
#             "00445", "00446", "00447", "00448", "00451", "00452", "29797",
#             "29798", "29805", "29806", "29813", "46005", "63786", "63787",
#             "90440", "99431", "99432")
#} else if(parameter == "pH"){
#    pCode <- c("00400") 
#} else if(parameter == "temp"){
#      pCode <- c("00010")
#}

#start.date <- "1900-01-01"; end.date <- "2021-12-31" 

# Search database by state using parameter code
wq_meta_data <- list()
for(j in coastal_states){
  wq_meta_data[[j]] <- whatWQPdata(stateCd=j,
                                   CharacteristicName = c("pH"),ProviderName = "NWIS")
}
# Place list into dataframe
wq_meta_data <- do.call("rbind", wq_meta_data)

# Count instances of provider and type
provider_name <- wq_meta_data %>% dplyr::count(ProviderName)
monitoring_type <- wq_meta_data %>% dplyr::count(MonitoringLocationTypeName)

# Remove all instances other than stream and river
wq_meta_data <- subset(wq_meta_data, MonitoringLocationTypeName %in% c("Stream", "River/Stream"))

# Only keep at least 50 instances
wq_meta_data <- subset(wq_meta_data, resultCount >= 50)

# Minimize distance to the coastline
# function to extract lat/long from shapefile
fortify.shape <- function(x){
  x@data$id <- rownames(x@data)
  x.f <- fortify(x, region = "id")
  x.join <- inner_join(x.f, x@data, by = "id")
}

subset.shape <- function(x, domain){
  x.subset <- filter(x, long > domain[1] & 
                       long < domain[2] & 
                       lat > domain[3] & 
                       lat < domain[4])
  x.subset
}

# high-resolution shapefile of the coast
coast.shp <- readOGR(dsn=".", layer="cb_2018_us_nation_20m") 

# Fortify the shapefile data using `fortify.shape()`:
dat.coast <- fortify.shape(coast.shp) 

# Set map domain (Alaska -> Florida)
domain <- c(-180, -65, 20, 75)

# Extract the coastline data for the desired domain using `subset.shape()`:
dat.coast.wc <- subset.shape(dat.coast, domain)
dat.coast.wc <- dat.coast.wc[,-c(3:12)]


# Pair down USGS dataframe and retrieve just lat/long
wq_latlon <- wq_meta_data[,-c(1,4:16)]
wq_latlon <- wq_latlon[,c(2,1)]

# Minimize distance to coastline
# This will take ~10 minutes!! 
for(i in 1:nrow(wq_latlon)){  
  distances<-geosphere::distGeo(wq_latlon[i,], dat.coast.wc)/1000
  ranking<-rank(distances, ties.method = "first")
  wq_latlon$shortest[i]<-which(ranking ==1)
  wq_latlon$shortestD[i]<-distances[wq_latlon$shortest[i]]
}

# Convert to separate dataframe (so we don't need to re-run the minimization code)
# add distances to the metadata dataframe
wq_sites <- cbind(wq_latlon, wq_meta_data)

D = 50 #distance in km

# Only keep stations within x-km of the coastline
wq_sites$shortestD[wq_sites$shortestD > D] <- NA #adjust distance here!
wq_sites <- na.omit(wq_sites)

#Remove sites with "RCE WRP" Organization Identifier (throws error when downloading data)
wq_sites$OrganizationIdentifier[wq_sites$OrganizationIdentifier == "RCE WRP"] <- NA 
wq_sites <- na.omit(wq_sites)
wq_sites_usgs <- wq_sites[wq_sites$ProviderName=="NWIS",] #only USGS sites for now


#Select sites closest to coastline for each HUC code


#Select stations with most data based on "activityCount" since we have already filtered to be <50km from coastline
huc_unique <- unique(wq_sites$HUCEightDigitCode)

huc_max <- list();
for(i in huc_unique){
  huc_temp <- wq_sites[wq_sites$HUCEightDigitCode==i,]
  huc_temp_max <- huc_temp[which.max(huc_temp$activityCount),]
  huc_max[[i]] <- huc_temp_max
}
huc_max_test <- do.call("rbind", huc_max)

leaflet() %>%
  addTiles() %>%
  addCircleMarkers(data=huc_max_test, lat=~lat, lng=~lon,
                   stroke=FALSE, color="#FFBF00", opacity=0.6, 
                   fillOpacity=0.6,radius=6,group="Filtered stations")

sites <- huc_max_test[,c(11)]

# Download data from selected sites for the selected variables
#Start with only pH data
pH_1 <- list();
for(j in sites){
  pH_1[[j]] <-  readWQPdata(siteNumber=j, 
                            CharacteristicName = c("Alkalinity", "Bicarbonate",
                                                   "pH", "Temperature, water"))
}

#Create a dataframe with all data
download_all <- do.call("rbind",pH_1)

#Create list with ph, alkalinity, temperature bound by date
all_datebind <- list()
for(j in sites){
  temp_all <- pH_1[[j]]
  pH_temp <- temp_all[temp_all$CharacteristicName=="pH",]
  alk_temp1 <- temp_all[temp_all$CharacteristicName=="Alkalinity",]
  alk_temp2 <- temp_all[temp_all$CharacteristicName=="Bicarbonate",]
  alk_temp <- rbind(alk_temp1,alk_temp2)
  temp_temp <- temp_all[temp_all$CharacteristicName=="Temperature, water",]
  df_join1 <- inner_join(pH_temp,alk_temp1,by = "ActivityStartDate")
  all_datebind[[j]] <- df_join1
  print(ggplot() + 
    geom_point(data=df_join1, aes(x=ActivityStartDate, y=ResultMeasureValue.x), color='green') + 
    geom_point(data=df_join1, aes(x=ActivityStartDate, y=ResultMeasureValue.y), color='red'))
  
}

plot <- ggplot() + 
  geom_point(data=df_join1, aes(x=ActivityStartDate, y=ResultMeasureValue.x), color='green') + 
  geom_point(data=df_join1, aes(x=ActivityStartDate, y=ResultMeasureValue.y), color='red')

print(plot)