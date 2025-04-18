## Uses final sites chosen 6/23/2022 plus adds additional sites from 8/25/22, downloads and organizes data
# Combines SRP and WF codes from neaa_final_download.R & USGS Download_220802.Rmd

# SRP updates to include additional sites as of 8/18/2023 review
# USGS-11303500 (San Joaquin River)

# 10/28/23 SRP Update: Only include observations with contemporaneous conductivity data

rm(list = ls())

library(dataRetrieval)
library(dplyr)
#library(rgeos)
library(leaflet)
#library(rgdal)
library(purrr); 
library(zoo); 
library(mapview); 
library(sp); 
#library(rgeos);
#library(maptools);
library(ggplot2);
library(tools);
library(lubridate);
library(readxl)
library(seacarb)
library(mblm)
library(leaflet)
library(wql)
library(birk)
library(EGRET);
library(tidyr)

#setwd("C:/Users/spacella/OneDrive - Environmental Protection Agency (EPA)/R Code May 2022 SRP Review")
setwd("C:/Users/spacella/OneDrive - Environmental Protection Agency (EPA)/National Estuary Acidification Assessment/NEAA Analysis Files 240814/R Code May 2022 SRP Review")

# Load in names of sites to use for final analysis  
final_sites <- read_excel("WQP_sites_manual_srp_update6.xlsx",col_names=FALSE) #reflects final stations as of 8/18/2023
final_sites_test <- data.frame(final_sites)
final_sites_test <- subset(final_sites_test, ...1!="USGS-02049500")
final_sites_names  <- c(final_sites_test[,1])
final_sites_names  <- c(final_sites_test[,1])
# Pull metadata for NEAA sites
siteNumbers <- sub("USGS-", "", final_sites_names)
siteINFO <- readNWISsite(siteNumbers)

download_neaa <- list()
for(j in final_sites_names){ 
  download_neaa[[j]] <-  readWQPdata(siteNumber=j, 
                                     characteristicName = c("Alkalinity",
                                                            "pH", 
                                                            "Temperature, water",
                                                            "Specific conductance"))
}

all_data_neaa <- do.call(rbind, download_neaa)

# Clean
all_data_neaa <- all_data_neaa[,c(1:9,21,22,27:29,32,34,35,46,63,64)]
all_data_neaa$ID <- paste(all_data_neaa$ActivityStartDate, all_data_neaa$MonitoringLocationIdentifier)
all_data_neaa <- dplyr::distinct(all_data_neaa)

# Specify date
all_data_neaa$Datetime <- ymd(all_data_neaa$ActivityStartDate)

# Only keep water media (sometimes biological tissue or sediment gets through)
media_type <- all_data_neaa %>% dplyr::count(ActivityMediaName)
all_data_neaa <- subset(all_data_neaa, ActivityMediaName %in% c("Water"))

# pH
pH <- all_data_neaa  %>%
  filter(CharacteristicName == "pH")
pH$ResultMeasureValue <- as.numeric(as.character(pH$ResultMeasureValue))
pH$ResultMeasureValue[pH$ResultMeasureValue < 5.5 | pH$ResultMeasureValue > 9.5] <- NA
pH <- na.omit(pH)

# Order by site ID and date
pH <- pH[order(pH$MonitoringLocationIdentifier,pH$ActivityStartDate),]

# Find duplicates samplings, and only keep field measurements when available
pH <- pH %>% group_by(ID) %>% filter(USGSPCode == min(USGSPCode))
# Keep first pH value for each day when there are multiple times
pH <- pH %>% group_by(ID) %>% filter(ActivityStartTime.Time == min(ActivityStartTime.Time))
pH <- dplyr::distinct(pH) #remove duplicate rows

# temperature
temperature <- all_data_neaa %>%
  filter(CharacteristicName == "Temperature, water")
temperature$ResultMeasureValue <- as.numeric(as.character(temperature$ResultMeasureValue))
temperature$ResultMeasureValue[temperature$ResultMeasureValue > 50 | temperature$ResultMeasureValue < 0] <- NA
temperature <- na.omit(temperature)
temperature <- temperature[order(temperature$MonitoringLocationIdentifier,temperature$ActivityStartDate),]
# Keep first temp value for each day when there are multiple times
temperature <- temperature %>% group_by(ID) %>% filter(ActivityStartTime.Time == min(ActivityStartTime.Time))
temperature <- dplyr::distinct(temperature) #remove duplicate rows


# alkalinity
alk <- all_data_neaa %>%
  filter(CharacteristicName == "Alkalinity")
alk$ResultMeasureValue <- as.numeric(as.character(alk$ResultMeasureValue))
alk$ResultMeasureValue <- as.numeric(alk$ResultMeasureValue*20)
alk$ResultMeasureValue[alk$ResultMeasureValue > 9000 | alk$ResultMeasureValue < 50] <- NA
alk <- na.omit(alk)
alk <- alk[order(alk$MonitoringLocationIdentifier,alk$ActivityStartDate),]
# Find duplicates samplings, and only keep field measurements when available
alk <- alk %>% group_by(ID) %>% filter(USGSPCode == min(USGSPCode))
# Keep first alk value for each day when there are multiple times
alk <- alk %>% group_by(ID) %>% filter(ActivityStartTime.Time == min(ActivityStartTime.Time))
alk <- dplyr::distinct(alk) #remove duplicate rows

# conductance
cond <- all_data_neaa %>%
  filter(CharacteristicName == "Specific conductance")
cond$ResultMeasureValue <- as.numeric(as.character(cond$ResultMeasureValue))
cond$ResultMeasureValue[cond$ResultMeasureValue > 5000 | cond$ResultMeasureValue < 0] <- NA
cond <- na.omit(cond)
cond <- cond[order(cond$MonitoringLocationIdentifier,cond$ActivityStartDate),]
# Keep first temp value for each day when there are multiple times
cond <- cond %>% group_by(ID) %>% filter(ActivityStartTime.Time == min(ActivityStartTime.Time))
cond <- dplyr::distinct(cond) #remove duplicate rows


#########

# Merge all temperature, conductivity, alkalinity, and pH data
# Merge together by ID
# join1 <- left_join(pH, temperature, by="ID")
# join2 <- na.omit(left_join(join1, alk, by="ID"))
# join3 <- na.omit(left_join(join2, cond, by="ID"))
# join3 <- join3 %>% filter(as.Date(ActivityStartDate.x) > "1950-01-01")

join1_test <- inner_join(pH,temperature, by="ID")
join2_test <- inner_join(join1_test, alk ,by="ID")
join3_test <- inner_join(join2_test, cond ,by="ID")
join3_test <- join3_test %>% filter(as.Date(ActivityStartDate.x) > "1950-01-01")

##put all data frames into list
#df_list_test <- list(pH, alk, temperature, cond)      

##merge all data frames together
#join_final <- df_list_test %>% reduce(inner_join, by='ID')
#join_final <- join_final %>% filter(as.Date(ActivityStartDate.x) > "1950-01-01")


# Reduce and rename columns
final_join <- join3_test[,c(21,1,7,11,16,18,20,38,40,42,59,61,63,80,82,84)]
names(final_join) <- c("ID", "OrganizationIdentifier",
                       "ActivityStartDate","MonitoringLocation",
                       "pH", "pHParameterCode", "pHTimeStamp", 
                       "Temperature", "TemperatureParameterCode", "TemperatureTimeStamp", 
                       "Alkalinity", "AlkalinityParameterCode",
                       "AlkalinityTimeStamp","Conductivity","ConductivityParameterCode","ConductivityTimeStamp")


# Find the number of stations remaining after harmonizing
site_unique_joined <- as.data.frame(unique(final_join$MonitoringLocation))


# Create list of data for each site
all_data_site <- list()
for(j in final_sites_names){ 
  data_site <- final_join %>% 
    filter(MonitoringLocation == j)
  
  # Detect and remove outliers using 3x IQR method
  alk_temp <- data_site$Alkalinity
  pH_temp <- data_site$pH
  temp_temp <- data_site$Temperature
  
  alk_summary <- summary(alk_temp)
  alk_IQR <- IQR(alk_temp)
  alk_min <- alk_summary[2]-(3*alk_IQR) 
  alk_max = alk_summary[5]+(3*alk_IQR) 
  
  pH_summary <- summary(pH_temp)
  pH_IQR <- IQR(pH_temp)
  pH_min <- pH_summary[2]-(3*pH_IQR) 
  pH_max = pH_summary[5]+(3*pH_IQR) 
  
  temp_summary <- summary(temp_temp)
  temp_IQR <- IQR(temp_temp)
  temp_min <- temp_summary[2]-(3*temp_IQR) 
  temp_max = temp_summary[5]+(3*temp_IQR) 
  
  #Remove rows with alkalinity outliers
  data_site_qa<- data_site[which(data_site$Alkalinity > alk_min & data_site$Alkalinity < alk_max), ]
  data_site_qa<- data_site_qa[which(data_site_qa$pH > pH_min & data_site_qa$pH < pH_max), ]
  data_site_qa<- data_site_qa[which(data_site_qa$Temperature > temp_min & data_site_qa$Temperature < temp_max), ]
  
  all_data_site[[j]]<-data_site_qa
}


#Add calculated salinity to dataframe
for(j in final_sites_names){
  for(i in 1:nrow(all_data_site[[j]])){
    all_data_site[[j]]$salinity[i] <- ec2pss(all_data_site[[j]]$Conductivity[i]/1000,all_data_site[[j]]$Temperature[i],p=0)
  }
}

salinity_means <- data.frame()
count = 0
for(j in final_sites_names){
  count <- count+1
  salinity_means[count,1] <- j
  if(is.na(mean(na.omit(all_data_site[[j]]$salinity))) == TRUE){
    salinity_means[count,2] <- 0
  }else{
    salinity_means[count,2] <- mean(na.omit(all_data_site[[j]]$salinity))
  } 
}

# # SRP manual salinity check on 10/29/2023 showed USGS-11042000 with average salinity = 1.2; removing this site from dataset
# high_sal_index <- which(salinity_means$V2 > 1)
# name_remove <- salinity_means$V1[high_sal_index]
# all_data_site[[name_remove]] = NULL
########################################################################################
# Subset data to only include observations with salinity < 1
all_data_site_fresh <- list()
for(i in final_sites_names){
  temp <- subset(all_data_site[[i]], salinity < 1)
  all_data_site_fresh[[i]] <- temp
}

# Subset data to only include sites with >5 observations after January 1, 1950 
#all_data_site_recent <- list()
all_data_site_final <- list()
for(i in final_sites_names){
  temp <- subset(all_data_site_fresh[[i]], ActivityStartDate > as.Date("1950-01-01"))
  if(nrow(temp)>5){ #greater than 5 observations
    all_data_site_final[[i]] <- temp
  }
}
names_final <- names(all_data_site_final)
# Pull metadata for NEAA sites
siteNumbers_final <- sub("USGS-", "", names_final)
siteINFO_final <- readNWISsite(siteNumbers_final)

# Subset data to only include sites with >10 observations after January 1, 1990 (following Stets et al 2017)
all_data_site_recent <- list()
for(i in final_sites_names){
  temp <- subset(all_data_site_fresh[[i]], ActivityStartDate > as.Date("1990-01-01"))
  if(nrow(temp)>10){ #greater than 10 observationws
    all_data_site_recent[[i]] <- temp
  }
}
names_recent <- names(all_data_site_recent)

# Pull metadata for NEAA sites
siteNumbers_recent <- sub("USGS-", "", names_recent)
siteINFO_recent <- readNWISsite(siteNumbers_recent)

# Subset data to only include sites with >10 observations after January 1, 2000
all_data_site_2000 <- list()
for(i in final_sites_names){
  temp <- subset(all_data_site_fresh[[i]], ActivityStartDate > as.Date("2000-01-01"))
  if(nrow(temp)>10){ #greater than 10 observationws
    all_data_site_2000[[i]] <- temp
  }
}

names_2000 <- names(all_data_site_2000)

# Pull metadata for NEAA sites
siteNumbers_2000 <- sub("USGS-", "", names_2000)
siteINFO_2000 <- readNWISsite(siteNumbers_2000)

leaflet(data=siteINFO_final) %>% 
  addProviderTiles("CartoDB.Positron") %>%
  addCircleMarkers(~dec_long_va,~dec_lat_va,
                   fillColor = "blue", 
                   radius=4, stroke=FALSE,
                   fillOpacity = 0.8, opacity = 0.8,
                   popup=~station_nm) %>%
  addCircleMarkers(siteINFO_recent$dec_long_va,siteINFO_recent$dec_lat_va,
                   fillColor = "orange", 
                   radius=4, stroke=FALSE,
                   fillOpacity = 0.8, opacity = 0.8,
                   popup=siteINFO_recent$station_nm) %>%
  addCircleMarkers(siteINFO_2000$dec_long_va,siteINFO_2000$dec_lat_va,
                   fillColor = "red", 
                   radius=4, stroke=FALSE,
                   fillOpacity = 0.8, opacity = 0.8,
                   popup=siteINFO_2000$station_nm) 



