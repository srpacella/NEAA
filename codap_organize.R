# Orgqnize relevant ocean data for each river station

library(EnvStats)
#2: Find closest CODAP ocean stations
source("close_ocean.R")
all_ocean_cc_calcs <- list()
for(i in final_sites_names){
  testsite <- i
  testsite <- as.character(testsite)
  ocean_ids <- close_ocean(testsite,200*1000) #200km radius; find ocean stations IDs within this radius
  #3.1: Calculate average Co2 system conditions using these stations (ALK, DIC, T, S)
  ocean_data <- CODAP_surface[match(ocean_ids,CODAP_surface$NEAA_ID), ] #create table with only data from these stations
  ocean_data <- ocean_data[!ocean_data$recommended_Salinity_flag==9, ] #remove data with bad salinity
  ocean_data <- ocean_data[!ocean_data$CTDTEMP_flag==9, ] #remove data with bad temperature
  ocean_data <- ocean_data[!ocean_data$DIC_flag==9, ] #remove data with bad DIC
  ocean_data <- ocean_data[!ocean_data$TALK_flag==9, ] #remove data with bad ALK


  ocean_alk_stats <- summaryStats(ocean_data$TALK, quartiles = TRUE)
  ocean_dic_stats <- summaryStats(ocean_data$DIC, quartiles = TRUE)
  ocean_temp_stats <- summaryStats(ocean_data$CTDTEMP_ITS90, quartiles = TRUE)
  ocean_sal_stats <- summaryStats(ocean_data$recommended_Salinity_PSS78, quartiles = TRUE)
  ocean_alk_mean <- ocean_alk_stats[2]/1000000 #mean ocean ALK in mol/kg
  ocean_dic_mean <- ocean_dic_stats[2]/1000000 #mean ocean DIC in mol/kg
  ocean_temp_mean <- ocean_temp_stats[2] #mean ocean temperature
  ocean_sal_mean <- ocean_sal_stats[2] #mean ocean temperature

  ocean_cc <- carbfull(15, ocean_data$TALK/1000000, ocean_data$DIC/1000000, S=ocean_data$recommended_Salinity_PSS78, T=ocean_data$CTDTEMP_ITS90, Patm=1, P=0, Pt=0, Sit=0,
                     k1k2="l", kf="x", ks="d", pHscale="T", b="u74", gas="potential",
                     NH4t=0, HSt=0)
  ocean_cc$Station_ID <- ocean_data$Station_ID
  ocean_cc_all <- merge(ocean_cc,ocean_data,by="Station_ID")
  
  all_ocean_cc_calcs[[i]] <-ocean_cc_all
}

all_data_ocean <- do.call(rbind, all_ocean_cc_calcs)
ocean_coords <- all_data_ocean[,c(1,61,62,119)]
ocean_coords <- unique(ocean_coords)
ocean_coords_all <- CODAP_surface[,c(7,16,17,74)]
ocean_coords_all <- unique(ocean_coords_all)
