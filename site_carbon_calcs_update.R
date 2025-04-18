#funtion to calculate river carbonate system, ocean carbonate system, and estuarine buffer factors 
# across salinity spectrum given a WQP Monitoring Location Identifier


#Steps:
#1: Calculate summary of CO2 system conditions for USGS NWIS station (pH, ALK, T, S=0)
#2: Find closest CODAP ocean stations
#3: Calculate average Co2 system conditions using these stations (ALK, DIC, T, S)
#4: Interpolate ALK, DIC, T, S values between ocean and river end members
#5: Calculate full CO2 system and buffer factors for all interpolated points and ocean & river endmembers
#6: Produce summary plots

#Set working directory
#model_directory <- "/Users/spacella/Desktop/R Code May 2022 SRP Review"
#setwd(model_directory)
library(EnvStats)

site_carbon_calcs_update <- function(siteNo) {
  testsite <- siteNo
  testsite <- as.character(testsite)
  data_site_qa <- all_data_site[[testsite]]
  
  # Calculate full carbonate system with QAed data
  pH_site <- data_site_qa$ResultMeasureValue.x
  alk_site <- data_site_qa$ResultMeasureValue/1000000
  temp_site <- data_site_qa$ResultMeasureValue.y
  sal_site <- data_site_qa$salinity
  
  riv_cc <- carbfull(8, pH_site, alk_site, S=sal_site, T=temp_site, Patm=1, P=0, Pt=0, Sit=0,
                     k1k2="m10", kf="x", ks="d", pHscale="T", b="u74", gas="potential",
                     NH4t=0, HSt=0)
  riv_cc_buff<- buffesm(8, pH_site, alk_site, S=sal_site, T=temp_site, Patm = 1, P = 0, Pt = 0, Sit = 0,
                        k1k2 = "m10", kf = "x", ks = "d", pHscale = "T", b = "u74", warn = "y",
                        eos = "eos80")
  
  #2: Find closest CODAP ocean stations
  source("close_ocean.R")
  ocean_ids <- close_ocean(testsite,100*1000) #100km radius; find ocean stations IDs within this radius
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
  
  ocean_cc_buff<- buffesm(15, ocean_data$TALK/1000000, ocean_data$DIC/1000000, S=ocean_data$recommended_Salinity_PSS78, T=ocean_data$CTDTEMP_ITS90, Patm = 1, P = 0, Pt = 0, Sit = 0,
                          k1k2 = "l", kf = "x", ks = "d", pHscale = "T", b = "u74", warn = "y",
                          eos = "eos80")
  
  #4: Interpolate ALK, DIC, T, S values between ocean and river end members
  riv_alk <- mean(riv_cc$ALK)
  riv_dic <- mean(riv_cc$DIC)
  riv_temp <- mean(riv_cc$T)
  f_mar <- seq(from = 0, to = 1, by = 0.05)
  est_dic <- f_mar*ocean_dic_mean + (1-f_mar)*riv_dic
  est_alk <- f_mar*ocean_alk_mean + (1-f_mar)*riv_alk
  est_temp <- f_mar*ocean_temp_mean + (1-f_mar)*riv_temp
  est_sal <- f_mar*ocean_sal_mean
  est_cc <- carbfull(15, est_alk, est_dic, S=est_sal, T=est_temp, Patm=1, P=0, Pt=0, Sit=0,
                     k1k2="m10", kf="x", ks="d", pHscale="T", b="u74", gas="potential",
                     NH4t=0, HSt=0)
  est_cc_buff<- buffesm(15, est_alk, est_dic, S=est_sal, T=est_temp, Patm = 1, P = 0, Pt = 0, Sit = 0,
                        k1k2 = "m10", kf = "x", ks = "d", pHscale = "T", b = "u74", warn = "y",
                        eos = "eos80")
  
  
  list(riv_cc,riv_cc_buff,ocean_cc,ocean_cc_buff,est_cc,est_cc_buff,data_site_qa)
  
  
}


