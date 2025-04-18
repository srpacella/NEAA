#funtion to calcualte buffer factors across estuarine salinity spectrum given a USGS NWIS station value


#Steps:
#1: Calculate summary of CO2 system conditions for USGS NWIS station (pH, ALK, T, S=0)
#2: Find closest CODAP ocean stations
#3: Calculate average Co2 system conditions using these stations (ALK, DIC, T, S)
#4: Interpolate ALK, DIC, T, S values between ocean and river end members
#5: Calculate full CO2 system and buffer factors for all interpolated points and ocean & river endmembers
#6: Produce summary plots

#Set working directory
model_directory <- "/Users/spacella/Desktop/R Code May 2022 SRP Review"
setwd(model_directory)
library(seacarb)
library(EnvStats)
source("close_ocean.R")
estuary_buffers_sm <- function(siteNo,alk,pH,temp) {
  #1. Calculate median and CO2 system conditions for USGS NWIS station (pH, ALK, T, S=0)
  #start.date <- "1900-01-01"
  #end.date <- "2020-01-01"
  #pCode_alk <- c("39086","00410","90410","29801","39036","00419","00417","29802","00440","95410") #alkalinity parameter codes
  #pCode_pH <- c("00400","00403","00408") #pH codes, 00434 = at 25C
  #pCode_temp <- c("00010") #temperature in Celsius
  # riv_alk <- readNWISqw(siteNumbers = siteNo, parameterCd = pCode_alk,
  #                       startDate = start.date, endDate = end.date) # download all ALK data from site
  # riv_pH <- readNWISqw(siteNumbers = siteNo, parameterCd = pCode_pH,
  #                      startDate = start.date, endDate = end.date) # download all pH data from site
  # riv_temp <- readNWISqw(siteNumbers = siteNo, parameterCd = pCode_temp,
  #                        startDate = start.date, endDate = end.date) # download all temperature data from site
  riv_alk_stats <- summaryStats(alk, quartiles = TRUE)
  riv_pH_stats <- summaryStats(pH, quartiles = TRUE)
  riv_temp_stats <- summaryStats(temp, quartiles = TRUE)
  riv_alk_mean <- riv_alk_stats[2]*0.02*1000/1000000 #mean river ALK in mol/kg
  riv_pH_mean <- riv_pH_stats[2] #mean river pH
  riv_temp_mean <- riv_temp_stats[2] #mean river temperature
  riv_cc <- carbfull(8, riv_pH_mean, riv_alk_mean, S=0, T=riv_temp_mean, Patm=1, P=0, Pt=0, Sit=0,
                     k1k2="m10", kf="x", ks="d", pHscale="T", b="u74", gas="potential",
                     NH4t=0, HSt=0)
  riv_cc_buff<- buffer(8, riv_pH_mean, riv_alk_mean, S=0, T = riv_temp_mean, Patm = 1, P = 0, Pt = 0, Sit = 0,
                       k1k2 = "m10", kf = "x", ks = "d", pHscale = "T", b = "u74", warn = "y",
                       eos = "eos80")
  
  #2: Find closest CODAP ocean stations
  ocean_ids <- close_ocean(siteNo,250*1000) #50km radius; find ocean stations IDs within this radius
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
  ocean_cc <- carbfull(15, ocean_alk_mean, ocean_dic_mean, S=ocean_sal_mean, T=ocean_temp_mean, Patm=1, P=0, Pt=0, Sit=0,
                       k1k2="l", kf="x", ks="d", pHscale="T", b="u74", gas="potential",
                       NH4t=0, HSt=0)
  
  #3.2 Calculate 25% and 75% quartiles for ocean and river, and interpolate through estuary
  riv_alk_25 <- riv_alk_stats[7]*0.02*1000/1000000 #mean river ALK in mol/kg
  riv_pH_25 <- riv_pH_stats[7] #mean river pH
  riv_temp_25 <- riv_temp_stats[7] #mean river temperature
  riv_alk_75 <- riv_alk_stats[8]*0.02*1000/1000000 #mean river ALK in mol/kg
  riv_pH_75 <- riv_pH_stats[8] #mean river pH
  riv_temp_75 <- riv_temp_stats[8] #mean river temperature
  riv_cc_25 <- carbfull(8, riv_pH_25, riv_alk_25, S=0, T=riv_temp_25, Patm=1, P=0, Pt=0, Sit=0,
                        k1k2="m10", kf="x", ks="d", pHscale="T", b="u74", gas="potential",
                        NH4t=0, HSt=0)
  riv_cc_buff_25<- buffesm(8, riv_pH_25, riv_alk_25, S=0, T=riv_temp_25, Patm = 1, P = 0, Pt = 0, Sit = 0,
                          k1k2 = "m10", kf = "x", ks = "d", pHscale = "T", b = "u74", warn = "y",
                          eos = "eos80")
  riv_cc_75 <- carbfull(8, riv_pH_75, riv_alk_75, S=0, T=riv_temp_75, Patm=1, P=0, Pt=0, Sit=0,
                        k1k2="m10", kf="x", ks="d", pHscale="T", b="u74", gas="potential",
                        NH4t=0, HSt=0)
  riv_cc_buff_75<- buffesm(8, riv_pH_75, riv_alk_75, S=0, T=riv_temp_75, Patm = 1, P = 0, Pt = 0, Sit = 0,
                          k1k2 = "m10", kf = "x", ks = "d", pHscale = "T", b = "u74", warn = "y",
                          eos = "eos80")
  
  ocean_alk_25 <- ocean_alk_stats[7]/1000000 #mean ocean ALK in mol/kg
  ocean_dic_25 <- ocean_dic_stats[7]/1000000 #mean ocean DIC in mol/kg
  ocean_temp_25 <- ocean_temp_stats[7] #mean ocean temperature
  ocean_sal_25 <- ocean_sal_stats[7] #mean ocean temperature
  ocean_cc_25 <- carbfull(15, ocean_alk_25, ocean_dic_25, S=ocean_sal_25, T=ocean_temp_25, Patm=1, P=0, Pt=0, Sit=0,
                          k1k2="l", kf="x", ks="d", pHscale="T", b="u74", gas="potential",
                          NH4t=0, HSt=0)
  ocean_alk_75 <- ocean_alk_stats[8]/1000000 #mean ocean ALK in mol/kg
  ocean_dic_75 <- ocean_dic_stats[8]/1000000 #mean ocean DIC in mol/kg
  ocean_temp_75 <- ocean_temp_stats[8] #mean ocean temperature
  ocean_sal_75 <- ocean_sal_stats[8] #mean ocean temperature
  ocean_cc_75 <- carbfull(15, ocean_alk_75, ocean_dic_75, S=ocean_sal_75, T=ocean_temp_75, Patm=1, P=0, Pt=0, Sit=0,
                          k1k2="l", kf="x", ks="d", pHscale="T", b="u74", gas="potential",
                          NH4t=0, HSt=0)
  
  #4: Interpolate ALK, DIC, T, S values between ocean and river end members
  f_mar <- seq(from = 0, to = 1, by = 0.05)
  est_dic <- f_mar*ocean_cc$DIC + (1-f_mar)*riv_cc$DIC
  est_alk <- f_mar*ocean_cc$ALK + (1-f_mar)*riv_cc$ALK
  est_temp <- f_mar*ocean_temp_mean + (1-f_mar)*riv_temp_mean
  est_sal <- f_mar*ocean_sal_mean
  est_cc <- carbfull(15, est_alk, est_dic, S=est_sal, T=est_temp, Patm=1, P=0, Pt=0, Sit=0,
                     k1k2="m10", kf="x", ks="d", pHscale="T", b="u74", gas="potential",
                     NH4t=0, HSt=0)
  est_cc_buff<- buffesm(15, est_alk, est_dic, S=est_sal, T=est_temp, Patm = 1, P = 0, Pt = 0, Sit = 0,
                       k1k2 = "m10", kf = "x", ks = "d", pHscale = "T", b = "u74", warn = "y",
                       eos = "eos80")
  #4.1: Interpolate ALK, DIC, T, S values between ocean and river end members
  est_dic_25 <- f_mar*ocean_cc_25$DIC + (1-f_mar)*riv_cc_25$DIC
  est_alk_25 <- f_mar*ocean_cc_25$ALK + (1-f_mar)*riv_cc_25$ALK
  est_temp_25 <- f_mar*ocean_temp_25 + (1-f_mar)*riv_temp_25
  est_sal_25 <- f_mar*ocean_sal_25
  est_cc_25 <- carbfull(15, est_alk_25, est_dic_25, S=est_sal_25, T=est_temp_25, Patm=1, P=0, Pt=0, Sit=0,
                        k1k2="m10", kf="x", ks="d", pHscale="T", b="u74", gas="potential",
                        NH4t=0, HSt=0)
  est_cc_buff_25<- buffesm(15, est_alk_25, est_dic_25, S=est_sal_25, T=est_temp_25, Patm = 1, P = 0, Pt = 0, Sit = 0,
                          k1k2 = "m10", kf = "x", ks = "d", pHscale = "T", b = "u74", warn = "y",
                          eos = "eos80")
  est_dic_75 <- f_mar*ocean_cc_75$DIC + (1-f_mar)*riv_cc_75$DIC
  est_alk_75 <- f_mar*ocean_cc_75$ALK + (1-f_mar)*riv_cc_75$ALK
  est_temp_75 <- f_mar*ocean_temp_75 + (1-f_mar)*riv_temp_75
  est_sal_75 <- f_mar*ocean_sal_75
  est_cc_75 <- carbfull(15, est_alk_75, est_dic_75, S=est_sal_75, T=est_temp_75, Patm=1, P=0, Pt=0, Sit=0,
                        k1k2="m10", kf="x", ks="d", pHscale="T", b="u74", gas="potential",
                        NH4t=0, HSt=0)
  est_cc_buff_75<- buffesm(15, est_alk_75, est_dic_75, S=est_sal_75, T=est_temp_75, Patm = 1, P = 0, Pt = 0, Sit = 0,
                          k1k2 = "m10", kf = "x", ks = "d", pHscale = "T", b = "u74", warn = "y",
                          eos = "eos80")
  
  revelle_plot <- ggplot(est_cc_buff, aes(x=est_cc$S,y=R)) +
    geom_line(color='black') +
    geom_line(aes(x=est_cc_25$S,y=est_cc_buff_25$R),color='blue') +
    geom_line(aes(x=est_cc_75$S,y=est_cc_buff_75$R),color='red') +
    theme_bw() +
    labs(title = paste("USGS Site#",siteNo),
         subtitle = paste(alk_sites$station_nm[which(alk_sites$site_no==siteNo)]),
         x = "Salinity",
         y = "Revelle Factor")
  list(est_cc,est_cc_buff)
  
  
}
