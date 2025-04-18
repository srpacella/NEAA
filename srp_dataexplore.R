# SRP data exlploration to work through how to pair pH and Alk observations
# Started 5/26/2022

rm(list=ls())

library(svDialogs); library(tidyverse); library(lubridate); library(xts)
library(dplyr); library(purrr); library(EGRET); library(dataRetrieval);
library(leaflet);library(ggplot2); library(maps); library(mapdata)

setwd("/Users/spacella/Desktop/R Code May 2022 SRP Review")


#Read in all alkalinty data written by Script 4
alk_all_txt <- read.delim("/Users/spacella/Desktop/R Code May 2022 SRP Review/NEAA_full_alk_data.txt",header = TRUE, sep = "\t", dec = ".")
alk_sites <- unique(alk_all_txt$.id)


#Read in all pH data written by Script 5
pH_all_txt <- read.delim("/Users/spacella/Desktop/R Code May 2022 SRP Review/NEAA_full_pH_data.txt",header = TRUE, sep = "\t", dec = ".")
pH_sites <- unique(pH_all_txt$.id)


common_sites <- intersect(alk_sites,pH_sites)
result_alk <- filter(alk_all_txt,alk_all_txt$.id == common_sites[3],)
result_pH <- filter(pH_all_txt,pH_all_txt$.id == common_sites[3],)

#Plot pH and Alk on same figure
ggplot() + 
  geom_point(data=result_alk, aes(x=sample_dt, y=result_va), color='green') + 
  geom_point(data=result_pH, aes(x=sample_dt, y=result_va), color='red')


p<-ggplot(result_alk, aes(x=sample_dt,y= result_va))+
  geom_point(color="firebrick") + 
  ylab("Alkalinity (ueq/kg)") +
  xlab("Year") +
  labs(title = paste("Selected USGS station")) +
  theme_minimal() +
  facet_wrap(~site_no)
print(p)
