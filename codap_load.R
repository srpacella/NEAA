# Read in CODAP data and remove flagged datapoints
setwd("C:/Users/spacella/OneDrive - Environmental Protection Agency (EPA)/National Estuary Acidification Assessment/NEAA Analysis Files 240814/R Code May 2022 SRP Review")

CODAP_NA_v2020 <- read_excel("CODAP_NA_v2020.xlsx", 
                             col_types = c("text", "text", "text", 
                                           "text", "text", "text", "text", "text", 
                                           "text", "text", "text", "text", "text", 
                                           "text", "date", "numeric", "numeric", "numeric", 
                                           "numeric", "numeric", "numeric", 
                                           "numeric", "text", "numeric", "text", 
                                           "numeric", "text", "numeric", "text", 
                                           "numeric", "text", "numeric", "text", 
                                           "numeric", "text", "numeric", "text", 
                                           "numeric", "text", "numeric", "text", 
                                           "numeric", "numeric", "text", "numeric", 
                                           "numeric", "numeric", "numeric", "text", 
                                           "numeric", "numeric", "numeric", "numeric", "text", 
                                           "numeric", "numeric", "numeric", "numeric", "numeric", 
                                           "numeric", "text", "numeric", "text", "numeric", 
                                           "text", "numeric", "text", "numeric", "text", 
                                           "numeric", "text", "numeric", "text"))

CODAP_work <- CODAP_NA_v2020[-c(1), ] #Working dataset with removed first row with unit metadata (mostly NA's)

CODAP_surface <- CODAP_work[!(CODAP_work$Depth>20), ] #pull data shallower than 20m
#Create unique cruise_station ID
CODAP_surface$NEAA_ID <- paste(CODAP_surface$Cruise_ID,CODAP_surface$Station_ID)

CODAP_surface <- CODAP_surface[!CODAP_surface$recommended_Salinity_flag==9, ] #remove data with missing salinity
CODAP_surface <- CODAP_surface[!CODAP_surface$CTDTEMP_flag==9, ] #remove data with missing temperature
CODAP_surface <- CODAP_surface[!CODAP_surface$DIC_flag==9, ] #remove data with missing DIC
CODAP_surface <- CODAP_surface[!CODAP_surface$TALK_flag==9, ] #remove data with missing ALK

CODAP_surface <- CODAP_surface[!CODAP_surface$recommended_Salinity_flag==3, ] #remove data with questionable salinity
CODAP_surface <- CODAP_surface[!CODAP_surface$CTDTEMP_flag==3, ] #remove data with questionable temperature
CODAP_surface <- CODAP_surface[!CODAP_surface$DIC_flag==3, ] #remove data with questionable DIC
CODAP_surface <- CODAP_surface[!CODAP_surface$TALK_flag==3, ] #remove data with questionable ALK

