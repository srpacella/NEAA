# writes neaa_final_download.R output to csv file for importing into matlab


setwd("C:/Users/spacella/OneDrive - Environmental Protection Agency (EPA)/National Estuary Acidification Assessment/NEAA Analysis Files 240814/R Code May 2022 SRP Review")

#load("~/Desktop/R Code May 2022 SRP Review/neaa_final_last.RData")
#load("C:/Users/spacella/OneDrive - Environmental Protection Agency (EPA)/R Code May 2022 SRP Review/neaa_download_220803.RData") #August 3 2022 update

#library(R.matlab)
# filename <- paste(tempfile(), ".mat", sep = "")
# 
# writeMat("~/Desktop/R Code May 2022 SRP Review/all_data_site.mat", all_date_site = all_data_site)
# 
# writeMat("~/Desktop/R Code May 2022 SRP Review/data_site.mat", data_site = data_site )

#Write all_data_site list to individual .csv files for each station
rm(i)
for(i in final_sites_names){
  path <- paste("C:/Users/spacella/OneDrive - Environmental Protection Agency (EPA)/National Estuary Acidification Assessment/NEAA Analysis Files 240814/R Code May 2022 SRP Review/neaa_data_exports/",i,".csv",sep="")
  write.csv(all_data_site_final[[i]], file = path, row.names = TRUE)
}

#Write calcium data list to individual .csv files for each station
rm(i)
for(i in final_sites_names){
  path <- paste("C:/Users/spacella/OneDrive - Environmental Protection Agency (EPA)/National Estuary Acidification Assessment/NEAA Analysis Files 240814/R Code May 2022 SRP Review/neaa_data_exports/",i,"_calcium.csv",sep="")
  write.csv(download_calcium2[[i]], file = path, row.names = TRUE)
}

#Write list of final stations
write.csv(names_final, file = paste("C:/Users/spacella/OneDrive - Environmental Protection Agency (EPA)/National Estuary Acidification Assessment/NEAA Analysis Files 240814/R Code May 2022 SRP Review/neaa_data_exports/","final_sites_names",".csv",sep=""), row.names = FALSE)

write.csv(siteINFO_final, file = paste("C:/Users/spacella/OneDrive - Environmental Protection Agency (EPA)/National Estuary Acidification Assessment/NEAA Analysis Files 240814/R Code May 2022 SRP Review/neaa_data_exports/","siteINFO",".csv",sep=""), row.names = FALSE)

write.csv(alkalinity_all, file = paste("C:/Users/spacella/OneDrive - Environmental Protection Agency (EPA)/National Estuary Acidification Assessment/NEAA Analysis Files 240814/R Code May 2022 SRP Review/neaa_data_exports/","alkalinity_all",".csv",sep=""), row.names = FALSE)

write.csv(trends_final, file = paste("C:/Users/spacella/OneDrive - Environmental Protection Agency (EPA)/National Estuary Acidification Assessment/NEAA Analysis Files 240814/R Code May 2022 SRP Review/neaa_data_exports/","ts_regressions",".csv",sep=""), row.names = FALSE)

write.csv(discharge_means, file = paste("C:/Users/spacella/OneDrive - Environmental Protection Agency (EPA)/National Estuary Acidification Assessment/NEAA Analysis Files 240814/R Code May 2022 SRP Review/neaa_data_exports/","discharge_means",".csv",sep=""), row.names = FALSE)
