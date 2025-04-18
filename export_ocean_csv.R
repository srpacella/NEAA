# writes codap data to csv file for importing into matlab

#Write all_data_site list to individual .csv files for each station
for(i in final_sites_names){
  path <- paste("C:/Users/spacella/OneDrive - Environmental Protection Agency (EPA)/National Estuary Acidification Assessment/NEAA Analysis Files 240814/R Code May 2022 SRP Review/neaa_data_exports/",i,"_ocean.csv",sep="")
  write.csv(all_ocean_cc_calcs[[i]], file = path, row.names = TRUE)
}
