download_calcium <- list(); # Download calcium data (filtered and unfiltered) for all rivers; 00910 units = mg/L as CaCO3; 00915 units = mg/L
download_calcium2 <- list(); # Download calcium data (filtered and unfiltered) for all rivers; 00910 units = mg/L as CaCO3; 00915 units = mg/L

for(j in names_final){ 
  download_calcium[[j]] <-  readWQPdata(siteNumber=j, 
                                          pCode = c("00910", #Calcium, water, unfiltered, milligrams per liter as calcium carbonate
                                                    "00915")) # Calcium, water, filtered, milligrams per liter
}
for(j in names_final){ 
  download_calcium2[[j]] <-  readWQPdata(siteNumber=j, 
                                        pCode = c("00915")) # Calcium, water, filtered, milligrams per liter
}

calcium_conversion = 1/1000*1/40.078

#Join calcium data with rest of dataset
calcium_join_test1 <- list()
for(j in names_final){ 
calcium_join_test1[[j]] <- inner_join(all_data_site_final[[j]],download_calcium2[[j]], by="ActivityStartDate")
}

#Calculate medians for each site and plot on map
calcium_medians <- data.frame()
count = 0
for(j in final_sites_names){
  count <- count+1
  calcium_medians[count,1] <- j
  if(is.na(mean(na.omit(calcium_join_test1[[j]]$ResultMeasureValue))) == TRUE){
    calcium_medians[count,2] <- 0
  }else{
    calcium_medians[count,2] <- mean(na.omit(calcium_join_test1[[j]]$ResultMeasureValue))
  } 
}


leaflet(data=siteINFO_final) %>% 
  addProviderTiles("CartoDB.Positron") %>%
  addCircleMarkers(~dec_long_va,~dec_lat_va,
                   fillColor = "blue", 
                   radius=calcium_medians$V2/2, stroke=FALSE,
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

