
#########################################################################################################
#Perform corrections to pH and Alkalinity data following Liu et al 2020 and re-calculate carbonate system
#########################################################################################################

### Alkalinity correction using DOC data #####

# Download DOC data available for each station
download_doc <- list();
for(j in names_final){ 
  download_doc[[j]] <-  readWQPdata(siteNumber=j, 
                                    pCode = c("00681"))
  download_doc[[j]]$umol <- download_doc[[j]]$ResultMeasureValue/1000/12.011*1000000 #convert from mg/L to umol/L
  download_doc[[j]]$org <- download_doc[[j]]$umol*0.8*0.1 #calculate upper constraint for organic alkalinity according to Liu et al 2020
}

# Download discharge data available for each station
download_discharge <- list();
for(j in names_final){ 
  download_discharge[[j]] <-  readWQPdata(siteNumber=j, 
                                    pCode = c("00061",
                                              "00060"))
}
for(j in names_final){ 
  download_discharge[[j]]<- download_discharge[[j]][which(download_discharge[[j]]$ResultMeasureValue > 0), ]
}

# Download conductivity data available for each station
download_conductivity <- list();
for(j in final_sites_names){ 
  download_conductivity[[j]] <-  readWQPdata(siteNumber=j, 
                                             characteristicName = "Specific conductance")
}

#Find which stations have at least 5 observations
download_doc <- download_doc[sapply(download_doc, nrow)>5]

# #Match up doc data with discharge
# doc_all <- do.call(rbind, download_doc_trim)
# discharge_all <- do.call(rbind, download_discharge)
# doc_discharge_match <- list()
# for(j in names_final){
#   download_doc[[j]] %>% drop_na(ResultMeasureValue)
#   download_discharge[[j]] %>% drop_na(ResultMeasureValue)
#   doc_discharge_match[[j]] <- left_join(download_doc[[j]], download_discharge[[j]], by="ActivityStartDate")
#   if(nrow(doc_discharge_match[[j]]) > 0){
#     plot(doc_discharge_match[[j]]$ResultMeasureValue.x,doc_discharge_match[[j]]$ResultMeasureValue.y)
#   }
# }



# #Find DOC data that is closest in time to observed pH,Alk,temp data, and append DOC data to full dataframe
# #Find index of DOC and Conductivity sampling date closest to full dataframe date.  Append this to full dataframe for each site
# library(birk)
# for(j in names_final){
#   for(i in 1:nrow(all_data_site_final[[j]])){
#     doc_index_closest <- which.closest(download_doc[[j]]$ActivityStartDate,all_data_site_final[[j]]$ActivityStartDate[i])
#     #conductivity_index_closest <- which.closest(download_conductivity[[j]]$ActivityStartDate,all_data_site_final[[j]]$ActivityStartDate[i])
#     all_data_site_final[[j]]$doc_umol[i] <- download_doc[[j]]$umol[doc_index_closest]
#     all_data_site_final[[j]]$org[i] <- download_doc[[j]]$org[doc_index_closest]
#     #all_data_site_final[[j]]$conductivity[i] <- download_conductivity[[j]]$ResultMeasureValue[conductivity_index_closest]
#   }
# }
#Plot organic alkalinity data

# Bind all org_alk data to dataframe
org_all <- do.call(rbind, download_doc)
conductivity_all <- do.call(rbind, download_conductivity)
discharge_all <- do.call(rbind, download_discharge)
summary_org <- summary(org_all$org)
org_all$ID <- paste(org_all$ActivityStartDate, org_all$MonitoringLocationIdentifier)
conductivity_all$ID <- paste(conductivity_all$ActivityStartDate, conductivity_all$MonitoringLocationIdentifier)
discharge_all$ID <- paste(discharge_all$ActivityStartDate, discharge_all$MonitoringLocationIdentifier)
join_org_discharge <- left_join(org_all, discharge_all, by="ID")
join_org_discharge <- join_org_discharge[which(is.na(join_org_discharge$ResultMeasureValue.y)==FALSE), ]

plot(join_org_discharge$ResultMeasureValue.y,join_org_discharge$ResultMeasureValue.x)

all_org_discharge <- list()
for(j in names_final){
  all_org_discharge[[j]] <- join_org_discharge %>% 
    filter(join_org_discharge$MonitoringLocationIdentifier.x == j)
  #plot(all_org_discharge[[j]]$ResultMeasureValue.y,all_org_discharge[[j]]$ResultMeasureValue.x)
}


#Calculate means for each site and compare with total alkalinity
orgalk_means <- list()
for(j in names_final){
  orgalk_means[[j]]$org <- mean(download_doc[[j]]$org)
  orgalk_means[[j]]$alk <- mean(all_data_site_final[[j]]$Alkalinity)
  #orgalk_means[j,1] <- download_doc[[j]]$MonitoringLocationIdentifier[1]
}
alkalinity_all <- do.call(rbind, orgalk_means)
alkalinity_all <- as.data.frame(alkalinity_all)
alkalinity_all$site <- sub("USGS-", "", names_final)
alkalinity_all$org <- as.numeric(alkalinity_all$org)
alkalinity_all$alk <- as.numeric(alkalinity_all$alk)
alkalinity_all$perc_org <- alkalinity_all$org/alkalinity_all$alk*100
plot(alkalinity_all$alk,alkalinity_all$org)
plot(alkalinity_all$alk,alkalinity_all$perc_org)

discharge_means <- data.frame()
count = 0
for(j in names_final){
  count <- count+1
  discharge_means[count,1] <- j
  if(is.na(mean(na.omit(download_discharge[[j]]$ResultMeasureValue))) == TRUE){
    discharge_means[count,2] <- 0
  }else{
  discharge_means[count,2] <- mean(na.omit(download_discharge[[j]]$ResultMeasureValue))
  } 
}

discharge_medians <- data.frame()
count = 0
for(j in names_final){
  count <- count+1
  discharge_medians[count,1] <- j
  if(is.na(mean(na.omit(download_discharge[[j]]$ResultMeasureValue))) == TRUE){
    discharge_medians[count,2] <- 0
  }else{
    discharge_medians[count,2] <- median(na.omit(download_discharge[[j]]$ResultMeasureValue))
  } 
}

conductivity_means <- data.frame()
count = 0
for(j in names_final){
  count <- count+1
  conductivity_means[count,1] <- j
  if(is.na(mean(na.omit(download_conductivity[[j]]$ResultMeasureValue))) == TRUE){
    conductivity_means[count,2] <- 0
  }else{
    conductivity_means[count,2] <- mean(na.omit(download_conductivity[[j]]$ResultMeasureValue))
  } 
}

# Correct pH and alkalinity data following Liu et al 2020 recommendation for quiescent waters
for(j in names_final){
  for(i in 1:nrow(all_data_site_final[[j]])){
    if("Conductivity" %in% colnames(all_data_site_final[[j]])){
      all_data_site_final[[j]]$ion[i] <- all_data_site_final[[j]]$Conductivity*1.3E-5 #calculate ionic strength from conductivity
      all_data_site_final[[j]]$pH_error[i] <- 0.03 + 0.05*log10(all_data_site_final[[j]]$ion[i])
      all_data_site_final[[j]]$pH_correct[i] <- all_data_site_final[[j]]$pH[i]-all_data_site_final[[j]]$pH_error[i]
    }
    if(is.na(alkalinity_all[j,1]) == "FALSE"){
      #all_data_site_final_qa[[j]]$carbalk[i] <- all_data_site_final[[j]]$ResultMeasureValue[i] - all_data_site_final[[j]]$org[i] #uses time matched values
      all_data_site_final[[j]]$carbalk[i] <- all_data_site_final[[j]]$Alkalinity[i] - orgalk_means[[j]]$org #uses mean orgalk for each site
    } else {
      all_data_site_final[[j]]$carbalk[i] <- "NaN"
    }
  }
}



# # Plot corrected alkalinity data
# for(i in names_final){
#   if("carbalk" %in% colnames(all_data_site_final[[i]])){
#     plottest <- ggplot(data = all_data_site_final[[i]], aes(x=ActivityStartDate, y=carbalk)) +
#       geom_point() +
#       ggtitle(i)
#     print(plottest)
#   }
# }    
# 
# for(i in 1:nrow(alkalinity_all)){
#   if(nrow(all_org_discharge[[i]]>1)){
#     plot(all_org_discharge[[i]]$ResultMeasureValue.y,all_org_discharge[[i]]$org,main=names(all_org_discharge[i]))
#   }
# }
# col_types <- c("darkblue","dodgerblue","green4","gold1","orange","brown","red")
# leg_vals <- unique(as.numeric(quantile(org_all$org, 
#                                        probs=c(0,0.01,0.1,0.25,0.5,0.75,0.9,.99,1), na.rm=TRUE)))
# pal = colorBin(col_types, org_all$org, bins = leg_vals)
# leaflet(data=all_data_site_final) %>% 
#   addProviderTiles("CartoDB.Positron") %>%
#   addCircleMarkers(all_conductivity_mean$V4,all_conductivity_mean$V3,
#                    color = pal(all_conductivity_mean$salinity), radius=4, stroke=FALSE,
#                    fillOpacity = 0.8, opacity = 0.8,
#                    popup=paste(all_conductivity_mean$V1)) %>%
#   addLegend(position = 'bottomleft',
#             pal=pal,
#             values=all_conductivity_mean$salinity,
#             opacity = 0.8,
#             labFormat = labelFormat(digits = 1), 
#             title = 'Salinity')