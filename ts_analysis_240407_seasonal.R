#### TS Slopes and estimated marginal means ######
# SRP Updated April 2024 to compare with results of rkt seasonal kendall

library(mblm);library(poorman);library(modelbased); library(leaflegend); library(rkt)

sen <- function(..., weights = NULL) {
  mblm::mblm(...)}

################################################################################
########################### Thiel-Sen Slope estimation #########################
# Calculate TS slopes for each individual station, output as text file
# and create a PDF of each site's pH and alk slope with scatterplot
# of data.

#Define origin date for DOY calculations
o <- as.Date("0000-01-01") 

# Create empty vectors
pH_ts_slope <- vector(mode="numeric")
pH_residuals <- vector(mode="numeric")
pH_ts_p_value <- vector(mode="numeric")
pH_lm_p_value <- vector(mode="numeric")
pH_t_value <- vector(mode="numeric")
pH_lm_slope <- vector(mode="numeric")
alk_ts_slope <- vector(mode="numeric")
alk_residuals <- vector(mode="numeric")
alk_ts_intercept_value <- vector(mode="numeric")
pH_ts_intercept_value <- vector(mode="numeric")
alk_p_value <- vector(mode="numeric")
alk_ts_p_value <- vector(mode="numeric")
alk_t_value <- vector(mode="numeric")
alk_lm_slope <- vector(mode="numeric")
site_number <- vector(mode="character")
max_date <- vector(mode="numeric")
min_date <- vector(mode="numeric")
temp_site <- vector(mode="numeric")
n <- vector(mode="numeric")
start_date <- vector(mode = "numeric")
end_date <- vector(mode = "numeric")

# Idenfity unique sites
uniq <- names_final

# Open a PDF and write to it
pdf("Station_ts_curves.pdf")
for (i in uniq){
  temp_site <- all_data_site_final[[i]]
  
  # Identify site, date, pH, and alk columns
  site_number_1 = i
  DOY <- as.numeric(difftime(as.Date(temp_site$ActivityStartDate), o, unit = "day"))
  pH <- 10^(-temp_site$pH_correct)
  alk <- temp_site$Alkalinity
  
  start_date_1 <- temp_site$ActivityStartDate[1]
  end_date_1 <- temp_site$ActivityStartDate[length(temp_site$ActivityStartDate)]
  
  # Do the regressions
  # Calculate pH fits
  lm_pH_1 <- lm(pH ~ DOY)
  ts_pH_1 <- mblm(pH ~ DOY)
  
  # Calculate alk fits
  lm_alk_1 <- lm(alk ~ DOY)
  ts_alk_1 <- mblm(alk ~ DOY)
  
  # # Do the regressions
  # # Calculate pH fits
  # lm_pH_1 <- lm(DOY ~ pH)
  # ts_pH_1 <- mblm(DOY ~ pH)
  # 
  # # Calculate alk fits
  # lm_alk_1 <- lm(DOY ~ alk)
  # ts_alk_1 <- mblm(DOY ~ alk)
  
  # Calculate coefficients for pH
  matrix_coef_ts <- summary(ts_pH_1)$coefficients
  matrix_coef_lm <- summary(lm_pH_1)$coefficients
  
  # Parse out coefficients
  pH_ts_slope_1 <- matrix_coef_ts[2,1] # theil-sen slope
  pH_residuals_1 <- median(ts_pH_1$residuals) # residuals from theil-sen
  pH_ts_p_value_1 <- matrix_coef_ts[2,4] # p-value should be less than 0.05
  pH_ts_intercept_value_1 <- matrix_coef_ts[1,1] # ts intercept
  pH_lm_p_value_1 <- matrix_coef_lm[2,4] # lm p-value
  pH_t_value_1 <- matrix_coef_lm[2,3] # t-value should be greater than 0
  pH_lm_slope_1 <- matrix_coef_lm[2,1] # slope from linear regression
  
  # Calculate coefficients for alk
  matrix_coef_ts <- summary(ts_alk_1)$coefficients
  matrix_coef_lm <- summary(lm_alk_1)$coefficients
  
  alk_ts_slope_1 <- matrix_coef_ts[2,1]
  alk_residuals_1 <- median(ts_alk_1$residuals)
  alk_p_value_1 <- matrix_coef_lm[2,4] # p-value should be less than 0.05
  alk_ts_intercept_value_1 <- matrix_coef_ts[1,1] # ts intercept
  alk_ts_p_value_1 <- matrix_coef_ts[2,4]
  alk_t_value_1 <- matrix_coef_lm[2,3] # t-value should be greater than 0
  alk_lm_slope_1 <- matrix_coef_lm[2,1]
  
  # Dates
  min_date_1 <- min(DOY)
  max_date_1 <- max(DOY)
  
  # Update vectors
  pH_ts_slope <- append(pH_ts_slope, pH_ts_slope_1)
  pH_residuals <- append(pH_residuals, pH_residuals_1)
  pH_ts_p_value <- append(pH_ts_p_value, pH_ts_p_value_1)
  pH_lm_p_value <- append(pH_lm_p_value, pH_lm_p_value_1)
  pH_t_value <- append(pH_t_value, pH_t_value_1)
  pH_lm_slope <- append(pH_lm_slope, pH_lm_slope_1)
  pH_ts_intercept_value <- append(pH_ts_intercept_value, pH_ts_intercept_value_1)
  
  # alk
  alk_ts_slope <- append(alk_ts_slope, alk_ts_slope_1)
  alk_residuals <- append(alk_residuals, alk_residuals_1)
  alk_p_value <- append(alk_p_value, alk_p_value_1)
  alk_ts_p_value <- append(alk_ts_p_value, alk_ts_p_value_1)
  alk_t_value <- append(alk_t_value, alk_t_value_1)
  alk_lm_slope <- append(alk_lm_slope, alk_lm_slope_1)
  alk_ts_intercept_value <- append(alk_ts_intercept_value, alk_ts_intercept_value_1)
  
  # dates
  max_date <- append(max_date, max_date_1)
  min_date <- append(min_date, min_date_1)
  site_number <- append(site_number, site_number_1)
  n_1 <- nrow(temp_site)
  n <- append(n, n_1)
  start_date <- append(start_date, start_date_1)
  end_date <- append(end_date, end_date_1)
  
  # Add plots
  #pdf(paste(unique(temp_site$geologicsite),"_pH",  ".pdf", sep=""), width=4, height=4)
  abc <- ggplot(temp_site, aes(ActivityStartDate, 10^-pH)) + 
    geom_point(alpha=.3, size=2) + ylab("[H+]")+ xlab("Year") +
    geom_smooth(method = loess, color="red") +
    geom_smooth(method = sen, color="blue") +
    ggtitle(label=paste(unique(temp_site$MonitoringLocation)," [H+] trend", sep=""),
            subtitle=paste(unique(temp_site$OrganizationIdentifier)))
  print(abc)
  #dev.off()
  
  #pdf(paste(unique(temp_site$geologicsite),"_alkalinity",  ".pdf", sep=""), width=4, height=4)
  cbd <- ggplot(temp_site, aes(ActivityStartDate, alk)) + 
    geom_point(alpha=0.3, size=2) + ylab("Alkalininty ueq/kg") + xlab("Year") +
    geom_smooth(method = loess, color="red") +
    geom_smooth(method = sen, color="blue")+
    ggtitle(paste(unique(temp_site$MonitoringLocation)," Alkalinity trend", sep=""),
            subtitle=paste(unique(temp_site$OrganizationIdentifier)))
  print(cbd)
  #dev.off()
}
dev.off()

# Put all of the TS slope infomation into a dataframe
DF1 <- data.frame(site_number,n, min_date, max_date, pH_ts_slope, 
                  pH_lm_slope, pH_ts_intercept_value, pH_residuals, pH_ts_p_value, pH_lm_p_value,
                  pH_t_value, alk_ts_slope, alk_ts_intercept_value, alk_lm_slope, alk_residuals, 
                  alk_p_value, alk_ts_p_value, alk_t_value, start_date, end_date)

DF1$duration <- DF1$max_date-DF1$min_date

trends_final <- DF1 %>%
  filter(duration > 3650) #10 year minimum
trends_final <- trends_final %>%
  filter(n > 24) #at least 25 observations
trends_final$duration_years <- trends_final$duration/365

#Indicate 

trends_final$siteNumbers <- sub("USGS-", "", trends_final$site_number)
count <- 0
for(i in trends_final$siteNumbers){
  count = count + 1
  trends_final$lat[count] <- siteINFO_final$dec_lat_va[which(siteINFO_final$site_no==i)]
  trends_final$lon[count] <- siteINFO_final$dec_long_va[which(siteINFO_final$site_no==i)]
  site_no_temp <- trends_final$site_number[count]
  trends_final$discharge[count] <- mean(download_discharge[[site_no_temp]]$ResultMeasureValue)
}

###### Produce scatterplot of pH and alkalinity trends
plot(trends_final$alk_ts_slope*365,trends_final$pH_ts_slope*365*1e9,xlab="Alkalinity trend (umol/kg/year)",ylab="H+ trend(nmol/kg/year)")




# Plot H+ trends with size by discharge
trends_sig <- trends_final %>%
  filter(pH_ts_p_value < 0.05)
trends_nonsig <- trends_final %>%
  filter(pH_ts_p_value > 0.05)
trends_sig_pos <- trends_sig %>%
  filter(pH_ts_slope > 0)
trends_sig_neg <- trends_sig %>%
  filter(pH_ts_slope < 0)

leaflet(data=trends_sig) %>%
  addProviderTiles("CartoDB.Positron") %>%
  addCircleMarkers(trends_sig_neg$lon,trends_sig_neg$lat,
                   color = "blue", radius=log10(abs(trends_sig_neg$pH_ts_slope*365*1e9))*7, stroke=FALSE,
                   fillOpacity = 0.8, opacity = 0.8,
                   popup=paste(trends_sig_neg$siteNumbers)) %>%
  addCircleMarkers(trends_sig_pos$lon,trends_sig_pos$lat,
                   color = "red", radius=log10(abs(trends_sig_pos$pH_ts_slope*365*1e9))*7, stroke=FALSE,
                   fillOpacity = 0.8, opacity = 0.8,
                   popup=paste(trends_sig_pos$siteNumbers)) %>%
  # addLegend(position = 'bottomleft',
  #           pal=pal,
  #           values=trends_sig_pos$pH_ts_slope*365*1e9,
  #           opacity = 0.8,
  #           labFormat = labelFormat(digits = 1),
  #           title = 'H+ Trend (nmol/year)')
  addLegendSize(
    values = log10(trends_sig_pos$pH_ts_slope*365*1e9)*7,
    color = 'red',
    fillColor = 'red',
    opacity = .5,
    title = 'H+ trend (nmol/kg/year)',
    shape = 'circle',
    orientation = 'vertical',
    breaks = 5)



# Plot H+ trends with size by discharge
col_types <- c("darkblue","dodgerblue","green4","gold1","orange","brown","red")
leg_vals <- unique(as.numeric(quantile(trends_sig$pH_ts_slope*365*1e9,
                                       probs=c(0,0.01,0.1,0.25,0.5,0.75,0.9,.99,1), na.rm=TRUE)))
pal = colorBin(col_types, trends_sig$pH_ts_slope*365*1e9, bins = leg_vals)
leaflet(data=trends_sig) %>%
  addProviderTiles("CartoDB.Positron") %>%
  addCircleMarkers(trends_nonsig$lon,trends_nonsig$lat,
                   color = "black", radius=4, stroke=FALSE,
                   fillOpacity = 0.8, opacity = 0.8,
                   popup=paste(trends_nonsig$site_number)) %>%
  addCircleMarkers(trends_sig$lon,trends_sig$lat,
                   color = pal(trends_sig$pH_ts_slope*365*1e9), radius=4, stroke=FALSE,
                   fillOpacity = 0.8, opacity = 0.8,
                   popup=paste(trends_sig$site_number)) %>%
  addLegend(position = 'bottomleft',
            pal=pal,
            values=trends_sig$pH_ts_slope*365*1e9,
            opacity = 0.8,
            labFormat = labelFormat(digits = 0),
            title = 'H+ Trend (nmol/year)')




# Plot alkalinity trends with size by discharge
trends_sig_alk <- trends_final %>%
  filter(alk_ts_p_value < 0.05)
trends_nonsig_alk <- trends_final %>%
  filter(alk_ts_p_value > 0.05)
trends_sig_pos_alk <- trends_sig %>%
  filter(alk_ts_slope > 0)
trends_sig_neg_alk <- trends_sig %>%
  filter(alk_ts_slope < 0)
# Plot Alkalinity trends with size by discharge
col_types <- c("darkblue","dodgerblue","green4","gold1","orange","brown","red")
leg_vals <- unique(as.numeric(quantile(trends_sig_alk$alk_ts_slope*365,
                                       probs=c(0,0.01,0.1,0.25,0.5,0.75,0.9,.99,1), na.rm=TRUE)))
pal = colorBin(col_types, trends_sig_alk$alk_ts_slope*365, bins = leg_vals)
leaflet(data=trends_sig_alk) %>%
  addProviderTiles("CartoDB.Positron") %>%
  addCircleMarkers(trends_nonsig_alk$lon,trends_nonsig_alk$lat,
                   color = "black", radius=4, stroke=FALSE,
                   fillOpacity = 0.8, opacity = 0.8,
                   popup=paste(trends_nonsig_alk$site_number)) %>%
  addCircleMarkers(trends_sig_alk$lon,trends_sig_alk$lat,
                   color = pal(trends_sig_alk$alk_ts_slope*365), radius=4, stroke=FALSE,
                   fillOpacity = 0.8, opacity = 0.8,
                   popup=paste(trends_sig_alk$site_number)) %>%
  addLegend(position = 'bottomleft',
            pal=pal,
            values=trends_sig_alk$alk_ts_slope*365,
            opacity = 0.8,
            labFormat = labelFormat(digits = 0),
            title = 'Alkalinity Trend (umol/year)')

######################################### Use TS analysis to calculate total changes in alkalinity & H+ through time

trends_final$ts_initial_alks <- trends_final$alk_ts_slope*trends_final$min_date + trends_final$alk_ts_intercept_value
trends_final$ts_final_alks <- trends_final$alk_ts_slope*trends_final$max_date + trends_final$alk_ts_intercept_value
trends_final$ts_delta_alks <- trends_final$ts_final_alks - trends_final$ts_initial_alks

trends_final$ts_initial_pH <- trends_final$pH_ts_slope*trends_final$min_date + trends_final$pH_ts_intercept_value
trends_final$ts_final_pH <- trends_final$pH_ts_slope*trends_final$max_date + trends_final$pH_ts_intercept_value
trends_final$ts_delta_pH <- trends_final$ts_final_pH - trends_final$ts_initial_pH

