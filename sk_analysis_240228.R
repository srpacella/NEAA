# Compute seasonal kendall trend analyses for river data intead of Thiel-Sen
# Uses package rkt
# rkt(date, y, block, cv, correct = F, rep = "e")

# date	
# a mandatory vector of numerical data representing dates, as years or years+decimal. If correction for intra-block correlation is required, dates will be truncated to the year, and no more than one value per block per year will be considered. If two equal dates (or truncated dates) are found, the behaviour of the program is determined by rep
# y	
# a mandatory vector of measured data. In this vector, missing data are allowed.
# block	
# an optional vector of positive integer numbers representing blocks, i.e. sites, seasons or months, or a code combining both sites and seasons/months. If no blocks are defined, the result will be the Mann-Kendall test.
# cv	
# an optional vector containing a covariable, such as river flow or deposition amount. In this vector, missing data are allowed, however all case for which covariable value is missing are deleted from the analysis. As a consequence, if a covariable with missing data is passed to this function, the Kendall score, tau and p-value will be different than without covariable.
# correct	
# a boolean value. If correct is FALSE, no correction for correlation between blocks is performed. If correct is TRUE, dates are truncated and the correction for correlation between blocks is performed. Note that the truncation is performed in any case, while the correction is performed only if there are more than one block, and more than nine years of data. Default value is FALSE.
# rep	
# a character value. If rep is set to "a", data sharing the same date (or truncated date if correct is TRUE) are averaged. If rep is set to "m", their median is used. For any other value of rep, an error is produced if two or more data share the same date (or truncated date if correct is TRUE). The latter is the default behaviour of the program.


# Compute seasonal kendall trend analysis, with monthly bins
kendall_test_month <- list()
#kendall_slopes_month <- vector(mode="numeric")
for(j in names_final){
  count = 1
  date <- decimal_date(all_data_site_final[[j]]$ActivityStartDate) #convert date to decimal year
  y <- 10^(-all_data_site_final[[j]]$pH_correct) #observed H+
  block <- month(all_data_site_final[[j]]$ActivityStartDate)
  kendall_test_month[[j]] <- rkt(date, y, block, correct = F, rep = "m")
}

kendall_test_month_all <- do.call(rbind.data.frame,kendall_test_month)

# Compute seasonal kendall trend analysis, with monthly bins
kendall_test <- list()
#kendall_slopes_month <- vector(mode="numeric")
for(j in names_final){
  count = 1
  date <- decimal_date(all_data_site_final[[j]]$ActivityStartDate) #convert date to decimal year
  y <- 10^(-all_data_site_final[[j]]$pH_correct) #observed H+
  kendall_test[[j]] <- rkt(date, y)
}

kendall_test_all <- do.call(rbind.data.frame,kendall_test)
#######################################################################3
# Load ts slopes for H+ from ts_analysis_220816.R
plot(DF1$pH_ts_slope*365,kendall_slopes_month2)
plot(kendall_slopes_none,kendall_slopes_month2)

# Add lat/lon to DF1
DF1$siteNumbers_edit <- DF1$site_number
DF1$siteNumbers_edit <- sub("USGS-", "", DF1$site_number)
count <- 0
for(i in DF1$siteNumbers_edit){
  count = count + 1
  DF1$lat[count] <- siteINFO_final$dec_lat_va[which(siteINFO_final$site_no==i)]
  DF1$lon[count] <- siteINFO_final$dec_long_va[which(siteINFO_final$site_no==i)]
  site_no_temp <- DF1$site_number[count]
  DF1$discharge[count] <- mean(download_discharge[[site_no_temp]]$ResultMeasureValue)
}

# Plot trends
leaflet(data=DF1) %>%
  addProviderTiles("CartoDB.Positron") %>%
  addCircleMarkers(DF1$lon,DF1$lat,
                   color = "blue", radius=DF1$pH_ts_slope*100000000000, stroke=FALSE,
                   fillOpacity = 0.8, opacity = 0.8,
                   popup=paste(DF1$siteNumbers)) 

leaflet(data=DF1) %>%
  addProviderTiles("CartoDB.Positron") %>%
  addCircleMarkers(DF1$lon,DF1$lat,
                   color = "blue", radius=unname(kendall_slopes_month2)*1000000000, stroke=FALSE,
                   fillOpacity = 0.8, opacity = 0.8,
                   popup=paste(DF1$siteNumbers)) 


# Plot trends
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







set.seed(1234)
x <- 1901:2000
y <- x+rnorm(100)
y[100] <- 200
fit_sk <- rkt(x, y)
fit <- mblm(y~x)
fit
fit_sk
summary(fit)
fit2 <- lm(y~x)
plot(x,y)
abline(fit)
plot(x,y)
abline(fit$coefficients[1],fit_sk$B, lty=2)
abline(fit2,lty=2)
plot(fit)
residuals(fit)
fitted(fit)
plot(density(fit$slopes))
plot(density(fit$intercepts))
anova(fit)
anova(fit2)
anova(fit,fit2)
confint(fit)
AIC(fit,fit2) 




fit_mblm <- mblm(all_data_site_final$`USGS-01021050`$pH~all_data_site_final$`USGS-01021050`$ActivityStartDate)

