
library(leaflet)
library(sf)
library(dplyr)
library(ggplot2)
library(cowplot)

sens_estuary_mean <- read_excel("sens_estuary_mean.xlsx",col_names=FALSE) 

# All final sites
leaflet(data=siteINFO) %>%
  addProviderTiles("CartoDB.Positron") %>%
  addCircleMarkers(~dec_long_va,~dec_lat_va,
                   color = "blue", radius=4, stroke=FALSE,
                   fillOpacity = 0.8, opacity = 0.8,
                   popup=siteINFO$station_nm)

# Sites colored & sized by mean sensitivity
col_types <- c("darkblue","dodgerblue","green4","gold1","orange","brown","red")
leg_vals <- unique(as.numeric(quantile(sens_estuary_mean$...1, 
                                       probs=c(0,0.01,0.1,0.25,0.5,0.75,0.9,.99,1), na.rm=TRUE)))
pal = colorBin(col_types, sens_estuary_mean$...1, bins = leg_vals)

leaflet(data=siteINFO) %>% 
  addProviderTiles("CartoDB.Positron") %>%
  addCircleMarkers(siteINFO$dec_long_va,siteINFO$dec_lat_va,
                   fillColor = pal(sens_estuary_mean$...1), 
                   radius=sens_estuary_mean$...1*10, stroke=FALSE,
                   fillOpacity = 0.8, opacity = 0.8,
                   popup=~station_nm) %>%
  addLegend(position = 'bottomleft',
            pal=pal,
            values=sens_estuary_mean,
            opacity = 0.8,
            labFormat = labelFormat(digits = 1), 
            title = '[H+]/[DIC]')

# Sites colored by mean sensitivity
col_types <- c("darkblue","dodgerblue","green4","gold1","orange","brown","red")
leg_vals <- unique(as.numeric(quantile(sens_estuary_mean$...1, 
                                       probs=c(0,0.01,0.1,0.25,0.5,0.75,0.9,.99,1), na.rm=TRUE)))
pal = colorBin(col_types, sens_estuary_mean$...1, bins = leg_vals)

leaflet(data=siteINFO) %>% 
  addProviderTiles("CartoDB.Positron") %>%
  addCircleMarkers(siteINFO$dec_long_va,siteINFO$dec_lat_va,
                   fillColor = pal(sens_estuary_mean$...1), 
                   radius=6, stroke=FALSE,
                   fillOpacity = 0.8, opacity = 0.8,
                   popup=~station_nm) %>%
  addLegend(position = 'bottomleft',
            pal=pal,
            values=sens_estuary_mean,
            opacity = 0.8,
            labFormat = labelFormat(digits = 1), 
            title = '[H+]/[DIC]')