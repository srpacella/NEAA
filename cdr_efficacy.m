% Calculate partial derivatives of stream carbonate systems with respect to
% to a change in alkalinity 

% This addresses the question "where will alkalinity addition be most
% effective at reducing CO2 emissions to the atmosphere".  This only looks
% at stream sites - important to think about mixing with ocean waters too.
%% Calculate partial derivatives w/r/t deltaALK for all river sites     
for i = 1:length(all_sites_names)
    if sum(index_highalk==i) > 0 %find sites with high alk and no required correction
        DIC = site_cc.(all_sites_names2{i})(:,2);
        ALK = all_data_sites.(all_sites_names2{i}).Alkalinity;
        temperature = all_data_sites.(all_sites_names2{i}).Temperature;
        salinity = all_data_sites.(all_sites_names2{i}).salinity;
        pressure = 0;
        SIL = 0;
        PO4 = 0;
        sens_sites_org_alk.(all_sites_names2{i}) = derivnum ('par1',ALK,DIC,1,2,salinity,temperature,temperature,pressure,pressure,... 
                                   SIL,PO4,0,0,...
                                   SCALE,K1K2,SO4,KF,BOR);
    elseif sum(index_lowalk_withorg==i) > 0 %find sites with low alk and correction possible
        DIC = site_cc.(all_sites_names2{i})(:,2);
        ALK = all_data_sites.(all_sites_names2{i}).Alkalinity.* (1-(alkalinity_all(i,4)./100));
        temperature = all_data_sites.(all_sites_names2{i}).Temperature;
        salinity = all_data_sites.(all_sites_names2{i}).salinity;
        pressure = 0;
        SIL = 0;
        PO4 = 0;
        sens_sites_org_alk.(all_sites_names2{i}) = derivnum ('par1',ALK,DIC,1,2,salinity,temperature,temperature,pressure,pressure,... 
                                   SIL,PO4,0,0,...
                                   SCALE,K1K2,SO4,KF,BOR);
    else
    end
    
end

%% Plot stream sensitivity factors

lat = str2double(all_sites_info(:,7));
lon = str2double(all_sites_info(:,8));
drainage = str2double(all_sites_info(:,30))*2.5899;

for i = 1:length(index_allcorrect)
    lat_orgstreams(i) = lat(index_allcorrect(i));
    lon_orgstreams(i) = lon(index_allcorrect(i));
    discharge_orgstreams(i) = discharge_means_numeric(index_allcorrect(i));
    pH_orgstreams_median(i) = nanmedian(site_cc_org.(names_allcorrect{i})(:,43));
    pco2_orgstreams_median(i) = nanmedian(site_cc_org.(names_allcorrect{i})(:,22));
    alk_orgstreams_median(i) = nanmedian(site_cc_org.(names_allcorrect{i})(:,1));
    sensco2_orgstreams_median(i) = nanmedian(sens_sites_org_alk.(names_allcorrect{i})(:,14)); %dpCO2/dALK
end

%% Figures

figure
geoscatter(lat_orgstreams,lon_orgstreams,-sensco2_orgstreams_median.*10,sensco2_orgstreams_median,'filled');
%gb.ColorVariable = sensH_orgstreams_median;
c = colorbar;
c.Label.String = "Median stream \deltapCO_2/\deltaALK";
geobasemap colorterrain
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';
cmocean('thermal')
%caxis([0 4]); % Truncates range to reduce effect of few outlier values - needs to be manually edited on figure

figure
%geoscatter(lat_orgstreams,lon_orgstreams,discharge_orgstreams./100,sensco2_orgstreams_median,'filled');
geoscatter(lat_orgstreams(high_discharge_index),lon_orgstreams(high_discharge_index),discharge_orgstreams(high_discharge_index).*conv_cfs/10,sensco2_orgstreams_median(high_discharge_index),'filled');
hold on
geoscatter(lat_orgstreams(low_discharge_index),lon_orgstreams(low_discharge_index),discharge_orgstreams(low_discharge_index).*conv_cfs/10,sensco2_orgstreams_median(low_discharge_index),'filled');
%gb.ColorVariable = sensH_orgstreams_median;
c = colorbar;
c.Label.String = "Median stream \deltapCO_2/\deltaALK";
geobasemap colorterrain
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';
cmocean('thermal')

figure
scatter(discharge_orgstreams,sensco2_orgstreams_median)



figure
scatter(discharge_orgstreams,sensco2_orgstreams_median)
xlabel('Stream discharge')
ylabel('\partialpCO_2/\partialALK')
box on
grid on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';

figure
scatter(alk_orgstreams_median,sensco2_orgstreams_median)
xlabel('Stream median alkalinity')
ylabel('\partialpCO_2/\partialALK')
box on
grid on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';

%%
% Find sites with large discharge, plot them before plotting the other
% sites to make the figure more legible
high_discharge_index = find(discharge_orgstreams > 1e5)
low_discharge_index = find(discharge_orgstreams < 1e5)
figure
subplot(2,2,1)
geoscatter(lat_orgstreams,lon_orgstreams,30,sensH_orgstreams_median,'filled');
%gb.ColorVariable = sensH_orgstreams_median;
c = colorbar;
c.Label.String = "Median stream H^+ sensitivity factor";
geobasemap landcover
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',24,'LineWidth',1)
fig = gcf
fig.Color='w';
cmocean('thermal')
caxis([0 3]);
subplot(2,2,2)
geoscatter(lat_orgstreams(high_discharge_index),lon_orgstreams(high_discharge_index),discharge_orgstreams(high_discharge_index).*conv_cfs/10,sensH_orgstreams_median(high_discharge_index),'filled');
hold on
geoscatter(lat_orgstreams(low_discharge_index),lon_orgstreams(low_discharge_index),discharge_orgstreams(low_discharge_index).*conv_cfs/10,sensH_orgstreams_median(low_discharge_index),'filled');
%gb.ColorVariable = sensH_orgstreams_median;
c = colorbar;
c.Label.String = "Median stream H^+ sensitivity factor";
geobasemap landcover
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',24,'LineWidth',1)
fig = gcf
fig.Color='w';
cmocean('thermal')
caxis([0 3]);

subplot(2,2,3)
geoscatter(lat_orgstreams,lon_orgstreams,30,sensco2_orgstreams_median,'filled');
%gb.ColorVariable = sensH_orgstreams_median;
c = colorbar;
c.Label.String = "Median stream \deltapCO_2/\deltaALK";
geobasemap landcover
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',24,'LineWidth',1)
fig = gcf
fig.Color='w';
cmocean('thermal')
subplot(2,2,4)
geoscatter(lat_orgstreams(high_discharge_index),lon_orgstreams(high_discharge_index),discharge_orgstreams(high_discharge_index).*conv_cfs/10,sensco2_orgstreams_median(high_discharge_index),'filled');
hold on
geoscatter(lat_orgstreams(low_discharge_index),lon_orgstreams(low_discharge_index),discharge_orgstreams(low_discharge_index).*conv_cfs/10,sensco2_orgstreams_median(low_discharge_index),'filled');
%gb.ColorVariable = sensH_orgstreams_median;
c = colorbar;
c.Label.String = "Median stream \deltapCO_2/\deltaALK";
geobasemap landcover
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',24,'LineWidth',1)
fig = gcf
fig.Color='w';
cmocean('thermal')


