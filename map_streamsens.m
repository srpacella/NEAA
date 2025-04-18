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
    sensH_orgstreams_median(i) = nanmedian(sens_sites_org.(names_allcorrect{i})(:,13));
end


% figure
% gb = geobubble(lat_orgstreams,lon_orgstreams,sensH_orgstreams_median);
% gb.SizeLegendTitle = 'Median H+ sensitivity factor';
% geobasemap colorterrain
% %geobasemap landcover
% 
% %%
figure
geoscatter(lat_orgstreams,lon_orgstreams,sensH_orgstreams_median.*50,sensH_orgstreams_median,'filled');
%gb.ColorVariable = sensH_orgstreams_median;
c = colorbar;
c.Label.String = "Median stream H^+ sensitivity factor";
geobasemap colorterrain
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';
caxis([0 4]); % Truncates range to reduce effect of few outlier values - needs to be manually edited on figure

figure
subplot(1,2,1)
geoscatter(lat_orgstreams,lon_orgstreams,30,sensH_orgstreams_median,'filled');
%gb.ColorVariable = sensH_orgstreams_median;
c = colorbar;
c.Label.String = "Median stream H^+ sensitivity factor";
geobasemap topographic
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';
cmocean('thermal')
caxis([0 3]);
subplot(1,2,2)
geoscatter(lat_orgstreams,lon_orgstreams,discharge_orgstreams.*conv_cfs/2,sensH_orgstreams_median,'filled');
%gb.ColorVariable = sensH_orgstreams_median;
c = colorbar;
c.Label.String = "Median stream H^+ sensitivity factor";
geobasemap topographic
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';
caxis([0 3]);
cmocean('thermal')
%%
% Find sites with large discharge, plot them before plotting the other
% sites to make the figure more legible
high_discharge_index = find(discharge_orgstreams > 1e5)
low_discharge_index = find(discharge_orgstreams < 1e5)
figure
subplot(1,2,1)
geoscatter(lat_orgstreams,lon_orgstreams,30,sensH_orgstreams_median,'filled');
%gb.ColorVariable = sensH_orgstreams_median;
c = colorbar;
c.Label.String = "Median stream H^+ sensitivity factor";
geobasemap satellite
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';
cmocean('thermal')
caxis([0 3]);
subplot(1,2,2)
geoscatter(lat_orgstreams(high_discharge_index),lon_orgstreams(high_discharge_index),discharge_orgstreams(high_discharge_index).*conv_cfs/10,sensH_orgstreams_median(high_discharge_index),'filled');
hold on
geoscatter(lat_orgstreams(low_discharge_index),lon_orgstreams(low_discharge_index),discharge_orgstreams(low_discharge_index).*conv_cfs/10,sensH_orgstreams_median(low_discharge_index),'filled');
%gb.ColorVariable = sensH_orgstreams_median;
c = colorbar;
c.Label.String = "Median stream H^+ sensitivity factor";
geobasemap landcover
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';
cmocean('thermal')
caxis([0 3]);

%% Make 6 panel plot with ph, alkalinity, H+ sensitivity
% Find sites with large discharge, plot them before plotting the other
% sites to make the figure more legible
high_discharge_index = find(discharge_orgstreams > 1e5)
low_discharge_index = find(discharge_orgstreams < 1e5)
figure
subplot(3,2,1)
geoscatter(lat_orgstreams,lon_orgstreams,30,pH_orgstreams_median,'filled');
%gb.ColorVariable = sensH_orgstreams_median;
c = colorbar;
c.Label.String = "Median stream pH_T";
geobasemap landcover
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';
cmocean('thermal')
caxis([5.7 8.7]);
charlbl =  compose("(%s)",('a':'z').'); 
text(0.025,0.95,charlbl{1},'Units','normalized','FontSize',12)
subplot(3,2,2)
geoscatter(lat_orgstreams(high_discharge_index),lon_orgstreams(high_discharge_index),discharge_orgstreams(high_discharge_index).*conv_cfs/10,pH_orgstreams_median(high_discharge_index),'filled');
hold on
geoscatter(lat_orgstreams(low_discharge_index),lon_orgstreams(low_discharge_index),discharge_orgstreams(low_discharge_index).*conv_cfs/10,pH_orgstreams_median(low_discharge_index),'filled');
%gb.ColorVariable = sensH_orgstreams_median;
c = colorbar;
c.Label.String = "Median stream pH_T";
geobasemap landcover
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';
cmocean('thermal')
caxis([5.7 8.7]);
charlbl =  compose("(%s)",('a':'z').'); 
text(0.025,0.95,charlbl{2},'Units','normalized','FontSize',12)

subplot(3,2,3)
geoscatter(lat_orgstreams,lon_orgstreams,30,alk_orgstreams_median,'filled');
%gb.ColorVariable = sensH_orgstreams_median;
c = colorbar;
c.Label.String = "Median stream alkalinity (\mueq L^-^1)";
geobasemap landcover
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';
cmocean('thermal')
caxis([0 6600]);
charlbl =  compose("(%s)",('a':'z').'); 
text(0.025,0.95,charlbl{3},'Units','normalized','FontSize',12)
subplot(3,2,4)
geoscatter(lat_orgstreams(high_discharge_index),lon_orgstreams(high_discharge_index),discharge_orgstreams(high_discharge_index).*conv_cfs/10,alk_orgstreams_median(high_discharge_index),'filled');
hold on
geoscatter(lat_orgstreams(low_discharge_index),lon_orgstreams(low_discharge_index),discharge_orgstreams(low_discharge_index).*conv_cfs/10,alk_orgstreams_median(low_discharge_index),'filled');
%gb.ColorVariable = sensH_orgstreams_median;
c = colorbar;
c.Label.String = "Median stream alkalinity (\mueq L^-^1)";
geobasemap landcover
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';
cmocean('thermal')
caxis([0 6600]);
charlbl =  compose("(%s)",('a':'z').'); 
text(0.025,0.95,charlbl{4},'Units','normalized','FontSize',12)

subplot(3,2,5)
geoscatter(lat_orgstreams,lon_orgstreams,30,sensH_orgstreams_median,'filled');
%gb.ColorVariable = sensH_orgstreams_median;
c = colorbar;
c.Label.String = "Median stream H^+ sensitivity factor";
geobasemap landcover
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';
cmocean('thermal')
caxis([0 3]);
charlbl =  compose("(%s)",('a':'z').'); 
text(0.025,0.95,charlbl{5},'Units','normalized','FontSize',12)
subplot(3,2,6)
geoscatter(lat_orgstreams(high_discharge_index),lon_orgstreams(high_discharge_index),discharge_orgstreams(high_discharge_index).*conv_cfs/10,sensH_orgstreams_median(high_discharge_index),'filled');
hold on
geoscatter(lat_orgstreams(low_discharge_index),lon_orgstreams(low_discharge_index),discharge_orgstreams(low_discharge_index).*conv_cfs/10,sensH_orgstreams_median(low_discharge_index),'filled');
%gb.ColorVariable = sensH_orgstreams_median;
c = colorbar;
c.Label.String = "Median stream H^+ sensitivity factor";
geobasemap landcover
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';
cmocean('thermal')
caxis([0 3]);
charlbl =  compose("(%s)",('a':'z').'); 
text(0.025,0.95,charlbl{6},'Units','normalized','FontSize',12)

%% Make 6 panel plot with ph, alkalinity, pCO2
% Find sites with large discharge, plot them before plotting the other
% sites to make the figure more legible
high_discharge_index = find(discharge_orgstreams > 1e5)
low_discharge_index = find(discharge_orgstreams < 1e5)
figure
subplot(3,2,1)
geoscatter(lat_orgstreams,lon_orgstreams,30,pH_orgstreams_median,'filled');
%gb.ColorVariable = sensH_orgstreams_median;
% c = colorbar;
% c.Label.String = "Median stream pH_T";
geobasemap landcover
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';
cmocean('thermal')
% caxis([5.7 8.7]);
charlbl =  compose("(%s)",('a':'z').'); 
text(0.025,0.95,charlbl{1},'Units','normalized','FontSize',12)
subplot(3,2,2)
geoscatter(lat_orgstreams(high_discharge_index),lon_orgstreams(high_discharge_index),discharge_orgstreams(high_discharge_index).*conv_cfs/10,pH_orgstreams_median(high_discharge_index),'filled');
hold on
geoscatter(lat_orgstreams(low_discharge_index),lon_orgstreams(low_discharge_index),discharge_orgstreams(low_discharge_index).*conv_cfs/10,pH_orgstreams_median(low_discharge_index),'filled');
%gb.ColorVariable = sensH_orgstreams_median;
c = colorbar;
c.Label.String = "Median stream pH_T";
geobasemap landcover
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';
cmocean('thermal')
caxis([5.7 8.7]);
charlbl =  compose("(%s)",('a':'z').'); 
text(0.025,0.95,charlbl{2},'Units','normalized','FontSize',12)

subplot(3,2,3)
geoscatter(lat_orgstreams,lon_orgstreams,30,alk_orgstreams_median,'filled');
%gb.ColorVariable = sensH_orgstreams_median;
% c = colorbar;
% c.Label.String = "Median stream alkalinity (\mueq L^-^1)";
geobasemap landcover
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';
cmocean('thermal')
% caxis([0 6600]);
charlbl =  compose("(%s)",('a':'z').'); 
text(0.025,0.95,charlbl{3},'Units','normalized','FontSize',12)
subplot(3,2,4)
geoscatter(lat_orgstreams(high_discharge_index),lon_orgstreams(high_discharge_index),discharge_orgstreams(high_discharge_index).*conv_cfs/10,alk_orgstreams_median(high_discharge_index),'filled');
hold on
geoscatter(lat_orgstreams(low_discharge_index),lon_orgstreams(low_discharge_index),discharge_orgstreams(low_discharge_index).*conv_cfs/10,alk_orgstreams_median(low_discharge_index),'filled');
%gb.ColorVariable = sensH_orgstreams_median;
c = colorbar;
c.Label.String = "Median stream alkalinity (\mueq L^-^1)";
geobasemap landcover
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';
cmocean('thermal')
caxis([0 6600]);
charlbl =  compose("(%s)",('a':'z').'); 
text(0.025,0.95,charlbl{4},'Units','normalized','FontSize',12)

subplot(3,2,5)
geoscatter(lat_orgstreams,lon_orgstreams,30,pco2_orgstreams_median,'filled');
%gb.ColorVariable = sensH_orgstreams_median;
% c = colorbar;
% c.Label.String = "Median stream pCO_2 (\muatm)";
geobasemap landcover
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';
cmocean('thermal')
% caxis([0 5000]);
charlbl =  compose("(%s)",('a':'z').'); 
text(0.025,0.95,charlbl{5},'Units','normalized','FontSize',12)
subplot(3,2,6)
geoscatter(lat_orgstreams(high_discharge_index),lon_orgstreams(high_discharge_index),discharge_orgstreams(high_discharge_index).*conv_cfs/10,pco2_orgstreams_median(high_discharge_index),'filled');
hold on
geoscatter(lat_orgstreams(low_discharge_index),lon_orgstreams(low_discharge_index),discharge_orgstreams(low_discharge_index).*conv_cfs/10,pco2_orgstreams_median(low_discharge_index),'filled');
%gb.ColorVariable = sensH_orgstreams_median;
c = colorbar;
c.Label.String = "Median stream pCO_2 (\muatm)";
geobasemap landcover
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';
cmocean('thermal')
caxis([0 5000]);
charlbl =  compose("(%s)",('a':'z').'); 
text(0.025,0.95,charlbl{6},'Units','normalized','FontSize',12)


