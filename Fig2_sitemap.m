%% Figure 2 - Map of all sites color coded by water resource region

% Create map of stream sites colored by water resource region
sites_01 = startsWith(hucs,'01'); % New England
sites_02 = startsWith(hucs,'02'); % Mid-Atlantic
sites_03 = startsWith(hucs,'03'); % South Atlantic-Gulf
sites_08 = startsWith(hucs,'08'); % Lower Mississippi
sites_12 = startsWith(hucs,'12'); % Texas-Gulf
sites_13 = startsWith(hucs,'13'); % Rio Grande
sites_18 = startsWith(hucs,'18'); % California
sites_17 = startsWith(hucs,'17'); % Pacific Northwest
sites_19 = startsWith(hucs,'19'); % Alaska

figure
for i = 1:length(all_sites_info)
    if sites_01(i) == 1
    lat=str2double(all_sites_info(i,7));
    lon = str2double(all_sites_info(i,8));
    geoscatter(lat,lon,50,[0, 0.4470, 0.7410],'filled','MarkerEdgeColor',[0 0 0]);
    end
    hold on
end
hold on
for i = 1:length(all_sites_info)
    if sites_02(i) == 1
    lat=str2double(all_sites_info(i,7));
    lon = str2double(all_sites_info(i,8));
    geoscatter(lat,lon,50,[0.8500, 0.3250, 0.0980],'filled','MarkerEdgeColor',[0 0 0]);
    end
    hold on
end
hold on
for i = 1:length(all_sites_info)
    if sites_03(i) == 1
    lat=str2double(all_sites_info(i,7));
    lon = str2double(all_sites_info(i,8));
    geoscatter(lat,lon,50,[0.9290, 0.6940, 0.1250],'filled','MarkerEdgeColor',[0 0 0]);
    end
    hold on
end
hold on
for i = 1:length(all_sites_info)
    if sites_08(i) == 1
    lat=str2double(all_sites_info(i,7));
    lon = str2double(all_sites_info(i,8));
    geoscatter(lat,lon,50,[0.4940, 0.1840, 0.5560],'filled','MarkerEdgeColor',[0 0 0]);
    end
    hold on
end
hold on
for i = 1:length(all_sites_info)
    if sites_12(i) == 1
    lat=str2double(all_sites_info(i,7));
    lon = str2double(all_sites_info(i,8));
    geoscatter(lat,lon,50,[0.4660, 0.6740, 0.1880],'filled','MarkerEdgeColor',[0 0 0]);
    end
    hold on
end
hold on
for i = 1:length(all_sites_info)
    if sites_13(i) == 1
    lat=str2double(all_sites_info(i,7));
    lon = str2double(all_sites_info(i,8));
    geoscatter(lat,lon,50,[0.3010, 0.7450, 0.9330],'filled','MarkerEdgeColor',[0 0 0]);
    end
    hold on
end
hold on
for i = 1:length(all_sites_info)
    if sites_18(i) == 1
    lat=str2double(all_sites_info(i,7));
    lon = str2double(all_sites_info(i,8));
    geoscatter(lat,lon,50,[0.6350, 0.0780, 0.1840],'filled','MarkerEdgeColor',[0 0 0]);
    end
    hold on
end
hold on
for i = 1:length(all_sites_info)
    if sites_17(i) == 1
    lat=str2double(all_sites_info(i,7));
    lon = str2double(all_sites_info(i,8));
    geoscatter(lat,lon,50,[0, 0.5, 0],'filled','MarkerEdgeColor',[0 0 0]);
    end
    hold on
end
hold on
for i = 1:length(all_sites_info)
    if sites_19(i) == 1
    lat=str2double(all_sites_info(i,7));
    lon = str2double(all_sites_info(i,8));
    geoscatter(lat,lon,50,[0.75, 0, 0.75],'filled','MarkerEdgeColor',[0 0 0]);
    end
    hold on
end
geobasemap topographic
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';