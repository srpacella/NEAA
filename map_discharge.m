%% Create summary statistics for discharge

%Important discharge data for all sites as of 8/26/22
discharge_means = readmatrix("C:\Users\spacella\OneDrive - Environmental Protection Agency (EPA)\National Estuary Acidification Assessment\NEAA Analysis Files 240814\R Code May 2022 SRP Review\neaa_data_exports/discharge_means.csv",'OutputType','char');
discharge_means = sortrows(discharge_means,1);
discharge_means_numeric = str2double(discharge_means(:,2));
discharge_means_names = discharge_means(:,1);
discharge_means_names=regexprep(discharge_means_names,'-','');
% discharge_means = readmatrix("/Users/spacella/Desktop/R Code May 2022 SRP Review/neaa_data_exports/discharge_means.csv");
% discharge_means = discharge_means(:,2);
conv_cfs = 1/35.3147 %converts cfs to m3/s

discharge_all = struct()
for i = 1:length(discharge_means)
    discharge_all.(discharge_means_names{i}) = discharge_means_numeric(i);
end


%Append discharge means to all_sites_info for each station
all_sites_info(:,43) = discharge_means(:,1);
all_sites_info(:,44) = discharge_means(:,2);

figure
lat=str2double(all_sites_info(:,7));
lon = str2double(all_sites_info(:,8));
gb = geobubble(lat,lon,discharge_means_numeric.*conv_cfs,'MapLayout','maximized','BubbleColorList',[0 0 1]);
gb.SizeLegendTitle = 'Mean discharge (m^3/s)';
%geobasemap colorterrain
geobasemap landcover

total_discharge = sum(discharge_means_numeric)
total_org_discharge = sum(discharge_means_numeric(index_allcorrect))
percent_org_discharge = total_org_discharge./total_discharge*100