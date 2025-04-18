clear
clear all

% Read in .csv files exported from R via "export_neaa_csv.R"

%addpath('/Users/spacella/Desktop/R Code May 2022 SRP Review/neaa_data_exports/')
addpath('C:\Users\spacella\OneDrive - Environmental Protection Agency (EPA)\National Estuary Acidification Assessment\NEAA Analysis Files 240814\R Code May 2022 SRP Review\neaa_data_exports')
addpath(genpath('C:\Users\spacella\OneDrive - Environmental Protection Agency (EPA)\National Estuary Acidification Assessment\NEAA Analysis Files 240814\R Code May 2022 SRP Review\CO2-System-Extd-master'))

%% Import names of all sites and remove header
all_sites_names = readmatrix("C:\Users\spacella\OneDrive - Environmental Protection Agency (EPA)\National Estuary Acidification Assessment\NEAA Analysis Files 240814\R Code May 2022 SRP Review\neaa_data_exports/final_sites_names.csv",'OutputType','char');
all_sites_names = all_sites_names(2:length(all_sites_names));
all_sites_names2=regexprep(all_sites_names,'-','');

%% Import site information
all_sites_info = readmatrix("C:\Users\spacella\OneDrive - Environmental Protection Agency (EPA)\National Estuary Acidification Assessment\NEAA Analysis Files 240814\R Code May 2022 SRP Review\neaa_data_exports/siteINFO.csv",'OutputType','char');
hucs = all_sites_info(:,24);
all_sites_names_test = strcat(all_sites_info(:,1),all_sites_info(:,2));

%% Harmonize order of site names "all_sites_names2" and site information "all_sites_info"
all_sites_names2 = sort(all_sites_names2);
all_sites_names = sort(all_sites_names);

%% Read in each site's data, rename columns, write to 
for i = 1:length(all_sites_names)
    filename = strcat("C:\Users\spacella\OneDrive - Environmental Protection Agency (EPA)\National Estuary Acidification Assessment\NEAA Analysis Files 240814\R Code May 2022 SRP Review\neaa_data_exports/",all_sites_names(i),".csv");
    test = readtable(filename);
    all_data_sites.(all_sites_names2{i}) = test;
end

%% Read in each site's coastal ocean data, rename columns, write to 
for i = 1:length(all_sites_names)
    filename = strcat("C:\Users\spacella\OneDrive - Environmental Protection Agency (EPA)\National Estuary Acidification Assessment\NEAA Analysis Files 240814\R Code May 2022 SRP Review\neaa_data_exports/",all_sites_names(i),"_ocean.csv");
    test = readtable(filename);
    all_ocean_cc.(all_sites_names2{i}) = test;
end

%% Read in each site's calcium data, rename columns, write to 
for i = 1:length(all_sites_names)
    filename = strcat("C:\Users\spacella\OneDrive - Environmental Protection Agency (EPA)\National Estuary Acidification Assessment\NEAA Analysis Files 240814\R Code May 2022 SRP Review\neaa_data_exports/",all_sites_names(i),"_calcium.csv");
    test = readtable(filename);
    all_calcium.(all_sites_names2{i}) = test;
end

%% Check for duplicate observations in stream all_data_sites
for i = 1:length(all_sites_names)
    [C,ia,ib] = unique(all_data_sites.(all_sites_names2{i}).ID) ; 
    data_unique = all_data_sites.(all_sites_names2{i})(ia,:);
    all_data_sites.(all_sites_names2{i}) = data_unique;
end

