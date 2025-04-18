% Calculate carbonate system of CODAP sites
% Requires output from neaa_import.m
%addpath('/Users/spacella/Desktop/R Code May 2022 SRP Review/neaa_data_exports/')
%addpath('/Users/spacella/Desktop/R Code May 2022 SRP Review/CO2-System-Extd-master/main/')


% Define constants for carbonate system calculations
SCALE  = 1; % Total scale
K1K2   = 10; %= Lueker et al, 2000                                  T:    2-35  S: 19-43. Total scale. Real seawater.
SO4    = 1; % Dickson (1990) KSO4
KF     = 2; % Perez & Fraga (1987) KF
BOR    = 2; % Lee et al (2010) TB



%% Calculate partial derivatves w/r/t deltaDIC for all river sites     
for i = 1:length(all_sites_names)
    DIC = all_ocean_cc.(all_sites_names2{i}).DIC_y;
    ALK = all_ocean_cc.(all_sites_names2{i}).TALK;
    temperature = all_ocean_cc.(all_sites_names2{i}).CTDTEMP_ITS90;
    salinity = all_ocean_cc.(all_sites_names2{i}).recommended_Salinity_PSS78;
    pressure = all_ocean_cc.(all_sites_names2{i}).CTDPRES;
    SIL = 0;
    PO4 = 0;
    sens_ocean.(all_sites_names2{i}) = derivnum ('par2',ALK,DIC,1,2,salinity,temperature,temperature,pressure,pressure,... 
                                   SIL,PO4,0,0,...
                                   SCALE,K1K2,SO4,KF,BOR);
end

%% Calculate propagated errors

ePAR1 = 0.002; %0.2% following Jiang et al
ePAR2 = 0.001; %0.1% following Jiang et al
eSAL = 0.01;
eTEMP = 0.01;
eSI = 0;
ePO4 = 0;
eNH4 = 0;
eH2S = 0;

for i = 1:length(all_sites_names)
    DIC = all_ocean_cc.(all_sites_names2{i}).DIC_y;
    ALK = all_ocean_cc.(all_sites_names2{i}).TALK;
    temperature = all_ocean_cc.(all_sites_names2{i}).CTDTEMP_ITS90;
    salinity = all_ocean_cc.(all_sites_names2{i}).recommended_Salinity_PSS78;
    pressure = all_ocean_cc.(all_sites_names2{i}).CTDPRES;
    SIL = 0;
    PO4 = 0;
    errors_ocean.(all_sites_names2{i}) = errors(ALK,DIC,1,2,salinity,temperature,temperature,pressure,pressure,SIL,PO4,0,0,ePAR1*nanmean(ALK),ePAR2*nanmean(DIC),eSAL,eTEMP,0,0,0,0,0,0,0,SCALE,K1K2,SO4,KF,BOR);
 
end


%% Subset ocean data to match Alkorg corrected data from neaa_carbcalc_alkorg_corrected_231102.m (140 stations)

for i = 1:length(names_allcorrect)
    sens_ocean_org.(names_allcorrect{i}) = sens_ocean.(names_allcorrect{i});
    errors_ocean_org.(names_allcorrect{i}) = errors_ocean.(names_allcorrect{i});
end

