% Calculate carbonate system of NEAA sites
% Requires output from neaa_import.m

% Does NOT correct for organic alkalinity
% Corrects omega outputs for river Ca using calcium_omega_correct


% Define constants for carbonate system calculations
SCALE  = 4; % NBS scale
K1K2   = 8; % Millero, 1979, FOR PURE WATER ONLY (i.e., Sal=0)    T:    0-50  S:     0. 
SO4    = 1; % Dickson (1990) KSO4
KF     = 2; % Perez & Fraga (1987) KF
BOR    = 2; % Lee et al (2010) TB


%% Calculate carbonate system for all river sites           
for i = 1:length(all_sites_names)
    pH = all_data_sites.(all_sites_names2{i}).pH_correct;
    alkalinity = all_data_sites.(all_sites_names2{i}).Alkalinity;
    temperature = all_data_sites.(all_sites_names2{i}).Temperature;
    %salinity = 0;
    salinity = all_data_sites.(all_sites_names2{i}).salinity;
    pressure = 0;
    SIL = 0;
    PO4 = 0;
    site_cc.(all_sites_names2{i}) = CO2SYS(alkalinity,pH,1,3,salinity,temperature,temperature,pressure,pressure,SIL,PO4,0,0,SCALE,K1K2,...
              SO4,KF,BOR);
end

%% Correct omega values for stream calcium concentrations
for i = 1:length(all_sites_names)
    calcium_median = nanmedian(all_calcium.(all_sites_names2{i}).ResultMeasureValue); %mg/L
    calcium_correct = calcium_median ./ 1000 ./ 40.078; %converts to mol/kg
    [arag_correct calc_correct] = calcium_omega_correct(site_cc.(all_sites_names2{i})(:,36),site_cc.(all_sites_names2{i})(:,35),calcium_correct,all_data_sites.(all_sites_names2{i}).salinity);
    site_cc.(all_sites_names2{i})(:,100) = arag_correct;
    site_cc.(all_sites_names2{i})(:,101) = calc_correct;
end

%% Calculate partial derivatives w/r/t deltaDIC for all river sites     
for i = 1:length(all_sites_names)
    DIC = site_cc.(all_sites_names2{i})(:,2);
    ALK = all_data_sites.(all_sites_names2{i}).Alkalinity;
    temperature = all_data_sites.(all_sites_names2{i}).Temperature;
    %salinity = 0;
    salinity = all_data_sites.(all_sites_names2{i}).salinity;
    pressure = 0;
    SIL = 0;
    PO4 = 0;
    sens_sites.(all_sites_names2{i}) = derivnum ('par2',ALK,DIC,1,2,salinity,temperature,temperature,pressure,pressure,... 
                                   SIL,PO4,0,0,...
                                   SCALE,K1K2,SO4,KF,BOR);
end
