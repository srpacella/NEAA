%% 9/1/2023 Updates by SRP

% 1. Compare data before and after 1990 Clean Air Act Amendments (Kaushal et
% al., 2013).
% 2. Changes means to medians for summarizing data within a site
% 3. Calculates discharge-normalized medians rather than overall mean
% amongst sites
% 4. Uses H+ instead of pH (code 33 in CO2SYS)
% 5. Only uses Alkorg corrected data
% 6. Adds in changes in aragonite saturation state, inlcuding corrections
% for stream Ca levels using calcium_omega_correct.m

%% Calculate median ocean end-member for each site
for i = 1:length(names_allcorrect)
    ocean_all.alkalinity(i) = nanmedian(all_ocean_cc.(names_allcorrect{i}).TALK);
    ocean_all.dic(i) = nanmedian(all_ocean_cc.(names_allcorrect{i}).DIC_y);
    ocean_all.salinity(i) = nanmedian(all_ocean_cc.(names_allcorrect{i}).S);
    ocean_all.temperature(i) = nanmedian(all_ocean_cc.(names_allcorrect{i}).T);
end


%% Pull data for "baseline" and "current" datasets
% Data from 1970 - 1989
river_baseline_cc = struct();
river_baseline_sens = struct();
for i = 1:length(names_allcorrect)
    index_temp = find(all_data_sites.(names_allcorrect{i}).ActivityStartDate > "1970-01-01" & all_data_sites.(names_allcorrect{i}).ActivityStartDate < "1989-12-31");%find indices of data ranging from 1980-1990
    if length(index_temp) < 10 %if fewer than 10 observations, continue with loop
        continue
    end
   cc_temp = site_cc_org.(names_allcorrect{i})(index_temp,:);%find calculated carbon system for this range
   sens_temp =  sens_sites_org.(names_allcorrect{i})(index_temp,:);%find calculated sensitivites for this range
   river_baseline_cc.(names_allcorrect{i})=cc_temp;% store carbonate system in structure
   river_baseline_sens.(names_allcorrect{i})=sens_temp;% store sensitivities in structure
end
names_baseline = fieldnames(river_baseline_cc);

%Calculate current river end-member pH from 1990-Present
river_current_cc = struct();
river_current_sens = struct();
for i = 1:length(names_baseline)
    index_temp = find(all_data_sites.(names_baseline{i}).ActivityStartDate > "1990-01-01" & all_data_sites.(names_baseline{i}).ActivityStartDate < "2023-01-01");%find indices of data ranging from 1980-1990
    if length(index_temp) < 10 %if fewer than 10 observations, continue with loop
        continue
    end
   cc_temp = site_cc_org.(names_baseline{i})(index_temp,:);%find calculated carbon system for this range
   sens_temp =  sens_sites_org.(names_baseline{i})(index_temp,:);%find calculated sensitivites for this range
   river_current_cc.(names_baseline{i})=cc_temp;% store carbonate system in structure
   river_current_sens.(names_baseline{i})=sens_temp;% store sensitivities in structure
end
names_current = fieldnames(river_current_cc);

%% Calculate mixing curves with baseline and current datasets. Quantify changes across salinity spectrum.

% Define constants for carbonate system calculations
SCALE  = 4; % NBS scale
K1K2_estuary   = 14; % 14 = Millero, 2010 T:    0-50  S:  1-50. Seaw. scale. Real seawater. 
SO4    = 1; % Dickson (1990) KSO4
KF     = 2; % Perez & Fraga (1987) KF
BOR    = 2; % Lee et al (2010) TB
pressure = 0;

river_baseline_cc_match = struct();
river_baseline_sens_match = struct();

%Subset baseline data to match sites with current data
for i = 1:length(names_current)
    idx = find(strcmp(names_baseline, names_current(i)));
    if idx > 0
        river_baseline_cc_match.(names_current{i}) = river_baseline_cc.(names_current{i});
        river_baseline_sens_match.(names_current{i}) = river_baseline_sens.(names_current{i});
    else
        continue
    end
end

%Calculate baseline mixing curves with present day ocean chemistry
for i = 1:length(names_current)
        %Calulate site-specific mixing line
        f_m = [0:0.02:1]; %ocean mixing fraction
        riv_dic = nanmedian(river_baseline_cc_match.(names_current{i})(:,2));
        riv_alk = nanmedian(river_baseline_cc_match.(names_current{i})(:,1));
        riv_salinity = nanmedian(river_baseline_cc_match.(names_current{i})(:,58));
        riv_temperature = nanmedian(river_baseline_cc_match.(names_current{i})(:,48));
        ocean_dic = nanmedian(all_ocean_cc.(names_current{i}).DIC_y);
        ocean_alk = nanmedian(all_ocean_cc.(names_current{i}).TALK);
        ocean_salinity = nanmedian(all_ocean_cc.(names_current{i}).recommended_Salinity_PSS78);
        ocean_temperature = nanmedian(all_ocean_cc.(names_current{i}).CTDTEMP_ITS90);
        salinity_mix = (f_m.*ocean_salinity) + ((1-f_m).*riv_salinity);
        dic_mix = (f_m.*ocean_dic) + ((1-f_m).*riv_dic);
        alk_mix = (f_m.*ocean_alk) + ((1-f_m).*riv_alk);
        temperature_mix = (f_m.*ocean_temperature) + ((1-f_m).*riv_temperature);
        baseline_estuary_cc.(names_current{i}) = CO2SYS(alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,temperature_mix,pressure,pressure,SIL,PO4,0,0,SCALE,K1K2_estuary,...
              SO4,KF,BOR);
        baseline_estuary_sens.(names_current{i}) = derivnum ('par2',alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,temperature_mix,pressure,pressure,... 
                                   SIL,PO4,0,0,...
                                   SCALE,K1K2_estuary,SO4,KF,BOR);
end
% Correct omega values for stream calcium concentrations
for i = 1:length(names_current)
    calcium_median = nanmedian(all_calcium.(names_current{i}).ResultMeasureValue); %mg/L
    calcium_correct = calcium_median ./ 1000 ./ 40.078; %converts to mol/kg
    Ca_fw = calcium_correct;
    S_sw = 35; % ocean end-member salinity
    Ca_sw = 0.010285; % ocean end-member calcium concentration in mol/kg
    S = baseline_estuary_cc.(names_current{i})(:,58);
    Ca_correct_beckwith = (Ca_fw./Ca_sw).*((S_sw./S)-1);
    baseline_estuary_cc.(names_current{i})(:,100) = baseline_estuary_cc.(names_current{i})(:,36) + (baseline_estuary_cc.(names_current{i})(:,36).*Ca_correct_beckwith);
    baseline_estuary_cc.(names_current{i})(:,101) = baseline_estuary_cc.(names_current{i})(:,35) + (baseline_estuary_cc.(names_current{i})(:,35).*Ca_correct_beckwith);
    if isnan(baseline_estuary_cc.(names_current{i})(1,100)) == 1
        baseline_estuary_cc.(names_current{i})(:,100) = baseline_estuary_cc.(names_current{i})(:,36)
        baseline_estuary_cc.(names_current{i})(:,101) = baseline_estuary_cc.(names_current{i})(:,35)
    end
end
%Calculate current mixing curves with present day ocean chemistry
for i = 1:length(names_current)
        %Calulate site-specific mixing line
        f_m = [0:0.02:1]; %ocean mixing fraction
        riv_dic = nanmedian(river_current_cc.(names_current{i})(:,2));
        riv_alk = nanmedian(river_current_cc.(names_current{i})(:,1));
        riv_salinity = nanmedian(river_current_cc.(names_current{i})(:,58));
        riv_temperature = nanmedian(river_current_cc.(names_current{i})(:,48));
        ocean_dic = nanmedian(all_ocean_cc.(names_current{i}).DIC_y);
        ocean_alk = nanmedian(all_ocean_cc.(names_current{i}).TALK);
        ocean_salinity = nanmedian(all_ocean_cc.(names_current{i}).recommended_Salinity_PSS78);
        ocean_temperature = nanmedian(all_ocean_cc.(names_current{i}).CTDTEMP_ITS90);
        salinity_mix = (f_m.*ocean_salinity) + ((1-f_m).*riv_salinity);
        dic_mix = (f_m.*ocean_dic) + ((1-f_m).*riv_dic);
        alk_mix = (f_m.*ocean_alk) + ((1-f_m).*riv_alk);
        temperature_mix = (f_m.*ocean_temperature) + ((1-f_m).*riv_temperature);
        current_estuary_cc.(names_current{i}) = CO2SYS(alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,temperature_mix,pressure,pressure,SIL,PO4,0,0,SCALE,K1K2_estuary,...
              SO4,KF,BOR);
        current_estuary_sens.(names_current{i}) = derivnum ('par2',alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,temperature_mix,pressure,pressure,... 
                                   SIL,PO4,0,0,...
                                   SCALE,K1K2_estuary,SO4,KF,BOR);
end

% Correct omega values for stream calcium concentrations
for i = 1:length(names_current)
    calcium_median = nanmedian(all_calcium.(names_current{i}).ResultMeasureValue); %mg/L
    calcium_correct = calcium_median ./ 1000 ./ 40.078; %converts to mol/kg
    %[omega_arag_correct omega_calc_correct] = calcium_omega_correct(current_estuary_cc.(names_current{i})(:,36),current_estuary_cc.(names_current{i})(:,35),calcium_correct,current_estuary_cc.(names_current{i})(58));
    Ca_fw = calcium_correct;
    S_sw = 35; % ocean end-member salinity
    Ca_sw = 0.010285; % ocean end-member calcium concentration in mol/kg
    S = current_estuary_cc.(names_current{i})(:,58);
    Ca_correct_beckwith = (Ca_fw./Ca_sw).*((S_sw./S)-1);
    current_estuary_cc.(names_current{i})(:,100) = current_estuary_cc.(names_current{i})(:,36) + (current_estuary_cc.(names_current{i})(:,36).*Ca_correct_beckwith);
    current_estuary_cc.(names_current{i})(:,101) = current_estuary_cc.(names_current{i})(:,35) + (current_estuary_cc.(names_current{i})(:,35).*Ca_correct_beckwith);
    if isnan(current_estuary_cc.(names_current{i})(1,100)) == 1
        current_estuary_cc.(names_current{i})(:,100) = current_estuary_cc.(names_current{i})(:,36)
        current_estuary_cc.(names_current{i})(:,101) = current_estuary_cc.(names_current{i})(:,35)
    end
end

%% Create stream discharge vector of sites used in pre/post 1990 analysis
discharge_means_names=regexprep(discharge_means(:,1),'-','');
discharge_change = []
count = 1
for i = 1:length(names_current)
    idx = find(strcmpi(discharge_means_names,names_current(i)));
    if idx > 0
        discharge_change(count) = discharge_means_numeric(idx);
    else
        continue
    end
    count = count+1;
end
discharge_weight = discharge_change./sum(discharge_change);
 %% Plots sites used in change analysis
 
% Map of all sites with sufficient data to look at change before/after 1990
% Make vector of 1=included and 2 = excluded
trend_category = []
for i = 1:length(names_allcorrect)
    if find(strcmpi(names_allcorrect(i),names_current)) > 0
        trend_category(i) = 1;
    else
        trend_category(i) = 0;
    end
end
trend_category = trend_category';
trend_category = categorical(trend_category);
figure
gb = geobubble(lat_orgstreams,lon_orgstreams,discharge_orgstreams*conv_cfs,'MapLayout','maximized','colordata',trend_category);
gb.SizeLegendTitle = 'Mean discharge (m^3/s)';
%geobasemap colorterrain
geobasemap landcover
%%
% %%%%%%%%%%%%%% Mixing diagrams

delta_estuary_H = [];
delta_estuary_pH = [];
delta_estuary_pCO2 = [];
delta_estuary_sensH = [];
delta_estuary_omega = [];
% Plot changes in estuary pH for each site
for i = 1:length(names_current)
    delta_estuary_H(:,i) = current_estuary_cc.(names_current{i})(:,33)-baseline_estuary_cc.(names_current{i})(:,33);
    delta_estuary_pH(:,i) = current_estuary_cc.(names_current{i})(:,43)-baseline_estuary_cc.(names_current{i})(:,43);
    delta_estuary_pCO2(:,i) = current_estuary_cc.(names_current{i})(:,22)-baseline_estuary_cc.(names_current{i})(:,22);
    delta_estuary_sensH(:,i) = current_estuary_sens.(names_current{i})(:,3)-baseline_estuary_sens.(names_current{i})(:,3);
    delta_estuary_omega(:,i) = current_estuary_cc.(names_current{i})(:,100)-baseline_estuary_cc.(names_current{i})(:,100);
    
end
%% Calculate changes in mixing from OA effects to ocean end-member
%Find atmospheric CO2 values for both periods
load CO_RCP_DATA
baseline_co2 = mean(rcp_co2(206:226,2)); %1970-1990
current_co2 = mean(rcp_co2(226:259,2)); % 1990-2023
%current_co2 = rcp_co2(259,2); % 2023
ocean_pco2_delta = current_co2-baseline_co2; % Change in pCO2 between current and baseline periods
   
%% Calculate mixing curves with baseline and current datasets. Quantify changes across salinity spectrum.

% Define constants for carbonate system calculations
SCALE  = 4; % NBS scale
K1K2_ocean   = 10 % = Lueker et al, 2000 T:    2-35  S: 19-43. Total scale. Real seawater.
SO4    = 1; % Dickson (1990) KSO4
KF     = 2; % Perez & Fraga (1987) KF
BOR    = 2; % Lee et al (2010) TB
pressure = 0; % Assumes surface water samples

ocean_baseline_cc_match = struct();
ocean_baseline_sens_match = struct();

%Subset current data to match sites with current river data
for i = 1:length(names_current)
    ocean_current_cc_match.(names_current{i}) = all_ocean_cc.(names_current{i});
end

%Calculate past mixing curves with altered ocean chemistry
for i = 1:length(names_current)
        %Calulate site-specific mixing line
        f_m = [0:0.02:1]; %ocean mixing fraction
        riv_dic = nanmedian(river_current_cc.(names_current{i})(:,2));
        riv_alk = nanmedian(river_current_cc.(names_current{i})(:,1));
        riv_salinity = nanmedian(river_current_cc.(names_current{i})(:,58));
        riv_temperature = nanmedian(river_current_cc.(names_current{i})(:,48));
        ocean_pCO2 = current_estuary_cc.(names_current{i})(length(f_m),22)-ocean_pco2_delta;
        %ocean_pCO2 = nanmedian(all_ocean_cc.(names_current{i}).pCO2-ocean_pco2_delta);
        ocean_alk = nanmedian(all_ocean_cc.(names_current{i}).TALK);
        ocean_salinity = nanmedian(all_ocean_cc.(names_current{i}).recommended_Salinity_PSS78);
        ocean_temperature = nanmedian(all_ocean_cc.(names_current{i}).CTDTEMP_ITS90);
        ocean_baseline_cc.(names_current{i})=CO2SYS(ocean_alk,ocean_pCO2,1,4,ocean_salinity,ocean_temperature,ocean_temperature,pressure,pressure,SIL,PO4,0,0,SCALE,K1K2_ocean,...
              SO4,KF,BOR);
        ocean_dic = ocean_baseline_cc.(names_current{i})(2);
        salinity_mix = (f_m.*ocean_salinity) + ((1-f_m).*riv_salinity);
        dic_mix = (f_m.*ocean_dic) + ((1-f_m).*riv_dic);
        alk_mix = (f_m.*ocean_alk) + ((1-f_m).*riv_alk);
        temperature_mix = (f_m.*ocean_temperature) + ((1-f_m).*riv_temperature);
        oceanbaseline_estuary_cc.(names_current{i}) = CO2SYS(alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,temperature_mix,pressure,pressure,SIL,PO4,0,0,SCALE,K1K2_estuary,...
              SO4,KF,BOR);
        oceanbaseline_estuary_sens.(names_current{i}) = derivnum ('par2',alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,temperature_mix,pressure,pressure,... 
                                   SIL,PO4,0,0,...
                                   SCALE,K1K2_estuary,SO4,KF,BOR);
end

% Correct omega values for stream calcium concentrations
for i = 1:length(names_current)
    calcium_median = nanmedian(all_calcium.(names_current{i}).ResultMeasureValue); %mg/L
    calcium_correct = calcium_median ./ 1000 ./ 40.078; %converts to mol/kg
    Ca_fw = calcium_correct;
    S_sw = 35; % ocean end-member salinity
    Ca_sw = 0.010285; % ocean end-member calcium concentration in mol/kg
    S = oceanbaseline_estuary_cc.(names_current{i})(:,58);
    Ca_correct_beckwith = (Ca_fw./Ca_sw).*((S_sw./S)-1);
    oceanbaseline_estuary_cc.(names_current{i})(:,100) = oceanbaseline_estuary_cc.(names_current{i})(:,36) + (oceanbaseline_estuary_cc.(names_current{i})(:,36).*Ca_correct_beckwith);
    oceanbaseline_estuary_cc.(names_current{i})(:,101) = oceanbaseline_estuary_cc.(names_current{i})(:,35) + (oceanbaseline_estuary_cc.(names_current{i})(:,35).*Ca_correct_beckwith);
    if isnan(oceanbaseline_estuary_cc.(names_current{i})(1,100)) == 1
        oceanbaseline_estuary_cc.(names_current{i})(:,100) = oceanbaseline_estuary_cc.(names_current{i})(:,36)
        oceanbaseline_estuary_cc.(names_current{i})(:,101) = oceanbaseline_estuary_cc.(names_current{i})(:,35)
    end
end


deltaocean_estuary_H = [];
deltaocean_estuary_pH = [];
deltaocean_estuary_pCO2 = [];
deltaocean_estuary_sensH = [];
deltaocean_estuary_omega = [];
% Plot changes in estuary pH for each site
for i = 1:length(names_current)
    deltaocean_estuary_H(:,i) = current_estuary_cc.(names_current{i})(:,33)-oceanbaseline_estuary_cc.(names_current{i})(:,33);
    deltaocean_estuary_pH(:,i) = current_estuary_cc.(names_current{i})(:,43)-oceanbaseline_estuary_cc.(names_current{i})(:,43);
    deltaocean_estuary_pCO2(:,i) = current_estuary_cc.(names_current{i})(:,22)-oceanbaseline_estuary_cc.(names_current{i})(:,22);
    deltaocean_estuary_sensH(:,i) = current_estuary_sens.(names_current{i})(:,3)-oceanbaseline_estuary_sens.(names_current{i})(:,3);
    deltaocean_estuary_omega(:,i) = current_estuary_cc.(names_current{i})(:,100)-oceanbaseline_estuary_cc.(names_current{i})(:,100);
end

%% Calculate mixing curves with past river & past ocean chemistry
for i = 1:length(names_current)
        %Calulate site-specific mixing line
        f_m = [0:0.02:1]; %ocean mixing fraction
        riv_dic = nanmedian(river_baseline_cc_match.(names_current{i})(:,2));
        riv_alk = nanmedian(river_baseline_cc_match.(names_current{i})(:,1));
        riv_salinity = nanmedian(river_baseline_cc_match.(names_current{i})(:,58));
        riv_temperature = nanmedian(river_baseline_cc_match.(names_current{i})(:,48));
        ocean_dic = nanmedian(ocean_baseline_cc.(names_current{i})(2));
        ocean_alk = nanmedian(all_ocean_cc.(names_current{i}).TALK);
        ocean_salinity = nanmedian(all_ocean_cc.(names_current{i}).recommended_Salinity_PSS78);
        ocean_temperature = nanmedian(all_ocean_cc.(names_current{i}).CTDTEMP_ITS90);
        salinity_mix = (f_m.*ocean_salinity) + ((1-f_m).*riv_salinity);
        dic_mix = (f_m.*ocean_dic) + ((1-f_m).*riv_dic);
        alk_mix = (f_m.*ocean_alk) + ((1-f_m).*riv_alk);
        temperature_mix = (f_m.*ocean_temperature) + ((1-f_m).*riv_temperature);
        allpast_estuary_cc.(names_current{i}) = CO2SYS(alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,temperature_mix,pressure,pressure,SIL,PO4,0,0,SCALE,K1K2_estuary,...
              SO4,KF,BOR);
        allpast_estuary_sens.(names_current{i}) = derivnum ('par2',alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,temperature_mix,pressure,pressure,... 
                                   SIL,PO4,0,0,...
                                   SCALE,K1K2_estuary,SO4,KF,BOR);
end
    
% Correct omega values for stream calcium concentrations
for i = 1:length(names_current)
    calcium_median = nanmedian(all_calcium.(names_current{i}).ResultMeasureValue); %mg/L
    calcium_correct = calcium_median ./ 1000 ./ 40.078; %converts to mol/kg
    Ca_fw = calcium_correct;
    S_sw = 35; % ocean end-member salinity
    Ca_sw = 0.010285; % ocean end-member calcium concentration in mol/kg
    S = allpast_estuary_cc.(names_current{i})(:,58);
    Ca_correct_beckwith = (Ca_fw./Ca_sw).*((S_sw./S)-1);
    allpast_estuary_cc.(names_current{i})(:,100) = allpast_estuary_cc.(names_current{i})(:,36) + (allpast_estuary_cc.(names_current{i})(:,36).*Ca_correct_beckwith);
    allpast_estuary_cc.(names_current{i})(:,101) = allpast_estuary_cc.(names_current{i})(:,35) + (allpast_estuary_cc.(names_current{i})(:,35).*Ca_correct_beckwith);
    if isnan(allpast_estuary_cc.(names_current{i})(1,100)) == 1
        allpast_estuary_cc.(names_current{i})(:,100) = allpast_estuary_cc.(names_current{i})(:,36)
        allpast_estuary_cc.(names_current{i})(:,101) = allpast_estuary_cc.(names_current{i})(:,35)
    end
end   
deltaallpast_estuary_H = [];
deltaallpast_estuary_pH = [];
deltaallpast_estuary_pCO2 = [];
deltaallpast_estuary_sensH = [];
deltaallpast_estuary_omega = [];
deltaallpast_estuary_omega_noca = [];
% Plot changes in estuary pH for each site
for i = 1:length(names_current)
    deltaallpast_estuary_H(:,i) = current_estuary_cc.(names_current{i})(:,33)-allpast_estuary_cc.(names_current{i})(:,33);
    deltaallpast_estuary_pH(:,i) = current_estuary_cc.(names_current{i})(:,43)-allpast_estuary_cc.(names_current{i})(:,43);
    deltaallpast_estuary_pCO2(:,i) = current_estuary_cc.(names_current{i})(:,22)-allpast_estuary_cc.(names_current{i})(:,22);
    deltaallpast_estuary_sensH(:,i) = current_estuary_sens.(names_current{i})(:,3)-allpast_estuary_sens.(names_current{i})(:,3);
    deltaallpast_estuary_omega(:,i) = current_estuary_cc.(names_current{i})(:,100)-allpast_estuary_cc.(names_current{i})(:,100); 
    deltaallpast_estuary_omega_noca(:,i) = current_estuary_cc.(names_current{i})(:,36)-allpast_estuary_cc.(names_current{i})(:,36); 
end

%% Calculate mixing curves with past river & past ocean chemistry while holding temperature at current levels
for i = 1:length(names_current)
        %Calulate site-specific mixing line
        f_m = [0:0.02:1]; %ocean mixing fraction
        riv_dic = nanmedian(river_baseline_cc_match.(names_current{i})(:,2));
        riv_alk = nanmedian(river_baseline_cc_match.(names_current{i})(:,1));
        riv_salinity = nanmedian(river_baseline_cc_match.(names_current{i})(:,58));
        riv_temperature = nanmedian(river_baseline_cc_match.(names_current{i})(:,48));
        ocean_dic = nanmedian(ocean_baseline_cc.(names_current{i})(2));
        ocean_alk = nanmedian(all_ocean_cc.(names_current{i}).TALK);
        ocean_salinity = nanmedian(all_ocean_cc.(names_current{i}).recommended_Salinity_PSS78);
        ocean_temperature = nanmedian(all_ocean_cc.(names_current{i}).CTDTEMP_ITS90);
        salinity_mix = (f_m.*ocean_salinity) + ((1-f_m).*riv_salinity);
        dic_mix = (f_m.*ocean_dic) + ((1-f_m).*riv_dic);
        alk_mix = (f_m.*ocean_alk) + ((1-f_m).*riv_alk);
        temperature_mix = (f_m.*ocean_temperature) + ((1-f_m).*riv_temperature);
        allpast_estuary_cc_temp.(names_current{i}) = CO2SYS(alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,current_estuary_cc.(names_current{i})(:,48),pressure,pressure,SIL,PO4,0,0,SCALE,K1K2_estuary,...
              SO4,KF,BOR);
        allpast_estuary_sens_temp.(names_current{i}) = derivnum ('par2',alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,current_estuary_cc.(names_current{i})(:,48),pressure,pressure,... 
                                   SIL,PO4,0,0,...
                                   SCALE,K1K2_estuary,SO4,KF,BOR);
end
    
% Correct omega values for stream calcium concentrations
for i = 1:length(names_current)
    calcium_median = nanmedian(all_calcium.(names_current{i}).ResultMeasureValue); %mg/L
    calcium_correct = calcium_median ./ 1000 ./ 40.078; %converts to mol/kg
    Ca_fw = calcium_correct;
    S_sw = 35; % ocean end-member salinity
    Ca_sw = 0.010285; % ocean end-member calcium concentration in mol/kg
    S = allpast_estuary_cc_temp.(names_current{i})(:,58);
    Ca_correct_beckwith = (Ca_fw./Ca_sw).*((S_sw./S)-1);
    allpast_estuary_cc_temp.(names_current{i})(:,100) = allpast_estuary_cc_temp.(names_current{i})(:,36) + (allpast_estuary_cc_temp.(names_current{i})(:,36).*Ca_correct_beckwith);
    allpast_estuary_cc_temp.(names_current{i})(:,101) = allpast_estuary_cc_temp.(names_current{i})(:,35) + (allpast_estuary_cc_temp.(names_current{i})(:,35).*Ca_correct_beckwith);
    if isnan(allpast_estuary_cc_temp.(names_current{i})(1,100)) == 1
        allpast_estuary_cc_temp.(names_current{i})(:,100) = allpast_estuary_cc_temp.(names_current{i})(:,36)
        allpast_estuary_cc_temp.(names_current{i})(:,101) = allpast_estuary_cc_temp.(names_current{i})(:,35)
    end
end      
deltaallpast_estuary_H_temp = [];
deltaallpast_estuary_pH_temp = [];
deltaallpast_estuary_pCO2_temp = [];
deltaallpast_estuary_sensH_temp = [];
deltaallpast_estuary_omega_temp = [];
% Plot changes in estuary pH for each site
for i = 1:length(names_current)
    deltaallpast_estuary_H_temp(:,i) = current_estuary_cc.(names_current{i})(:,33)-allpast_estuary_cc_temp.(names_current{i})(:,33);
    deltaallpast_estuary_pH_temp(:,i) = current_estuary_cc.(names_current{i})(:,43)-allpast_estuary_cc_temp.(names_current{i})(:,43);
    deltaallpast_estuary_pCO2_temp(:,i) = current_estuary_cc.(names_current{i})(:,22)-allpast_estuary_cc_temp.(names_current{i})(:,22);
    deltaallpast_estuary_sensH_temp(:,i) = current_estuary_sens.(names_current{i})(:,3)-allpast_estuary_sens_temp.(names_current{i})(:,3);
    deltaallpast_estuary_omega_temp(:,i) = current_estuary_cc.(names_current{i})(:,100)-allpast_estuary_cc_temp.(names_current{i})(:,100); 
end
%% Multi-panel figure with pco2, H+, and H+ sensitivity
figure
subplot(3,3,1)
X = f_m;
Y_noisy_pCO2_deltas = deltaocean_estuary_pCO2';
plot_distribution_prctile(X,Y_noisy_pCO2_deltas,'Prctile',[25 50 75 90]);
ylabel('Ocean \Delta{\itp}CO_2 (\muatm)')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,2)
Y_noisy_pCO2_deltas = delta_estuary_pCO2';
plot_distribution_prctile(X,Y_noisy_pCO2_deltas,'Prctile',[25 50 75 90]);
ylabel('Stream \Delta{\itp}CO_2 (\muatm)')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,3)
Y_noisy_pCO2_deltas = deltaallpast_estuary_pCO2';
plot_distribution_prctile(X,Y_noisy_pCO2_deltas,'Prctile',[25 50 75 90]);
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
ylabel('Combined \Delta{\itp}CO_2 (\muatm)')


subplot(3,3,4)
X = f_m;
Y_noisy_H_deltas = deltaocean_estuary_H';
plot_distribution_prctile(X,Y_noisy_H_deltas,'Prctile',[25 50 75 90]);
ylabel('Ocean \Delta[H^+]')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,5)
Y_noisy_H_deltas = delta_estuary_H';
plot_distribution_prctile(X,Y_noisy_H_deltas,'Prctile',[25 50 75 90]);
ylabel('Stream \Delta[H^+]')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,6)
Y_noisy_H_deltas = deltaallpast_estuary_H';
plot_distribution_prctile(X,Y_noisy_H_deltas,'Prctile',[25 50 75 90]);
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
ylabel('Combined \Delta[H^+]')


subplot(3,3,7)
X = f_m;
Y_noisy_sensH_deltas = deltaocean_estuary_sensH';
plot_distribution_prctile(X,Y_noisy_sensH_deltas,'Prctile',[25 50 75 90]);
ylabel('Ocean \Delta \DeltaH^+/\DeltaDIC')
xlabel('Fraction ocean water')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,8)
Y_noisy_sensH_deltas = delta_estuary_sensH';
plot_distribution_prctile(X,Y_noisy_sensH_deltas,'Prctile',[25 50 75 90]);
ylabel('Stream \Delta \DeltaH^+/\DeltaDIC')
xlabel('Fraction ocean water')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,9)
Y_noisy_sensH_deltas = deltaallpast_estuary_sensH';
plot_distribution_prctile(X,Y_noisy_sensH_deltas,'Prctile',[25 50 75 90]);
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
ylabel('Combined \Delta \DeltaH^+/\DeltaDIC')
xlabel('Fraction ocean water')

%% Multi-panel figure with pco2, ph, and H+ sensitivity
figure
subplot(3,3,1)
X = f_m;
Y_noisy_pCO2_deltas = deltaocean_estuary_pCO2';
plot_distribution_prctile(X,Y_noisy_pCO2_deltas,'Prctile',[25 50 75 90]);
ylabel('Ocean \Delta{\itp}CO_2 (\muatm)')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,2)
Y_noisy_pCO2_deltas = delta_estuary_pCO2';
plot_distribution_prctile(X,Y_noisy_pCO2_deltas,'Prctile',[25 50 75 90]);
ylabel('Stream \Delta{\itp}CO_2 (\muatm)')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,3)
Y_noisy_pCO2_deltas = deltaallpast_estuary_pCO2';
plot_distribution_prctile(X,Y_noisy_pCO2_deltas,'Prctile',[25 50 75 90]);
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
ylabel('Combined \Delta{\itp}CO_2 (\muatm)')


subplot(3,3,4)
X = f_m;
Y_noisy_pH_deltas = deltaocean_estuary_pH';
plot_distribution_prctile(X,Y_noisy_pH_deltas,'Prctile',[25 50 75 90]);
ylabel('Ocean \DeltapH_T')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,5)
Y_noisy_pH_deltas = delta_estuary_pH';
plot_distribution_prctile(X,Y_noisy_pH_deltas,'Prctile',[25 50 75 90]);
ylabel('Stream \DeltapH_T')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,6)
Y_noisy_pH_deltas = deltaallpast_estuary_pH';
plot_distribution_prctile(X,Y_noisy_pH_deltas,'Prctile',[25 50 75 90]);
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
ylabel('Combined \DeltapH_T')


subplot(3,3,7)
X = f_m;
Y_noisy_sensH_deltas = deltaocean_estuary_sensH';
plot_distribution_prctile(X,Y_noisy_sensH_deltas,'Prctile',[25 50 75 90]);
ylabel('Ocean \Delta \DeltaH^+/\DeltaDIC')
xlabel('Fraction ocean water')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,8)
Y_noisy_sensH_deltas = delta_estuary_sensH';
plot_distribution_prctile(X,Y_noisy_sensH_deltas,'Prctile',[25 50 75 90]);
ylabel('Stream \Delta \DeltaH^+/\DeltaDIC')
xlabel('Fraction ocean water')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,9)
Y_noisy_sensH_deltas = deltaallpast_estuary_sensH';
plot_distribution_prctile(X,Y_noisy_sensH_deltas,'Prctile',[25 50 75 90]);
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
ylabel('Combined \Delta \DeltaH^+/\DeltaDIC')
xlabel('Fraction ocean water')

%% Multi-panel figure with pco2, ph, and H+ sensitivity AND flow-weighted average change
figure
subplot(3,3,1)
X = f_m;
Y_noisy_pCO2_deltas = deltaocean_estuary_pCO2';
plot_distribution_prctile(X,Y_noisy_pCO2_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaocean_estuary_pCO2
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Ocean \Delta{\itp}CO_2 (\muatm)')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,2)
Y_noisy_pCO2_deltas = delta_estuary_pCO2';
plot_distribution_prctile(X,Y_noisy_pCO2_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*delta_estuary_pCO2
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Stream \Delta{\itp}CO_2 (\muatm)')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,3)
Y_noisy_pCO2_deltas = deltaallpast_estuary_pCO2';
plot_distribution_prctile(X,Y_noisy_pCO2_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaallpast_estuary_pCO2
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
ylabel('Combined \Delta{\itp}CO_2 (\muatm)')


subplot(3,3,4)
X = f_m;
Y_noisy_pH_deltas = deltaocean_estuary_pH';
plot_distribution_prctile(X,Y_noisy_pH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaocean_estuary_pH
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Ocean \DeltapH_T')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,5)
Y_noisy_pH_deltas = delta_estuary_pH';
plot_distribution_prctile(X,Y_noisy_pH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*delta_estuary_pH
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Stream \DeltapH_T')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,6)
Y_noisy_pH_deltas = deltaallpast_estuary_pH';
plot_distribution_prctile(X,Y_noisy_pH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaallpast_estuary_pH
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
ylabel('Combined \DeltapH_T')


subplot(3,3,7)
X = f_m;
Y_noisy_sensH_deltas = deltaocean_estuary_sensH';
plot_distribution_prctile(X,Y_noisy_sensH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaocean_estuary_sensH
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Ocean \Delta \partialH^+/\partialDIC')
xlabel('Fraction ocean water')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,8)
Y_noisy_sensH_deltas = delta_estuary_sensH';
plot_distribution_prctile(X,Y_noisy_sensH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*delta_estuary_sensH
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Stream \Delta \partialH^+/\partialDIC')
xlabel('Fraction ocean water')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,9)
Y_noisy_sensH_deltas = deltaallpast_estuary_sensH';
plot_distribution_prctile(X,Y_noisy_sensH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaallpast_estuary_sensH
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
ylabel('Combined \Delta \partialH^+/\partialDIC')
xlabel('Fraction ocean water')

%% Multi-panel figure with pco2, ph, and H+ sensitivity AND flow-weighted average change and omega
figure
subplot(4,3,1)
X = f_m;
Y_noisy_pCO2_deltas = deltaocean_estuary_pCO2';
plot_distribution_prctile(X,Y_noisy_pCO2_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaocean_estuary_pCO2
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Ocean \Delta{\itp}CO_2 (\muatm)')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(4,3,2)
Y_noisy_pCO2_deltas = delta_estuary_pCO2';
plot_distribution_prctile(X,Y_noisy_pCO2_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*delta_estuary_pCO2
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Stream \Delta{\itp}CO_2 (\muatm)')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(4,3,3)
Y_noisy_pCO2_deltas = deltaallpast_estuary_pCO2';
plot_distribution_prctile(X,Y_noisy_pCO2_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaallpast_estuary_pCO2
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
ylabel('Combined \Delta{\itp}CO_2 (\muatm)')


subplot(4,3,4)
X = f_m;
Y_noisy_pH_deltas = deltaocean_estuary_pH';
plot_distribution_prctile(X,Y_noisy_pH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaocean_estuary_pH
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Ocean \DeltapH_T')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(4,3,5)
Y_noisy_pH_deltas = delta_estuary_pH';
plot_distribution_prctile(X,Y_noisy_pH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*delta_estuary_pH
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Stream \DeltapH_T')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(4,3,6)
Y_noisy_pH_deltas = deltaallpast_estuary_pH';
plot_distribution_prctile(X,Y_noisy_pH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaallpast_estuary_pH
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
ylabel('Combined \DeltapH_T')


subplot(4,3,7)
X = f_m;
Y_noisy_sensH_deltas = deltaocean_estuary_sensH';
plot_distribution_prctile(X,Y_noisy_sensH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaocean_estuary_sensH
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Ocean \Delta \partialH^+/\partialDIC')
xlabel('Fraction ocean water')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(4,3,8)
Y_noisy_sensH_deltas = delta_estuary_sensH';
plot_distribution_prctile(X,Y_noisy_sensH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*delta_estuary_sensH
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Stream \Delta \partialH^+/\partialDIC')
xlabel('Fraction ocean water')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(4,3,9)
Y_noisy_sensH_deltas = deltaallpast_estuary_sensH';
plot_distribution_prctile(X,Y_noisy_sensH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaallpast_estuary_sensH
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
ylabel('Combined \Delta \partialH^+/\partialDIC')
xlabel('Fraction ocean water')
% Aragonite saturation state plots
subplot(4,3,10)
X = f_m;
Y_noisy_omega_deltas = deltaocean_estuary_omega';
plot_distribution_prctile(X,Y_noisy_omega_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaocean_estuary_omega
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Ocean \Delta\Omega_a_r_a_g')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(4,3,11)
Y_noisy_omega_deltas = delta_estuary_omega';
plot_distribution_prctile(X,Y_noisy_omega_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*delta_estuary_omega
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Stream \Delta\Omega_a_r_a_g')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(4,3,12)
Y_noisy_omega_deltas = deltaallpast_estuary_omega';
plot_distribution_prctile(X,Y_noisy_omega_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaallpast_estuary_omega
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
ylabel('Combined \Delta\Omega_a_r_a_g')

%% Multi-panel figure with pco2, ph, and H+ sensitivity AND flow-weighted average change and omega
% NO TEMPERATURE CHANGE
figure
title('Constant temperature')
subplot(4,3,1)
X = f_m;
Y_noisy_pCO2_deltas = deltaocean_estuary_pCO2';
plot_distribution_prctile(X,Y_noisy_pCO2_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaocean_estuary_pCO2
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Ocean \Delta{\itp}CO_2 (\muatm)')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(4,3,2)
Y_noisy_pCO2_deltas = delta_estuary_pCO2';
plot_distribution_prctile(X,Y_noisy_pCO2_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*delta_estuary_pCO2
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Stream \Delta{\itp}CO_2 (\muatm)')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(4,3,3)
Y_noisy_pCO2_deltas = deltaallpast_estuary_pCO2_temp';
plot_distribution_prctile(X,Y_noisy_pCO2_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaallpast_estuary_pCO2_temp
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
ylabel('Combined \Delta{\itp}CO_2 (\muatm)')


subplot(4,3,4)
X = f_m;
Y_noisy_pH_deltas = deltaocean_estuary_pH';
plot_distribution_prctile(X,Y_noisy_pH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaocean_estuary_pH
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Ocean \DeltapH_T')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(4,3,5)
Y_noisy_pH_deltas = delta_estuary_pH';
plot_distribution_prctile(X,Y_noisy_pH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*delta_estuary_pH
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Stream \DeltapH_T')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(4,3,6)
Y_noisy_pH_deltas = deltaallpast_estuary_pH_temp';
plot_distribution_prctile(X,Y_noisy_pH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaallpast_estuary_pH_temp
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
ylabel('Combined \DeltapH_T')


subplot(4,3,7)
X = f_m;
Y_noisy_sensH_deltas = deltaocean_estuary_sensH';
plot_distribution_prctile(X,Y_noisy_sensH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaocean_estuary_sensH
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Ocean \Delta \partialH^+/\partialDIC')
xlabel('Fraction ocean water')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(4,3,8)
Y_noisy_sensH_deltas = delta_estuary_sensH';
plot_distribution_prctile(X,Y_noisy_sensH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*delta_estuary_sensH
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Stream \Delta \partialH^+/\partialDIC')
xlabel('Fraction ocean water')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(4,3,9)
Y_noisy_sensH_deltas = deltaallpast_estuary_sensH_temp';
plot_distribution_prctile(X,Y_noisy_sensH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaallpast_estuary_sensH_temp
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
ylabel('Combined \Delta \partialH^+/\partialDIC')
xlabel('Fraction ocean water')
% Aragonite saturation state plots
subplot(4,3,10)
X = f_m;
Y_noisy_omega_deltas = deltaocean_estuary_omega';
plot_distribution_prctile(X,Y_noisy_omega_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaocean_estuary_omega
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Ocean \Delta\Omega_a_r_a_g')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(4,3,11)
Y_noisy_omega_deltas = delta_estuary_omega';
plot_distribution_prctile(X,Y_noisy_omega_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*delta_estuary_omega
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Stream \Delta\Omega_a_r_a_g')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(4,3,12)
Y_noisy_omega_deltas = deltaallpast_estuary_omega_temp';
plot_distribution_prctile(X,Y_noisy_omega_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaallpast_estuary_omega_temp
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
ylabel('Combined \Delta\Omega_a_r_a_g')
title('Constant temp')

%% Multi-panel figure with pco2, ph, and H+ sensitivity AND flow-weighted average change and omega
% add line for constant temp
figure
subplot(4,3,1)
X = f_m;
Y_noisy_pCO2_deltas = deltaocean_estuary_pCO2';
plot_distribution_prctile(X,Y_noisy_pCO2_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaocean_estuary_pCO2
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Ocean \Delta{\itp}CO_2 (\muatm)')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
xticks([0:0.2:1])
subplot(4,3,2)
Y_noisy_pCO2_deltas = delta_estuary_pCO2';
plot_distribution_prctile(X,Y_noisy_pCO2_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*delta_estuary_pCO2
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Stream \Delta{\itp}CO_2 (\muatm)')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
xticks([0:0.2:1])
subplot(4,3,3)
Y_noisy_pCO2_deltas = deltaallpast_estuary_pCO2';
plot_distribution_prctile(X,Y_noisy_pCO2_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaallpast_estuary_pCO2
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
test4 = discharge_weight.*deltaallpast_estuary_pCO2_temp
test5 = sum(test4,2)
plot(f_m,test5,'.r','LineWidth',3)
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
ylabel('Combined \Delta{\itp}CO_2 (\muatm)')
xticks([0:0.2:1])


subplot(4,3,4)
X = f_m;
Y_noisy_pH_deltas = deltaocean_estuary_pH';
plot_distribution_prctile(X,Y_noisy_pH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaocean_estuary_pH
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Ocean \DeltapH_T')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
xticks([0:0.2:1])
subplot(4,3,5)
Y_noisy_pH_deltas = delta_estuary_pH';
plot_distribution_prctile(X,Y_noisy_pH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*delta_estuary_pH
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Stream \DeltapH_T')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
xticks([0:0.2:1])
subplot(4,3,6)
Y_noisy_pH_deltas = deltaallpast_estuary_pH';
plot_distribution_prctile(X,Y_noisy_pH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaallpast_estuary_pH
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
test4 = discharge_weight.*deltaallpast_estuary_pH_temp
test5 = sum(test4,2)
plot(f_m,test5,'.r','LineWidth',3)
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
ylabel('Combined \DeltapH_T')
xticks([0:0.2:1])

subplot(4,3,7)
X = f_m;
Y_noisy_sensH_deltas = deltaocean_estuary_sensH';
plot_distribution_prctile(X,Y_noisy_sensH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaocean_estuary_sensH
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Ocean \Delta \partialH^+/\partialDIC')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
xticks([0:0.2:1])
subplot(4,3,8)
Y_noisy_sensH_deltas = delta_estuary_sensH';
plot_distribution_prctile(X,Y_noisy_sensH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*delta_estuary_sensH
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Stream \Delta \partialH^+/\partialDIC')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
xticks([0:0.2:1])
subplot(4,3,9)
Y_noisy_sensH_deltas = deltaallpast_estuary_sensH';
plot_distribution_prctile(X,Y_noisy_sensH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaallpast_estuary_sensH
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
test4 = discharge_weight.*deltaallpast_estuary_sensH_temp
test5 = sum(test4,2)
plot(f_m,test5,'.r','LineWidth',3)
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
ylabel('Combined \Delta \partialH^+/\partialDIC')
xticks([0:0.2:1])
% Aragonite saturation state plots
subplot(4,3,10)
X = f_m;
Y_noisy_omega_deltas = deltaocean_estuary_omega';
plot_distribution_prctile(X,Y_noisy_omega_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaocean_estuary_omega
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Ocean \Delta\Omega_a_r_a_g')
xlabel('Fraction ocean water')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
xlabel('Fraction ocean water')
xticks([0:0.2:1])
subplot(4,3,11)
Y_noisy_omega_deltas = delta_estuary_omega';
plot_distribution_prctile(X,Y_noisy_omega_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*delta_estuary_omega
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
ylabel('Stream \Delta\Omega_a_r_a_g')
xlabel('Fraction ocean water')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
xlabel('Fraction ocean water')
xticks([0:0.2:1])
subplot(4,3,12)
Y_noisy_omega_deltas = deltaallpast_estuary_omega';
plot_distribution_prctile(X,Y_noisy_omega_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaallpast_estuary_omega
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
test4 = discharge_weight.*deltaallpast_estuary_omega_temp
test5 = sum(test4,2)
plot(f_m,test5,'.r','LineWidth',3)
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
ylabel('Combined \Delta\Omega_a_r_a_g')
xlabel('Fraction ocean water')
xticks([0:0.2:1])

%% Multi-panel figure with pco2, ph, and H+ sensitivity AND flow-weighted average change and omega
% add line for constant temp
figure
subplot(4,3,1)
X = f_m;
Y_noisy_pCO2_deltas = deltaocean_estuary_pCO2';
plot_distribution_prctile(X,Y_noisy_pCO2_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaocean_estuary_pCO2
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
xticks([0:0.2:1])
subplot(4,3,2)
Y_noisy_pCO2_deltas = delta_estuary_pCO2';
plot_distribution_prctile(X,Y_noisy_pCO2_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*delta_estuary_pCO2
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
xticks([0:0.2:1])
subplot(4,3,3)
Y_noisy_pCO2_deltas = deltaallpast_estuary_pCO2';
plot_distribution_prctile(X,Y_noisy_pCO2_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaallpast_estuary_pCO2
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
test4 = discharge_weight.*deltaallpast_estuary_pCO2_temp
test5 = sum(test4,2)
plot(f_m,test5,'.r','LineWidth',3)
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
xticks([0:0.2:1])


subplot(4,3,4)
X = f_m;
Y_noisy_pH_deltas = deltaocean_estuary_pH';
plot_distribution_prctile(X,Y_noisy_pH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaocean_estuary_pH
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
xticks([0:0.2:1])
subplot(4,3,5)
Y_noisy_pH_deltas = delta_estuary_pH';
plot_distribution_prctile(X,Y_noisy_pH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*delta_estuary_pH
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
xticks([0:0.2:1])
subplot(4,3,6)
Y_noisy_pH_deltas = deltaallpast_estuary_pH';
plot_distribution_prctile(X,Y_noisy_pH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaallpast_estuary_pH
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
test4 = discharge_weight.*deltaallpast_estuary_pH_temp
test5 = sum(test4,2)
plot(f_m,test5,'.r','LineWidth',3)
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
xticks([0:0.2:1])

subplot(4,3,7)
X = f_m;
Y_noisy_sensH_deltas = deltaocean_estuary_sensH';
plot_distribution_prctile(X,Y_noisy_sensH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaocean_estuary_sensH
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
xticks([0:0.2:1])
subplot(4,3,8)
Y_noisy_sensH_deltas = delta_estuary_sensH';
plot_distribution_prctile(X,Y_noisy_sensH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*delta_estuary_sensH
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
xticks([0:0.2:1])
subplot(4,3,9)
Y_noisy_sensH_deltas = deltaallpast_estuary_sensH';
plot_distribution_prctile(X,Y_noisy_sensH_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaallpast_estuary_sensH
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
test4 = discharge_weight.*deltaallpast_estuary_sensH_temp
test5 = sum(test4,2)
plot(f_m,test5,'.r','LineWidth',3)
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
xticks([0:0.2:1])
% Aragonite saturation state plots
subplot(4,3,10)
X = f_m;
Y_noisy_omega_deltas = deltaocean_estuary_omega';
plot_distribution_prctile(X,Y_noisy_omega_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaocean_estuary_omega
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
xlabel('Fraction ocean water')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
xlabel('Fraction ocean water')
xticks([0:0.2:1])
subplot(4,3,11)
Y_noisy_omega_deltas = delta_estuary_omega';
plot_distribution_prctile(X,Y_noisy_omega_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*delta_estuary_omega
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
xlabel('Fraction ocean water')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
xlabel('Fraction ocean water')
xticks([0:0.2:1])
subplot(4,3,12)
Y_noisy_omega_deltas = deltaallpast_estuary_omega';
plot_distribution_prctile(X,Y_noisy_omega_deltas,'Prctile',[25 50 75 90]);
hold on
test2 = discharge_weight.*deltaallpast_estuary_omega
test3 = sum(test2,2)
plot(f_m,test3,'-k','LineWidth',3)
test4 = discharge_weight.*deltaallpast_estuary_omega_temp
test5 = sum(test4,2)
plot(f_m,test5,'.r','LineWidth',3)
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
xlabel('Fraction ocean water')
xticks([0:0.2:1])

