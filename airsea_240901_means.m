%% 240901 Combined script to assess 1). changes in sea-air pCO2 differentials, 2). changes in associated DICdiseq

%% Calculate mean coastal stream sea-air pCO2 disequilibria for Pre- and Post-1990 (year-specific)

% Create structure of sea-air pCO2 differentials for all stream
% observations (with alk and pH corrections)
for i = 1:length(names_allcorrect)
    year_atm = year(all_data_sites.(names_allcorrect{i}).ActivityStartDate); %Extract vector of years
    for z = 1:length(year_atm)
        atm_co2(z) = rcp_co2(find(rcp_co2(:,1) == year_atm(z)),2);
    end
    airsea_observed.(names_allcorrect{i})(:,3) = site_cc_org.(names_allcorrect{i})(:,22) - atm_co2';
    airsea_observed.(names_allcorrect{i})(:,2) = atm_co2';
    airsea_observed.(names_allcorrect{i})(:,1) = year_atm;
    clear atm_co2 year_atm
end

% Use Wilcoxon test for significant changes sea-air CO2 at comparison sites
for i = 1:length(names_current)
    x_find = find(airsea_observed.(names_current{i})(:,1) < 1990);
    y_find = find(airsea_observed.(names_current{i})(:,1) >= 1990);
    x = airsea_observed.(names_current{i})(x_find,3);
    y = airsea_observed.(names_current{i})(y_find,3);
    [p,h,stats] = ranksum(x,y);
    wilcoxon_airsea(i,1) = p;
    wilcoxon_airsea(i,2) = h;
end

% Calculate means sea-air disequilibria for Pre-and Post-1990 at each
% site
for i = 1:length(names_current)
    pre_find = find(airsea_observed.(names_current{i})(:,1) < 1990);
    post_find = find(airsea_observed.(names_current{i})(:,1) >= 1990);
    airsea_means(i,1) = nanmean(airsea_observed.(names_current{i})(pre_find,3));
    airsea_means(i,2) = nanmean(airsea_observed.(names_current{i})(post_find,3));
    airsea_means(i,3) = airsea_means(i,2) - airsea_means(i,1);
end

% Calculated discharged-weighted mean of sea-air disequilibrium change
airsea_fwa_pre1990= sum(discharge_weight'.*airsea_means(:,1))
airsea_fwa_post1990= sum(discharge_weight'.*airsea_means(:,2))
airsea_fwa_diseq= sum(discharge_weight'.*airsea_means(:,3))

airsea_fwa_all = [airsea_fwa_pre1990;airsea_fwa_post1990;airsea_fwa_diseq];

change_in_airseadeltapco2 = airsea_fwa_post1990-airsea_fwa_pre1990
change_in_airseadeltapco2_percent = (airsea_fwa_post1990-airsea_fwa_pre1990)./airsea_fwa_pre1990.*100
%% Calculate disequilibirum DIC inventory in river waters pre and post 1990 (does NOT use year-specific atm co2)
% 1. Find atmospheric CO2 value specific to year of stream observation
% 2. Calculate equilibrium stream CO2 system with year-specific atm CO2
% values
% 3. Calculate DIC disequilibrium concentration by difference of
% equilibrium DIC and observed DIC concentrations
% 4. Calculate mean DICdiseq pre and post 1990 for each site


% Define constants for carbonate system calculations
SCALE  = 4; % NBS scale
K1K2   = 8; % Millero, 1979, FOR PURE WATER ONLY (i.e., Sal=0)    T:    0-50  S:     0. 
SO4    = 1; % Dickson (1990) KSO4
KF     = 2; % Perez & Fraga (1987) KF
BOR    = 2; % Lee et al (2010) TB

site_cc_org_eq = [];
dic_diseq = [];
% Plot changes in estuary pH for each site
for i = 1:length(names_current)
     year_atm = year(all_data_sites.(names_current{i}).ActivityStartDate); %Extract vector of years
    for z = 1:length(year_atm)
        atm_co2(z) = rcp_co2(find(rcp_co2(:,1) == year_atm(z)),2);
    end
    alkalinity = site_cc_org.(names_current{i})(:,1);
    temperature = site_cc_org.(names_current{i})(:,48);
    salinity = site_cc_org.(names_current{i})(:,58);
    pressure = 0;
    SIL = 0;
    PO4 = 0;
    site_cc_org_eq.(names_current{i}) = CO2SYS(alkalinity,atm_co2,1,4,salinity,temperature,temperature,pressure,pressure,SIL,PO4,0,0,SCALE,K1K2,...
              SO4,KF,BOR);
    dic_diseq.(names_current{i}) = site_cc_org.(names_current{i})(:,2)-site_cc_org_eq.(names_current{i})(:,2);
    co2_diseq.(names_current{i}) = site_cc_org.(names_current{i})(:,22)-site_cc_org_eq.(names_current{i})(:,22);
    
    index_pre1990 = find(year_atm < 1990);
    index_post1990 = find(year_atm >= 1990);
    
    dic_diseq_delta.(names_current{i}) = mean(dic_diseq.(names_current{i})(index_post1990))-mean(dic_diseq.(names_current{i})(index_pre1990));
    co2_diseq_delta.(names_current{i}) = mean(co2_diseq.(names_current{i})(index_post1990))-mean(co2_diseq.(names_current{i})(index_pre1990));

    dic_diseq_delta_list(i) = dic_diseq_delta.(names_current{i});
    co2_diseq_delta_list(i) = co2_diseq_delta.(names_current{i});
    
    pre1990_DIC_diseq(i) = mean(dic_diseq.(names_current{i})(index_pre1990));
    post1990_DIC_diseq(i) = mean(dic_diseq.(names_current{i})(index_post1990));
    pre1990_co2_diseq(i) = mean(co2_diseq.(names_current{i})(index_pre1990));
    post1990_co2_diseq(i) = mean(co2_diseq.(names_current{i})(index_post1990));
    discharge_current(i) = discharge_all.(names_current{i});
    load_dicdiseq_pre1990(i) = pre1990_DIC_diseq(i) ./1e6 .* 12.01 .* 1000 .* discharge_current(i) .* conv_cfs .* 31536000; %g/year
    load_dicdiseq_post1990(i) = post1990_DIC_diseq(i) ./1e6 .* 12.01 .* 1000 .* discharge_current(i) .* conv_cfs .* 31536000; %g/year
    load_dictotal_pre1990(i) = allpast_estuary_cc.(names_current{i})(1,2)./1e6 .* 12.01 .* 1000 .* discharge_current(i) .* conv_cfs .* 31536000; %g/year;
    load_dictotal_post1990(i) = current_estuary_cc.(names_current{i})(1,2)./1e6 .* 12.01 .* 1000 .* discharge_current(i) .* conv_cfs .* 31536000; %g/year;
clear atm_co2 alkalinity temperature salinity
end

dic_diseq_delta_dischargeweight = sum(discharge_weight .* dic_diseq_delta_list)
co2_diseq_delta_dischargeweight = sum(discharge_weight .* co2_diseq_delta_list)
%%
delta_dicdiseq_total_grams = sum(load_dicdiseq_post1990)-sum(load_dicdiseq_pre1990) %grams/year
delta_dicdiseq_total_Mt = delta_dicdiseq_total_grams ./ 1e12 %Mt C per year
delta_dicdiseq_total_MtCO2 = delta_dicdiseq_total_grams ./ 1e12 .* 44.01 ./ 12.01 % Mt CO2 per year
total_load_dicdiseq_pre1990_MtCO2 = sum(load_dicdiseq_pre1990)./ 1e12 .* 44.01 ./ 12.01 % Mt CO2 per year
total_load_dicdiseq_post1990_MtCO2 = sum(load_dicdiseq_post1990)./ 1e12 .* 44.01 ./ 12.01 % Mt CO2 per year
total_load_dictotal_pre1990_MtCO2 = sum(load_dictotal_pre1990)./ 1e12 .* 44.01 ./ 12.01 % Mt CO2 per year
total_load_dictotal_post1990_MtCO2 = sum(load_dictotal_post1990)./ 1e12 .* 44.01 ./ 12.01 % Mt CO2 per year
total_load_dictotal_pre1990_Mt = sum(load_dictotal_pre1990)./ 1e12 % Mt C per year
total_load_dictotal_post1990_Mt = sum(load_dictotal_post1990)./ 1e12 % Mt C per year

percent_reduction_emissions_of_totalload = delta_dicdiseq_total_MtCO2 / total_load_dictotal_post1990_MtCO2 *100
percent_reduction_emissions_of_baselineemissions  = delta_dicdiseq_total_MtCO2 / total_load_dicdiseq_pre1990_MtCO2 * 100

