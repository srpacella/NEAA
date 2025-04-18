%% This code uses averaged atmospheric values for two time bins (pre/post 1990), not year-specific values

%Calculate and plot: 1). changes in DIC levels in river discharge, 2).
% changes in air-sea disequilibrium for pCO2 and DIC, 3). Disequilibrium
% DIC inventory in each estuary and how this has changed 

% Calculate and plot how river air-sea pCO2 has changed

%% DIC delivered by rivers will end up being transferred to the 1). atmosphere via gas exchange, 
% or 2).the coastal ocean.  Air-river gasex will be controlled by relative
% rate of change of deltapco2 between river and air, and can be estiamted
% by looking at DIC disequilibrium inventory in river through time.  The
% coastal ocean acts as a buffer, reacting with riverine CO2 to form HCO3.
% OA will reduce this role of the coastal ocean, resulting in more coastal
% ocean acidification from acidic rivers AND more transport of this CO2 to
% the atmosphere due to reduced buffering of transported riverine CO2.  We
% can look at riverince DICdiseq through time to look at how carbon flux
% has been altered.  

% think about partial derivatives - how is coastal ocean attenuating river
% co2 for a given unit of river co2 - deltaoceanpCO2 / deltariverpCO2

%% Calculate changes in river DIC disequilibrium, and calculate how much this has changed pre/post 1990


airsea_allpast_estuary_pCO2 = [];
airsea_current_estuary_pCO2 = [];
% Plot changes in estuary pCO2 for each site
for i = 1:length(names_current)
     airsea_allpast_estuary_pCO2(:,i) = allpast_estuary_cc.(names_current{i})(:,22)-baseline_co2;
     airsea_current_estuary_pCO2(:,i) = current_estuary_cc.(names_current{i})(:,22)-current_co2;
end

figure
X = f_m;
Y_noisy_pCO2_airsea = airsea_allpast_estuary_pCO2';
plot_distribution_prctile(X,Y_noisy_pCO2_airsea,'Prctile',[25 50 75 90]);
hold on
Y_noisy_pCO2_airsea = airsea_current_estuary_pCO2';
plot_distribution_prctile(X,Y_noisy_pCO2_airsea,'Prctile',[25 50 75 90]);
plot_distribution_prctile(X,Y_noisy_pCO2_airsea,'Prctile',[25 50 75 90]);
ylabel('Sea-air \Delta{\itp}CO_2 (\muatm)')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
legend('Pre-1990','Post-1990')

figure
X = f_m;
Y_noisy_pCO2_airsea_delta = airsea_current_estuary_pCO2'./airsea_allpast_estuary_pCO2';
plot_distribution_prctile(X,Y_noisy_pCO2_airsea_delta,'Prctile',[25 50]);
ylabel('Post-1990:Pre-1900 ratio of \Delta{\itp}CO_2_,_s_e_a_-_a_i_r (\muatm)')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';

%Weight by discharge
figure
X = f_m;
Y_noisy_pCO2_airsea_delta = airsea_current_estuary_pCO2'./airsea_allpast_estuary_pCO2';
for i = 1:length(f_m)
    pCO2_airsea_delta_dischargeweighted(i) = sum(Y_noisy_pCO2_airsea_delta(:,i).* discharge_weight');
end

plot(X,pCO2_airsea_delta_dischargeweighted');
ylabel('Post-1990:Pre-1900 ratio of \Delta{\itp}CO_2_,_s_e_a_-_a_i_r (\muatm)')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
%% Calculate disequilibirum DIC inventory in river waters pre and post 1990
% Define constants for carbonate system calculations
SCALE  = 4; % NBS scale
K1K2   = 8; % Millero, 1979, FOR PURE WATER ONLY (i.e., Sal=0)    T:    0-50  S:     0. 
SO4    = 1; % Dickson (1990) KSO4
KF     = 2; % Perez & Fraga (1987) KF
BOR    = 2; % Lee et al (2010) TB

pre1990_DIC_diseq = [];
post1990_DIC_diseq = [];
% Plot changes in estuary pH for each site
for i = 1:length(names_current)
    alkalinity_old = allpast_estuary_cc.(names_current{i})(1,1);
    temperature_old = allpast_estuary_cc.(names_current{i})(1,48);
    alkalinity_new = current_estuary_cc.(names_current{i})(1,1);
    temperature_new = current_estuary_cc.(names_current{i})(1,48);
    salinity = 0;
    pressure = 0;
    SIL = 0;
    PO4 = 0;
    cc_dic_eq_pre1990 = CO2SYS(alkalinity_old,baseline_co2,1,4,salinity,temperature_old,temperature_old,pressure,pressure,SIL,PO4,0,0,SCALE,K1K2,...
              SO4,KF,BOR);
    cc_dic_eq_post1990 = CO2SYS(alkalinity_new,current_co2,1,4,salinity,temperature_old,temperature_new,pressure,pressure,SIL,PO4,0,0,SCALE,K1K2,...
              SO4,KF,BOR);
    pre1990_DIC_diseq(i) = allpast_estuary_cc.(names_current{i})(1,2) - cc_dic_eq_pre1990(2);
    post1990_DIC_diseq(i) = current_estuary_cc.(names_current{i})(1,2) - cc_dic_eq_post1990(2);
    discharge_current(i) = discharge_all.(names_current{i});
    load_dicdiseq_pre1990(i) = pre1990_DIC_diseq(i) ./1e6 .* 12.01 .* 1000 .* discharge_current(i) .* conv_cfs .* 31536000; %g/year
    load_dicdiseq_post1990(i) = post1990_DIC_diseq(i) ./1e6 .* 12.01 .* 1000 .* discharge_current(i) .* conv_cfs .* 31536000; %g/year
    load_dictotal_pre1990(i) = allpast_estuary_cc.(names_current{i})(1,2)./1e6 .* 12.01 .* 1000 .* discharge_current(i) .* conv_cfs .* 31536000; %g/year;
    load_dictotal_post1990(i) = current_estuary_cc.(names_current{i})(1,2)./1e6 .* 12.01 .* 1000 .* discharge_current(i) .* conv_cfs .* 31536000; %g/year;
end

delta_dicdiseq_total_grams = sum(load_dicdiseq_post1990)-sum(load_dicdiseq_pre1990) %grams/year
delta_dicdiseq_total_Mt = delta_dicdiseq_total_grams ./ 1e12 %Mt C per year
delta_dicdiseq_total_MtCO2 = delta_dicdiseq_total_grams ./ 1e12 .* 44.01 ./ 12.01 % Mt CO2 per year
total_load_dicdiseq_pre1990_MtCO2 = sum(load_dicdiseq_pre1990)./ 1e12 .* 44.01 ./ 12.01 % Mt CO2 per year
total_load_dicdiseq_post1990_MtCO2 = sum(load_dicdiseq_post1990)./ 1e12 .* 44.01 ./ 12.01 % Mt CO2 per year
total_load_dictotal_pre1990_MtCO2 = sum(load_dictotal_pre1990)./ 1e12 .* 44.01 ./ 12.01 % Mt CO2 per year
total_load_dictotal_post1990_MtCO2 = sum(load_dictotal_post1990)./ 1e12 .* 44.01 ./ 12.01 % Mt CO2 per year
total_load_dictotal_pre1990_Mt = sum(load_dictotal_pre1990)./ 1e12 % Mt C per year
total_load_dictotal_post1990_Mt = sum(load_dictotal_post1990)./ 1e12 % Mt C per year

percent_reduction_emissions_of_totalload = delta_dicdiseq_total_MtCO2 / total_load_dictotal_post1990_Mt *100
percent_reduction_emissions_of_baselineemissions  = delta_dicdiseq_total_MtCO2 / total_load_dicdiseq_pre1990_MtCO2 * 100
%%
% Calculate discharge-weighted average of median changes in stream
% alkalinity concentration
delta_alk_weightedsum = sum(delta_stream_alk.*discharge_current)/sum(discharge_current); %umol/kg
delta_dic_weightedsum = sum(delta_stream_dic.*discharge_current)/sum(discharge_current); %umol/kg
delta_alk_load = sum(load_alk_post1990-load_alk_pre1990); %mol/year
delta_dic_load = sum(load_dic_post1990-load_dic_pre1990); %mol/year
delta_dic_load_grams = delta_dic_load .* 12.01; %g/year
total_dic_load_post1990 = sum(load_dic_post1990) .*12.01; %g/year

% All sites
for i = 1:length(all_sites_names2)
    median_alk_conc(i) = nanmedian(site_cc.(all_sites_names2{i})(:,1)); %umol/kg
    median_dic_conc(i) = nanmedian(site_cc.(all_sites_names2{i})(:,2)); %umol/kg
    total_alk_load_all(i) = discharge_means_numeric(i) .* nanmedian(site_cc.(all_sites_names2{i})(:,1)); %mol/year
    total_dic_load_all(i) = discharge_means_numeric(i) .* nanmedian(site_cc.(all_sites_names2{i})(:,2)); %mol/year
end

national_dic_load = sum(total_dic_load_all).* 12.01; %g/year
national_dic_conc = sum(median_dic_conc'.*(discharge_means_numeric./sum(discharge_means_numeric))); %umol/kg
national_dic_conc_mgL = national_dic_conc ./ 1e6 .* 12.01 .* 1000 %mg/kg

median(site_cc.(all_sites_names2{118})(:,2))
median(site_cc.(all_sites_names2{118})(:,2))./ 1e6 .* 12.01 .* 1000 %mg/kg


