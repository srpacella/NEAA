% Calculate carbonate system of NEAA sites
% Requires output from neaa_import.m

% CORRECTS for organic alkalinity
% CORRECTS omega values for stream calcium concentrations with
% calcium_omega_correct.m

% Uses % error for alkalinity error propagation calcs


% Define constants for carbonate system calculations
SCALE  = 4; % NBS scale
K1K2   = 8; % Millero, 1979, FOR PURE WATER ONLY (i.e., Sal=0)    T:    0-50  S:     0. 
SO4    = 1; % Dickson (1990) KSO4
KF     = 2; % Perez & Fraga (1987) KF
BOR    = 2; % Lee et al (2010) TB

%% Import mean alkalinity data
alkalinity_all = readmatrix("C:\Users\spacella\OneDrive - Environmental Protection Agency (EPA)\National Estuary Acidification Assessment\NEAA Analysis Files 240814\R Code May 2022 SRP Review\neaa_data_exports/alkalinity_all.csv");
alkalinity_all = sortrows(alkalinity_all,3);

%% Find sites with mean alkalinity > 1000 uM ( no need to correct), and sites with <1000 uM and Alkorg available (correctable).
index_highalk = find(alkalinity_all(:,2) > 1000);
index_lowalk_withorg = find(alkalinity_all(:,2) < 1000 & alkalinity_all(:,4) > 0);
index_allcorrect = [index_highalk;index_lowalk_withorg];
index_allcorrect = sort(index_allcorrect); % indices of sites in corrected dataset
names_allcorrect = all_sites_names2(index_allcorrect);


%% Calculate carbonate system for all river sites with corrected ALKorg           
for i = 1:length(all_sites_names)
    if sum(index_highalk==i) > 0 %find sites with high alk and no required correction
        pH = all_data_sites.(all_sites_names2{i}).pH_correct;
        alkalinity = all_data_sites.(all_sites_names2{i}).Alkalinity;
        temperature = all_data_sites.(all_sites_names2{i}).Temperature;
        salinity = all_data_sites.(all_sites_names2{i}).salinity;
        pressure = 0;
        SIL = 0;
        PO4 = 0;
        site_cc_org.(all_sites_names2{i}) = CO2SYS(alkalinity,pH,1,3,salinity,temperature,temperature,pressure,pressure,SIL,PO4,0,0,SCALE,K1K2,...
              SO4,KF,BOR);
    elseif sum(index_lowalk_withorg==i) > 0 %find sites with low alk and correction possible
        pH = all_data_sites.(all_sites_names2{i}).pH_correct;
        alkalinity = all_data_sites.(all_sites_names2{i}).Alkalinity.* (1-(alkalinity_all(i,4)./100));
        temperature = all_data_sites.(all_sites_names2{i}).Temperature;
        salinity = all_data_sites.(all_sites_names2{i}).salinity;
        pressure = 0;
        SIL = 0;
        PO4 = 0;
        site_cc_org.(all_sites_names2{i}) = CO2SYS(alkalinity,pH,1,3,salinity,temperature,temperature,pressure,pressure,SIL,PO4,0,0,SCALE,K1K2,...
              SO4,KF,BOR);
    else
    end
    
end

%% Correct omega values for stream calcium concentrations (no change in streams because of
for i = 1:length(names_allcorrect)
    calcium_median = nanmedian(all_calcium.(all_sites_names2{i}).ResultMeasureValue); %mg/L
    calcium_correct = calcium_median ./ 1000 ./ 40.078; %converts to mol/kg
    Ca_fw = calcium_correct;
    S_sw = 35; % ocean end-member salinity
    Ca_sw = 0.010285; % ocean end-member calcium concentration in mol/kg
    S = site_cc_org.(names_allcorrect{i})(:,58);
    Ca_correct_beckwith = (Ca_fw./Ca_sw).*((S_sw./S)-1);
    site_cc_org.(names_allcorrect{i})(:,100) = site_cc_org.(names_allcorrect{i})(:,36) + (site_cc_org.(names_allcorrect{i})(:,36).*Ca_correct_beckwith);
    site_cc_org.(names_allcorrect{i})(:,101) = site_cc_org.(names_allcorrect{i})(:,35) + (site_cc_org.(names_allcorrect{i})(:,35).*Ca_correct_beckwith);
    if isnan(site_cc_org.(names_allcorrect{i})(1,100)) == 1
        site_cc_org.(names_allcorrect{i})(:,100) = site_cc_org.(names_allcorrect{i})(:,36)
        site_cc_org.(names_allcorrect{i})(:,101) = site_cc_org.(names_allcorrect{i})(:,35)
    end
end 

%Find instances of sites without available calcium correction data
calcium_names = fieldnames(all_calcium);
for i = 1:length(names_allcorrect)
    ind = contains(calcium_names,names_allcorrect(i))
    if sum(ind) > 0
        index = find(ind,1)
        truth = isempty(all_calcium.(calcium_names{i}))
        if truth == 1
            calcium_corrections(i) = 0
        else
            calcium_corrections(i) = 1
        end
    else 
        calcium_corrections(i) = 0
    end
end

sum(calcium_corrections)
        
        
    

%% Calculate partial derivatives w/r/t deltaDIC for all river sites     
for i = 1:length(all_sites_names)
    if sum(index_highalk==i) > 0 %find sites with high alk and no required correction
        DIC = site_cc.(all_sites_names2{i})(:,2);
        ALK = all_data_sites.(all_sites_names2{i}).Alkalinity;
        temperature = all_data_sites.(all_sites_names2{i}).Temperature;
        salinity = all_data_sites.(all_sites_names2{i}).salinity;
        pressure = 0;
        SIL = 0;
        PO4 = 0;
        sens_sites_org.(all_sites_names2{i}) = derivnum ('par2',ALK,DIC,1,2,salinity,temperature,temperature,pressure,pressure,... 
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
        sens_sites_org.(all_sites_names2{i}) = derivnum ('par2',ALK,DIC,1,2,salinity,temperature,temperature,pressure,pressure,... 
                                   SIL,PO4,0,0,...
                                   SCALE,K1K2,SO4,KF,BOR);
    else
    end
    
end

%% Calculate errors associated with carbonate calculations

% **** SYNTAX:
%  [err, headers, units] = errors(PAR1,PAR2,PAR1TYPE,PAR2TYPE,...
%                                 SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,...
%                                 NH4,H2S,ePAR1,ePAR2,eSAL,eTEMP,eSI,ePO4,...
%                                 eNH4,eH2S,epK,eBt,r,pHSCALEIN,K1K2CONSTANTS,...
%                                 KSO4CONSTANT,KFCONSTANT,BORON,eCAL(optional))

ePAR1_perc = .02; % 2% measurement error from https://pubs.usgs.gov/of/2002/ofr-02-0223/E04PresetEndpoint_M.pdf
ePAR2 = 0.14; % Bias of -0.14 unites from n=122 freshwater samples from Liu et al. 2020
eSAL = 0.1;
eTEMP = 0.1;
eSI = 0;
ePO4 = 0;
eNH4 = 0;
eH2S = 0;

for i = 1:length(all_sites_names)
    if sum(index_highalk==i) > 0 %find sites with high alk and no required correction
        pH = all_data_sites.(all_sites_names2{i}).pH_correct;
        alkalinity = all_data_sites.(all_sites_names2{i}).Alkalinity;
        temperature = all_data_sites.(all_sites_names2{i}).Temperature;
        salinity = all_data_sites.(all_sites_names2{i}).salinity;
        pressure = 0;
        SIL = 0;
        PO4 = 0;
        ePAR1 = ePAR1_perc .* nanmedian(alkalinity);
        %errors_org.(all_sites_names2{i}) = errors(alkalinity,pH,1,3,salinity,temperature,temperature,pressure,pressure,SIL,PO4,0,0,ePAR1,ePAR2,eSAL,eTEMP,eSI,ePO4,eNH4,eH2S,'','',0,SCALE,K1K2,...
        %      SO4,KF,BOR);
        errors_org.(all_sites_names2{i}) = errors(alkalinity,pH,1,3,salinity,temperature,temperature,pressure,pressure,SIL,PO4,0,0,ePAR1,ePAR2,eSAL,eTEMP,0,0,0,0,0,0,0,SCALE,K1K2,SO4,KF,BOR);
   
    elseif sum(index_lowalk_withorg==i) > 0 %find sites with low alk and correction possible
        pH = all_data_sites.(all_sites_names2{i}).pH_correct;
        alkalinity = all_data_sites.(all_sites_names2{i}).Alkalinity.* (1-(alkalinity_all(i,4)./100));
        temperature = all_data_sites.(all_sites_names2{i}).Temperature;
        salinity = all_data_sites.(all_sites_names2{i}).salinity;
        pressure = 0;
        SIL = 0;
        PO4 = 0;
        ePAR1 = ePAR1_perc .* nanmedian(alkalinity);
        %errors_org.(all_sites_names2{i}) = errors(alkalinity,pH,1,3,salinity,temperature,temperature,pressure,pressure,SIL,PO4,0,0,ePAR1,ePAR2,eSAL,eTEMP,eSI,ePO4,eNH4,eH2S,'','',0,SCALE,K1K2,...
        %      SO4,KF,BOR);
        errors_org.(all_sites_names2{i}) = errors(alkalinity,pH,1,3,salinity,temperature,temperature,pressure,pressure,SIL,PO4,0,0,ePAR1,ePAR2,eSAL,eTEMP,0,0,0,0,0,0,0,SCALE,K1K2,SO4,KF,BOR);
    else
    end
    
end
%%
%Revelle error
figure
scatter(site_cc_org.(names_allcorrect{1})(:,34),errors_org.(names_allcorrect{1})(:,19)./site_cc_org.(names_allcorrect{1})(:,34).*100)

%pCO2 error
figure
scatter(site_cc_org.(names_allcorrect{1})(:,22),errors_org.(names_allcorrect{1})(:,14)./site_cc_org.(names_allcorrect{1})(:,22).*100)

%DIC error
figure
scatter(site_cc_org.(names_allcorrect{1})(:,2),errors_org.(names_allcorrect{1})(:,2)./site_cc_org.(names_allcorrect{1})(:,2).*100)


%Revelle error
figure
for i = 1:length(names_allcorrect)
    scatter(nanmedian(site_cc_org.(names_allcorrect{i})(:,34)),nanmedian(errors_org.(names_allcorrect{i})(:,19)./site_cc_org.(names_allcorrect{i})(:,34).*100))
    hold on
end
title('Revelle Factor Error')
ylabel('Median stream Revelle Factor error (%)')
xlabel('Median stream Revelle Factor')

figure
for i = 1:length(names_allcorrect)
    scatter(nanmedian(site_cc_org.(names_allcorrect{i})(:,2)),nanmedian(errors_org.(names_allcorrect{i})(:,2)./site_cc_org.(names_allcorrect{i})(:,2).*100))
    hold on
end
title('DIC Error')
ylabel('Median stream DIC error (%)')
xlabel('Median stream DIC (\mumol kg^-^1)')


figure
for i = 1:length(names_allcorrect)
    scatter(nanmedian(site_cc_org.(names_allcorrect{i})(:,22)),nanmedian(errors_org.(names_allcorrect{i})(:,14)./site_cc_org.(names_allcorrect{i})(:,22).*100))
    hold on
end
title('Revelle Factor Error')
ylabel('Median stream pCO2 error (%)')
xlabel('Median stream pCO2 (\muatm)')

figure
for i = 1:length(names_allcorrect)
    scatter(nanmedian(site_cc_org.(names_allcorrect{i})(:,22)),nanmedian(errors_org.(names_allcorrect{i})(:,14)))
    hold on
end
title('Revelle Factor Error')
ylabel('Median stream pCO2 error (\muatm)')
xlabel('Median stream pCO2 (\muatm)')

%% Collect median error information for all sites

for i = 1:length(names_allcorrect)
    errors_org_median.(names_allcorrect{i}) = nanmedian(errors_org.(names_allcorrect{i}),1); %calculate medians of all error metrics for each site
    errors_org_median_matrix(i,:) = errors_org_median.(names_allcorrect{i}); % Place same median error metrics in a matrix for convenient summary stats
    errors_org_median_matrix_percent(i,1) = errors_org_median_matrix(i,13)./1000./nanmedian(site_cc_org.(names_allcorrect{i})(:,33)).*100; %calculate % H+ error for each site
    errors_org_median_matrix_percent(i,2) = errors_org_median_matrix(i,14)./nanmedian(site_cc_org.(names_allcorrect{i})(:,22)).*100; %calculate % pco2 error for each site
    errors_org_median_matrix_percent(i,3) = errors_org_median_matrix(i,2)./nanmedian(site_cc_org.(names_allcorrect{i})(:,2)).*100; %calculate % dic error for each site
end

error_H_mean_allsites = mean(errors_org_median_matrix_percent(:,1)); % average H+ error 
error_pco2_mean_allsites = mean(errors_org_median_matrix_percent(:,2)); % average H+ error 
error_dic_mean_allsites = mean(errors_org_median_matrix_percent(:,3)); % average H+ error 
    
figure
xgroupdata1 = ones(length(errors_org_median_matrix_percent(:,1)),1);
xgroupdata2 = ones(length(errors_org_median_matrix_percent(:,1)),1).*2;
xgroupdata3 = ones(length(errors_org_median_matrix_percent(:,1)),1).*3;
ydata1 = errors_org_median_matrix_percent(:,1);
ydata2 = errors_org_median_matrix_percent(:,2);
ydata3 = errors_org_median_matrix_percent(:,3);
boxchart([xgroupdata1;xgroupdata2;xgroupdata3],[ydata1;ydata2;ydata3])
%legend('H+ % error','pCO2 % error','DIC % error')





