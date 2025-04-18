%% Calculate mixing curves for all Alk and pH corrected sites using end-member median values
% Only uses sites with Alkorg correction (n=140 sites)
% Includes calcium correction for rivers

% Define constants for carbonate system calculations
SCALE  = 1; % Total scale
K1K2   = 14; % Millero, 2010  T:    0-50  S:  1-50. Seaw. scale. Real seawater.
SO4    = 1; % Dickson (1990) KSO4
KF     = 2; % Perez & Fraga (1987) KF
BOR    = 2; % Lee et al (2010) TB
pressure = 0;
SIL = 0;
PO4 = 0;

for i = 1:length(names_allcorrect);
    river_data = site_cc_org.(names_allcorrect{i});
    
            %Define ocean end-member
            ocean_dic = nanmedian(all_ocean_cc.(names_allcorrect{i}).DIC_y);
            ocean_alk = nanmedian(all_ocean_cc.(names_allcorrect{i}).TALK);
            ocean_sal = nanmedian(all_ocean_cc.(names_allcorrect{i}).S);
            ocean_temp = nanmedian(all_ocean_cc.(names_allcorrect{i}).T);
        
            %Define river end-member
            river_dic = nanmedian(site_cc_org.(names_allcorrect{i})(:,2));
            river_alk = nanmedian(site_cc_org.(names_allcorrect{i})(:,1));
            river_sal = nanmedian(site_cc_org.(names_allcorrect{i})(:,58));
            river_temp = nanmedian(site_cc_org.(names_allcorrect{i})(:,48));

            f_m = [0:0.02:1]; %ocean mixing fraction

            salinity_mix = (f_m.*ocean_sal) + ((1-f_m).*river_sal);
            dic_mix = (f_m.*ocean_dic) + ((1-f_m).*river_dic);
            alk_mix = (f_m.*ocean_alk) + ((1-f_m).*river_alk);
            temperature_mix = (f_m.*ocean_temp) + ((1-f_m).*river_temp);
            cc_estuary_org.(names_allcorrect{i}) = CO2SYS(alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,temperature_mix,pressure,pressure,SIL,PO4,0,0,SCALE,K1K2,...
                  SO4,KF,BOR);
            sens_estuary_org.(names_allcorrect{i}) = derivnum ('par2',alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,temperature_mix,pressure,pressure,... 
                                       SIL,PO4,0,0,...
                                       SCALE,K1K2,SO4,KF,BOR);
            sens_estuary_org_alk.(names_allcorrect{i}) = derivnum ('par1',alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,temperature_mix,pressure,pressure,... 
                                       SIL,PO4,0,0,...
                                       SCALE,K1K2,SO4,KF,BOR);                       
            
 
end

% Correct omega values for stream calcium concentrations
for i = 1:length(names_allcorrect)
    calcium_median = nanmedian(all_calcium.(names_allcorrect{i}).ResultMeasureValue); %mg/L
    calcium_correct = calcium_median ./ 1000 ./ 40.078; %converts to mol/kg
    %[omega_arag_correct omega_calc_correct] = calcium_omega_correct(current_estuary_cc.(names_current{i})(:,36),current_estuary_cc.(names_current{i})(:,35),calcium_correct,current_estuary_cc.(names_current{i})(58));
    Ca_fw = calcium_correct;
    S_sw = 35; % ocean end-member salinity
    Ca_sw = 0.010285; % ocean end-member calcium concentration in mol/kg
    S = cc_estuary_org.(names_allcorrect{i})(:,58);
    Ca_correct_beckwith = (Ca_fw./Ca_sw).*((S_sw./S)-1);
    cc_estuary_org.(names_allcorrect{i})(:,100) = cc_estuary_org.(names_allcorrect{i})(:,36) + (cc_estuary_org.(names_allcorrect{i})(:,36).*Ca_correct_beckwith);
    cc_estuary_org.(names_allcorrect{i})(:,101) = cc_estuary_org.(names_allcorrect{i})(:,35) + (cc_estuary_org.(names_allcorrect{i})(:,35).*Ca_correct_beckwith);
    if isnan(cc_estuary_org.(names_allcorrect{i})(1,100)) == 1
        cc_estuary_org.(names_allcorrect{i})(:,100) = cc_estuary_org.(names_allcorrect{i})(:,36)
        cc_estuary_org.(names_allcorrect{i})(:,101) = cc_estuary_org.(names_allcorrect{i})(:,35)
    end
end

%% Map of mean estuary sensitivity factors

for i = 1:length(names_allcorrect);
    estuary_sensH_means_org(i) = nanmean(sens_estuary_org.(names_allcorrect{i})(:,13)); %H+ sensitivity
    estuary_sensco2_means_org(i) = nanmean(sens_estuary_org_alk.(names_allcorrect{i})(:,14)); %deltapco2/deltaALK
end


figure
geoscatter(lat_orgstreams,lon_orgstreams,estuary_sensH_means_org.*200,estuary_sensH_means_org,'filled');
%gb.ColorVariable = sensH_orgstreams_median;
c = colorbar;
c.Label.String = "Mean estuary H^+ sensitivity factor";
geobasemap colorterrain
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';
caxis([0 0.8]); % Truncates range to reduce effect of few outlier values - needs to be manually edited on figure

figure
geoscatter(lat_orgstreams,lon_orgstreams,-estuary_sensco2_means_org.*20,-estuary_sensco2_means_org,'filled');
%gb.ColorVariable = sensH_orgstreams_median;
c = colorbar;
c.Label.String = "Mean estuary \partialpCO_2/\partialALK";
geobasemap colorterrain
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';
%caxis([0 0.8]); % Truncates range to reduce effect of few outlier values - needs to be manually edited on figure

figure
scatter(estuary_sensH_means_org,estuary_sensco2_means_org)
xlabel('Mean estuary \partialH^+/\partialDIC')
ylabel('Mean estuary \partialpCO_2/\partialALK')
box on
grid on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';