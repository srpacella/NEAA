%% Create table with summary information for stations and write to file

% USGS Station code
% Station name
% Latitude
% Longitude
% HUC code
% Water resource region
% Drainage area
% Average discharge
% Median pHcorrect
% Median alkalinity (uncorrect)
% ALKorg
% Median pCO2 (calculated)
% Median DIC (calculated)
% Median stream dH/dDIC
% Median estuary dH/dDIC
% Median pH error (%)
% Median pCO2 error (%)
% Median DIC error (%)


% Create each table column

table_sitecode = all_sites_info(:,2);
table_sitename = all_sites_info(:,3);
table_lats = all_sites_info(:,7);
table_lons = all_sites_info(:,8);
table_hucs = hucs;
table_drainage = drainage;
table_discharge = discharge_means_numeric;
table_alkalinity = alkalinity_all(:,2);
table_alkorg = alkalinity_all(:,1);
 
%% Run loops to calculate summary stats and write to columns for table

for i = 1:length(all_sites_names2)
    table_pH(i) = nanmedian(site_cc.(all_sites_names2{i})(:,43));
    if sum(strcmp(names_allcorrect,all_sites_names2(i))) > 0
        table_streamsens(i) = nanmedian(sens_sites_org.(all_sites_names2{i})(:,13));
        table_estuarysens(i) = nanmean(sens_estuary_org.(all_sites_names2{i})(:,13));
        table_streamsens_alk(i) = nanmedian(sens_sites_org_alk.(all_sites_names2{i})(:,14));
        table_estuarysens_alk(i) = nanmean(sens_estuary_org_alk.(all_sites_names2{i})(:,14));
        table_pco2(i) = nanmedian(site_cc_org.(all_sites_names2{i})(:,22));
        table_dic(i) = nanmedian(site_cc_org.(all_sites_names2{i})(:,2));
        table_error_H(i) = nanmedian(errors_org.(all_sites_names2{i})(:,13))./1000./nanmedian(site_cc_org.(all_sites_names2{i})(:,33)).*100;
        table_error_pco2(i) = nanmedian(errors_org.(all_sites_names2{i})(:,14))./nanmedian(site_cc_org.(all_sites_names2{i})(:,22)).*100;
        table_error_dic(i) = nanmedian(errors_org.(all_sites_names2{i})(:,2))./nanmedian(site_cc_org.(all_sites_names2{i})(:,2)).*100;
    else
        table_streamsens(i) = NaN;
        table_estuarysens(i) = NaN;
        table_streamsens_alk(i) = NaN;
        table_estuarysens_alk(i) = NaN;
        table_pco2(i) = NaN;
        table_dic(i) = NaN
        table_error_H(i) = NaN;
        table_error_pco2(i) = NaN;
        table_error_dic(i) = NaN;
    end
end
table_pH = table_pH';
table_streamsens = table_streamsens';
table_estuarysens = table_estuarysens';
table_streamsens_alk = table_streamsens_alk';
table_estuarysens_alk = table_estuarysens_alk';
table_pco2 = table_pco2';
table_dic = table_dic';
table_error_H = table_error_H';
table_error_pco2 = table_error_pco2';
table_error_dic = table_error_dic';
    
%% Write columns to table

table_master = table(table_sitecode,table_sitename,table_lats,table_lons,table_drainage,table_discharge,table_pH,table_alkalinity,table_alkorg,table_streamsens,table_estuarysens,table_streamsens_alk,table_estuarysens_alk,table_pco2,table_dic,table_error_H,table_error_pco2,table_error_dic)

table_master_cut = table_master(table_master.table_alkalinity >= 200, :);

figure
geoscatter(str2double(table_master.table_lats),str2double(table_master.table_lons),table_master.table_estuarysens.*200,table_master.table_estuarysens,'filled');
%gb.ColorVariable = sensH_orgstreams_median;
c = colorbar;
c.Label.String = "Mean estuary H^+ sensitivity factor";
geobasemap colorterrain
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';

figure
geoscatter(str2double(table_master_cut.table_lats),str2double(table_master_cut.table_lons),table_master_cut.table_estuarysens.*200,table_master_cut.table_estuarysens,'filled');
%gb.ColorVariable = sensH_orgstreams_median;
c = colorbar;
c.Label.String = "Mean estuary H^+ sensitivity factor";
geobasemap colorterrain
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';

%% Write table to excel file

filename = 'neaa_table_master.xlsx';
writetable(table_master,filename,'Sheet',1)

median_ph_error = nanmedian(table_error_H)
median_pco2_error = nanmedian(table_error_pco2)
median_dic_error = nanmedian(table_error_dic)

%% Summary figures from median calcs

figure
scatter(table_alkalinity,table_pH,[],table_pco2,'filled')
c = colorbar;
c.Label.String = "Median stream pCO_2 (\muatm)";
axis tight
cmocean('thermal')
caxis([0 5000]);
ylabel('Stream pH_T')
xlabel('Stream alkalinity (\mumol kg^-^1)')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';


