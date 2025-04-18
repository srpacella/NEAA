%% Calculate changes in stream alkalinity 

delta_stream_alk = [];
delta_stream_dic = [];

for i = 1:length(names_current)
    delta_stream_alk(i) = current_estuary_cc.(names_current{i})(1,1)-allpast_estuary_cc.(names_current{i})(1,1);
    delta_stream_dic(i) = current_estuary_cc.(names_current{i})(1,2)-allpast_estuary_cc.(names_current{i})(1,2);
    discharge_current(i) = discharge_all.(names_current{i});
    load_alk_pre1990(i) = allpast_estuary_cc.(names_current{i})(1,1) .* discharge_current(i) .* conv_cfs ./ 1000 .* 31536000 ./ 1000000; %mol/year;
    load_alk_post1990(i) = current_estuary_cc.(names_current{i})(1,1) .* discharge_current(i) .* conv_cfs ./ 1000 .* 31536000 ./ 1000000; %mol/year;
    load_dic_pre1990(i) = allpast_estuary_cc.(names_current{i})(1,2) .* discharge_current(i) .* conv_cfs ./ 1000 .* 31536000 ./ 1000000; %mol/year;
    load_dic_post1990(i) = current_estuary_cc.(names_current{i})(1,2) .* discharge_current(i) .* conv_cfs ./ 1000 .* 31536000 ./ 1000000; %mol/year;

end


% Calculate discharge-weighted average of median changes in stream
% alkalinity concentration
delta_alk_weightedsum = sum(delta_stream_alk.*discharge_current)/sum(discharge_current); %umol/kg
delta_dic_weightedsum = sum(delta_stream_dic.*discharge_current)/sum(discharge_current); %umol/kg
delta_alk_load = sum(load_alk_post1990-load_alk_pre1990); %mol/year
delta_dic_load = sum(load_dic_post1990-load_dic_pre1990); %mol/year
delta_dic_load_grams = delta_dic_load .* 12.01; %g/year
total_dic_load_post1990 = sum(load_dic_post1990) .*12.01; %g/year
total_alk_load_post1990 = sum(load_alk_post1990); %mol/year
total_alk_load_pre1990 = sum(load_alk_pre1990); %mol/year

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

