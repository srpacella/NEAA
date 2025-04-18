%% Add atmospheric CO2 values to observed stream datasets based on year and calculate the sea-air pCO2 differential

% Create strcuture of sea-air pCO2 differentials for all stream
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

% Calculate medians sea-air disequilibria for Pre-and Post-1990 at each
% site
for i = 1:length(names_current)
    pre_find = find(airsea_observed.(names_current{i})(:,1) < 1990);
    post_find = find(airsea_observed.(names_current{i})(:,1) >= 1990);
    airsea_medians(i,1) = nanmedian(airsea_observed.(names_current{i})(pre_find,3));
    airsea_medians(i,2) = nanmedian(airsea_observed.(names_current{i})(post_find,3));
    airsea_medians(i,3) = airsea_medians(i,2) - airsea_medians(i,1);
end

% Calculated discharged-weighted mean of sea-air disequilibrium change
airsea_fwa_pre1990= sum(discharge_weight'.*airsea_medians(:,1))
airsea_fwa_post1990= sum(discharge_weight'.*airsea_medians(:,2))
airsea_fwa_diseq= sum(discharge_weight'.*airsea_medians(:,3))

airsea_fwa_all = [airsea_fwa_pre1990;airsea_fwa_post1990;airsea_fwa_diseq];

%% Make boxplots of Pre/Post 1990 sea-air CO2 disequilibria
figure
boxchart(airsea_medians,'MarkerStyle','none')
hold on
plot(airsea_fwa_all,'-o')
grid on
ylabel('Stream \DeltapCO_2 disequilibrium (\muatm)')
xlabel('Pre-1990, Post-1990, Post-Pre')

%% Make scatterplots of Pre/Post 1990 sea-air CO2 disequilibria
figure
scatter(airsea_medians(:,1),airsea_medians(:,2),discharge_current/1000+10,discharge_current,'filled')
hold on
plot([0 10000],[0 10000],'--k','LineWidth',2)
xlabel('Pre-1990 Stream \DeltapCO_2 disequilibrium (\muatm)')
ylabel('Post-1990 Stream \DeltapCO_2 disequilibrium (\muatm)')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';
% %%
% groupname = {'A','B','C'};
% colors = 'rgb';
% colors = colors(findgroups(groupname));
% figure
% hold on
% %boxchart(airsea_medians(:,1),'MarkerStyle','none','BoxFaceColor',colors(1),'XData',1*ones(length(airsea_medians),1))
% boxchart(airsea_medians(:,2),'MarkerStyle','none','BoxFaceColor',colors(2),'XData',2*ones(length(airsea_medians),1))
% boxchart(airsea_medians(:,3),'MarkerStyle','none','BoxFaceColor',colors(3),'XData',3*ones(length(airsea_medians),1))
% 

