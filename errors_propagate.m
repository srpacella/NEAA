%% Propagates errors from stream and ocean end-members through the full estuarine salinity spectrum

% Test error propogation using 1 site

clear u_stream_pco2
clear u_ocean_pco2
clear u_mix_pco2
figure
u_stream_pco2 = [];
u_ocean_pco2 = [];
u_mix_pco2 = [];
for i = 1:length(names_allcorrect)
    u_stream_pco2.(names_allcorrect{i}) = nanmean(errors_org.(names_allcorrect{i})(:,14));
    u_ocean_pco2.(names_allcorrect{i}) = nanmean(errors_ocean_org.(names_allcorrect{i})(:,14));
    u_mix_pco2.(names_allcorrect{i}) = sqrt(((1-f_m).*u_stream_pco2.(names_allcorrect{i})).^2 + ((f_m).*u_ocean_pco2.(names_allcorrect{i})).^2);
    plot(f_m,u_mix_pco2.(names_allcorrect{i}))
    hold on
end
ylabel('pCO2 error')
xlabel('f_m')

%% Calculate discharge-weighted mixing curves and show propagated uncertainties - 140 sites

discharge_weight_orgstreams = discharge_orgstreams'./sum(discharge_orgstreams);

% % First test with one site
% figure
% x=f_m';
% y = cc_estuary_org.(names_allcorrect{i})(:,22);
% err = u_mix_pco2.(names_allcorrect{i})';
% plotUnc(x,y,err,err,'EdgeColor','none','FaceAlpha',0.2);
% hold on
% plot(x,y,'k')
% 
for i = 1:length(names_allcorrect)
    medians_estuary_pco2(i,:) = cc_estuary_org.(names_allcorrect{i})(:,22);
    errors_estuary_pco2(i,:) = u_mix_pco2.(names_allcorrect{i});
end

for i = 1:length(f_m)
    fwm_estuary_pco2(i) = sum(discharge_weight_orgstreams.*medians_estuary_pco2(:,i));
    fwm_errors_pco2(i) = sum(discharge_weight_orgstreams.*errors_estuary_pco2(:,i));
    %fwm_std_pco2(i) = nanstd(
end
figure
x=f_m;
y = fwm_estuary_pco2;
err = fwm_errors_pco2;
plotUnc(x,y,err',err','EdgeColor','none','FaceAlpha',0.2);
hold on
plot(x,y,'k')


%% Calculate for Pre- and Post-1990 dataset and plot together, 91 sites
%Uses median pCo2 for pre/post 1990
for i = 1:length(f_m)
    fwm_airsea_allpast_pco2(:,i) = sum(airsea_allpast_estuary_pCO2(i,:).*discharge_weight);
    fwm_airsea_current_pco2(:,i) = sum(airsea_current_estuary_pCO2(i,:).*discharge_weight);
end




figure
x=f_m;
y1 = fwm_airsea_allpast_pco2;
y2 = fwm_airsea_current_pco2;
err = fwm_errors_pco2;
plotUnc(x,y1,err',err','EdgeColor','none','FaceAlpha',0.2);
hold on
plotUnc(x,y2,err',err','EdgeColor','none','FaceAlpha',0.2);
hold on
plot(x,y1,'--k','LineWidth',2)
plot(x,y2,'k','LineWidth',2)
ylabel('Sea-air \Delta{\itp}CO_2 (\muatm)')
xlabel('Fraction ocean water')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';
legend('Pre-1990','Post-1990')

% %Use updated calcs from airsea_240901.m
% figure
% x=f_m;
% y1 = airsea_fwa_pre1990;
% y2 = airsea_fwa_post1990;
% err = fwm_errors_pco2;
% plotUnc(x,y1,err',err','EdgeColor','none','FaceAlpha',0.2);
% hold on
% plotUnc(x,y2,err',err','EdgeColor','none','FaceAlpha',0.2);
% hold on
% plot(x,y1,'--k','LineWidth',2)
% plot(x,y2,'k','LineWidth',2)
% ylabel('Sea-air \Delta{\itp}CO_2 (\muatm)')
% xlabel('Fraction ocean water')
% grid on
% box on
% set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
% fig = gcf
% fig.Color='w';
% legend('Pre-1990','Post-1990')

% figure
% x=f_m;
% y1 = fwm_airsea_allpast_pco2;
% y2 = fwm_airsea_current_pco2;
% err = fwm_errors_pco2;
% plotUnc(x,y1,err',err','EdgeColor','none','FaceAlpha',0.2);
% hold on
% plotUnc(x,y2,err',err','EdgeColor','none','FaceAlpha',0.2);
% hold on
% plot(x,y1,'--k','LineWidth',2)
% plot(x,y2,'k','LineWidth',2)
% ylabel('Sea-air \Delta{\itp}CO_2 (\muatm)')
% xlabel('Fraction ocean water')
% grid on
% box on
% set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
% fig = gcf
% fig.Color='w';
% legend('Pre-1990','Post-1990')

