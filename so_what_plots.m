% Take the present day observations and add +50umol kg-1 DIC to mixing
% curve and calculate changes

%cc_estuary_org = mixing curves for each estuary using median end-member
%values

%site_cc_org = stream carbonate calculations

% add 50 umol kg-1 DIC to each mixing curve and recalculate
DIC_added = 50; %amount of DIC added to estuary in umol kg-1
ALK_added = 50;
SCALE  = 1; % Total scale
K1K2   = 14; % Millero, 2010  T:    0-50  S:  1-50. Seaw. scale. Real seawater.
SO4    = 1; % Dickson (1990) KSO4
KF     = 2; % Perez & Fraga (1987) KF
BOR    = 2; % Lee et al (2010) TB
pressure = 0;
SIL = 0;
PO4 = 0;
for i = 1:length(names_allcorrect)
    DIC = cc_estuary_org.(names_allcorrect{i})(:,2);
    DIC_eutro = DIC + DIC_added;
    alk = cc_estuary_org.(names_allcorrect{i})(:,1);
    salinity = cc_estuary_org.(names_allcorrect{i})(:,58);
    temperature = cc_estuary_org.(names_allcorrect{i})(:,48);
    cc_estuary_org_eutro.(names_allcorrect{i}) = CO2SYS(alk,DIC_eutro,1,2,salinity,temperature,temperature,pressure,pressure,SIL,PO4,0,0,SCALE,K1K2,...
                  SO4,KF,BOR);
    cc_estuary_org_eutro2.(names_allcorrect{i}) = CO2SYS(alk+ALK_added,DIC_eutro,1,2,salinity,temperature,temperature,pressure,pressure,SIL,PO4,0,0,SCALE,K1K2,...
                  SO4,KF,BOR);
end

for i = 1:length(names_allcorrect);
    estuary_deltapH_eutro(i) = nanmean(cc_estuary_org_eutro.(names_allcorrect{i})(:,43)) - nanmean(cc_estuary_org.(names_allcorrect{i})(:,43)); %change in pH
    estuary_deltaH_eutro(i) = nanmean(cc_estuary_org_eutro.(names_allcorrect{i})(:,33)) - nanmean(cc_estuary_org.(names_allcorrect{i})(:,33)); %change in H+
    ocean_sensH_org_median(i) = nanmedian(sens_ocean_org.(names_allcorrect{i})(:,13));
    percent_delta_H(i) = estuary_deltaH_eutro(i)./nanmean(cc_estuary_org.(names_allcorrect{i})(:,33)).*100;
    percent_delta_pH(i) = estuary_deltapH_eutro(i)./nanmean(cc_estuary_org.(names_allcorrect{i})(:,43)).*100;
    %Alk added
    estuary_deltapH_eutro2(i) = nanmean(cc_estuary_org_eutro2.(names_allcorrect{i})(:,43)) - nanmean(cc_estuary_org.(names_allcorrect{i})(:,43)); %change in pH
    estuary_deltaH_eutro2(i) = nanmean(cc_estuary_org_eutro2.(names_allcorrect{i})(:,33)) - nanmean(cc_estuary_org.(names_allcorrect{i})(:,33)); %change in H+
    ocean_sensH_org_median(i) = nanmedian(sens_ocean_org.(names_allcorrect{i})(:,13));
    percent_delta_H2(i) = estuary_deltaH_eutro2(i)./nanmean(cc_estuary_org.(names_allcorrect{i})(:,33)).*100;
    percent_delta_pH2(i) = estuary_deltapH_eutro2(i)./nanmean(cc_estuary_org.(names_allcorrect{i})(:,43)).*100;
end


figure
geoscatter(lat_orgstreams,lon_orgstreams,abs(estuary_deltapH_eutro)*1000,estuary_deltapH_eutro,'filled');
%gb.ColorVariable = sensH_orgstreams_median;
c = colorbar;
c.Label.String = "Mean \DeltapH_t with +50\mumol kg^-^1 DIC";
geobasemap colorterrain
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';
cmocean('thermal')
%caxis([0 0.8]); % Truncates range to reduce effect of few outlier values - needs to be manually edited on figure

figure
geoscatter(lat_orgstreams,lon_orgstreams,abs(percent_delta_H),percent_delta_H,'filled');
%gb.ColorVariable = sensH_orgstreams_median;
c = colorbar;
c.Label.String = "Mean estuary %\DeltaH^+ with +50\mumol kg^-^1 DIC";
geobasemap topographic
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';
cmocean('thermal')
%caxis([0 0.8]); % Truncates range to reduce effect of few outlier values - needs to be manually edited on figure

figure
scatter(sensH_orgstreams_median,estuary_deltapH_eutro)
xlabel('Median stream H^+ sensitivity factor')
ylabel('Mean estuary \DeltapH_T with +50\mumol kg^-^1 DIC')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';

figure
scatter(sensH_orgstreams_median,estuary_deltaH_eutro)
xlabel('Median stream H^+ sensitivity factor')
ylabel('Mean estuary \DeltaH^+ with +50\mumol kg^-^1 DIC')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';

figure
scatter(ocean_sensH_org_median,estuary_deltaH_eutro)
xlabel('Median ocean H^+ sensitivity factor')
ylabel('Mean estuary \DeltaH^+ with +50\mumol kg^-^1 DIC')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';

%Percent change in acidity
figure
scatter(ocean_sensH_org_median,percent_delta_H)
xlabel('Median ocean H^+ sensitivity factor')
ylabel('Mean % estuary \DeltaH^+ with +50\mumol kg^-^1 DIC')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';

figure
scatter(sensH_orgstreams_median,percent_delta_pH)
xlabel('Median stream H^+ sensitivity factor')
ylabel('Mean % estuary \DeltapH with +50\mumol kg^-^1 DIC')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';

figure
scatter(estuary_sensH_means_org,estuary_deltaH_eutro)
xlabel('Mean estuary H^+ sensitivity factor')
ylabel('Mean estuary \DeltaH^+ with +50\mumol kg^-^1 DIC')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';

%%
figure

t = tiledlayout(1,1);
ax1 = axes(t);
%p1 = plot(ax1,sensH_orgstreams_median,estuary_deltaH_eutro,'k*');
p1 = scatter(ax1,sensH_orgstreams_median,estuary_deltaH_eutro,'MarkerFaceColor','g','MarkerEdgeColor','k');
xlabel(ax1,'Median stream H^+ sensitivity factor');
ax1.XAxis.Exponent = 0;
ax1.Box = 'off';
ylabel('Mean estuary \DeltaH^+ with +50\mumol kg^-^1 DIC')
ax2 = axes(t);
%p2 = plot(ax2,ocean_sensH_org_median,estuary_deltaH_eutro,'rx');
p2 = scatter(ax2,ocean_sensH_org_median,estuary_deltaH_eutro,'Marker','square','MarkerFaceColor','b','MarkerEdgeColor','k');
xlabel(ax2, 'Median ocean H^+ sensitivity factor');
ax2.XAxisLocation = 'top';
ax2.XAxis.Exponent = 0;
ax2.Color = 'none';
ax2.Box = 'off';
grid on
legend([p1, p2], {'Stream', 'Ocean'}, 'Location', 'north')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';

figure
t = tiledlayout(1,1);
ax1 = axes(t);
%p1 = plot(ax1,sensH_orgstreams_median,estuary_deltapH_eutro,'k*');
p1 = scatter(ax1,sensH_orgstreams_median,estuary_deltapH_eutro,'MarkerFaceColor','g','MarkerEdgeColor','k');
xlabel(ax1,'Median stream H^+ sensitivity factor');
ax1.XAxis.Exponent = 0;
ax1.Box = 'off';
ylabel('Mean estuary \DeltapH_T with +50\mumol kg^-^1 DIC')
ax2 = axes(t);
%p2 = plot(ax2,ocean_sensH_org_median,estuary_deltapH_eutro,'rx');
p2 = scatter(ax2,ocean_sensH_org_median,estuary_deltapH_eutro,'Marker','square','MarkerFaceColor','b','MarkerEdgeColor','k');
xlabel(ax2, 'Median ocean H^+ sensitivity factor');
ax2.XAxisLocation = 'top';
ax2.XAxis.Exponent = 0;
ax2.Color = 'none';
ax2.Box = 'off';
grid on
legend([p1, p2], {'Stream', 'Ocean'}, 'Location', 'north')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';


%%
figure

t = tiledlayout(1,1);
ax1 = axes(t);
p1 = plot(ax1,sensH_orgstreams_median,percent_delta_H,'k*');
xlabel(ax1,'Median stream H^+ sensitivity factor');
ax1.XAxis.Exponent = 0;
ax1.Box = 'off';
ylabel('Mean % estuary \DeltaH^+ with +50\mumol kg^-^1 DIC')
ax2 = axes(t);
p2 = plot(ax2,ocean_sensH_org_median,percent_delta_H,'rx');
xlabel(ax2, 'Median ocean H^+ sensitivity factor');
ax2.XAxisLocation = 'top';
ax2.XAxis.Exponent = 0;
ax2.Color = 'none';
ax2.Box = 'off';
grid on
legend([p1, p2], {'Stream', 'Ocean'}, 'Location', 'north')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';

figure
t = tiledlayout(1,1);
ax1 = axes(t);
p1 = plot(ax1,sensH_orgstreams_median,estuary_deltapH_eutro,'k*');
xlabel(ax1,'Median stream H^+ sensitivity factor');
ax1.XAxis.Exponent = 0;
ax1.Box = 'off';
ylabel('Mean estuary \DeltapH_T with +50\mumol kg^-^1 DIC')
ax2 = axes(t);
p2 = plot(ax2,ocean_sensH_org_median,estuary_deltapH_eutro,'rx');
xlabel(ax2, 'Median ocean H^+ sensitivity factor');
ax2.XAxisLocation = 'top';
ax2.XAxis.Exponent = 0;
ax2.Color = 'none';
ax2.Box = 'off';
grid on
legend([p1, p2], {'Stream', 'Ocean'}, 'Location', 'north')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';

%%
figure

t = tiledlayout(1,1);
ax1 = axes(t);
p1 = plot(ax1,sensH_orgstreams_median,percent_delta_H2,'k*');
xlabel(ax1,'Median stream H^+ sensitivity factor');
ax1.XAxis.Exponent = 0;
ax1.Box = 'off';
ylabel('Mean % estuary \DeltaH^+ with +50\mumol kg^-^1 DIC')
ax2 = axes(t);
p2 = plot(ax2,ocean_sensH_org_median,percent_delta_H2,'rx');
xlabel(ax2, 'Median ocean H^+ sensitivity factor');
ax2.XAxisLocation = 'top';
ax2.XAxis.Exponent = 0;
ax2.Color = 'none';
ax2.Box = 'off';
grid on
legend([p1, p2], {'Stream', 'Ocean'}, 'Location', 'north')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';

figure
t = tiledlayout(1,1);
ax1 = axes(t);
p1 = plot(ax1,sensH_orgstreams_median,estuary_deltapH_eutro2,'k*');
xlabel(ax1,'Median stream H^+ sensitivity factor');
ax1.XAxis.Exponent = 0;
ax1.Box = 'off';
ylabel('Mean estuary \DeltapH_T with +50\mumol kg^-^1 DIC')
ax2 = axes(t);
p2 = plot(ax2,ocean_sensH_org_median,estuary_deltapH_eutro2,'rx');
xlabel(ax2, 'Median ocean H^+ sensitivity factor');
ax2.XAxisLocation = 'top';
ax2.XAxis.Exponent = 0;
ax2.Color = 'none';
ax2.Box = 'off';
grid on
legend([p1, p2], {'Stream', 'Ocean'}, 'Location', 'north')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',16,'LineWidth',1)
fig = gcf
fig.Color='w';

