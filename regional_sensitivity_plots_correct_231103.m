%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Calculate sensitivities by region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Correct HUC data to only include sites corrected for ALKorg
hucs_org = hucs(index_allcorrect);
%% 01 = New England
%Find which stations in dataset are in the Pacific Northwest water resource
%region
sites_01 = startsWith(hucs_org,'01');
%pH_correct observations
figure
for i = 1:length(names_allcorrect)
    if sites_01(i) == 1
        scatter(all_data_sites.(names_allcorrect{i}).ActivityStartDate,all_data_sites.(names_allcorrect{i}).pH_correct)
        hold on
    end
end

%sensitivity factor observations
figure
for i = 1:length(names_allcorrect)
    if sites_01(i) == 1
        scatter(all_data_sites.(names_allcorrect{i}).ActivityStartDate,sens_sites_org.(names_allcorrect{i})(:,3))
        hold on
    end
end
ylabel('\Delta[H^+]/\Delta[DIC] (nmol/\mumol)')
grid on
xlabel('Salinity')
title('New England')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';

% Define constants for carbonate system calculations
SCALE  = 1; % Total scale
K1K2   = 14; % Millero, 2010  T:    0-50  S:  1-50. Seaw. scale. Real seawater.
SO4    = 1; % Dickson (1990) KSO4
KF     = 2; % Perez & Fraga (1987) KF
BOR    = 2; % Lee et al (2010) TB
pressure = 0;
SIL = 0;
PO4 = 0;
%Mixing lines
count = 0;
for i = 1:length(names_allcorrect)
    if sites_01(i) == 1
        count = count + 1;
        %Calulate site-specific mixing line
        f_m = [0:0.02:1]; %ocean mixing fraction
        riv_dic = nanmean(site_cc_org.(names_allcorrect{i})(:,2));
        riv_alk = nanmean(site_cc_org.(names_allcorrect{i})(:,1));
        riv_salinity = nanmean(site_cc_org.(names_allcorrect{i})(:,58));
        riv_temperature = nanmean(site_cc_org.(names_allcorrect{i})(:,48));
        ocean_dic = ocean_all.dic(i);
        ocean_alk = ocean_all.alkalinity(i);
        ocean_salinity = ocean_all.salinity(i);
        ocean_temperature = ocean_all.temperature(i);
        salinity_mix = (f_m.*ocean_salinity) + ((1-f_m).*riv_salinity);
        dic_mix = (f_m.*ocean_dic) + ((1-f_m).*riv_dic);
        alk_mix = (f_m.*ocean_alk) + ((1-f_m).*riv_alk);
        temperature_mix = (f_m.*ocean_temperature) + ((1-f_m).*riv_temperature);
        cc_estuary_all.(names_allcorrect{i}) = CO2SYS(alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,temperature_mix,pressure,pressure,SIL,PO4,0,0,SCALE,K1K2,...
              SO4,KF,BOR);
        sens_estuary_all.(names_allcorrect{i}) = derivnum ('par2',alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,temperature_mix,pressure,pressure,... 
                                   SIL,PO4,0,0,...
                                   SCALE,K1K2,SO4,KF,BOR);
        sites_01_sens_all(count,:) = sens_estuary_all.(names_allcorrect{i})(:,3);         
        sites_01_h_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,33);
        sites_01_sal_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,58);
        sites_01_pH_correct_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,43);
        sites_01_alk_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,1);
    end
end

figure
for i = 1:length(names_allcorrect)
    if sites_01(i) == 1
        plot(cc_estuary_all.(names_allcorrect{i})(:,58),sens_estuary_all.(names_allcorrect{i})(:,3))
        hold on
    end
end
ylabel('\Delta[H^+]/\Delta[DIC] (nmol/\mumol)')
grid on
xlabel('Salinity')
title('New England')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';

%% 02 = Mid-Atlantic
%Find which stations in dataset are in the Pacific Northwest water resource
%region
sites_02 = startsWith(hucs_org,'02');
%pH_correct observations
figure
for i = 1:length(names_allcorrect)
    if sites_02(i) == 1
        scatter(all_data_sites.(names_allcorrect{i}).ActivityStartDate,all_data_sites.(names_allcorrect{i}).pH_correct)
        hold on
    end
end
ylabel('USGS uncorrected pH_correct')
grid on
xlabel('Salinity')
title('Mid-Atlantic')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';


%sensitivity factor observations
figure
for i = 1:length(names_allcorrect)
    if sites_02(i) == 1
        scatter(all_data_sites.(names_allcorrect{i}).ActivityStartDate,sens_sites_org.(names_allcorrect{i})(:,3))
        hold on
    end
end
ylabel('\Delta[H^+]/\Delta[DIC] (nmol/\mumol)')
grid on
xlabel('Salinity')
title('Mid-Atlantic')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';

% Define constants for carbonate system calculations
SCALE  = 1; % Total scale
K1K2   = 14; % Millero, 2010  T:    0-50  S:  1-50. Seaw. scale. Real seawater.
SO4    = 1; % Dickson (1990) KSO4
KF     = 2; % Perez & Fraga (1987) KF
BOR    = 2; % Lee et al (2010) TB
pressure = 0;
SIL = 0;
PO4 = 0;
%Mixing lines
count = 0;
for i = 1:length(names_allcorrect)
    if sites_02(i) == 1
        count = count + 1;
        %Calulate site-specific mixing line
        f_m = [0:0.02:1]; %ocean mixing fraction
        riv_dic = nanmean(site_cc_org.(names_allcorrect{i})(:,2));
        riv_alk = nanmean(site_cc_org.(names_allcorrect{i})(:,1));
        riv_salinity = nanmean(site_cc_org.(names_allcorrect{i})(:,58));
        riv_temperature = nanmean(site_cc_org.(names_allcorrect{i})(:,48));
        ocean_dic = ocean_all.dic(i);
        ocean_alk = ocean_all.alkalinity(i);
        ocean_salinity = ocean_all.salinity(i);
        ocean_temperature = ocean_all.temperature(i);
        salinity_mix = (f_m.*ocean_salinity) + ((1-f_m).*riv_salinity);
        dic_mix = (f_m.*ocean_dic) + ((1-f_m).*riv_dic);
        alk_mix = (f_m.*ocean_alk) + ((1-f_m).*riv_alk);
        temperature_mix = (f_m.*ocean_temperature) + ((1-f_m).*riv_temperature);
        cc_estuary_all.(names_allcorrect{i}) = CO2SYS(alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,temperature_mix,pressure,pressure,SIL,PO4,0,0,SCALE,K1K2,...
              SO4,KF,BOR);
        sens_estuary_all.(names_allcorrect{i}) = derivnum ('par2',alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,temperature_mix,pressure,pressure,... 
                                   SIL,PO4,0,0,...
                                   SCALE,K1K2,SO4,KF,BOR);
        sites_02_sens_all(count,:) = sens_estuary_all.(names_allcorrect{i})(:,3);  
        sites_02_h_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,33);
        sites_02_sal_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,58);
        sites_02_pH_correct_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,43);
        sites_02_alk_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,1);
    end
end

figure
for i = 1:length(names_allcorrect)
    if sites_02(i) == 1
        plot(cc_estuary_all.(names_allcorrect{i})(:,58),sens_estuary_all.(names_allcorrect{i})(:,3))
        hold on
    end
end
ylabel('\Delta[H^+]/\Delta[DIC] (nmol/\mumol)')
grid on
xlabel('Salinity')
title('Mid-Atlantic')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';

%% 03 = South Atlantic-Gulf
%Find which stations in dataset are in the Pacific Northwest water resource
%region
sites_03 = startsWith(hucs_org,'03');
%pH_correct observations
figure
for i = 1:length(names_allcorrect)
    if sites_03(i) == 1
        scatter(all_data_sites.(names_allcorrect{i}).ActivityStartDate,all_data_sites.(names_allcorrect{i}).pH_correct)
        hold on
    end
end
ylabel('USGS uncorrected pH_correct')
grid on
xlabel('Salinity')
title('South Atlantic-Gulf')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';


%sensitivity factor observations
figure
for i = 1:length(names_allcorrect)
    if sites_03(i) == 1
        scatter(all_data_sites.(names_allcorrect{i}).ActivityStartDate,sens_sites_org.(names_allcorrect{i})(:,3))
        hold on
    end
end
ylabel('\Delta[H^+]/\Delta[DIC] (nmol/\mumol)')
grid on
xlabel('Salinity')
title('South Atlantic-Gulf')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';

% Define constants for carbonate system calculations
SCALE  = 1; % Total scale
K1K2   = 14; % Millero, 2010  T:    0-50  S:  1-50. Seaw. scale. Real seawater.
SO4    = 1; % Dickson (1990) KSO4
KF     = 2; % Perez & Fraga (1987) KF
BOR    = 2; % Lee et al (2010) TB
pressure = 0;
SIL = 0;
PO4 = 0;
%Mixing lines
count = 0;
for i = 1:length(names_allcorrect)
    if sites_03(i) == 1
        count = count + 1;
        %Calulate site-specific mixing line
        f_m = [0:0.02:1]; %ocean mixing fraction
        riv_dic = nanmean(site_cc_org.(names_allcorrect{i})(:,2));
        riv_alk = nanmean(site_cc_org.(names_allcorrect{i})(:,1));
        riv_salinity = nanmean(site_cc_org.(names_allcorrect{i})(:,58));
        riv_temperature = nanmean(site_cc_org.(names_allcorrect{i})(:,48));
        ocean_dic = ocean_all.dic(i);
        ocean_alk = ocean_all.alkalinity(i);
        ocean_salinity = ocean_all.salinity(i);
        ocean_temperature = ocean_all.temperature(i);
        salinity_mix = (f_m.*ocean_salinity) + ((1-f_m).*riv_salinity);
        dic_mix = (f_m.*ocean_dic) + ((1-f_m).*riv_dic);
        alk_mix = (f_m.*ocean_alk) + ((1-f_m).*riv_alk);
        temperature_mix = (f_m.*ocean_temperature) + ((1-f_m).*riv_temperature);
        cc_estuary_all.(names_allcorrect{i}) = CO2SYS(alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,temperature_mix,pressure,pressure,SIL,PO4,0,0,SCALE,K1K2,...
              SO4,KF,BOR);
        sens_estuary_all.(names_allcorrect{i}) = derivnum ('par2',alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,temperature_mix,pressure,pressure,... 
                                   SIL,PO4,0,0,...
                                   SCALE,K1K2,SO4,KF,BOR);
        sites_03_sens_all(count,:) = sens_estuary_all.(names_allcorrect{i})(:,3);    
        sites_03_h_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,33);
        sites_03_sal_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,58);
        sites_03_pH_correct_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,43);
        sites_03_alk_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,1);
    end
end

figure
for i = 1:length(names_allcorrect)
    if sites_03(i) == 1
        plot(cc_estuary_all.(names_allcorrect{i})(:,58),sens_estuary_all.(names_allcorrect{i})(:,3))
        hold on
    end
end
ylabel('\Delta[H^+]/\Delta[DIC] (nmol/\mumol)')
grid on
xlabel('Salinity')
title('South Atlantic-Gulf')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';

%% 08 = Lower Mississippi
%Find which stations in dataset are in the Pacific Northwest water resource
%region
sites_08 = startsWith(hucs_org,'08');
%pH_correct observations
figure
for i = 1:length(names_allcorrect)
    if sites_08(i) == 1
        scatter(all_data_sites.(names_allcorrect{i}).ActivityStartDate,all_data_sites.(names_allcorrect{i}).pH_correct)
        hold on
    end
end
ylabel('USGS uncorrected pH_correct')
grid on
xlabel('Salinity')
title('Lower Mississippi')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';


%sensitivity factor observations
figure
for i = 1:length(names_allcorrect)
    if sites_08(i) == 1
        scatter(all_data_sites.(names_allcorrect{i}).ActivityStartDate,sens_sites_org.(names_allcorrect{i})(:,3))
        hold on
    end
end
ylabel('\Delta[H^+]/\Delta[DIC] (nmol/\mumol)')
grid on
xlabel('Salinity')
title('Lower Mississippi')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';

% Define constants for carbonate system calculations
SCALE  = 1; % Total scale
K1K2   = 14; % Millero, 2010  T:    0-50  S:  1-50. Seaw. scale. Real seawater.
SO4    = 1; % Dickson (1990) KSO4
KF     = 2; % Perez & Fraga (1987) KF
BOR    = 2; % Lee et al (2010) TB
pressure = 0;
SIL = 0;
PO4 = 0;
%Mixing lines
count = 0;
for i = 1:length(names_allcorrect)
    if sites_08(i) == 1
        count = count + 1;
        %Calulate site-specific mixing line
        f_m = [0:0.02:1]; %ocean mixing fraction
        riv_dic = nanmean(site_cc_org.(names_allcorrect{i})(:,2));
        riv_alk = nanmean(site_cc_org.(names_allcorrect{i})(:,1));
        riv_salinity = nanmean(site_cc_org.(names_allcorrect{i})(:,58));
        riv_temperature = nanmean(site_cc_org.(names_allcorrect{i})(:,48));
        ocean_dic = ocean_all.dic(i);
        ocean_alk = ocean_all.alkalinity(i);
        ocean_salinity = ocean_all.salinity(i);
        ocean_temperature = ocean_all.temperature(i);
        salinity_mix = (f_m.*ocean_salinity) + ((1-f_m).*riv_salinity);
        dic_mix = (f_m.*ocean_dic) + ((1-f_m).*riv_dic);
        alk_mix = (f_m.*ocean_alk) + ((1-f_m).*riv_alk);
        temperature_mix = (f_m.*ocean_temperature) + ((1-f_m).*riv_temperature);
        cc_estuary_all.(names_allcorrect{i}) = CO2SYS(alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,temperature_mix,pressure,pressure,SIL,PO4,0,0,SCALE,K1K2,...
              SO4,KF,BOR);
        sens_estuary_all.(names_allcorrect{i}) = derivnum ('par2',alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,temperature_mix,pressure,pressure,... 
                                   SIL,PO4,0,0,...
                                   SCALE,K1K2,SO4,KF,BOR);
        sites_08_sens_all(count,:) = sens_estuary_all.(names_allcorrect{i})(:,3);          
        sites_08_h_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,33);
        sites_08_sal_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,58);
        sites_08_pH_correct_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,43);
        sites_08_alk_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,1);
    end
end

figure
for i = 1:length(names_allcorrect)
    if sites_08(i) == 1
        plot(cc_estuary_all.(names_allcorrect{i})(:,58),sens_estuary_all.(names_allcorrect{i})(:,3))
        hold on
    end
end
ylabel('\Delta[H^+]/\Delta[DIC] (nmol/\mumol)')
grid on
xlabel('Salinity')
title('Lower Mississippi')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';

%% 12 = Texas-Gulf
%Find which stations in dataset are in the Pacific Northwest water resource
%region
sites_12 = startsWith(hucs_org,'12');
%pH_correct observations
figure
for i = 1:length(names_allcorrect)
    if sites_12(i) == 1
        scatter(all_data_sites.(names_allcorrect{i}).ActivityStartDate,all_data_sites.(names_allcorrect{i}).pH_correct)
        hold on
    end
end

%sensitivity factor observations
figure
for i = 1:length(names_allcorrect)
    if sites_12(i) == 1
        scatter(all_data_sites.(names_allcorrect{i}).ActivityStartDate,sens_sites_org.(names_allcorrect{i})(:,3))
        hold on
    end
end
ylabel('\Delta[H^+]/\Delta[DIC] (nmol/\mumol)')
grid on
xlabel('Salinity')
title('Texas-Gulf')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';

% Define constants for carbonate system calculations
SCALE  = 1; % Total scale
K1K2   = 14; % Millero, 2010  T:    0-50  S:  1-50. Seaw. scale. Real seawater.
SO4    = 1; % Dickson (1990) KSO4
KF     = 2; % Perez & Fraga (1987) KF
BOR    = 2; % Lee et al (2010) TB
pressure = 0;
SIL = 0;
PO4 = 0;
%Mixing lines
count = 0;
for i = 1:length(names_allcorrect)
    if sites_12(i) == 1
        count = count + 1;
        %Calulate site-specific mixing line
        f_m = [0:0.02:1]; %ocean mixing fraction
        riv_dic = nanmean(site_cc_org.(names_allcorrect{i})(:,2));
        riv_alk = nanmean(site_cc_org.(names_allcorrect{i})(:,1));
        riv_salinity = nanmean(site_cc_org.(names_allcorrect{i})(:,58));
        riv_temperature = nanmean(site_cc_org.(names_allcorrect{i})(:,48));
        ocean_dic = ocean_all.dic(i);
        ocean_alk = ocean_all.alkalinity(i);
        ocean_salinity = ocean_all.salinity(i);
        ocean_temperature = ocean_all.temperature(i);
        salinity_mix = (f_m.*ocean_salinity) + ((1-f_m).*riv_salinity);
        dic_mix = (f_m.*ocean_dic) + ((1-f_m).*riv_dic);
        alk_mix = (f_m.*ocean_alk) + ((1-f_m).*riv_alk);
        temperature_mix = (f_m.*ocean_temperature) + ((1-f_m).*riv_temperature);
        cc_estuary_all.(names_allcorrect{i}) = CO2SYS(alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,temperature_mix,pressure,pressure,SIL,PO4,0,0,SCALE,K1K2,...
              SO4,KF,BOR);
        sens_estuary_all.(names_allcorrect{i}) = derivnum ('par2',alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,temperature_mix,pressure,pressure,... 
                                   SIL,PO4,0,0,...
                                   SCALE,K1K2,SO4,KF,BOR);
        sites_12_sens_all(count,:) = sens_estuary_all.(names_allcorrect{i})(:,3);    
        sites_12_h_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,33);
        sites_12_sal_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,58);
        sites_12_pH_correct_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,43);
        sites_12_alk_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,1);
    end
end

figure
for i = 1:length(names_allcorrect)
    if sites_12(i) == 1
        plot(cc_estuary_all.(names_allcorrect{i})(:,58),sens_estuary_all.(names_allcorrect{i})(:,3))
        hold on
    end
end
ylabel('\Delta[H^+]/\Delta[DIC] (nmol/\mumol)')
grid on
xlabel('Salinity')
title('Texas-Gulf')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';

%% 13 = Rio Grande
%Find which stations in dataset are in the Pacific Northwest water resource
%region
sites_13 = startsWith(hucs_org,'13');
%pH_correct observations
figure
for i = 1:length(names_allcorrect)
    if sites_13(i) == 1
        scatter(all_data_sites.(names_allcorrect{i}).ActivityStartDate,all_data_sites.(names_allcorrect{i}).pH_correct)
        hold on
    end
end

%sensitivity factor observations
figure
for i = 1:length(names_allcorrect)
    if sites_13(i) == 1
        scatter(all_data_sites.(names_allcorrect{i}).ActivityStartDate,sens_sites_org.(names_allcorrect{i})(:,3))
        hold on
    end
end
ylabel('\Delta[H^+]/\Delta[DIC] (nmol/\mumol)')
grid on
xlabel('Salinity')
title('Rio Grande')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';

% Define constants for carbonate system calculations
SCALE  = 1; % Total scale
K1K2   = 14; % Millero, 2010  T:    0-50  S:  1-50. Seaw. scale. Real seawater.
SO4    = 1; % Dickson (1990) KSO4
KF     = 2; % Perez & Fraga (1987) KF
BOR    = 2; % Lee et al (2010) TB
pressure = 0;
SIL = 0;
PO4 = 0;
%Mixing lines
count = 0;
for i = 1:length(names_allcorrect)
    if sites_13(i) == 1
        count = count + 1;
        %Calulate site-specific mixing line
        f_m = [0:0.02:1]; %ocean mixing fraction
        riv_dic = nanmean(site_cc_org.(names_allcorrect{i})(:,2));
        riv_alk = nanmean(site_cc_org.(names_allcorrect{i})(:,1));
        riv_salinity = nanmean(site_cc_org.(names_allcorrect{i})(:,58));
        riv_temperature = nanmean(site_cc_org.(names_allcorrect{i})(:,48));
        ocean_dic = ocean_all.dic(i);
        ocean_alk = ocean_all.alkalinity(i);
        ocean_salinity = ocean_all.salinity(i);
        ocean_temperature = ocean_all.temperature(i);
        salinity_mix = (f_m.*ocean_salinity) + ((1-f_m).*riv_salinity);
        dic_mix = (f_m.*ocean_dic) + ((1-f_m).*riv_dic);
        alk_mix = (f_m.*ocean_alk) + ((1-f_m).*riv_alk);
        temperature_mix = (f_m.*ocean_temperature) + ((1-f_m).*riv_temperature);
        cc_estuary_all.(names_allcorrect{i}) = CO2SYS(alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,temperature_mix,pressure,pressure,SIL,PO4,0,0,SCALE,K1K2,...
              SO4,KF,BOR);
        sens_estuary_all.(names_allcorrect{i}) = derivnum ('par2',alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,temperature_mix,pressure,pressure,... 
                                   SIL,PO4,0,0,...
                                   SCALE,K1K2,SO4,KF,BOR);
        sites_13_sens_all(count,:) = sens_estuary_all.(names_allcorrect{i})(:,3);    
        sites_13_h_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,33);
        sites_13_sal_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,58);
        sites_13_pH_correct_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,43);
        sites_13_alk_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,1);
    end
end

figure
for i = 1:length(names_allcorrect)
    if sites_13(i) == 1
        plot(cc_estuary_all.(names_allcorrect{i})(:,58),sens_estuary_all.(names_allcorrect{i})(:,3))
        hold on
    end
end
ylabel('\Delta[H^+]/\Delta[DIC] (nmol/\mumol)')
grid on
xlabel('Salinity')
title('Rio Grande')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';

%% 18 = California
%Find which stations in dataset are in the Pacific Northwest water resource
%region
sites_18 = startsWith(hucs_org,'18');
%pH_correct observations
figure
for i = 1:length(names_allcorrect)
    if sites_18(i) == 1
        scatter(all_data_sites.(names_allcorrect{i}).ActivityStartDate,all_data_sites.(names_allcorrect{i}).pH_correct)
        hold on
    end
end
ylabel('USGS uncorrected pH_correct')
grid on
xlabel('Salinity')
title('California')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';


%sensitivity factor observations
figure
for i = 1:length(names_allcorrect)
    if sites_18(i) == 1
        scatter(all_data_sites.(names_allcorrect{i}).ActivityStartDate,sens_sites_org.(names_allcorrect{i})(:,3))
        hold on
    end
end
ylabel('\Delta[H^+]/\Delta[DIC] (nmol/\mumol)')
grid on
xlabel('Salinity')
title('California')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';

% Define constants for carbonate system calculations
SCALE  = 1; % Total scale
K1K2   = 14; % Millero, 2010  T:    0-50  S:  1-50. Seaw. scale. Real seawater.
SO4    = 1; % Dickson (1990) KSO4
KF     = 2; % Perez & Fraga (1987) KF
BOR    = 2; % Lee et al (2010) TB
pressure = 0;
SIL = 0;
PO4 = 0;
%Mixing lines
count = 0;
for i = 1:length(names_allcorrect)
    if sites_18(i) == 1
        count = count + 1;
        %Calulate site-specific mixing line
        f_m = [0:0.02:1]; %ocean mixing fraction
        riv_dic = nanmean(site_cc_org.(names_allcorrect{i})(:,2));
        riv_alk = nanmean(site_cc_org.(names_allcorrect{i})(:,1));
        riv_salinity = nanmean(site_cc_org.(names_allcorrect{i})(:,58));
        riv_temperature = nanmean(site_cc_org.(names_allcorrect{i})(:,48));
        ocean_dic = ocean_all.dic(i);
        ocean_alk = ocean_all.alkalinity(i);
        ocean_salinity = ocean_all.salinity(i);
        ocean_temperature = ocean_all.temperature(i);
        salinity_mix = (f_m.*ocean_salinity) + ((1-f_m).*riv_salinity);
        dic_mix = (f_m.*ocean_dic) + ((1-f_m).*riv_dic);
        alk_mix = (f_m.*ocean_alk) + ((1-f_m).*riv_alk);
        temperature_mix = (f_m.*ocean_temperature) + ((1-f_m).*riv_temperature);
        cc_estuary_all.(names_allcorrect{i}) = CO2SYS(alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,temperature_mix,pressure,pressure,SIL,PO4,0,0,SCALE,K1K2,...
              SO4,KF,BOR);
        sens_estuary_all.(names_allcorrect{i}) = derivnum ('par2',alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,temperature_mix,pressure,pressure,... 
                                   SIL,PO4,0,0,...
                                   SCALE,K1K2,SO4,KF,BOR);
        sites_18_sens_all(count,:) = sens_estuary_all.(names_allcorrect{i})(:,3); 
        sites_18_h_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,33);
        sites_18_sal_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,58);
        sites_18_pH_correct_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,43);
        sites_18_alk_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,1);
    end
end

figure
for i = 1:length(names_allcorrect)
    if sites_18(i) == 1
        plot(cc_estuary_all.(names_allcorrect{i})(:,58),sens_estuary_all.(names_allcorrect{i})(:,3))
        hold on
    end
end
ylabel('\Delta[H^+]/\Delta[DIC] (nmol/\mumol)')
grid on
xlabel('Salinity')
title('California')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';

%% 17 = Pacific Northwest
%Find which stations in dataset are in the Pacific Northwest water resource
%region
sites_17 = startsWith(hucs_org,'17');
%pH_correct observations
figure
for i = 1:length(names_allcorrect)
    if sites_17(i) == 1
        scatter(all_data_sites.(names_allcorrect{i}).ActivityStartDate,all_data_sites.(names_allcorrect{i}).pH_correct)
        hold on
    end
end

%sensitivity factor observations
figure
for i = 1:length(names_allcorrect)
    if sites_17(i) == 1
        scatter(all_data_sites.(names_allcorrect{i}).ActivityStartDate,sens_sites_org.(names_allcorrect{i})(:,3))
        hold on
    end
end
ylabel('\Delta[H^+]/\Delta[DIC] (nmol/\mumol)')
grid on
xlabel('Salinity')
title('Pacific Northwest')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';

% Define constants for carbonate system calculations
SCALE  = 1; % Total scale
K1K2   = 14; % Millero, 2010  T:    0-50  S:  1-50. Seaw. scale. Real seawater.
SO4    = 1; % Dickson (1990) KSO4
KF     = 2; % Perez & Fraga (1987) KF
BOR    = 2; % Lee et al (2010) TB
pressure = 0;
SIL = 0;
PO4 = 0;
%Mixing lines
count = 0;
for i = 1:length(names_allcorrect)
    if sites_17(i) == 1
        count = count + 1;
        %Calulate site-specific mixing line
        f_m = [0:0.02:1]; %ocean mixing fraction
        riv_dic = nanmean(site_cc_org.(names_allcorrect{i})(:,2));
        riv_alk = nanmean(site_cc_org.(names_allcorrect{i})(:,1));
        riv_salinity = nanmean(site_cc_org.(names_allcorrect{i})(:,58));
        riv_temperature = nanmean(site_cc_org.(names_allcorrect{i})(:,48));
        ocean_dic = ocean_all.dic(i);
        ocean_alk = ocean_all.alkalinity(i);
        ocean_salinity = ocean_all.salinity(i);
        ocean_temperature = ocean_all.temperature(i);
        salinity_mix = (f_m.*ocean_salinity) + ((1-f_m).*riv_salinity);
        dic_mix = (f_m.*ocean_dic) + ((1-f_m).*riv_dic);
        alk_mix = (f_m.*ocean_alk) + ((1-f_m).*riv_alk);
        temperature_mix = (f_m.*ocean_temperature) + ((1-f_m).*riv_temperature);
        cc_estuary_all.(names_allcorrect{i}) = CO2SYS(alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,temperature_mix,pressure,pressure,SIL,PO4,0,0,SCALE,K1K2,...
              SO4,KF,BOR);
        sens_estuary_all.(names_allcorrect{i}) = derivnum ('par2',alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,temperature_mix,pressure,pressure,... 
                                   SIL,PO4,0,0,...
                                   SCALE,K1K2,SO4,KF,BOR);
        sites_17_sens_all(count,:) = sens_estuary_all.(names_allcorrect{i})(:,3);
        sites_17_h_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,33);
        sites_17_sal_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,58);
        sites_17_pH_correct_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,43);
        sites_17_alk_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,1);
    end
end

figure
for i = 1:length(names_allcorrect)
    if sites_17(i) == 1
        plot(cc_estuary_all.(names_allcorrect{i})(:,58),sens_estuary_all.(names_allcorrect{i})(:,3))
        hold on
    end
end
plot(cc_estuary_all.(names_allcorrect{1})(:,58),mean(sites_17_sens_all,1),'k','LineWidth',3)
ylabel('\Delta[H^+]/\Delta[DIC] (nmol/\mumol)')
grid on
xlabel('Salinity')
title('Pacific Northwest')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';

%% 19 = Alaska
%Find which stations in dataset are in the Pacific Northwest water resource
%region
sites_19 = startsWith(hucs_org,'19');
%pH_correct observations
figure
for i = 1:length(names_allcorrect)
    if sites_19(i) == 1
        scatter(all_data_sites.(names_allcorrect{i}).ActivityStartDate,all_data_sites.(names_allcorrect{i}).pH_correct)
        hold on
    end
end

%sensitivity factor observations
figure
for i = 1:length(names_allcorrect)
    if sites_19(i) == 1
        scatter(all_data_sites.(names_allcorrect{i}).ActivityStartDate,sens_sites_org.(names_allcorrect{i})(:,3))
        hold on
    end
end
ylabel('\Delta[H^+]/\Delta[DIC] (nmol/\mumol)')
grid on
xlabel('Salinity')
title('Alaska')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';

% Define constants for carbonate system calculations
SCALE  = 1; % Total scale
K1K2   = 14; % Millero, 2010  T:    0-50  S:  1-50. Seaw. scale. Real seawater.
SO4    = 1; % Dickson (1990) KSO4
KF     = 2; % Perez & Fraga (1987) KF
BOR    = 2; % Lee et al (2010) TB
pressure = 0;
SIL = 0;
PO4 = 0;
%Mixing lines
count = 0;
for i = 1:length(names_allcorrect)
    if sites_19(i) == 1
        count = count + 1;
        %Calulate site-specific mixing line
        f_m = [0:0.02:1]; %ocean mixing fraction
        riv_dic = nanmean(site_cc_org.(names_allcorrect{i})(:,2));
        riv_alk = nanmean(site_cc_org.(names_allcorrect{i})(:,1));
        riv_salinity = nanmean(site_cc_org.(names_allcorrect{i})(:,58));
        riv_temperature = nanmean(site_cc_org.(names_allcorrect{i})(:,48));
        ocean_dic = ocean_all.dic(i);
        ocean_alk = ocean_all.alkalinity(i);
        ocean_salinity = ocean_all.salinity(i);
        ocean_temperature = ocean_all.temperature(i);
        salinity_mix = (f_m.*ocean_salinity) + ((1-f_m).*riv_salinity);
        dic_mix = (f_m.*ocean_dic) + ((1-f_m).*riv_dic);
        alk_mix = (f_m.*ocean_alk) + ((1-f_m).*riv_alk);
        temperature_mix = (f_m.*ocean_temperature) + ((1-f_m).*riv_temperature);
        cc_estuary_all.(names_allcorrect{i}) = CO2SYS(alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,temperature_mix,pressure,pressure,SIL,PO4,0,0,SCALE,K1K2,...
              SO4,KF,BOR);
        sens_estuary_all.(names_allcorrect{i}) = derivnum ('par2',alk_mix,dic_mix,1,2,salinity_mix,temperature_mix,temperature_mix,pressure,pressure,... 
                                   SIL,PO4,0,0,...
                                   SCALE,K1K2,SO4,KF,BOR);
        sites_19_sens_all(count,:) = sens_estuary_all.(names_allcorrect{i})(:,3);
        sites_19_h_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,33);
        sites_19_sal_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,58);
        sites_19_pH_correct_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,43);
        sites_19_alk_all(count,:) = cc_estuary_all.(names_allcorrect{i})(:,1);
    end
end

figure
for i = 1:length(names_allcorrect)
    if sites_19(i) == 1
        plot(cc_estuary_all.(names_allcorrect{i})(:,58),sens_estuary_all.(names_allcorrect{i})(:,3))
        hold on
    end
end
ylabel('\Delta[H^+]/\Delta[DIC] (nmol/\mumol)')
grid on
xlabel('Salinity')
title('Alaska')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';

%% Compare regional average sensitivities with national average

%Sensitivities %%%%%NEED TO CHANGE SALINITY ON X AXIS FOR EACH
%REGION!!!!!!!!!!!!!!
figure
%plot(cc_nat_avg(:,58),sens_nat_avg(:,3),'k','LineWidth',3)
%hold on
plot(mean(sites_01_sal_all,1),mean(sites_01_sens_all,1),'LineWidth',1)
hold on
plot(mean(sites_02_sal_all,1),mean(sites_02_sens_all,1),'LineWidth',1)
plot(mean(sites_03_sal_all,1),mean(sites_03_sens_all,1),'LineWidth',1)
plot(mean(sites_08_sal_all,1),mean(sites_08_sens_all,1),'LineWidth',1)
plot(mean(sites_12_sal_all,1),mean(sites_12_sens_all,1),'LineWidth',1)
plot(mean(sites_13_sal_all,1),mean(sites_13_sens_all,1),'LineWidth',1)
plot(mean(sites_18_sal_all,1),mean(sites_18_sens_all,1),'LineWidth',1)
plot(mean(sites_17_sal_all,1),mean(sites_17_sens_all,1),'LineWidth',1)
plot(mean(sites_19_sal_all,1),mean(sites_19_sens_all,1),'LineWidth',1)
ylabel('\Delta[H^+]/\Delta[DIC] (nmol/\mumol)')
grid on
xlabel('Salinity')
legend('New England','Mid-Atlantic','South Atlantic-Gulf','Lower Mississippi','Texas-Gulf','Rio Grande','California','Pacific Northwest','Alaska')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
%%
figure
plot(mean(sites_01_sal_all,1),mean(sites_01_sens_all,1),'LineWidth',1)
hold on
plot(mean(sites_02_sal_all,1),mean(sites_02_sens_all,1),'LineWidth',1)
plot(mean(sites_03_sal_all,1),mean(sites_03_sens_all,1),'LineWidth',1)
plot(mean(sites_08_sal_all,1),mean(sites_08_sens_all,1),'LineWidth',1)
plot(mean(sites_12_sal_all,1),mean(sites_12_sens_all,1),'LineWidth',1)
plot(mean(sites_13_sal_all,1),mean(sites_13_sens_all,1),'LineWidth',1)
plot(mean(sites_18_sal_all,1),mean(sites_18_sens_all,1),'LineWidth',1)
plot(mean(sites_17_sal_all,1),mean(sites_17_sens_all,1),'LineWidth',1)
plot(mean(sites_19_sal_all,1),mean(sites_19_sens_all,1),'LineWidth',1)
ylabel('\Delta[H^+]/\Delta[DIC] (nmol/\mumol)')
grid on
xlabel('Salinity')
legend('New England','Mid-Atlantic','South Atlantic-Gulf','Lower Mississippi','Texas-Gulf','Rio Grande','California','Pacific Northwest','Alaska')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
%%
figure
X = mean(sites_01_sal_all,1);
Y_noisy = sites_01_sens_all;
plot_distribution_prctile(X,Y_noisy,'Prctile',[75]);
hold on
X = mean(sites_02_sal_all,1);
Y_noisy = sites_02_sens_all;
plot_distribution_prctile(X,Y_noisy,'Prctile',[75]);
X = mean(sites_03_sal_all,1);
Y_noisy = sites_03_sens_all;
plot_distribution_prctile(X,Y_noisy,'Prctile',[75]);
X = mean(sites_08_sal_all,1);
Y_noisy = sites_08_sens_all;
plot_distribution_prctile(X,Y_noisy,'Prctile',[75]);
X = mean(sites_12_sal_all,1);
Y_noisy = sites_12_sens_all;
plot_distribution_prctile(X,Y_noisy,'Prctile',[75]);
hold on
plot(mean(sites_13_sal_all,1),mean(sites_13_sens_all,1),'LineWidth',1)
hold on
X = mean(sites_18_sal_all,1);
Y_noisy = sites_18_sens_all;
plot_distribution_prctile(X,Y_noisy,'Prctile',[75]);
X = mean(sites_17_sal_all,1);
Y_noisy = sites_17_sens_all;
plot_distribution_prctile(X,Y_noisy,'Prctile',[75]);
X = mean(sites_19_sal_all,1);
Y_noisy = sites_19_sens_all;
plot_distribution_prctile(X,Y_noisy,'Prctile',[75]);
ylabel('\Delta[H^+]/\Delta[DIC] (nmol/\mumol)')
grid on
box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';

%% Multi-Panel figure of median + ranges for each region

figure
subplot(3,3,1)
X = mean(sites_01_sal_all,1);
Y_noisy = sites_01_sens_all;
plot_distribution_prctile(X,Y_noisy,'Prctile',[25 50 75 90]);
title('New England')
axis tight
grid on
%box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,2)
X = mean(sites_02_sal_all,1);
Y_noisy = sites_02_sens_all;
plot_distribution_prctile(X,Y_noisy,'Prctile',[25 50 75 90]);
title('Mid-Atlantic')
axis tight
grid on
%box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,3)
X = mean(sites_03_sal_all,1);
Y_noisy = sites_03_sens_all;
plot_distribution_prctile(X,Y_noisy,'Prctile',[25 50 75 90]);
title('South Atlantic-Gulf')
axis tight
grid on
%box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,4)
X = mean(sites_08_sal_all,1);
Y_noisy = sites_08_sens_all;
plot_distribution_prctile(X,Y_noisy,'Prctile',[25 50 75 90]);
title('Lower Mississippi')
axis tight
grid on
%box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,5)
X = mean(sites_12_sal_all,1);
Y_noisy = sites_12_sens_all;
plot_distribution_prctile(X,Y_noisy,'Prctile',[25 50 75 90]);
title('Texas-Gulf')
axis tight
grid on
%box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,6)
plot(mean(sites_13_sal_all,1),mean(sites_13_sens_all,1),'LineWidth',1)
title('Rio Grande')
axis tight
grid on
%box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,7)
X = mean(sites_18_sal_all,1);
Y_noisy = sites_18_sens_all;
plot_distribution_prctile(X,Y_noisy,'Prctile',[25 50 75 90]);
title('California')
axis tight
grid on
%box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,8)
X = mean(sites_17_sal_all,1);
Y_noisy = sites_17_sens_all;
plot_distribution_prctile(X,Y_noisy,'Prctile',[25 50 75 90]);
title('Pacific Northwest')
axis tight
grid on
%box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,9)
X = mean(sites_19_sal_all,1);
Y_noisy = sites_19_sens_all;
plot_distribution_prctile(X,Y_noisy,'Prctile',[25 50 75 90]);
title('Alaska')
axis tight
grid on
%box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';

%% Multi-Panel figure of median + ranges for each region; make x and y axes all the same

figure
subplot(3,3,1)
X = mean(sites_01_sal_all,1);
Y_noisy = sites_01_sens_all;
plot_distribution_prctile(X,Y_noisy,'Prctile',[25 50 75 90]);
title('New England')
xlim([0 35]);
ylim([0 3]);
grid on
%box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,2)
X = mean(sites_02_sal_all,1);
Y_noisy = sites_02_sens_all;
plot_distribution_prctile(X,Y_noisy,'Prctile',[25 50 75 90]);
title('Mid-Atlantic')
xlim([0 35]);
ylim([0 3]);
grid on
%box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,3)
X = mean(sites_03_sal_all,1);
Y_noisy = sites_03_sens_all;
plot_distribution_prctile(X,Y_noisy,'Prctile',[25 50 75 90]);
title('South Atlantic-Gulf')
xlim([0 35]);
ylim([0 3]);
grid on
%box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,4)
X = mean(sites_08_sal_all,1);
Y_noisy = sites_08_sens_all;
plot_distribution_prctile(X,Y_noisy,'Prctile',[25 50 75 90]);
title('Lower Mississippi')
xlim([0 35]);
ylim([0 3]);
grid on
%box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,5)
X = mean(sites_12_sal_all,1);
Y_noisy = sites_12_sens_all;
plot_distribution_prctile(X,Y_noisy,'Prctile',[25 50 75 90]);
title('Texas-Gulf')
xlim([0 35]);
ylim([0 3]);
grid on
%box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,6)
plot(mean(sites_13_sal_all,1),mean(sites_13_sens_all,1),'LineWidth',1)
title('Rio Grande')
xlim([0 35]);
ylim([0 3]);
grid on
%box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,7)
X = mean(sites_18_sal_all,1);
Y_noisy = sites_18_sens_all;
plot_distribution_prctile(X,Y_noisy,'Prctile',[25 50 75 90]);
title('California')
xlim([0 35]);
ylim([0 3]);
grid on
%box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,8)
X = mean(sites_17_sal_all,1);
Y_noisy = sites_17_sens_all;
plot_distribution_prctile(X,Y_noisy,'Prctile',[25 50 75 90]);
title('Pacific Northwest')
xlim([0 35]);
ylim([0 3]);
grid on
%box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
subplot(3,3,9)
X = mean(sites_19_sal_all,1);
Y_noisy = sites_19_sens_all;
plot_distribution_prctile(X,Y_noisy,'Prctile',[25 50 75 90]);
title('Alaska')
xlim([0 35]);
ylim([0 3]);
grid on
%box on
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';

%% Old figures
% %%
% %1/Beta DIC
% %%%%%NEED TO CHANGE SALINITY ON X AXIS FOR EACH
% %REGION!!!!!!!!!!!!!!
% figure
% plot(cc_estuary_all.(names_allcorrect{1})(:,58),mean(sites_01_sens_all,1)./mean(sites_01_h_all,1)-(sens_nat_avg(:,3)./cc_nat_avg(:,33)),'LineWidth',1)
% hold on
% plot(cc_estuary_all.(names_allcorrect{1})(:,58),mean(sites_02_sens_all,1)./mean(sites_02_h_all,1)-(sens_nat_avg(:,3)./cc_nat_avg(:,33)),'LineWidth',1)
% plot(cc_estuary_all.(names_allcorrect{1})(:,58),mean(sites_03_sens_all,1)./mean(sites_03_h_all,1)-(sens_nat_avg(:,3)./cc_nat_avg(:,33)),'LineWidth',1)
% plot(cc_estuary_all.(names_allcorrect{1})(:,58),mean(sites_08_sens_all,1)./mean(sites_08_h_all,1)-(sens_nat_avg(:,3)./cc_nat_avg(:,33)),'LineWidth',1)
% plot(cc_estuary_all.(names_allcorrect{1})(:,58),mean(sites_12_sens_all,1)./mean(sites_12_h_all,1)-(sens_nat_avg(:,3)./cc_nat_avg(:,33)),'LineWidth',1)
% plot(cc_estuary_all.(names_allcorrect{1})(:,58),mean(sites_13_sens_all,1)./mean(sites_13_h_all,1)-(sens_nat_avg(:,3)./cc_nat_avg(:,33)),'LineWidth',1)
% plot(cc_estuary_all.(names_allcorrect{1})(:,58),mean(sites_18_sens_all,1)./mean(sites_18_h_all,1)-(sens_nat_avg(:,3)./cc_nat_avg(:,33)),'LineWidth',1)
% plot(cc_estuary_all.(names_allcorrect{1})(:,58),mean(sites_17_sens_all,1)./mean(sites_17_h_all,1)-(sens_nat_avg(:,3)./cc_nat_avg(:,33)),'LineWidth',1)
% plot(cc_estuary_all.(names_allcorrect{1})(:,58),mean(sites_19_sens_all,1)./mean(sites_19_h_all,1)-(sens_nat_avg(:,3)./cc_nat_avg(:,33)),'LineWidth',1)
% ylabel('1/\beta_D_I_C')
% grid on
% xlabel('Salinity')
% legend('United States','New England','Mid-Atlantic','South Atlantic-Gulf','Lower Mississippi','Texas-Gulf','Rio Grande','California','Pacific Northwest','Alaska')
% set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
% fig = gcf
% fig.Color='w';
% %%
% figure
% plot(cc_nat_avg(:,58),sens_nat_avg(:,3)-sens_nat_avg(:,3),'k','LineWidth',3)
% hold on
% plot(mean(sites_01_sal_all,1),mean(sites_01_sens_all,1)-sens_nat_avg(:,3),'LineWidth',1)
% plot(mean(sites_02_sal_all,1),mean(sites_02_sens_all,1)-sens_nat_avg(:,3),'LineWidth',1)
% plot(mean(sites_03_sal_all,1),mean(sites_03_sens_all,1)-sens_nat_avg(:,3),'LineWidth',1)
% plot(mean(sites_08_sal_all,1),mean(sites_08_sens_all,1)-sens_nat_avg(:,3),'LineWidth',1)
% plot(mean(sites_12_sal_all,1),mean(sites_12_sens_all,1)-sens_nat_avg(:,3),'LineWidth',1)
% plot(mean(sites_13_sal_all,1),mean(sites_13_sens_all,1)-sens_nat_avg(:,3),'LineWidth',1)
% plot(mean(sites_18_sal_all,1),mean(sites_18_sens_all,1)-sens_nat_avg(:,3),'LineWidth',1)
% plot(mean(sites_17_sal_all,1),mean(sites_17_sens_all,1)-sens_nat_avg(:,3),'LineWidth',1)
% plot(mean(sites_19_sal_all,1),mean(sites_19_sens_all,1)-sens_nat_avg(:,3),'LineWidth',1)
% ylabel('\Delta[H^+]/\Delta[DIC] (nmol/\mumol)')
% grid on
% xlabel('Salinity')
% legend('United States','New England','Mid-Atlantic','South Atlantic-Gulf','Lower Mississippi','Texas-Gulf','Rio Grande','California','Pacific Northwest','Alaska')
% set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
% fig = gcf
% fig.Color='w';

%% Calculate mean sensitivies for all river and ocean end-members; produce scatterplot

sens_river_mean = [];
sens_ocean_mean = [];
for i = 1:length(names_allcorrect)
    sens_ocean_mean(i) = mean(sens_ocean.(names_allcorrect{i})(:,3));
    sens_river_mean(i) = mean(sens_sites_org.(names_allcorrect{i})(:,3));
end

figure
scatter(sens_ocean_mean,sens_river_mean)
ylabel('River \DeltaH^+/\DeltaDIC')
xlabel('Ocean \DeltaH^+/\DeltaDIC')

%%
% Scatter H+ vs sens, by region
figure
scatter(mean(sites_01_h_all,2),mean(sites_01_sens_all,2),'LineWidth',4)
hold on
scatter(mean(sites_02_h_all,2),mean(sites_02_sens_all,2),'LineWidth',2)
scatter(mean(sites_03_h_all,2),mean(sites_03_sens_all,2),'LineWidth',2)
scatter(mean(sites_08_h_all,2),mean(sites_08_sens_all,2),'LineWidth',2)
scatter(mean(sites_12_h_all,2),mean(sites_12_sens_all,2),'LineWidth',2)
scatter(mean(sites_13_h_all,2),mean(sites_13_sens_all,2),'LineWidth',2)
scatter(mean(sites_18_h_all,2),mean(sites_18_sens_all,2),'LineWidth',2)
scatter(mean(sites_17_h_all,2),mean(sites_17_sens_all,2),'LineWidth',2)
scatter(mean(sites_19_h_all,2),mean(sites_19_sens_all,2),'LineWidth',2)
ylabel('Mean estuary sensitivity (\Delta[H^+]/\Delta[DIC] (nmol/\mumol))')
grid on
xlabel('Mean estuary acidity ([H+])')
legend('New England','Mid-Atlantic','South Atlantic-Gulf','Lower Mississippi','Texas-Gulf','Rio Grande','California','Pacific Northwest','Alaska')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';

%%
% Scatter pH_correct vs sens, by region
figure
scatter(mean(sites_01_pH_correct_all,2),mean(sites_01_sens_all,2),'LineWidth',2)
hold on
scatter(mean(sites_02_pH_correct_all,2),mean(sites_02_sens_all,2),'LineWidth',2)
scatter(mean(sites_03_pH_correct_all,2),mean(sites_03_sens_all,2),'LineWidth',2)
scatter(mean(sites_08_pH_correct_all,2),mean(sites_08_sens_all,2),'LineWidth',2)
scatter(mean(sites_12_pH_correct_all,2),mean(sites_12_sens_all,2),'LineWidth',2)
scatter(mean(sites_13_pH_correct_all,2),mean(sites_13_sens_all,2),'LineWidth',2)
scatter(mean(sites_18_pH_correct_all,2),mean(sites_18_sens_all,2),'LineWidth',2)
scatter(mean(sites_17_pH_correct_all,2),mean(sites_17_sens_all,2),'LineWidth',2)
scatter(mean(sites_19_pH_correct_all,2),mean(sites_19_sens_all,2),'LineWidth',2)
ylabel('Mean estuary sensitivity (\Delta[H^+]/\Delta[DIC] (nmol/\mumol))')
grid on
xlabel('Mean estuary pH')
legend('New England','Mid-Atlantic','South Atlantic-Gulf','Lower Mississippi','Texas-Gulf','Rio Grande','California','Pacific Northwest','Alaska')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
%%
% Scatter coastal stream pH_correct vs sens, by region
figure
scatter(sites_01_pH_correct_all(:,1),sites_01_sens_all(:,1),'LineWidth',2)
hold on
scatter(sites_02_pH_correct_all(:,1),sites_02_sens_all(:,1),'LineWidth',2)
scatter(sites_03_pH_correct_all(:,1),sites_03_sens_all(:,1),'LineWidth',2)
scatter(sites_08_pH_correct_all(:,1),sites_08_sens_all(:,1),'LineWidth',2)
scatter(sites_12_pH_correct_all(:,1),sites_12_sens_all(:,1),'LineWidth',2)
scatter(sites_13_pH_correct_all(:,1),sites_13_sens_all(:,1),'LineWidth',2)
scatter(sites_18_pH_correct_all(:,1),sites_18_sens_all(:,1),'LineWidth',2)
scatter(sites_17_pH_correct_all(:,1),sites_17_sens_all(:,1),'LineWidth',2)
scatter(sites_19_pH_correct_all(:,1),sites_19_sens_all(:,1),'LineWidth',2)
ylabel('Stream sensitivity (\Delta[H^+]/\Delta[DIC] (nmol/\mumol))')
grid on
xlabel('Stream pH')
legend('New England','Mid-Atlantic','South Atlantic-Gulf','Lower Mississippi','Texas-Gulf','Rio Grande','California','Pacific Northwest','Alaska')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';

figure
scatter(sites_01_alk_all(:,1),sites_01_sens_all(:,1),'LineWidth',2)
hold on
scatter(sites_02_alk_all(:,1),sites_02_sens_all(:,1),'LineWidth',2)
scatter(sites_03_alk_all(:,1),sites_03_sens_all(:,1),'LineWidth',2)
scatter(sites_08_alk_all(:,1),sites_08_sens_all(:,1),'LineWidth',2)
scatter(sites_12_alk_all(:,1),sites_12_sens_all(:,1),'LineWidth',2)
scatter(sites_13_alk_all(:,1),sites_13_sens_all(:,1),'LineWidth',2)
scatter(sites_18_alk_all(:,1),sites_18_sens_all(:,1),'LineWidth',2)
scatter(sites_17_alk_all(:,1),sites_17_sens_all(:,1),'LineWidth',2)
scatter(sites_19_alk_all(:,1),sites_19_sens_all(:,1),'LineWidth',2)
ylabel('Stream H^ sensitivity(nmol/\mumol))')
grid on
xlabel('Stream alkalinity (\mueq kg^-^1)')
legend('New England','Mid-Atlantic','South Atlantic-Gulf','Lower Mississippi','Texas-Gulf','Rio Grande','California','Pacific Northwest','Alaska')
set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
fig = gcf
fig.Color='w';
% %%
% % Scatter pH_correct vs sens, anomalies from national average
% figure
% scatter(mean(sites_01_pH_correct_all,2)-mean(cc_nat_avg(:,43)),mean(sites_01_sens_all,2)-mean(sens_nat_avg(:,3)),'LineWidth',2)
% hold on
% scatter(mean(sites_02_pH_correct_all,2)-mean(cc_nat_avg(:,43)),mean(sites_02_sens_all,2)-mean(sens_nat_avg(:,3)),'LineWidth',2)
% scatter(mean(sites_03_pH_correct_all,2)-mean(cc_nat_avg(:,43)),mean(sites_03_sens_all,2)-mean(sens_nat_avg(:,3)),'LineWidth',2)
% scatter(mean(sites_08_pH_correct_all,2)-mean(cc_nat_avg(:,43)),mean(sites_08_sens_all,2)-mean(sens_nat_avg(:,3)),'LineWidth',2)
% scatter(mean(sites_12_pH_correct_all,2)-mean(cc_nat_avg(:,43)),mean(sites_12_sens_all,2)-mean(sens_nat_avg(:,3)),'LineWidth',2)
% scatter(mean(sites_13_pH_correct_all,2)-mean(cc_nat_avg(:,43)),mean(sites_13_sens_all,2)-mean(sens_nat_avg(:,3)),'LineWidth',2)
% scatter(mean(sites_18_pH_correct_all,2)-mean(cc_nat_avg(:,43)),mean(sites_18_sens_all,2)-mean(sens_nat_avg(:,3)),'LineWidth',2)
% scatter(mean(sites_17_pH_correct_all,2)-mean(cc_nat_avg(:,43)),mean(sites_17_sens_all,2)-mean(sens_nat_avg(:,3)),'LineWidth',2)
% scatter(mean(sites_19_pH_correct_all,2)-mean(cc_nat_avg(:,43)),mean(sites_19_sens_all,2)-mean(sens_nat_avg(:,3)),'LineWidth',2)
% ylabel('Mean estuarine sensitivity anomaly (\Delta[H^+]/\Delta[DIC] (nmol/\mumol))')
% grid on
% xlabel('Mean estuarine pH anomaly')
% legend('New England','Mid-Atlantic','South Atlantic-Gulf','Lower Mississippi','Texas-Gulf','Rio Grande','California','Pacific Northwest','Alaska')
% set(findobj(gcf,'type','axes'),'FontName','Times','FontSize',12,'LineWidth',1)
% fig = gcf
% fig.Color='w';

    
    