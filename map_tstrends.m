%% Import T-S H+ trend data (only sites with >10 years of data and >25 observations); alkalinity not corrected
ts_regressions = readmatrix("C:\Users\spacella\OneDrive - Environmental Protection Agency (EPA)\National Estuary Acidification Assessment\NEAA Analysis Files 240814\R Code May 2022 SRP Review\neaa_data_exports/ts_regressions.csv",'OutputType','char');
ts_regressions(:,1)=regexprep(ts_regressions(:,1),'-','');
% Make map of trends (sig +, sig -, sig no, non-sig

regressions_categories = [];
regressions_discharge = [];
% H+ first
for i = 1:height(all_sites_names2)
    if isempty(find(strcmpi(all_sites_names2(i),ts_regressions(:,1)))) == 1 %was the site tested for a trend?
        regressions_categories(i) = 1; % 1 = not tested
        regressions_discharge(i) = discharge_means_numeric(i);
        continue
    end
    index = find(strcmpi(ts_regressions(:,1),all_sites_names2(i))); %find index of site in ts_regressions
    if str2double(ts_regressions(index,9)) > 0.01
        regressions_categories(i) = 2; % 2 = not signficant
        regressions_discharge(i) = discharge_means_numeric(i);
        continue
    end
    if str2double(ts_regressions(index,5)) > 0
        regressions_categories(i) = 3; % 3 = signficant +
        regressions_discharge(i) = discharge_means_numeric(i);
    else
        regressions_categories(i) = 4; % 4 = signficant -
        regressions_discharge(i) = discharge_means_numeric(i);
    end
end
regressions_categories = regressions_categories';
regressions_categories = categorical(regressions_categories);
categories(regressions_categories)
regressions_categories = renamecats(regressions_categories,{'Insufficient Data','Not significant','Increasing acidity','Decreasing acidity'})
regressions_discharge = regressions_discharge';

index_ts = [];
counter = 1;
for i = 1:height(ts_regressions)
    index_temp = find(strcmpi(ts_regressions(i,1),all_sites_names2));
    if index_temp > 0
        index_ts(counter) = index_temp;
        counter = counter+1;
    end
end
index_ts = index_ts';
discharge_ts = discharge_means_numeric(index_ts);
    

figure
gb = geobubble(lat(index_ts),lon(index_ts),discharge_means_numeric(index_ts)*conv_cfs,'MapLayout','maximized','colordata',regressions_categories(index_ts));
gb.SizeLegendTitle = 'Mean discharge (m^3/s)';
%gb.BubbleColorList = [.7 .7 .7; 0 0.4470 0.7410; 0.8350 0.0780 0.1840; 0.4660 0.8740 0.1880];
geobasemap colorterrain

figure
gb = geobubble(lat,lon,discharge_means_numeric*conv_cfs,'MapLayout','maximized','colordata',regressions_categories);
gb.SizeLegendTitle = 'Mean discharge (m^3/s)';
%gb.BubbleColorList = [.7 .7 .7; 0 0.4470 0.7410; 0.8350 0.0780 0.1840; 0.4660 0.8740 0.1880];
geobasemap colorterrain



%calculate % total discharge for each category of trend
index_increasing = find(regressions_categories == {'Increasing'});
index_decreasing = find(regressions_categories == {'Decreasing'});
index_notsig = find(regressions_categories == {'Not significant'});
index_insufficient = find(regressions_categories == {'Insufficient Data'});
discharge_increasing = sum(regressions_discharge(index_increasing));
discharge_decreasing = sum(regressions_discharge(index_decreasing));
discharge_notsig = sum(regressions_discharge(index_notsig));
discharge_insufficient = sum(regressions_discharge(index_insufficient));
discharge_total = sum(discharge_means_numeric);
discharge_percent_increasing = discharge_increasing/discharge_total*100;
discharge_percent_decreasing = discharge_decreasing/discharge_total*100;
discharge_percent_notsig = discharge_notsig/discharge_total*100;
discharge_percent_insufficient = discharge_insufficient/discharge_total*100;


