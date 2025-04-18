%% Master script aas of 11/25/2024

neaa_import_241218 %Imports data from R output files
neaa_carbcalc_corrected_240418 % SRP QC 11/25/24 calculates full carbonate systems & sensitivity factors of river stations - uses corrrect pH following Liu et al 2020 but does not correct or ALKorg;
neaa_carbcalc_alkorg_corrected_240530 % SRP QC 11/25/24; calculates full carbonate systems & sensitivity factors of river stations - uses corrrect pH and Alkorg following Liu et al 2020; excludes statins with <1000uM Alk and no DOC data
ocean_carbcalc_231102 % SRP QC 11/25/24; calculates ocean sensitivity factors (DIC and ALK directly measured in CODAP-NA dataset) and subsets with Alkorg correct stations
Fig2_sitemap % SRP QC 11/25/24; Map of all sites color coded by water resource region
map_discharge % Map of all sites sized by stream discharge
map_streamsens % Map of all stream sites by H+ sensitivity
cdr_efficacy % Calculate sensitivity factors with alkalinity change
map_tstrends % Map of all stream site H+ trends
mixing_deltas_240412_omega % Calculate changes in estuarine baseline chemistry Pre/Post 1990
alk_loads_1990 % Calculates changes in alk and dic delivery by coastal streams
mixing_curves_medians % Calculates baseline estuarine chemistry and summary stats, including mean sensitivities
neaa_tables_alksens %Write summaries for each site to an excel file as a table
airsea_240901_means % calculate pre/post 1990 air-sea pCO2 and DIC disequilibria
so_what_plots %looks at scenarios of added DIC and Alk