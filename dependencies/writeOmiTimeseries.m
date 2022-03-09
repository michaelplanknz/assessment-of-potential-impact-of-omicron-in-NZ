function [] = writeOmiTimeseries(savefolder, scenario_letter, daily_inf, daily_cases, daily_hosp, cdeaths, hosp_beds, transSc, boostprop, boost_scenario)

% Example output filename = scenarioA_dailyinf_lowTrans_70boost1Jan.csv

filename_dailyinf = append('timeseries/', savefolder, '/scenario', scenario_letter, '_dailyinf_', transSc, 'Trans_', boostprop ,'boost', boost_scenario, '.csv');
writematrix(daily_inf, filename_dailyinf)

filename_dailycases = append('timeseries/', savefolder, '/scenario', scenario_letter, '_dailycases_', transSc, 'Trans_', boostprop ,'boost', boost_scenario, '.csv');
writematrix(daily_cases, filename_dailycases)

filename_dailyhosp = append('timeseries/', savefolder, '/scenario', scenario_letter, '_dailyhosp_', transSc, 'Trans_', boostprop ,'boost', boost_scenario, '.csv');
writematrix(daily_hosp, filename_dailyhosp)

filename_cdeaths = append('timeseries/', savefolder, '/scenario', scenario_letter, '_cumuldeaths_', transSc, 'Trans_', boostprop ,'boost', boost_scenario, '.csv');
writematrix(cdeaths, filename_cdeaths)

filename_hospbeds = append('timeseries/', savefolder, '/scenario', scenario_letter, '_hospbeds_', transSc, 'Trans_', boostprop ,'boost', boost_scenario, '.csv');
writematrix(hosp_beds, filename_hospbeds)


end