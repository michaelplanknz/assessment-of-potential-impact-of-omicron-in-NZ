function [dailyinf, dailycases, dailyhosp, cdeaths, hospbeds] = readOmiTimeseries(savefolder, scenario_letter, transSc, boostprop, boost_scenario)


filename_dailyinf = append('timeseries/', savefolder, '/scenario', scenario_letter, '_dailyinf_', transSc, 'Trans_', boostprop ,'boost', boost_scenario, '.csv');
dailyinf = csvread(filename_dailyinf);

filename_dailycases = append('timeseries/', savefolder, '/scenario', scenario_letter, '_dailycases_', transSc, 'Trans_', boostprop ,'boost', boost_scenario, '.csv');
dailycases = csvread(filename_dailycases);

filename_dailyhosp = append('timeseries/', savefolder, '/scenario', scenario_letter, '_dailyhosp_', transSc, 'Trans_', boostprop ,'boost', boost_scenario, '.csv');
dailyhosp = csvread(filename_dailyhosp);

filename_cdeaths = append('timeseries/', savefolder, '/scenario', scenario_letter, '_cumuldeaths_', transSc, 'Trans_', boostprop ,'boost', boost_scenario, '.csv');
cdeaths = csvread(filename_cdeaths);

filename_hospbeds = append('timeseries/', savefolder, '/scenario', scenario_letter, '_hospbeds_', transSc, 'Trans_', boostprop ,'boost', boost_scenario, '.csv');
hospbeds = csvread(filename_hospbeds);



end