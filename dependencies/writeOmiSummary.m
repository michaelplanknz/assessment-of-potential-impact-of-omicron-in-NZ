function [] = writeOmiSummary(scenarioname, scenarios, genScenario)

nscenarios = size(scenarios, 2);

% Summary matrix initialissation
[date0, scenario] = deal(strings(nscenarios, 1));
[R0scenario, boostCov, infections, cases, hosp, peakBeds, deaths, ...
    maxPrev10, maxAvgCases7, PCdate, PBdate] = deal(zeros(nscenarios, 1));
summary = table(scenario, R0scenario, boostCov, date0, infections, ...
    cases, hosp, peakBeds, deaths, maxPrev10, maxAvgCases7, PCdate, PBdate);

for jj = 1:nscenarios
    
    scenario_letter = char(ceil(jj/(length(genScenario) * length(unique(scenarios(1, :))))) + 64);
    
    genSc = ["high", "med", "low"];
    R0s = [4.3, 3.3, 2.7];

    idate0 = scenarios(1, jj);
    iboostProp = scenarios(2, jj);
    igenScGI = scenarios(5, jj);
    
    [dailyinfTS, dailycasesTS, dailyhospTS, cdeathsTS, hospbedsTS] = ...
        readOmiTimeseries(scenarioname, scenario_letter, genSc(igenScGI), int2str(100*iboostProp), datestr(idate0, "ddmmm"));
    
     [~, pcmaxdate] = max(median(dailycasesTS, 1));
     [~, pbmaxdate] = max(median(hospbedsTS, 1));
    
%     t1 = 14;
%     t2 = 28;
%     [pf, ~] = polyfit(t1:t2, log(mean(dailyinfTS(:, t1:t2), 1)), 1);
%     dt = log(2) / pf(1);
    summary(jj, :) = {scenario_letter, R0s(igenScGI), iboostProp, ...
        datestr(idate0, "ddmmm"),...
        round(median(sum(dailyinfTS, 2))), ...
        round(median(sum(dailycasesTS, 2))), ...
        round(median(sum(dailyhospTS, 2))), ...
        round(median(max(hospbedsTS, [], 2))), ...
        round(median(cdeathsTS(:, end))), ...
        round(median(max(movsum(dailyinfTS, [0, 9], 2), [], 2)) / 5112280, 2), ...
        round(median(max(movsum(dailycasesTS, [0, 6], 2)./7, [], 2))), ...
        pcmaxdate, pbmaxdate};
    writetable(summary, append("summaries/", scenarioname, "_summary_", datestr(now,'yyyymmdd'), ".csv"))
end

end