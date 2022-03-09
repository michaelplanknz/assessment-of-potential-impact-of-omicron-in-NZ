function [] = plotOmiTS(scenario, transScen, boostprop, scenario_letter, saveplot)


newcolors = [0 0.4470 0.7410; 0.4660 0.6740 0.1880; 0.9290 0.6940 0.1250; 0.8500 0.3250 0.0980; 0.9 0 0.9];
ncolors = {'#009E73', '#000000', '#D55E00', '#E69F00', '#56B4E9'};
for c=1:5
    col = ncolors{c};
    newcolors(c, :) = sscanf(col(2:end),'%2x%2x%2x',[1 3])/255;
end

% Booster scenarios
b_scenarios = ["01Jan", "01Feb", "01Mar", "01Apr", "01May"];
b_scenarios_leg = ["1 Jan", "1 Feb", "1 Mar", "1 Apr", "1 May"];

f = figure;
f.Position = [300 300 800 500];
tl = tiledlayout(2,2);
title(tl, append(transScen, " transmission"));

%%%% Plot daily reported cases
nexttile
title("(a)")
hold on
for bs = 1:length(b_scenarios)
    for transSc = transScen
        [~, dailycases, ~, ~, ~] = readOmiTimeseries(scenario, scenario_letter, transSc, boostprop, b_scenarios(bs));
        t = 1:size(dailycases, 2);
        plot(t, median(dailycases, 1))
        [~, maxi] = max(median(dailycases, 1));
        fprintf("%s trans, %s%% boost cov, %s scenario, cases peak on %s, after %s days\n", transSc, boostprop, b_scenarios(bs), datetime(b_scenarios(bs), 'Format', 'ddMMM') + maxi, int2str(maxi))
    end
end
hold off
leg = legend(b_scenarios_leg, 'Location', 'northeast');
title(leg,'Start date')
xlabel("time since outbreak start (days)")
ylabel("daily reported cases")

%%%% Plot daily hospitalisations
nexttile
title("(b)")
hold on
for bs = 1:length(b_scenarios)
    for transSc = transScen
        [~, ~, dailyhosp, ~, ~] = readOmiTimeseries(scenario, scenario_letter, transSc, boostprop, b_scenarios(bs));
        t = 1:size(dailycases, 2);
        plot(t, median(dailyhosp, 1))
    end
end
hold off
xlabel("time since outbreak start (days)")
ylabel("new daily hospitalisations")

%%%% Plot cumulative deaths
nexttile
title("(c)")
hold on
for bs = 1:length(b_scenarios)
    for transSc = transScen
        [~, ~, ~, cdeaths, ~] = readOmiTimeseries(scenario, scenario_letter, transSc, boostprop, b_scenarios(bs));
        t = 1:size(dailycases, 2); 
        plot(t, median(cdeaths, 1))
    end
end
hold off
xlabel("time since outbreak start (days)")
ylabel("cumulative deaths")

%%%% Plot hospital beds occupied
nexttile
title("(d)")
hold on
for bs = 1:length(b_scenarios)
    for transSc = transScen
        [~, ~, ~, ~, hospbeds] = readOmiTimeseries(scenario, scenario_letter, transSc, boostprop, b_scenarios(bs));
        t = 1:size(dailycases, 2); 
        plot(t, median(hospbeds, 1))
    end
end
hold off
xlabel("time since outbreak start (days)")
ylabel("hospital beds occupied")


colororder(newcolors)

if saveplot == 1
    saveas(f,append('plots/', scenario, '_', transSc, 'Transmission_',  boostprop ,'boost.png'));
end


end