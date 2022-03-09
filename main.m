%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   BRANCHING PROCESS MODEL
%      OMICRON MODEL WITH VAX BOOSTERS AND WANING IMMUNITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc


%_________________________ SIMULATION PARAMETERS __________________________

nReps = 10;  % Number of simulation repetitions

% Paths for dependencies
addpath('data');
addpath('dependencies');

% Folder where timeseries will be saved (make sure to create a "timeseries"
% folder, and within it a folder with the same name as savefolder
savefolder = "baseline";

writesum = 1; % If =1, a summary table will be saved in the "summary" folder
plots = 1; % If = 1, plots will appear at the end of the simulation
saveplots = 1; % If = 1, plots will be saved in the "plots" folder

% Simulation start dates
dates = ["01JAN2022", "01FEB2022", "01MAR2022", "01APR2022", "01MAY2022"];
startDate = datenum(dates);
nChange = length(startDate);

% Simulation duration (days)
itEnd = 200; %200

% Proportion of eligible population getting boosted
boostProp = [0.9]; %[0.9, 0.7];

% Number of weekly border cases
wkBorderSeeds = [0]; %[65, 650];

% Border cases transmission reduction (from isolation at arrival)
BC_transReduc = [0]; %[0.8, 0.25];

% Transmission reduction corresponding to the effect of control measures
transReduc = [0.2]; % Baseline = 0.2

% Vaccination efficacy ["original", "low"];
vaxEffn = [1]; % Baseline = 1 (original)

% Transmission level tested (1 = high, 2 = med, 3 = low)
genScenario = [3, 2, 1]; % Baseline [1, 2, 3]

% Community cases isolation effectiveness
isolEff = [0.8]; % Baseline = 0.8

% Probability of tracing contacts of confirmed cases
pTrace = [0.25]; % Baseline = 0.25

% Contacts isolation effectiveness
isolEffCT = [0.5]; % Baseline = 0.5

% Probability of detecting clinical/symptomatic cases
pTestClin = [0.3]; % Baseline = 0.3

% If =1, contact matrix gets adjusted to reflect younger age groups interacting more
cmAdjustBool = [0]; % Baseline = 0

% Scenarios combinations
scenarios_combinations = combvec(startDate', boostProp, transReduc, vaxEffn, ...
    genScenario, wkBorderSeeds, BC_transReduc, isolEff, pTrace, isolEffCT, pTestClin, cmAdjustBool);
[nvar, nscenarios] = size(scenarios_combinations);


earlyReject.tData = [];
earlyReject.nCasesData = [];
earlyReject.threshold = inf;
earlyReject.rejectElimFlag = 1;


%___________________________ SCENARIOS LOOP _______________________________
parfor jj = 1:nscenarios
    
    scenarios = scenarios_combinations;
    scenario_letter = char(ceil(jj/(length(genScenario) * size(dates, 2))) + 64);
    
    vaxEff = ["original", "low"];
    transSc = ["high", "med", "low"];
    
    idate0 = scenarios(1, jj);
    iboostProp = scenarios(2, jj);
    itransReduc = scenarios(3, jj);
    ivaxEff = vaxEff(scenarios(4, jj));
    igenScenario = scenarios(5, jj);
    iborderSeeds = scenarios(6, jj);
    iborderTransRed = scenarios(7, jj);
    iisolEff = scenarios(8, jj);
    ipTrace = scenarios(9, jj);
    iisolEffCT = scenarios(10, jj);
    ipTestClin = scenarios(11, jj);
    icmAdjustBool = scenarios(12, jj);
    
    % Timeseries matrices initialisation
    [dailyinfTS, dailycasesTS, dailyhospTS, cdeathsTS, hospbedsTS] = ...
        deal(zeros(nReps, itEnd));
    
    [dailycasesTS_0d, dailycasesTS_1d, dailycasesTS_2d, dailycasesTS_3d, ...
        dailyhospTS_0d, dailyhospTS_1d, dailyhospTS_2d, dailyhospTS_3d] = ...
        deal(zeros(itEnd, 16, nReps));
    
    % Get sim parameters
    par = getParOmiWane(idate0, itransReduc, iboostProp, ivaxEff, igenScenario, itEnd, iborderSeeds, ...
        iborderTransRed, iisolEff, ipTrace, iisolEffCT, ipTestClin, icmAdjustBool);
    
    % Sim reps loop:
    for iRep = 1:nReps
        fprintf("Parameter combination %i of %i (scenario %s, %s transmission), rep %i of %i\n", jj, nscenarios, scenario_letter, transSc(igenScenario), iRep, nReps)
        % Run simulation
        [cases] = runSimWaning(par, earlyReject);
        % Process cases data
        [nInf, nIsol, cases_0dose, cases_1dose, cases_2dose, cases_3dose, ...
            nHosp, dailyHosp_0dose, dailyHosp_1dose, dailyHosp_2dose, dailyHosp_3dose, ...
            nDisc, nICUIn, nICUOut, nDeaths] = postProcess(cases, par);
        
        dailyinfTS(iRep, :) = sum(nInf(1:par.tEnd, :), 2);
        dailycasesTS(iRep, :) = sum(nIsol(1:par.tEnd, :), 2);
        dailyhospTS(iRep, :) = sum(nHosp(1:par.tEnd, :), 2);
        cdeathsTS(iRep, :) = cumsum(sum(nDeaths(1:par.tEnd, :), 2), 1);
        hospbedsTS(iRep, :) = cumsum(sum(nHosp(1:par.tEnd, :)-nDisc(1:par.tEnd, :), 2), 1);
        
        dailycasesTS_0d(:, :, iRep) = cases_0dose(1:par.tEnd, :);
        dailycasesTS_1d(:, :, iRep) = cases_1dose(1:par.tEnd, :);
        dailycasesTS_2d(:, :, iRep) = cases_2dose(1:par.tEnd, :);
        dailycasesTS_3d(:, :, iRep) = cases_3dose(1:par.tEnd, :);
        dailyhospTS_0d(:, :, iRep) = dailyHosp_0dose(1:par.tEnd, :);
        dailyhospTS_1d(:, :, iRep) = dailyHosp_1dose(1:par.tEnd, :);
        dailyhospTS_2d(:, :, iRep) = dailyHosp_2dose(1:par.tEnd, :);
        dailyhospTS_3d(:, :, iRep) = dailyHosp_3dose(1:par.tEnd, :);
    end
    
    % Write timeseries results, one line per repetition, one column per day:
    writeOmiTimeseries(savefolder, scenario_letter, dailyinfTS, dailycasesTS, dailyhospTS, cdeathsTS, ...
        hospbedsTS, transSc(igenScenario), int2str(100*iboostProp), datestr(idate0, "ddmmm"))
   
    
end

%__________________________ Summary writing _______________________________

if writesum == 1
    writeOmiSummary(savefolder, scenarios_combinations, genScenario)
end



%_______________________________ Plots ____________________________________

if plots == 1
    plot_scenario = savefolder;
    trans_scenarios_toplot = ["low", "med", "high"];
    boostProp_toplot = ["90"];
    scenario_letter = "A";
    for plot_transSc = trans_scenarios_toplot
        for plot_boostprop = boostProp_toplot
            plotOmiTS(plot_scenario, plot_transSc, plot_boostprop, scenario_letter, 1)
        end
    end
end

