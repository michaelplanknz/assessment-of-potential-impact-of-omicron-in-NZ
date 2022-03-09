function par = getParOmiWane(istartDate, iTransReduc, iboostProp, ...
    ivaxEff, igenScenario, itEnd, iborderSeeds, iborderTransRed, iisolEff, ipTrace, iisolEffCT, ipTestClin, icmAdjustBool)

% Function that outputs a table of parameters that can then be fed to
% "runSimWaning"

% INPUTS:
% - istartDate: first date of simulation, in datenum format
% - iTransScenario: relative transmission scenario (1, 2, or 3), each corresponding to a different control measure setting
% - iboostProp: maximum booster coverage reached
% - ivaxEff: vaccine efficacy setting (1 - baseline, 2 - low)
% - itEnd: number of days simulated
% - iborderSeeds: number of weekly border cases arriving from "borderOpenDate"
% - iborderTransRed: transmission reduction from isolation of border cases
% - iisolEff: transmission reduction for isolated community cases
% - ipTrace: proportion of infected contacts found via contact tracing
% - iisolEffCT: transmission reduction for isolated contacts
% - ipTestClin: probability of testing for clinical (symptomatic) cases

% OUTPUT:
% - par: table of parameter values


%------------- Seed and control Parameters --------------

par.tEnd = itEnd; % 200
par.date0 = istartDate;
par.borderOpenDate = max(par.date0 + 1, datenum("28FEB2022")); % Only used if border seeds ("wkBorderSeeds" in the main) != 0

% Relative transmission (R multiplier corresponding to control measures)
par.relTransBaseAL = (1 - iTransReduc) * ones(1, par.tEnd+1);

% Generation time distribution, 2 by default, see below for details
R0 = [4.3, 3.3, 2.7];
par.R0 = R0(igenScenario);


par.minDetectTime = 0;    % time of outbreak detection
par.followUpTime = 7;    % cases with an isolation time prior to detection are distributed over this time period post-detection

% Community seed cases
par.relInfSeedCases = 1;
par.seedPeriod = 7;
par.meanSeedsPerDay = 500/7;
par.genSeedVaxStatus = @genSeedVaxStatus_cat2d;

% Border seed cases
par.borderSeedsPerDay = iborderSeeds/7;
par.borderSeedPeriod = par.tEnd - (par.borderOpenDate - par.date0);
par.relInfBorderCases = 1 - iborderTransRed;



par.maxActiveCases = inf;     % 500
par.minActiveCases = -1;
par.alarmActiveCases = inf;
par.maxHospBeds = inf;        % 500
par.minHospBeds = -1;
par.alarmHospBeds = inf;
par.tierUp = 4;
par.tierDown = 1;


%------------- Branching Process Parameters --------------

par.cSub = 0.5; % Relative infectiousness of subclinicals
par.ssk = 0.5; % Overdispersion/superspreading parameter k
par.maxInfectTime = 14; % Maximum infectious perod (for computational efficiency)

if igenScenario == 1 % long GI scenario, high transmission
    par.genA = 5.665; par.genB = 2.826; % Generation time distribution parameters
    par.incA = 5.8; par.incB = 0.95; % Exposure to Onset distribution parameters
elseif igenScenario == 2 % baseline GI scenario, medium transmission
    par.genA = 3.7016; par.genB = 2.826; % Generation time distribution parameters
    par.incA = 5.8; par.incB = 0.62; % Exposure to Onset distribution parameters
elseif igenScenario == 3 % short GI scenario, low transmission
    par.genA = 3; par.genB = 2.826; % Generation time distribution parameters
    par.incA = 5.8; par.incB = 0.5; % Exposure to Onset distribution parameters
end

par.isolA = 1; par.isolB = 4; % Onset to isolation distribution parameters
par.traceA = 3; par.traceB = 3/3; % Parent isolation to quarantine distribution parameters
par.hospA = 1; par.hospB = 5; % Onset to hospitalisation distribution parameters
par.losA = 1; par.losB = 4; % Hospital LOS distribution parameters

par.pTestClin = ipTestClin; % 0.3; % Probability of detecting symptomatic cases
par.pTestSub = 0;  % Probability of detecting subclinical cases
par.pTrace = ipTrace; %0.25; % Probability of detecting a case by contact tracing
par.traceCapacity = inf;        % 1000/10   (now tracing capacity for average daily cases over last 7 days)

par.cIsol = 1 - iisolEff; %0.2;
par.cQuar = 1 - iisolEffCT; %0.5;


%------------- Vaccine Effectiveness Parameters --------------
% after 0, 1, 2a, ... doses
if ivaxEff == "original"
    % Vaccine effectiveness from United Kingdom Health Security Agency, 2022 #326;Keeling, 2022 #323United Kingdom Health Security Agency, 2022 #326;Keeling, 2022 #323
    % [V0    V1    V2a V2b V2c V2d    V3a V3b V3c V3d]
    par.VEi = [0    0    0.62 0.55 0.4 0.28    0.64 0.57 0.47 0.4]';
    par.VEs = [0    0    0 0 0 0               0 0 0 0]';
    par.VEt = [0    0    0 0 0 0               0 0 0 0]';
    par.VEd = [0    0    0.55 0.56 0.58 0.61   0.78 0.72 0.68 0.65]';
    par.VEm = [0    0    0.78 0.78 0.79 0.81   0.89 0.86 0.84 0.83]';
elseif ivaxEff == "low"
    % Low vaccine effectiveness from Nick Golding's model
    par.VEi = [0    0    0.37	0.24	0.18	0.14	0.58	0.43	0.35	0.28]';
    par.VEs = [0    0    0.10	0.07	0.05	0.04	0.15	0.11	0.09	0.07]';
    par.VEt = [0    0    0.12	0.06	0.04	0.03	0.26	0.16	0.11	0.08]';
    par.VEd = [0    0    0.67	0.56	0.50	0.43	0.78	0.70	0.65	0.59]';
    par.VEm = [0    0    0.65	0.55	0.48	0.42	0.77	0.69	0.64	0.58]';
end


%------------- Disease Rate Data --------------
par.nAgeGroups = 16;
par.IDR = [0.5440, 0.5550, 0.5770, 0.5985, 0.6195, 0.6395, 0.6585, 0.6770, 0.6950, 0.7117, 0.7272, 0.7418, 0.7552, 0.7680, 0.7800, 0.8008]'; % Fraser group
par.ui = [0.4000, 0.3950, 0.3850, 0.4825, 0.6875, 0.8075, 0.8425, 0.8450, 0.8150, 0.8050, 0.8150, 0.8350, 0.8650, 0.8450, 0.7750, 0.7400 ];    % Davies relative suscewptibility
par.pSub = 1-par.IDR;

[par.IHR, par.pICU, par.IFR] = getHerreraRatesOmi();



%------------- Specify Population Structure --------------
fs = 'nzpopdist.xlsx';
popSizeData = readmatrix(fs); % Load NZ population structure from data folder

par.popCount = zeros(par.nAgeGroups, 1); % Create popDist vector
par.popCount(1:par.nAgeGroups-1) = popSizeData(1:par.nAgeGroups-1, 2); % Fill entries with population distribution
par.popCount(par.nAgeGroups) = sum(popSizeData(par.nAgeGroups:end, 2)); % Aggregate 75+ age-groups
par.totalPopSize = sum(par.popCount);
par.popDist = par.popCount/sum(par.popCount);

par.age = (2.5:5:77.5)'; % Define final age groups (matching contact matrix)


%------------- Get next generation matrix ---------

% Parameters to adjusts contact matrix to reflect younger age groups 
% interacting more
par.contactPar.adjustContacts = icmAdjustBool;
par.contactPar.diagBlockIndices = {[1:3], [4:5], [6:7], [8:12], [13:16]};
par.contactPar.numDiagBlocks = length(par.contactPar.diagBlockIndices);
par.contactPar.diagBlockMultipliers = {1.4, 1.35, 1.1, 0.15, 0.1};
par.contactPar.offBlockMultipliers = {[0.9,0.8,0.5,0.3]; [0.8,0.5,0.3]; [0.5,0.3]; [0.3]}; % from to ; from to

[par.NGM, par.NGMclin] = getNGM(par);


%------------- Define vaccine coverage either static or time-dependent ---------
par.boostProp = iboostProp;
par.vaccImmDelay = 14;     % delay from vaccination to immunity
fNamePart = "_national";

% Creates matrices for 1st, 2nd and 3rd dose coverage by time and age group, starting
% on date0
[~, par.cov1, par.cov2a, par.cov2b, par.cov2c, par.cov2d, par.cov3a, par.cov3b, par.cov3c, par.cov3d] = getSegmentedVaccineData(fNamePart, par);


end




