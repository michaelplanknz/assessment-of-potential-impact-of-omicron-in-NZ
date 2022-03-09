function [cases, dist, relTransCurrentAL, ReffEmp, Rvt] = runSimWaning(par, earlyReject)

% transmission process will stop if this many people have been infected
maxCases = 5100000;

% Set up time array
t = par.date0 + (0:1:par.tEnd);
nSteps = length(t)-1;

% Set the relative transmission rate due to alert level settings, as a function of time
% (this is the base setting and can be later modified by dynamic control rules)
relTransCurrentAL = par.relTransBaseAL;

% calculate the area under the curve of the generation time distribution in 1-day time
% steps - this gives the relative amount of transmitting each person does
% on each day of their infection
tArr = 0:1:par.maxInfectTime+1;
C = wblcdf(tArr, par.genA, par.genB);
auc = diff(C).';
auc = auc/sum(auc); % Renormalise auc so that it definitely sums to 1 :)


% initialise variables for the case table
caseID = (1:maxCases).';
parentID = nan(maxCases, 1);
gen = nan(maxCases, 1);
nOff = zeros(maxCases, 1);
ageGroup = nan(maxCases, 1);
Rimult = genRi(maxCases, 1, par);   % pre generate each individual's values of Yi
subclinFlag = nan(maxCases, 1);
vaccDoses = zeros(maxCases, 1);
tInfect = nan(maxCases, 1);
tOnset = genOnsetDelay(maxCases, 1, par);   % pre-generate each individuals incubation period
tIsol = nan(maxCases, 1);   % testing and isolation times assumed simulteneous
tQuar = nan(maxCases, 1);   % testing and isolation times assumed simulteneous
tHosp = nan(maxCases, 1);
tDisc = nan(maxCases, 1);
icuFlag = zeros(maxCases, 1);
diedFlag = zeros(maxCases, 1);

% initialise case table
cases = table(caseID, parentID, gen, nOff, ageGroup, Rimult, subclinFlag, vaccDoses, tInfect, tOnset, tQuar, tIsol, tHosp, tDisc, icuFlag, diedFlag);

nSeed0 = poissrnd(par.meanSeedsPerDay * par.seedPeriod);

% Draw total number of border seed cases over the simulated phase
nBorderCases = poissrnd(par.borderSeedsPerDay * par.borderSeedPeriod);

% Total number of seed cases is the existing seed cases in first 7 days + 
% border cases arriving during simulated period
nSeedCases = nSeed0 + nBorderCases;
    
% Set parentID of all seed cases to 0
cases.parentID(1:nSeedCases) = 0;

% Set generation of seed cases to 1 
cases.gen(1:nSeedCases) = 1;

% Set ages of seed cases to be distributed the same as the overall popn
cases.ageGroup(1:nSeedCases) = discRand2(par.popDist, nSeedCases, 1); 

% Reduce infectiousness of seed cases by multiplying Rimult by 
% relInfSeedCases <1, to represent home isolation of seed cases
cases.Rimult(1:nSeed0) = par.relInfSeedCases * cases.Rimult(1:nSeed0);
cases.Rimult(nSeed0+1:nSeedCases) = par.relInfBorderCases * cases.Rimult(nSeed0+1:nSeedCases);

% Set vaccination status of seed cases
cases.vaccDoses(1:nSeedCases) = par.genSeedVaxStatus(cases, nSeedCases, par);

% Set subclinical status for seed cases from probability of being
% subclinical for each age group
cases.subclinFlag(1:nSeedCases) = rand(nSeedCases, 1) < par.VEs(1+cases.vaccDoses(1:nSeedCases)) + (1-par.VEs(1+cases.vaccDoses(1:nSeedCases))).*par.pSub(cases.ageGroup(1:nSeedCases));

% Set time of appearance of each seed case, randomly drawn from the
% simulated time period
cases.tInfect(1:nSeed0) = t(1) + floor(par.seedPeriod*rand(nSeed0, 1));
cases.tInfect(nSeed0+1:nSeedCases) = t(par.borderOpenDate - par.date0) + floor(par.borderSeedPeriod*rand(nBorderCases, 1));

% Set wait till symptom onset period of seed cases from randomly drawn inc periods
cases.tOnset(1:nSeedCases) = cases.tOnset(1:nSeedCases) + cases.tInfect(1:nSeedCases);

% Set probability that seed cases go into isolation from the probabilities 
% of symptomatic and subclinical cases getting tested
pIsol = par.pTestClin*(cases.subclinFlag(1:nSeedCases) == 0 ) + par.pTestSub*(cases.subclinFlag(1:nSeedCases) == 1 );

% Set isolation status of seed cases 
isolFlag = rand(nSeedCases, 1) < pIsol;
nIsol = sum(isolFlag);

% Set time when seed cases went into isolation
tIsol = cases.tOnset(isolFlag) + genIsolDelay(nIsol, 1, par);

% Find which seed cases started isolating before the minimum detection time
ind = tIsol < t(1)+par.minDetectTime;
% Make sure all seed cases start isolating after min detection time
tIsol(ind) = t(1)+par.minDetectTime + rand(sum(ind), 1)*par.followUpTime;
cases.tIsol(isolFlag) = tIsol;

% Fractions of each age group that is in the susceptible class and has had
% 0, 1 or 2 doses of vaccines. These will initially sum to 1 but will
% become depleted as people get infected
susFrac0 = zeros(nSteps, par.nAgeGroups);
susFrac1 = zeros(nSteps, par.nAgeGroups);
susFrac2a = zeros(nSteps, par.nAgeGroups);
susFrac2b = zeros(nSteps, par.nAgeGroups);
susFrac2c = zeros(nSteps, par.nAgeGroups);
susFrac2d = zeros(nSteps, par.nAgeGroups);
susFrac3a = zeros(nSteps, par.nAgeGroups);
susFrac3b = zeros(nSteps, par.nAgeGroups);
susFrac3c = zeros(nSteps, par.nAgeGroups);
susFrac3d = zeros(nSteps, par.nAgeGroups);

cumCov1 = par.cov1+par.cov2a+par.cov2b+par.cov2c+par.cov2d+par.cov3a+par.cov3b+par.cov3c+par.cov3d;
cumCov2a = par.cov2a+par.cov2b+par.cov2c+par.cov2d+par.cov3a+par.cov3b+par.cov3c+par.cov3d;
cumCov2b = par.cov2b+par.cov2c+par.cov2d+par.cov3a+par.cov3b+par.cov3c+par.cov3d;
cumCov2c = par.cov2c+par.cov2d+par.cov3a+par.cov3b+par.cov3c+par.cov3d;
cumCov2d = par.cov2d+par.cov3a+par.cov3b+par.cov3c+par.cov3d;
cumCov3a = par.cov3a+par.cov3b+par.cov3c+par.cov3d;
cumCov3b = par.cov3b+par.cov3c+par.cov3d;
cumCov3c = par.cov3c+par.cov3d;
cumCov3d = par.cov3d;

susFrac0(1, :) = 1-cumCov1(1, :);
susFrac1(1, :) = par.cov1(1, :);
susFrac2a(1, :) = par.cov2a(1, :);
susFrac2b(1, :) = par.cov2b(1, :);
susFrac2c(1, :) = par.cov2c(1, :);
susFrac2d(1, :) = par.cov2d(1, :);
susFrac3a(1, :) = par.cov3a(1, :);
susFrac3b(1, :) = par.cov3b(1, :);
susFrac3c(1, :) = par.cov3c(1, :);
susFrac3d(1, :) = par.cov3d(1, :);

% Fractions moving between dose compartments each time step
q01 = min(1, max(0, (cumCov1(2:end, :)-cumCov1(1:end-1, :))./(1-cumCov1(1:end-1, :))));
q12a = min(1, max(0, (cumCov2a(2:end, :)-cumCov2a(1:end-1, :))./par.cov1(1:end-1, :)));
q2ab = min(1, max(0, (cumCov2b(2:end, :)-cumCov2b(1:end-1, :))./par.cov2a(1:end-1, :)));
q2bc = min(1, max(0, (cumCov2c(2:end, :)-cumCov2c(1:end-1, :))./par.cov2b(1:end-1, :)));
q2cd = min(1, max(0, (cumCov2d(2:end, :)-cumCov2d(1:end-1, :))./par.cov2c(1:end-1, :)));
q2d3a = min(1, max(0, (cumCov3a(2:end, :)-cumCov3a(1:end-1, :))./par.cov2d(1:end-1, :)));
q3ab = min(1, max(0, (cumCov3b(2:end, :)-cumCov3b(1:end-1, :))./par.cov3a(1:end-1, :)));
q3bc = min(1, max(0, (cumCov3c(2:end, :)-cumCov3c(1:end-1, :))./par.cov3b(1:end-1, :)));
q3cd = min(1, max(0, (cumCov3d(2:end, :)-cumCov3d(1:end-1, :))./par.cov3c(1:end-1, :)));


% Relative transmisison rates as a result of vaccination and
% quarantine/isolation
relTransIsol = [1; par.cQuar; par.cIsol];

% Make an array (propList) whose coluumns are the possible combinations of: (1) age group of
% offspring, (2) dose category of offspring
% ageGroupList = repmat((1:par.nAgeGroups)', 10, 1);
% vaxCatList = [zeros(par.nAgeGroups, 1); ones(par.nAgeGroups, 1); 2*ones(par.nAgeGroups, 1); 3*ones(par.nAgeGroups, 1); 4*ones(par.nAgeGroups, 1); 5*ones(par.nAgeGroups, 1); 6*ones(par.nAgeGroups, 1); 7*ones(par.nAgeGroups, 1); 8*ones(par.nAgeGroups, 1); 9*ones(par.nAgeGroups, 1)] ;
% To get ordering right, do combinations in differnet order
ageGroupList = repelem((1:par.nAgeGroups)', 10);
vaxCatList = repmat( (0:9)', par.nAgeGroups, 1);

nCases0 = histcounts(cases.ageGroup(1:nSeedCases), 1:par.nAgeGroups+1);
nCases = nCases0;
ReffEmp = zeros(size(t));
Rvt = zeros(size(t));
dailyAvg = 0;
nActive = 0;
nFuture = nSeedCases;
iStep = 1;
dist = 0;
while iStep < nSteps && sum(nCases) < maxCases && (nActive + nFuture > 0 || relTransCurrentAL(iStep+1) < 1) && dist < earlyReject.threshold
    
    iStep = iStep+1;
    
    % number of seed cases whose infeciotn time is still in the future:
    nFuture = sum( t(iStep) <= cases.tInfect );
    
    % Get IDs of active cases at current time step:
    activeID = cases.caseID(t(iStep) > cases.tInfect & t(iStep) <= cases.tInfect+par.maxInfectTime);
    nActive = length(activeID);
    
    if nActive > 0
        % Simulate new secondary cases on day i from each currently active
        % case:
        
        % Area under the curve of the transmission rate for each active case at current time step
        auci = auc(t(iStep)-cases.tInfect(activeID));
        
        % isolation status of active cases (0 for nothing, 1 for
        % quarantined, 2 for isolated)
        isolStatus = (t(iStep) > cases.tQuar(activeID) & ~(t(iStep) > cases.tIsol(activeID))) + 2*(t(iStep) > cases.tIsol(activeID));
        
        % Calculate the expected number of offspring in each age group during the current
        % time step from each active case, assuming a fully susceptible population.
        % This defines a matrix expOff whose i,j element is the expected number of
        % offspring from parent case i in age group j
        % This is the product of the following factors for each case:
        % - relative transmission due to current alert level
        % - relative transmission due to vaccination status of active case i
        % - relative reproduction number (typically gamma distributed) of case i
        % - relative transmission rate for clinical/subclinical status of case i
        % - relative transmission due to isolation/quarantine status of case i
        % - time-dependent transmission rate for the parent case i, quantified by auci
        % - jth column of the NGM for an unvaccinated clinical individual in the age group of the parent case
        expOff = relTransCurrentAL(iStep) * (1-par.VEt(1+cases.vaccDoses(activeID))) .* ...
            relTransIsol(1+isolStatus) .* cases.Rimult(activeID) .* ...
            (1 - (1-par.cSub)*cases.subclinFlag(activeID)) .* auci .* (par.NGMclin(:, cases.ageGroup(activeID))).';
        
        % split expOff into expected offspring who have had 0, 1 or 2 doses
        % In doing this, some putative infections are prevented by the
        % vaccine (and others by immunity from prior infection is susFrac0,
        % susFrac1 and susFrac2 sum to less than 1).
        expOff0 = expOff .* susFrac0(iStep-1, :);
        expOff1 = (1-par.VEi(2)) * expOff .* susFrac1(iStep-1, :);
        expOff2a = (1-par.VEi(3)) * expOff .* susFrac2a(iStep-1, :);
        expOff2b = (1-par.VEi(4)) * expOff .* susFrac2b(iStep-1, :);
        expOff2c = (1-par.VEi(5)) * expOff .* susFrac2c(iStep-1, :);
        expOff2d = (1-par.VEi(6)) * expOff .* susFrac2d(iStep-1, :);
        expOff3a = (1-par.VEi(7)) * expOff .* susFrac3a(iStep-1, :);
        expOff3b = (1-par.VEi(8)) * expOff .* susFrac3b(iStep-1, :);
        expOff3c = (1-par.VEi(9)) * expOff .* susFrac3c(iStep-1, :);
        expOff3d = (1-par.VEi(10)) * expOff .* susFrac3d(iStep-1, :);
        expOff_withImmunity = expOff0 + expOff1 + expOff2a  + expOff2b + expOff2c + expOff2d + expOff3a + expOff3b + expOff3c + expOff3d;
        ReffEmp(iStep) = sum(sum(expOff0+expOff1+expOff2a+expOff2b+expOff2c+expOff2d+expOff3a+expOff3b+expOff3c+expOff3d))/sum(auci);
        
        % Generate nOff actual number of offspring from parent case i in age group j, thinned due to immunity in each age group as a result of infection prevention in vaccinated individuals and prior infection prevention in vaccinated offspring
        expOff_allParents_byDoseCat = [sum(expOff0, 1); sum(expOff1, 1); sum(expOff2a, 1); sum(expOff2b, 1);  sum(expOff2c, 1); sum(expOff2d, 1); sum(expOff3a, 1); sum(expOff3b, 1); sum(expOff3c, 1); sum(expOff3d, 1)];
        nOff_byDoseCat = poissrnd(expOff_allParents_byDoseCat);
        
        % Total offspring in each age group
        nOffAll = sum(nOff_byDoseCat, 1);
        
        % Total number of offspring summed across all parent cases and all
        % age groups:
        nOffTot = sum(nOffAll);
        
        if nOffTot > 0
            
            secIDs = (sum(nCases)+1:sum(nCases)+nOffTot).';       % IDs for today's newly infected cases
            
            % assign age groups, parent ID, clinical status and vaccination status for new
            % cases based on nOff matrices
            
            % The frequency (number of cases with each age gorup, parent ID
            % and dose combination) of each row of propList is given by
            % the elements in the nOff matrices. Generate a sample X whose
            % columns are the values of these three properties for each
            % secondary case:
            % NOTE: transpose (original version) if propList is [1  0]
            %                                                   [2  0]
            %                                                   [.  .]
            %                                                   [1  1]
            %                                                   ,etc.
            % Don't transpose if it is [1 0]
            %                          [1 1], etc.
            %freqMat = [nOff0; nOff1; nOff2a; nOff2b; nOff2c; nOff2d; nOff3a; nOff3b; nOff3c; nOff3d]';
            %freqMat = [nOff0; nOff1; nOff2a; nOff2b; nOff2c; nOff2d; nOff3a; nOff3b; nOff3c; nOff3d];
            
            X = repelem( [ageGroupList, vaxCatList], nOff_byDoseCat(:) , 1);
            
            % Note cases by secIDs are in age group order - important for
            % efficiently assigning parent cases below
            cases.ageGroup(secIDs) = X(:, 1);
            cases.vaccDoses(secIDs) = X(:, 2);
            
            pParentByOffspringAge = expOff_withImmunity./sum(expOff_withImmunity);
            xm = mnrnd(nOffAll', pParentByOffspringAge')';
            parentList = repmat(1:nActive, par.nAgeGroups, 1)';
            % Parent IDs are also in order of the age group of the
            % offspring to match with the ordering of offspring (see above)
            pid = repelem( parentList(:), xm(:));
            cases.parentID(secIDs) = activeID(pid);
            parentIsolFlag = ~isnan(cases.tIsol(cases.parentID(secIDs)));
            cases.gen(secIDs) = cases.gen(cases.parentID(secIDs))+1;    % Generation of each new cases is generation of parent + 1
            
            pOffSubclin = par.VEs(1+cases.vaccDoses(secIDs)) + (1-par.VEs(1+cases.vaccDoses(secIDs))).*par.pSub(cases.ageGroup(secIDs));
            cases.subclinFlag(secIDs) = rand(nOffTot, 1) < pOffSubclin;
            cases.tInfect(secIDs) = t(iStep);                           % Infection time for each new cases is today
            cases.tOnset(secIDs) = t(iStep) + cases.tOnset(secIDs);     % For efficiency infection to onset delay is pre-stored in cases.tOnset
            cases.nOff(activeID) = cases.nOff(activeID)+sum(nOffAll, 2);
            
            % simulate case testing and isolation effects for new cases
            pIsol = par.pTestClin*(cases.subclinFlag(secIDs) == 0 ) + par.pTestSub*(cases.subclinFlag(secIDs) == 1 ) ;
            isolFlag = rand(nOffTot, 1) < pIsol;
            nIsol = sum(isolFlag);
            tIsol = cases.tOnset(secIDs(isolFlag)) + genIsolDelay(nIsol, 1, par);
            ind = tIsol < t(1)+par.minDetectTime;
            tIsol(ind) = t(1)+par.minDetectTime + rand(sum(ind), 1)*par.followUpTime;
            cases.tIsol(secIDs(isolFlag)) = tIsol;
            
            % simulate contact tracing effects for new cases
            pTrace = par.pTrace * (~(dailyAvg > par.traceCapacity)) * parentIsolFlag;
            traceFlag = rand(nOffTot, 1) < pTrace;
            nTrace = sum(traceFlag);
            cases.tQuar(secIDs(traceFlag)) = cases.tIsol(cases.parentID(secIDs(traceFlag))) + genTraceDelay(nTrace, 1, par);
            % optional: individuals traced prior to onset go into full isolation (as opposed to quarantine) on symptom onset:
            % this is applied to all individuals including subclinical on
            % assumption that asymtoma%tic contacts get tested. This allows
            % offspring of subclinicals to be traced
            cases.tIsol(secIDs(traceFlag)) = min(cases.tIsol(secIDs(traceFlag)), max(cases.tQuar(secIDs(traceFlag)), cases.tOnset(secIDs(traceFlag))));
            
            % simulate clinical outcomes (hsopitalisation, ICU and death) for new cases
            pHospClin = (1-par.VEd(1+cases.vaccDoses(secIDs)))./(1-par.VEs(1+cases.vaccDoses(secIDs))).* par.IHR(cases.ageGroup(secIDs))./par.IDR(cases.ageGroup(secIDs));
            hospFlag = (cases.subclinFlag(secIDs) == 0) & (rand(nOffTot, 1) < pHospClin);
            cases.tHosp(secIDs(hospFlag)) = cases.tOnset(secIDs(hospFlag)) + genHospDelay(sum(hospFlag), 1, par);
            cases.tDisc(secIDs(hospFlag)) = cases.tHosp(secIDs(hospFlag)) + genHospLOS(sum(hospFlag), 1, par);
            % cases are detected and isolated once hospitalised
            cases.tIsol(secIDs(hospFlag)) = min(cases.tIsol(secIDs(hospFlag)), cases.tHosp(secIDs(hospFlag)));
            
            % Update cumulative infections to date (in each age group)
            nCases = nCases+sum(nOffAll, 1);
            %toc
        end
    else
        nOff_byDoseCat = zeros(10, par.nAgeGroups);
    end
    %tic
    % Update susceptible fractions according to vaccinations given
    % and new infections this time step
    susFrac0(iStep, :) =  max(0, (susFrac0(iStep-1, :)   - sum(nOff_byDoseCat(1, :), 1)./(par.popCount')).*(1-q01(iStep-1, :)));
    susFrac1(iStep, :) =  max(0, (susFrac1(iStep-1, :)   - sum(nOff_byDoseCat(2, :), 1)./(par.popCount')).*(1-q12a(iStep-1, :)) + (susFrac0(iStep-1, :) - sum(nOff_byDoseCat(1, :), 1)./(par.popCount')).*q01(iStep-1, :) ) ;
    susFrac2a(iStep, :) = max(0, (susFrac2a(iStep-1, :)  - sum(nOff_byDoseCat(3, :), 1)./(par.popCount')).*(1-q2ab(iStep-1, :)) + (susFrac1(iStep-1, :) - sum(nOff_byDoseCat(2, :), 1)./(par.popCount')).*q12a(iStep-1, :) );
    susFrac2b(iStep, :) = max(0, (susFrac2b(iStep-1, :)  - sum(nOff_byDoseCat(4, :), 1)./(par.popCount')).*(1-q2bc(iStep-1, :)) + (susFrac2a(iStep-1, :)- sum(nOff_byDoseCat(3, :), 1)./(par.popCount')).*q2ab(iStep-1, :) );
    susFrac2c(iStep, :) = max(0, (susFrac2c(iStep-1, :)  - sum(nOff_byDoseCat(5, :), 1)./(par.popCount')).*(1-q2cd(iStep-1, :)) + (susFrac2b(iStep-1, :)- sum(nOff_byDoseCat(4, :), 1)./(par.popCount')).*q2bc(iStep-1, :) );
    susFrac2d(iStep, :) = max(0, (susFrac2d(iStep-1, :)  - sum(nOff_byDoseCat(6, :), 1)./(par.popCount')).*(1-q2d3a(iStep-1, :))+ (susFrac2c(iStep-1, :)- sum(nOff_byDoseCat(5, :), 1)./(par.popCount')).*q2cd(iStep-1, :) );
    susFrac3a(iStep, :) = max(0, (susFrac3a(iStep-1, :)  - sum(nOff_byDoseCat(7, :), 1)./(par.popCount')).*(1-q3ab(iStep-1, :)) + (susFrac2d(iStep-1, :)- sum(nOff_byDoseCat(6, :), 1)./(par.popCount')).*q2d3a(iStep-1, :) );
    susFrac3b(iStep, :) = max(0, (susFrac3b(iStep-1, :)  - sum(nOff_byDoseCat(8, :), 1)./(par.popCount')).*(1-q3bc(iStep-1, :)) + (susFrac3a(iStep-1, :)- sum(nOff_byDoseCat(7, :), 1)./(par.popCount')).*q3ab(iStep-1, :) );
    susFrac3c(iStep, :) = max(0, (susFrac3c(iStep-1, :)  - sum(nOff_byDoseCat(9, :), 1)./(par.popCount')).*(1-q3cd(iStep-1, :)) + (susFrac3b(iStep-1, :)- sum(nOff_byDoseCat(8, :), 1)./(par.popCount')).*q3bc(iStep-1, :) );
    susFrac3d(iStep, :) = max(0, (susFrac3d(iStep-1, :)  - sum(nOff_byDoseCat(10, :), 1)./(par.popCount'))                      + (susFrac3c(iStep-1, :)- sum(nOff_byDoseCat(9, :), 1)./(par.popCount')).*q3cd(iStep-1, :) );
    %     figure(50)
    %     plot(1:16, susFrac0(iStep, :), 1:16, susFrac1(iStep, :),1:16, susFrac2a(iStep, :),1:16, susFrac2b(iStep, :),1:16, susFrac2c(iStep, :),1:16, susFrac2d(iStep, :),1:16, susFrac3a(iStep, :),1:16, susFrac3b(iStep, :),1:16, susFrac3c(iStep, :),1:16, susFrac3d(iStep, :))
    %     hold on
    %     plot(1:16, susFrac0(iStep, :)+susFrac1(iStep, :)+susFrac2a(iStep, :)+susFrac2b(iStep, :)+susFrac2c(iStep, :)+susFrac2d(iStep, :)+susFrac3a(iStep, :)+susFrac3b(iStep, :)+susFrac3c(iStep, :)+susFrac3d(iStep, :), '--')
    %     legend('0', '1', '2a', '2b', '2c', '2d', '3a', '3b', '3c', '3d', 'all')
    %     title(sprintf('%i', iStep))
    %     hold off
    %     drawnow
    
    
    
    
    % Check triggers and update current alert level accordingly
    dailyAvg = sum( t(iStep) > min(cases.tIsol, cases.tQuar) & t(iStep) <= min(cases.tIsol, cases.tQuar)+7 )/7;
    %     nHospBeds = sum( t(iStep) > cases.tHosp & t(iStep) <= cases.tDisc);
    %     if relTransCurrentAL(iStep) > par.relTransAL(par.tierUp) & (dailyAvg > par.maxActiveCases | nHospBeds > par.maxHospBeds)
    %         relTransCurrentAL(iStep+1:end) = par.relTransAL(par.tierUp)
    %     elseif dailyAvg > par.alarmActiveCases | nHospBeds > par.alarmHospBeds   % go up one extra AL if exceed trigger by a factor of 3
    %         relTransCurrentAL(iStep+1:end) = par.relTransAL(4)
    %     elseif ~(dailyAvg > par.minActiveCases) & ~(nHospBeds > par.minHospBeds)
    %        relTransCurrentAL(iStep+1:end) = par.relTransAL(par.tierDown)
    %     end
    
    % calculate distance metric to allow early rejection of simulation
    % during fitting process
    nIsolTemp = histcounts(cases.tIsol, [t(1):t(iStep)+1] );
    dist = calcError(earlyReject.tData, earlyReject.nCasesData, t(1):t(iStep), nIsolTemp);
    %toc
    %fprintf('\n')
end

if iStep < nSteps
    dist = 2*earlyReject.threshold;
end

totCases = sum(nCases);
cases = cases(1:totCases, :);
cases.icuFlag = ~isnan(cases.tHosp) & rand(totCases, 1) < par.pICU(cases.ageGroup)  ;
cases.diedFlag = ~isnan(cases.tHosp) & rand(totCases, 1) <  (1-par.VEm(1+cases.vaccDoses))./(1-par.VEd(1+cases.vaccDoses)) .* par.IFR(cases.ageGroup)./par.IHR(cases.ageGroup);


infBlock = 1-par.cov1(1, :) - ...
    par.cov2a(1, :)-par.cov2b(1, :)-par.cov2c(1, :)-par.cov2d(1, :) - ...
    par.cov3a(1, :)-par.cov3b(1, :)-par.cov3c(1, :)-par.cov3d(1, :) + ...
    (1-par.VEi(2))*par.cov1(1, :) + ...
    (1-par.VEi(3))*par.cov2a(1, :) + (1-par.VEi(4))*par.cov2b(1, :) + (1-par.VEi(5))*par.cov2c(1, :) + (1-par.VEi(5))*par.cov2d(1, :) + ...
    (1-par.VEi(6))*par.cov3a(1, :) + (1-par.VEi(7))*par.cov3b(1, :) + (1-par.VEi(8))*par.cov3c(1, :) + (1-par.VEi(9))*par.cov3d(1, :);    % row
transBlockNum0 = (par.IDR + par.cSub*(1-par.IDR)) .* (1-par.cov1(1, :) - ...
    par.cov2a(1, :)-par.cov2b(1, :)-par.cov2c(1, :)-par.cov2d(1, :) - ...
    par.cov3a(1, :)-par.cov3b(1, :)-par.cov3c(1, :)-par.cov3d(1, :))';
transBlockNum1 = (1-par.VEi(2)) * (1-par.VEt(2)) * ((1-par.VEs(2))*par.IDR + par.cSub*(1-(1-par.VEs(2))*par.IDR)) .* par.cov1(1, :)';
transBlockNum2a = (1-par.VEi(3)) * (1-par.VEt(3)) * ((1-par.VEs(3))*par.IDR + par.cSub*(1-(1-par.VEs(3))*par.IDR)) .* par.cov2a(1, :)';
transBlockNum2b = (1-par.VEi(4)) * (1-par.VEt(4)) * ((1-par.VEs(4))*par.IDR + par.cSub*(1-(1-par.VEs(4))*par.IDR)) .* par.cov2b(1, :)';
transBlockNum2c = (1-par.VEi(5)) * (1-par.VEt(5)) * ((1-par.VEs(5))*par.IDR + par.cSub*(1-(1-par.VEs(5))*par.IDR)) .* par.cov2c(1, :)';
transBlockNum2d = (1-par.VEi(6)) * (1-par.VEt(6)) * ((1-par.VEs(6))*par.IDR + par.cSub*(1-(1-par.VEs(6))*par.IDR)) .* par.cov2d(1, :)';
transBlockNum3a = (1-par.VEi(7)) * (1-par.VEt(7)) * ((1-par.VEs(7))*par.IDR + par.cSub*(1-(1-par.VEs(7))*par.IDR)) .* par.cov3a(1, :)';
transBlockNum3b = (1-par.VEi(8)) * (1-par.VEt(8)) * ((1-par.VEs(8))*par.IDR + par.cSub*(1-(1-par.VEs(8))*par.IDR)) .* par.cov3b(1, :)';
transBlockNum3c = (1-par.VEi(9)) * (1-par.VEt(9)) * ((1-par.VEs(9))*par.IDR + par.cSub*(1-(1-par.VEs(9))*par.IDR)) .* par.cov3c(1, :)';
transBlockNum3d = (1-par.VEi(10)) * (1-par.VEt(10)) * ((1-par.VEs(10))*par.IDR + par.cSub*(1-(1-par.VEs(10))*par.IDR)) .* par.cov3d(1, :)';
transBlock = (transBlockNum0+transBlockNum1+...
    transBlockNum2a+transBlockNum2b+transBlockNum2c+transBlockNum2d+...
    transBlockNum3a+transBlockNum3b+transBlockNum3c+transBlockNum3d)./infBlock';                                     % column
NGMv = par.NGMclin .* infBlock .* transBlock;
l = eigs(NGMv);
Rv = abs(l(1));

l0 = eigs(par.NGM);
Rv0 = abs(l0(1));

vaxEff = (1 - Rv / Rv0);
%     disp(vaxEff)

% calculate TTIQ effect overall and over time
wI = 1-wblcdf(floor(cases.tIsol)-cases.tInfect, par.genA, par.genB); % amount of transmission prevented if isolation was perfect
wQ = 1-wblcdf(floor(min(cases.tQuar, cases.tIsol))-cases.tInfect, par.genA, par.genB); % amount of transmission prevented if quarantine was perfect    % wQ >= wI
wI(isnan(wI)) = 0;
wQ(isnan(wQ)) = 0;
ReffReduc = (1-par.cQuar)*(wQ-wI) + (1-par.cIsol)*(wI);
infAfterDetectFlag = cases.tInfect >= t(1)+par.minDetectTime;
indSub = infAfterDetectFlag & cases.subclinFlag == 0;
indClin = infAfterDetectFlag & cases.subclinFlag == 1;
TTIQeff = 1 - (sum(1-ReffReduc(indClin)) + par.cSub*sum(1-ReffReduc(indSub))) / (sum(indClin) + par.cSub*sum(indSub));
%     disp(TTIQeff)

Reff = par.R0 * par.relTransBaseAL(1) * (1 - TTIQeff) * (1 - vaxEff);
%     disp(Reff)
GR = calcGRfromReff(Reff, par.genA, par.genB);
% disp(log(2)/GR)

end



