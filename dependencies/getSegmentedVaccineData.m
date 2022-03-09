function [dates, V1, V2a, V2b, V2c, V2d, V3a, V3b, V3c, V3d] = getSegmentedVaccineData(fNamePart, par)

% combine data on vaccine doses with simple model for booster and 5-11 rollout

% booster rollout model
boostGap = 134;     % 120 days to become eligible plus 7 days (on average)
boostProp = par.boostProp;
boostElig = [0 0 0 2/5 ones(1, 12)];

%5-11 vax model
startDate = datenum('17JAN2022') + par.vaccImmDelay;
childGap = 8*7;
maxChildCov = 0.75;


[dates, V1, V2, V3] = getVaccineData(fNamePart);    % read time-dependent vaccination data
dates = dates+par.vaccImmDelay;     % shift vaccination dates to allow for delay

iDate = find(datenum(dates) == par.date0);

nPad = par.tEnd - (length(dates) - iDate);
if nPad > 0
    dates = [dates, dates(end)+1:dates(end)+nPad];
    V1 = [V1(1:end, :); repmat(V1(end, :), nPad, 1)];
    V2 = [V2(1:end, :); repmat(V2(end, :), nPad, 1)];
    V3 = [V3(1:end, :); repmat(V3(end, :), nPad, 1)];
end

V123 = V1+V2+V3;
V23 = V2+V3;

% Assume everyone single dosed gets their second dose 35 days later if they
% haven't already had it
V23(1+35:end, :) = max(V23(1+35:end, :), V123(1:end-35, :));

pChild1 = maxChildCov * ( (datenum(dates') >= startDate ) .*          min(1, (datenum(dates')+1-startDate)/childGap ) );
pChild2 = maxChildCov * ( (datenum(dates') >= startDate+childGap ) .* min(1, (datenum(dates')+1-startDate-childGap)/childGap) );
V123(:, 2) = V123(:, 2) + pChild1;
V123(:, 3) = V123(:, 3) + 2/5*pChild1;
V23(:, 2) = V23(:, 2) + pChild2;
V23(:, 3) = V23(:, 3) + 2/5*pChild2;

V1 = max(0, V123-V23);


% optional assume x% of adults who are double-dosed by time t are boosted by time
% (t + boostGap)
V3(1+boostGap:end, :) = max(V3(1+boostGap:end, :), boostProp*V23(1:end-boostGap, :).*boostElig);


% Segment second and third dose compartments by time since last dose:

% Because data is already time shifted by 2 weeks (to allow for time for
% immune response), this segments into 2-5 weeks, 5-10 weeks, 10-15 weeks
% and 15+ weeks after dose
cutPoints = [21, 56, 91];

V2a = [ nan(cutPoints(1), par.nAgeGroups); V23(1+cutPoints(1):end, :) - V23(1:end-cutPoints(1), :) ];
V2b = [ nan(cutPoints(2), par.nAgeGroups); V23(1+cutPoints(2)-cutPoints(1):end-cutPoints(1), :) - V23(1:end-cutPoints(2), :) ];
V2c = [ nan(cutPoints(3), par.nAgeGroups); V23(1+cutPoints(3)-cutPoints(2):end-cutPoints(2), :) - V23(1:end-cutPoints(3), :) ];
V2d = [ nan(cutPoints(3), par.nAgeGroups); V23(1:end-cutPoints(3), :) - V3(1+cutPoints(3):end, :) ];
V3a = [ nan(cutPoints(1), par.nAgeGroups); V3(1+cutPoints(1):end, :) - V3(1:end-cutPoints(1), :) ];
V3b = [ nan(cutPoints(2), par.nAgeGroups); V3(1+cutPoints(2)-cutPoints(1):end-cutPoints(1), :) - V3(1:end-cutPoints(2), :) ];
V3c = [ nan(cutPoints(3), par.nAgeGroups); V3(1+cutPoints(3)-cutPoints(2):end-cutPoints(2), :) - V3(1:end-cutPoints(3), :) ];
V3d = [ nan(cutPoints(3), par.nAgeGroups); V3(1:end-cutPoints(3), :) ];

dates = dates(iDate:end);
V1 = V1(iDate:end, :);
V23 = V23(iDate:end, :);
V2a = V2a(iDate:end, :);
V2b = V2b(iDate:end, :);
V2c = V2c(iDate:end, :);
V2d = V2d(iDate:end, :);
V3 = V3(iDate:end, :);
V3a = V3a(iDate:end, :);
V3b = V3b(iDate:end, :);
V3c = V3c(iDate:end, :);
V3d = V3d(iDate:end, :);



V1_all = sum( V1.*par.popDist', 2); 
V23_all = sum( V23.*par.popDist', 2);
V2a_all = sum( V2a.*par.popDist', 2); 
V2b_all = sum( V2b.*par.popDist', 2);
V2c_all = sum( V2c.*par.popDist', 2);
V2d_all = sum( V2d.*par.popDist', 2);
V3a_all = sum( V3a.*par.popDist', 2); 
V3b_all = sum( V3b.*par.popDist', 2);
V3c_all = sum( V3c.*par.popDist', 2);
V3d_all = sum( V3d.*par.popDist', 2);
V3_all = sum( V3.*par.popDist', 2);

% figure;
% plot(dates, V1_all, dates, V2a_all, dates, V2b_all, dates, V2c_all, dates, V2d_all, dates, V3a_all, dates, V3b_all, dates, V3c_all, dates, V3d_all, dates, V1_all+V2a_all+V2b_all+V2c_all+V2d_all+V3a_all+V3b_all+V3c_all+V3d_all, '--', dates, V2a_all+V2b_all+V2c_all+V2d_all+V3a_all+V3b_all+V3c_all+V3d_all, '--', dates, V3a_all+V3b_all+V3c_all+V3d_all, '--')
% legend("1", "2a", "2b", "2c", "2d", "3a", "3b", "3c", "3d", "1+", "2+", "3")
% drawnow
% 
% fprintf("%.2f %.2f %.2f %.2f %.2f\n", V3_all(1), V3_all(1+31), V3_all(1+31+28), V3_all(1+31+28+31), V3_all(1+31+28+31+30))
% figure;
% V3_all = V3a_all+V3b_all+V3c_all+V3d_all;
% plot(dates, V1_all+V2a_all+V2b_all+V2c_all+V2d_all+V3a_all+V3b_all+V3c_all+V3d_all, dates, V2a_all+V2b_all+V2c_all+V2d_all+V3a_all+V3b_all+V3c_all+V3d_all, dates, V3_all)
% ylabel('proportion of total population')
% legend("1+ doses", "2+ doses", "3 doses")
% drawnow
% ylim([0 1])
% pause
end
