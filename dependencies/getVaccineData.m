function [dates, V1, V2, V3] = getVaccineData(fNamePart)

fNameDates = "data/dates.csv";
fName1 = "data/firstDoseProp" + fNamePart + ".csv";
fName2 = "data/secondDoseProp" + fNamePart + ".csv";
fName3 = "data/thirdDoseProp" + fNamePart + ".csv";

% fprintf('   Loading vaccination data:    %s\n                                %s\n                                %s\n', fNameDates, fName1, fName2)
dates = readmatrix(fNameDates, 'OutputType', 'datetime');
V1 = readmatrix(fName1).';
V2 = readmatrix(fName2).';
V3 = readmatrix(fName3).';

end


