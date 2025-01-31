function [average] = getAverageF107()

rawdata = readmatrix('historicalF107DATA.txt');
f107data = rawdata((end-81:end), 26);
average = mean(f107data);
end