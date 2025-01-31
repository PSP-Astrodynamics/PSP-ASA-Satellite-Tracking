function [average] = getAverageF107()

historicalF107DATA = websave('historicalF107DATA.txt','https://www-app3.gfz-potsdam.de/kp_index/Kp_ap_Ap_SN_F107_since_1932.txt');
rawdata = readmatrix('historicalF107Data.txt','CommentStyle','#');
f107data = rawdata((end-81:end), 26);
average = mean(f107data);

end