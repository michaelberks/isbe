%
% Script for experimenting with parameters for subtracting background
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Michael Berks
% Experiment date 11/09/2006

if 0
    
%%%%%% keep n1 = 0 n2 = 20 spacing = 5, vary sigma %%%%%%%%%%%%%%%%%
[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 0, 20, 5, 3, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_0_20_5_3 errors;
er_0_20_5_var(:,1) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 0, 20, 5, 5, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_0_20_5_5 errors;
er_0_20_5_var(:,2) = errors;
clear d erros;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 0, 20, 5, 7, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_0_20_5_7 errors;
er_0_20_5_var(:,3) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 0, 20, 5, 9, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_0_20_5_9 errors;
er_0_20_5_var(:,4) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 0, 20, 5, 11, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_0_20_5_11 errors;
er_0_20_5_var(:,5) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 0, 20, 5, 13, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_0_20_5_13 errors;
er_0_20_5_var(:,6) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 0, 20, 5, 15, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_0_20_5_15 errors;
er_0_20_5_var(:,7) = errors;
save C:\isbe\mammograms\test\sub\tps\er_0_20_5_var er_0_20_5_var; 
clear er_0_20_5_var d errors;
pack;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% keep n1 = 10 n2 = 20 spacing = 5, vary sigma %%%%%%%%%%%%%%%%%
[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 10, 20, 5, 3, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_10_20_5_3 errors;
er_10_20_5_var(:,1) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 10, 20, 5, 5, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_10_20_5_5 errors;
er_10_20_5_var(:,2) = errors;
clear d erros;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 10, 20, 5, 7, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_10_20_5_7 errors;
er_10_20_5_var(:,3) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 10, 20, 5, 9, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_10_20_5_9 errors;
er_10_20_5_var(:,4) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 10, 20, 5, 11, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_10_20_5_11 errors;
er_10_20_5_var(:,5) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 10, 20, 5, 13, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_10_20_5_13 errors;
er_10_20_5_var(:,6) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 10, 20, 5, 15, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_10_20_5_15 errors;
er_10_20_5_var(:,7) = errors;
save C:\isbe\mammograms\test\sub\tps\er_10_20_5_var er_10_20_5_var; 
clear er_10_20_5_var d errors;
pack;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% keep n1 = 30 n2 = 20 spacing = 5, vary sigma %%%%%%%%%%%%%%%%%
[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 30, 20, 5, 3, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_30_20_5_3 errors;
er_30_20_5_var(:,1) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 30, 20, 5, 5, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_30_20_5_5 errors;
er_30_20_5_var(:,2) = errors;
clear d erros;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 30, 20, 5, 7, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_30_20_5_7 errors;
er_30_20_5_var(:,3) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 30, 20, 5, 9, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_30_20_5_9 errors;
er_30_20_5_var(:,4) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 30, 20, 5, 11, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_30_20_5_11 errors;
er_30_20_5_var(:,5) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 30, 20, 5, 13, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_30_20_5_13 errors;
er_30_20_5_var(:,6) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 30, 20, 5, 15, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_30_20_5_15 errors;
er_30_20_5_var(:,7) = errors;
save C:\isbe\mammograms\test\sub\tps\er_30_20_5_var er_30_20_5_var; 
clear er_30_20_5_var d errors;
pack;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% keep n1 = 40 n2 = 20 spacing = 5, vary sigma %%%%%%%%%%%%%%%%%
[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 40, 20, 5, 3, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_40_20_5_3 errors;
er_40_20_5_var(:,1) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 40, 20, 5, 5, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_40_20_5_5 errors;
er_40_20_5_var(:,2) = errors;
clear d erros;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 40, 20, 5, 7, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_40_20_5_7 errors;
er_40_20_5_var(:,3) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 40, 20, 5, 9, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_40_20_5_9 errors;
er_40_20_5_var(:,4) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 40, 20, 5, 11, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_40_20_5_11 errors;
er_40_20_5_var(:,5) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 40, 20, 5, 13, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_40_20_5_13 errors;
er_40_20_5_var(:,6) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 40, 20, 5, 15, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_40_20_5_15 errors;
er_40_20_5_var(:,7) = errors;
save C:\isbe\mammograms\test\sub\tps\er_40_20_5_var er_40_20_5_var; 
clear er_40_20_5_var d errors;
pack;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% keep n1 = 50 n2 = 20 spacing = 5, vary sigma %%%%%%%%%%%%%%%%%
[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 50, 20, 5, 3, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_50_20_5_3 errors;
er_50_20_5_var(:,1) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 50, 20, 5, 5, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_50_20_5_5 errors;
er_50_20_5_var(:,2) = errors;
clear d erros;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 50, 20, 5, 7, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_50_20_5_7 errors;
er_50_20_5_var(:,3) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 50, 20, 5, 9, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_50_20_5_9 errors;
er_50_20_5_var(:,4) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 50, 20, 5, 11, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_50_20_5_11 errors;
er_50_20_5_var(:,5) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 50, 20, 5, 13, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_50_20_5_13 errors;
er_50_20_5_var(:,6) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 50, 20, 5, 15, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\e_50_20_5_15 errors;
er_50_20_5_var(:,7) = errors;
save C:\isbe\mammograms\test\sub\tps\er_50_20_5_var er_50_20_5_var; 
clear er_50_20_5_var d errors;
pack;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%% keep n1 = 20 spacing = 5 sigma = 5, vary n2 %%%%%%%%%%%%%%%%%
[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 20, 10, 5, 5, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\er_20_10_5_5 errors;
er_20_var_5_5(:,1) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 20, 20, 5, 5, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\er_20_20_5_5 errors;
er_20_var_5_5(:,2) = errors;
clear d erros;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 20, 30, 5, 5, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\er_20_30_5_5 errors;
er_20_var_5_5(:,3) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 20, 40, 5, 5, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\er_20_40_5_5 errors;
er_20_var_5_5(:,4) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 20, 50, 5, 5, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\er_20_50_5_5 errors;
er_20_var_5_5(:,5) = errors;
save C:\isbe\mammograms\test\sub\tps\er_20_var_5_5 er_20_var_5_5;
clear d errors er_20_var_5_5;
pack;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% keep n1 = 20 spacing = 5 sigma = 5, vary n2 %%%%%%%%%%%%%%%%%
[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 20, 20, 3, 5, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\er_20_20_3_5 errors;
er_20_20_var_5(:,1) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 20, 20, 5, 5, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\er_20_20_5_5 errors;
er_20_20_var_5(:,2) = errors;
clear d erros;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 20, 20, 7, 5, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\er_20_20_7_5 errors;
er_20_20_var_5(:,3) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 20, 20, 9, 5, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\er_20_20_9_5 errors;
er_20_20_var_5(:,4) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 20, 20, 11, 5, 0, 1, 'biharmTPS');
save C:\isbe\mammograms\test\sub\tps\er_20_20_11_5 errors;
er_20_20_var_5(:,5) = errors;
save C:\isbe\mammograms\test\sub\tps\er_20_20_var_5 er_20_20_var_5;
clear d errors er_20_20_var_5;
pack;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USE CLAMPED-PLATE SPLINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% keep n1 = 0 n2 = 20 spacing = 5, vary sigma %%%%%%%%%%%%%%%%%
[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 0, 20, 5, 3, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_0_20_5_3 errors;
er_0_20_5_var(:,1) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 0, 20, 5, 5, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_0_20_5_5 errors;
er_0_20_5_var(:,2) = errors;
clear d erros;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 0, 20, 5, 7, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_0_20_5_7 errors;
er_0_20_5_var(:,3) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 0, 20, 5, 9, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_0_20_5_9 errors;
er_0_20_5_var(:,4) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 0, 20, 5, 11, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_0_20_5_11 errors;
er_0_20_5_var(:,5) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 0, 20, 5, 13, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_0_20_5_13 errors;
er_0_20_5_var(:,6) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 0, 20, 5, 15, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_0_20_5_15 errors;
er_0_20_5_var(:,7) = errors;
save C:\isbe\mammograms\test\sub\cps\er_0_20_5_var er_0_20_5_var; 
clear er_0_20_5_var d errors;
pack;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% keep n1 = 10 n2 = 20 spacing = 5, vary sigma %%%%%%%%%%%%%%%%%
[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 10, 20, 5, 3, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_10_20_5_3 errors;
er_10_20_5_var(:,1) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 10, 20, 5, 5, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_10_20_5_5 errors;
er_10_20_5_var(:,2) = errors;
clear d erros;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 10, 20, 5, 7, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_10_20_5_7 errors;
er_10_20_5_var(:,3) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 10, 20, 5, 9, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_10_20_5_9 errors;
er_10_20_5_var(:,4) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 10, 20, 5, 11, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_10_20_5_11 errors;
er_10_20_5_var(:,5) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 10, 20, 5, 13, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_10_20_5_13 errors;
er_10_20_5_var(:,6) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 10, 20, 5, 15, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_10_20_5_15 errors;
er_10_20_5_var(:,7) = errors;
save C:\isbe\mammograms\test\sub\cps\er_10_20_5_var er_10_20_5_var; 
clear er_10_20_5_var d errors;
pack;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% keep n1 = 20 n2 = 20 spacing = 5, vary sigma %%%%%%%%%%%%%%%%%
[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 20, 20, 5, 3, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_20_20_5_3 errors;
er_20_20_5_var(:,1) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 20, 20, 5, 5, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_20_20_5_5 errors;
er_20_20_5_var(:,2) = errors;
clear d erros;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 20, 20, 5, 7, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_20_20_5_7 errors;
er_20_20_5_var(:,3) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 20, 20, 5, 9, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_20_20_5_9 errors;
er_20_20_5_var(:,4) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 20, 20, 5, 11, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_20_20_5_11 errors;
er_20_20_5_var(:,5) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 20, 20, 5, 13, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_20_20_5_13 errors;
er_20_20_5_var(:,6) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 20, 20, 5, 15, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_20_20_5_15 errors;
er_20_20_5_var(:,7) = errors;
save C:\isbe\mammograms\test\sub\cps\er_20_20_5_var er_20_20_5_var; 
clear er_20_20_5_var d errors;
pack;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% keep n1 = 30 n2 = 20 spacing = 5, vary sigma %%%%%%%%%%%%%%%%%
[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 30, 20, 5, 3, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_30_20_5_3 errors;
er_30_20_5_var(:,1) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 30, 20, 5, 5, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_30_20_5_5 errors;
er_30_20_5_var(:,2) = errors;
clear d erros;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 30, 20, 5, 7, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_30_20_5_7 errors;
er_30_20_5_var(:,3) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 30, 20, 5, 9, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_30_20_5_9 errors;
er_30_20_5_var(:,4) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 30, 20, 5, 11, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_30_20_5_11 errors;
er_30_20_5_var(:,5) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 30, 20, 5, 13, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_30_20_5_13 errors;
er_30_20_5_var(:,6) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 30, 20, 5, 15, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_30_20_5_15 errors;
er_30_20_5_var(:,7) = errors;
save C:\isbe\mammograms\test\sub\cps\er_30_20_5_var er_30_20_5_var; 
clear er_30_20_5_var d errors;
pack;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% keep n1 = 40 n2 = 20 spacing = 5, vary sigma %%%%%%%%%%%%%%%%%
[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 40, 20, 5, 3, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_40_20_5_3 errors;
er_40_20_5_var(:,1) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 40, 20, 5, 5, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_40_20_5_5 errors;
er_40_20_5_var(:,2) = errors;
clear d erros;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 40, 20, 5, 7, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_40_20_5_7 errors;
er_40_20_5_var(:,3) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 40, 20, 5, 9, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_40_20_5_9 errors;
er_40_20_5_var(:,4) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 40, 20, 5, 11, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_40_20_5_11 errors;
er_40_20_5_var(:,5) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 40, 20, 5, 13, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_40_20_5_13 errors;
er_40_20_5_var(:,6) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 40, 20, 5, 15, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_40_20_5_15 errors;
er_40_20_5_var(:,7) = errors;
save C:\isbe\mammograms\test\sub\cps\er_40_20_5_var er_40_20_5_var; 
clear er_40_20_5_var d errors;
pack;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% keep n1 = 50 n2 = 20 spacing = 5, vary sigma %%%%%%%%%%%%%%%%%
[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 50, 20, 5, 3, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_50_20_5_3 errors;
er_50_20_5_var(:,1) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 50, 20, 5, 5, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_50_20_5_5 errors;
er_50_20_5_var(:,2) = errors;
clear d erros;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 50, 20, 5, 7, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_50_20_5_7 errors;
er_50_20_5_var(:,3) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 50, 20, 5, 9, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_50_20_5_9 errors;
er_50_20_5_var(:,4) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 50, 20, 5, 11, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_50_20_5_11 errors;
er_50_20_5_var(:,5) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 50, 20, 5, 13, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_50_20_5_13 errors;
er_50_20_5_var(:,6) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 50, 20, 5, 15, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\e_50_20_5_15 errors;
er_50_20_5_var(:,7) = errors;
save C:\isbe\mammograms\test\sub\cps\er_50_20_5_var er_50_20_5_var; 
clear er_50_20_5_var d errors;
pack;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% keep n1 = 20 spacing = 5 sigma = 5, vary n2 %%%%%%%%%%%%%%%%%
[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 20, 10, 5, 5, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\er_20_10_5_5 errors;
er_20_var_5_5(:,1) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 20, 20, 5, 5, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\er_20_20_5_5 errors;
er_20_var_5_5(:,2) = errors;
clear d erros;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 20, 30, 5, 5, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\er_20_30_5_5 errors;
er_20_var_5_5(:,3) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 20, 40, 5, 5, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\er_20_40_5_5 errors;
er_20_var_5_5(:,4) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 20, 50, 5, 5, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\er_20_50_5_5 errors;
er_20_var_5_5(:,5) = errors;
save C:\isbe\mammograms\test\sub\cps\er_20_var_5_5 er_20_var_5_5;
clear d errors er_20_var_5_5;
pack;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% keep n1 = 20 n2 = 2 sigma = 5, vary spacing %%%%%%%%%%%%%%%%%
[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 20, 20, 3, 5, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\er_20_20_3_5 errors;
er_20_20_var_5(:,1) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 20, 20, 5, 5, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\er_20_20_5_5 errors;
er_20_20_var_5(:,2) = errors;
clear d erros;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 20, 20, 7, 5, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\er_20_20_7_5 errors;
er_20_20_var_5(:,3) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 20, 20, 9, 5, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\er_20_20_9_5 errors;
er_20_20_var_5(:,4) = errors;
clear d errors;
pack;

[d errors] = subtract_mass_bg(test_data, 'C:\isbe\mammograms\test\temp_in', 'C:\isbe\mammograms\test\temp_out', 20, 20, 11, 5, 0, 1, 'biharmCPS');
save C:\isbe\mammograms\test\sub\cps\er_20_20_11_5 errors;
er_20_20_var_5(:,5) = errors;
save C:\isbe\mammograms\test\sub\cps\er_20_20_var_5 er_20_20_var_5;
clear d errors er_20_20_var_5;
pack;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%