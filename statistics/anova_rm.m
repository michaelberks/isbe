function [p, table] = anova_rm(X, displayopt)
%   [p, table] = anova_rm(X, displayopt)
%   Single factor, tepeated measures ANOVA.
%
%   [p, table] = anova_rm(X, displayopt) performs a repeated measures ANOVA
%   for comparing the means of two or more columns (time) in one or more
%   samples(groups). Unbalanced samples (i.e. different number of subjects 
%   per group) is supported though the number of columns (followups)should 
%   be the same. 
%
%   DISPLAYOPT can be 'on' (the default) to display the ANOVA table, or 
%   'off' to skip the display. For a design with only one group and two or 
%   more follow-ups, X should be a matrix with one row for each subject. 
%   In a design with multiple groups, X should be a cell array of matrixes.
% 
%   Example: Gait-Cycle-times of a group of 7 PD patients have been
%   measured 3 times, in one baseline and two follow-ups:
%
%   patients = [
%    1.1015    1.0675    1.1264
%    0.9850    1.0061    1.0230
%    1.2253    1.2021    1.1248
%    1.0231    1.0573    1.0529
%    1.0612    1.0055    1.0600
%    1.0389    1.0219    1.0793
%    1.0869    1.1619    1.0827 ];
%
%   more over, a group of 8 controls has been measured in the same protocol:
%
%   controls = [
%     0.9646    0.9821    0.9709
%     0.9768    0.9735    0.9576
%     1.0140    0.9689    0.9328
%     0.9391    0.9532    0.9237
%     1.0207    1.0306    0.9482
%     0.9684    0.9398    0.9501
%     1.0692    1.0601    1.0766
%     1.0187    1.0534    1.0802 ];
%
%   We are interested to see if the performance of the patients for the
%   followups were the same or not:
%  
%   p = anova_rm(patients);
%
%   By considering the both groups, we can also check to see if the 
%   follow-ups were significantly different and also check two see that the
%   two groups had a different performance:
%
%   p = anova_rm({patients controls});
%
%
%   ref: Statistical Methods for the Analysis of Repeated Measurements, 
%     C. S. Daivs, Springer, 2002
%
%   Copyright 2008, Arash Salarian
%   mailto://arash.salarian@ieee.org
%

if nargin < 2
    displayopt = 'on';
end

if ~iscell(X)
    X = {X};
end

%number of groups
g = size(X,2);  

%subjects per group 
n_g = zeros(g, 1);
for i_g=1:g
    n_g(i_g) = size(X{i_g}, 1);    
end
n = sum(n_g);

%number of follow-ups
t = size(X{1},2);   

% overall mean
y = 0;
for i_g=1:g
    y = y + sum(sum(X{i_g}));
end
y = y / (n * t);

% allocate means
y_g = zeros(g,1); %Means of groups
y_t = zeros(t,1); %Means of time repeats
y_gt = zeros(g,t); %Means of each group (over subjects) at each time
y_gs = cell(g,1); %Means of each subject (over time) in each group
for i_g=1:g
    y_gs{i_g} = zeros(n_g(i_g),1);
end

% Group means
for i_g=1:g
    y_g(i_g) = sum(sum(X{i_g})) / (n_g(i_g) * t);
end

% Means of time repeats
for i_t=1:t
    y_t(i_t) = 0;
    for i_g=1:g
        y_t(i_t) = y_t(i_t) + sum(X{i_g}(:,i_t));
    end
    y_t(i_t) = y_t(i_t) / n;
end

% Means of each group (over subjects) at each time
for i_g=1:g
    for i_t=1:t
        y_gt(i_g,i_t) = sum(X{i_g}(:,i_t) / n_g(i_g));
    end
end

% Means of each subject (over time) in each group
for i_g = 1:g
    for i_s = 1:n_g(i_g)
        y_gs{i_g}(i_s) = sum(X{i_g}(i_s,:)) / t;
    end
end

% calculate the sum of squares
ssG = 0; %SS of groups
ssSG = 0; %SS of each of subjects within each group
ssT = 0; %SS of time repeats
ssGT = 0; %SS of time repeats for each group
ssR = 0; %Residual SS (to compute unexplained error)

for i_g = 1:g %Loop through groups
    for i_s = 1:n_g(i_g) %Loop through subject in each group
        for i_t = 1:t %Loop through each time point
            ssG  = ssG  + (y_g(i_g) - y)^2;
            ssSG = ssSG + (y_gs{i_g}(i_s) - y_g(i_g))^2;
            ssT  = ssT  + (y_t(i_t) - y)^2;
            ssGT = ssGT + (y_gt(i_g,i_t) - y_g(i_g) - y_t(i_t) + y)^2;
            ssR  = ssR  + (X{i_g}(i_s,i_t) - y_gt(i_g,i_t) - y_gs{i_g}(i_s) + y_g(i_g))^2;
        end
    end
end

% calculate mean square errors
if g > 1
    msG  = ssG  / (g-1);
    msGT = ssGT / ((g-1)*(t-1));
end
msSG = ssSG / (n-g);
msT  = ssT  / (t-1);
msR  = ssR  / ((n-g)*(t-1));


% calculate the F-statistics
if g > 1
    FG  = msG  / msSG;
    FGT = msGT / msR;
end
FT  = msT  / msR;
FSG = msSG / msR;


% single or multiple sample designs?
if g > 1
    % case for multiple samples
    pG  = 1 - fcdf(FG, g-1, n-g);
    pT  = 1 - fcdf(FT, t-1, (n-g)*(t-1));
    pGT = 1 - fcdf(FGT, (g-1)*(t-1), (n-g)*(t-1));
    pSG = 1 - fcdf(FSG, n-g, (n-g)*(t-1));

    p = [pT, pG, pSG, pGT];

    table = { 'Source' 'SS' 'df' 'MS' 'F' 'Prob>F'
        'Time'  ssT t-1 msT FT pT
        'Group' ssG g-1 msG FG pG
        'Ineratcion' ssGT (g-1)*(t-1) msGT FGT pGT
        'Subjects (matching)' ssSG n-g msSG FSG pSG
        'Error' ssR (n-g)*(t-1) msR  [] []
        'Total' [] [] [] [] []
        };
    table{end, 2} = sum([table{2:end-1,2}]);
    table{end, 3} = sum([table{2:end-1,3}]);

    if (isequal(displayopt, 'on'))
        digits = [-1 -1 0 -1 2 4];
        statdisptable(table, 'multi-sample repeated measures ANOVA', 'ANOVA Table', '', digits);
    end
else
    % case for only one sample
    pT  = 1 - fcdf(FT, t-1, (n-g)*(t-1));
    pSG = 1 - fcdf(FSG, n-g, (n-g)*(t-1));

    p = [pT, pSG];

    table = { 'Source' 'SS' 'df' 'MS' 'F' 'Prob>F'
        'Time'  ssT t-1 msT FT pT
        'Subjects (matching)' ssSG n-g msSG FSG pSG
        'Error' ssR (n-g)*(t-1) msR  [] []
        'Total' [] [] [] [] []
        };
    table{end, 2} = sum([table{2:end-1,2}]);
    table{end, 3} = sum([table{2:end-1,3}]);

    if (isequal(displayopt, 'on'))
        digits = [-1 -1 0 -1 2 4];
        statdisptable(table, 'repeated measures ANOVA', 'ANOVA Table', '', digits);
    end
end
