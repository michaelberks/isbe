%--------------------------------------------------------------------------
%---- Script to create training and test sets of bars ---------------------
%--------------------------------------------------------------------------

%%
% Create a single bar using randomly selected parameters
%--------------------------------------------------------------------------

%get directory listing of backgrounds
bg_list = dir('*mam_dir*\chen\data\normal_smooth128\*.mat');

%randomly select background - will load a matrix 'bg'
bg_idx = ceil(length(bg_list)*rand);
load(['*mam_dir*\chen\data\normal_smooth128\', bg_list(bg_idx).name]);

%randomly select width between 4 and 32
width = 4 + 28*rand;

%randomly select contrast between 8 and 16
contrast = 8 + 8*rand;

%randomly select orientation between 0 and 180 degrees
orientation = 180*rand;

%randomly select bar type
g_bar = rand > .5;

% make bar in 128x128 image

if g_bar
    %make gaussian bar in 128x128 image
    bar = []; %use your method here
    truth_file = []; %use your method here
else
    %make rectangular bar
    bar = []; %use your method here
    truth_file = []; %use your method here
end

bar = bar + bg;
%--------------------------------------------------------------------------
%%