function [] = spicule_classification_crossfold(fold, n_fold, level, win_size, all_or_max, decomp)
%
%
%
% Inputs:
%      fold- *Insert description of input variable here*
%
%      n_fold- *Insert description of input variable here*
%
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 03-Feb-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester
display(['Spicule crossfold win size = ' num2str(win_size), ' levels = ', num2str(level)', ' subband = ', all_or_max ', fold ', zerostr(fold,2), ' is running']);

%Load in mass and norm data
load([mberksroot, 'classification/data/spicules/exp_samples_' decomp '_' num2str(win_size) '_' all_or_max '_' num2str(level), '.mat']);

%This loads in mass_data, norm_data, norm_pts_cumsum, spic_pts_cumsum

%Get number of mass and normal points
n_mass = size(mass_data,1); %#ok
n_norm = size(norm_data,1); %#ok

%3. Workout start and end indices for fold
if fold == 1
    start_m = 1; %#ok
    start_m_i = 1;
else
    start_m = find((spic_pts_cumsum / n_mass)>= (fold-1)/n_fold, 1) + 1;
    start_m_i = spic_pts_cumsum(start_m-1) + 1;
end
end_m = find((spic_pts_cumsum / n_mass) >= fold/n_fold, 1);
end_m_i = spic_pts_cumsum(end_m);

if fold == 1
    start_n = 1; %#ok
    start_n_i = 1;
else
    start_n = find((norm_pts_cumsum / n_mass)>= (fold-1)/n_fold, 1) + 1;
    start_n_i = norm_pts_cumsum(start_n-1) + 1;
end
end_n = find((norm_pts_cumsum / n_mass) >= fold/n_fold, 1);
end_n_i = norm_pts_cumsum(end_n);

%4. Extract test data for this fold
test_data = [mass_data(start_m_i:end_m_i,:); norm_data(start_n_i:end_n_i,:)]; %#ok
test_labels = [true(end_m_i-start_m_i+1,1); false(end_n_i-start_n_i+1,1)]; %#ok

%5. Extract training data for this fold
forest_args.X = [mass_data([1:start_m_i-1 end_m_i+1:end],:); norm_data([1:start_n_i-1 end_n_i+1:end],:)]; %#ok
forest_args.y = [true(n_mass-end_m_i+start_m_i-1,1); false(n_norm-end_n_i+start_n_i-1,1)];
clear mass_data norm_data

%6. Set remaining forest args
forest_args.n_trees = 200;
forest_args.do_oob = 0;
forest_args.do_oob_trace = 0;

%7. Build forest
random_forest = mb_random_forest_class_train_boot(forest_args);

%9. Save forest and results of test data
forest_name = [mberksroot, 'classification/data/spicules/rf_spic_' decomp '_', num2str(win_size), '_', num2str(level)', '_', all_or_max '_', zerostr(fold,2), '.mat'];
save(forest_name, 'random_forest', 'start_m', 'start_n', 'end_m', 'end_n');

%8. Run forest on test data
[test_predictions test_votes] = mb_random_forest_class_predict_boot(random_forest, test_data); %#ok

results_name = [mberksroot, 'classification/data/spicules/rf_spic_' decomp '_', num2str(win_size), '_', num2str(level)', '_', all_or_max '_', zerostr(fold,2), '_results.mat'];
save(results_name, 'test_predictions', 'test_votes', 'test_labels');

display(results_name);
display('Finished!');

