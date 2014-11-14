for win_size = 3
    for level = 5
        for ii = 2:10
            spicule_classification_crossfold(ii, 10, level, win_size, 'all', 'g2d');
            pack;
        end
    end
end
%%
f = figure; 
hold all; axis equal; axis([0 1 0 1]);
%%
pooled_test_labels = [];
pooled_test_scores = [];
fold_aucs = zeros(10,1);

for ii = 1:10 %loop through different folds of cross-fold evaluation

    %load in data
    load(...
        ['C:\isbe\asymmetry_project\experiments\spicule_classification\' ...
         'rf_spic_g2d_1_5_all_', zerostr(ii,2) '_results.mat']);
    fold_accuracy = mean(test_labels == str2double(test_predictions));
    display(['accuracy of fold ', num2str(ii), ' = ', num2str(fold_accuracy)]);

    test_scores = test_votes(:,2) / sum(test_votes(1,:));

    %Pool labels and scores for this fold
    pooled_test_labels = [pooled_test_labels; test_labels]; %#ok
    pooled_test_scores = [pooled_test_scores; test_scores]; %#ok

end 

% Comptue ROC statistics for pooled data
[pooled_roc_pts pooled_auc dum dum auc_se] = calculate_roc_curve(pooled_test_scores,pooled_test_labels,linspace(-0.0001,1.0001,100));
ci_95 = 1.96*auc_se;

%Plot curve
plot(pooled_roc_pts(:,1), pooled_roc_pts(:,2), '-', 'linewidth', 2);
%%
test_im = u_load('D:\isbe\dev\image_data\masses512x512\mass028.mat');
[line_strength, orientation_map, scale_map] = gaussian_2nd_derivative_line(test_im, [1 2 4 8]);
[grad_strength, grad_orientation] = gaussian_1st_derivative_gradient(test_im, 10);

discard_map = ...
    ((abs(mb_mod(orientation_map - grad_orientation + pi/2, pi)) < pi/6) &...
    (grad_strength > 25)) | ...
    (line_strength > 0);

line_map = abs(line_strength);
line_map(discard_map) = 0;
line_map = line_map ./ sum(line_map(:));
%%
circle_area = pi*(20*50/9)^2;

rf_f1_ab = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\rfs_froc_abnormals_f1.mat');
rf_f1_norm = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\rfs_froc_normals_f1.mat');

wrf_f1_ab = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\wrfs_froc_abnormals_f1.mat');
wrf_f1_norm = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\wrfs_froc_normals_f1.mat');

g2_f1_ab = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\g2s_froc_abnormals_f1.mat');
g2_f1_norm = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\g2s_froc_normals_f1.mat');

figure; hold all;
plot(mean(rf_f1_norm.fp_pixels) / circle_area, mean(rf_f1_ab.tp > 0), 'x');
plot(mean(wrf_f1_norm.fp_pixels) / circle_area, mean(wrf_f1_ab.tp > 0), 'x');
plot(mean(g2_f1_norm.fp_pixels) / circle_area, mean(g2_f1_ab.tp > 0), 'x');
legend({'RF maps', 'weighted RF maps', 'gaussian maps'});
%
rf_f2_ab = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\rfs_froc_abnormals_f2.mat');
rf_f2_norm = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\rfs_froc_normals_f2.mat');

wrf_f2_ab = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\wrfs_froc_abnormals_f2.mat');
wrf_f2_norm = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\wrfs_froc_normals_f2.mat');

g2_f2_ab = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\g2s_froc_abnormals_f2.mat');
g2_f2_norm = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\g2s_froc_normals_f2.mat');

figure; hold all;
plot(mean(rf_f2_norm.fp_pixels) / circle_area, mean(rf_f2_ab.tp > 0), 'x');
plot(mean(wrf_f2_norm.fp_pixels) / circle_area, mean(wrf_f2_ab.tp > 0), 'x');
plot(mean(g2_f2_norm.fp_pixels) / circle_area, mean(g2_f2_ab.tp > 0), 'x');
legend({'RF maps', 'weighted RF maps', 'gaussian maps'});
%%
circle_area = pi*(20*50/9)^2;

rf_f1_ab = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\rfs_froc_abnormals_f1.mat');
rf_f1_norm = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\rfs_froc_normals_f1.mat');

wrf_f1_ab = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\wrfs_froc_abnormals_f1.mat');
wrf_f1_norm = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\wrfs_froc_normals_f1.mat');

g2_f1_ab = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\g2s_froc_abnormals_f1.mat');
g2_f1_norm = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\g2s_froc_normals_f1.mat');

figure; hold all;
plot(mean(rf_f1_norm.fp), mean(rf_f1_ab.tp > 0), 'x');
plot(mean(wrf_f1_norm.fp), mean(wrf_f1_ab.tp > 0), 'x');
plot(mean(g2_f1_norm.fp), mean(g2_f1_ab.tp > 0), 'x');
legend({'RF maps', 'weighted RF maps', 'gaussian maps'});
%
rf_f2_ab = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\rfs_froc_abnormals_f2.mat');
rf_f2_norm = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\rfs_froc_normals_f2.mat');

wrf_f2_ab = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\wrfs_froc_abnormals_f2.mat');
wrf_f2_norm = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\wrfs_froc_normals_f2.mat');

g2_f2_ab = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\g2s_froc_abnormals_f2.mat');
g2_f2_norm = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\g2s_froc_normals_f2.mat');

figure; hold all;
plot(mean(rf_f2_norm.fp), mean(rf_f2_ab.tp > 0), 'x');
plot(mean(wrf_f2_norm.fp), mean(wrf_f2_ab.tp > 0), 'x');
plot(mean(g2_f2_norm.fp), mean(g2_f2_ab.tp > 0), 'x');
legend({'RF maps', 'weighted RF maps', 'gaussian maps'});