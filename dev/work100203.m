generate_bar_training_data(...
    'num_images', 100,...
    'image_dir', 'C:\isbe\dev\classification\data\bg+bar_128_test\',...
    'num_levels', 4,...
    'compute_dt', 1);
generate_bar_training_data(...
    'num_images', 160,...
    'image_dir', 'C:\isbe\dev\classification\data\bg+bar_128\',...
    'num_levels', 4,...
    'compute_dt', 1);
%%
sampling_method_args.num_samples = 1e5;
sampling_method_args.image_dir = [mberksroot, 'classification\data\bg+bar_128\'];
sampling_method_args.total_samples = 128*128*160;

forest_args.sampling_method = 'sample_image_training_data';
forest_args.sampling_method_args = sampling_method_args;
forest_args.d = 10;
forest_args.n_trees = 2;

for ii = 1:3

    %forest_args.tree_dir = [mberksroot, 'classification\rf\rf_', num2str(ii), '_trees\'];
    %forest_args.save_path = [mberksroot, 'classification\rf\rf_', num2str(ii), '.mat'];
    forest_args.tree_dir = ['\\isbe-san1\mberks\dev\classification\rf\rf_', num2str(ii), '_trees\'];
    forest_args.save_path = ['\\isbe-san1\mberks\dev\classification\rf\rf_', num2str(ii), '.mat'];
    mb_random_forest_class_train(forest_args);
end
%
random_forest = mb_combine_rfs('rf_dir', '\\isbe-san1\mberks\dev\classification\rf\',...
    'tree_dir', 'C:\isbe\dev\classification\rf\rf_all\', 'delete_trees', 0,...
    'save_path', 'C:\isbe\dev\classification\rf\rf_all.mat');
%%
for ii = 1:100
    load(['C:\isbe\dev\classification\data\bg+bar_128_test\bar', zerostr(ii,3), '.mat']);
    [probability_image] = classify_image(image_out, random_forest);
    save(['C:\isbe\dev\classification\data\bg+bar_128_test_results_isbe\bar', zerostr(ii,3), '_results.mat'], 'probability_image');
%     figure; 
%     subplot(1,2,1); imagesc(image_out); axis image;
%     subplot(1,2,2); imagesc(probability_image); axis image;
end
%%
random_forest = mb_combine_rfs('rf_dir', '\\isbe-san1\mberks\dev\classification\rf\',...
    'tree_dir', 'C:\isbe\dev\classification\rf\rf_bar_128\', 'delete_trees', 0,...
    'save_path', 'C:\isbe\dev\classification\rf\rf_bar_128.mat');
%%
for ii = 1:20
    pb_brei = u_load([mberksroot, 'classification/data/bg+bar_128_test_results_brei/bar', zerostr(ii,3), '_results.mat']);
    pb_isbe = u_load([mberksroot, 'classification/data/bg+bar_128_test_results_isbe/bar', zerostr(ii,3), '_results.mat']);
    load(['C:\isbe\dev\classification\data\bg+bar_128_test\bar', zerostr(ii,3), '.mat']);
    figure;
    subplot(1,3,1); imagesc(image_out); axis image;
    subplot(1,3,2); imagesc(pb_brei); axis image;
    subplot(1,3,3); imagesc(pb_isbe); axis image;
end
%%
for ii = 1:10
    load(['\\isbe-san1\mberks\dev\classification\rf\rf_', num2str(ii)]);
    random_forest.tree_dir = ['classification/rf/rf_', num2str(ii), '_trees/'];
    random_forest.tree_root = '\\isbe-san1\mberks\dev\';
    save(['\\isbe-san1\mberks\dev\classification\rf\rf_', num2str(ii)], 'random_forest');
end
%%
tp_counts_brei = 0;
fp_counts_brei = 0;
tp_counts_isbe = 0;
fp_counts_isbe = 0;
t_counts = 0;
f_counts = 0;

for ii = 1:100
    pb_brei = u_load([mberksroot, 'classification/data/bg+bar_128_test_results_brei/bar', zerostr(ii,3), '_results.mat']);
    pb_isbe = u_load([mberksroot, 'classification/data/bg+bar_128_test_results_isbe/bar', zerostr(ii,3), '_results.mat']);
    load(['C:\isbe\dev\classification\data\bg+bar_128_test\bar', zerostr(ii,3), '.mat']);
    
    [d d tp_count fp_count] = calculate_roc_curve(1-pb_brei(:),label(:),(-1:100)/100);
    tp_counts_brei = tp_counts_brei + tp_count;
    fp_counts_brei = fp_counts_brei + fp_count;
    
    [d d tp_count fp_count] = calculate_roc_curve(1-pb_isbe(:),label(:),(-1:100)/100);
    tp_counts_isbe = tp_counts_isbe + tp_count;
    fp_counts_isbe = fp_counts_isbe + fp_count;
    
    t_counts = t_counts + sum(label(:));
    f_counts = f_counts + sum(~label(:));
end
roc_pts_brei = [fp_counts_brei / f_counts tp_counts_brei / t_counts];
roc_pts_isbe = [fp_counts_isbe / f_counts tp_counts_isbe / t_counts];

figure; 
hold on;
plot(roc_pts_brei(:,1), roc_pts_brei(:,2));
plot(roc_pts_isbe(:,1), roc_pts_isbe(:,2), 'g');

auc_brie = sum( (roc_pts_brei(2:end,1)-roc_pts_brei(1:end-1,1)) .* roc_pts_brei(1:end-1,2)) + ...
        0.5*sum( (roc_pts_brei(2:end,1)-roc_pts_brei(1:end-1,1)) .* (roc_pts_brei(2:end,2)-roc_pts_brei(1:end-1,2)) );
    
auc_isbe = sum( (roc_pts_isbe(2:end,1)-roc_pts_isbe(1:end-1,1)) .* roc_pts_isbe(1:end-1,2)) + ...
        0.5*sum( (roc_pts_isbe(2:end,1)-roc_pts_isbe(1:end-1,1)) .* (roc_pts_isbe(2:end,2)-roc_pts_isbe(1:end-1,2)) );
    
%%
tic;
[training_data training_labels parameters] = sample_bar_training_data('num_samples', 2e5);
toc;

%%
%profile on;
sampling_method_args.num_samples = 2e5;
sampling_method_args.bg_dir = [mberksroot, 'classification\data\normal_smooth128\'];

forest_args.sampling_method = 'sample_bar_training_data';
forest_args.sampling_method_args = sampling_method_args;
forest_args.d = 10;
forest_args.n_trees = 100;
forest_args.tree_dir = [mberksroot, 'classification\rf\rf_detect_', zerostr(3,2), '_trees\'];
forest_args.save_path = [mberksroot, 'classification\rf\rf_detect_', zerostr(3,2), '.mat'];
random_forest = mb_random_forest_class_train(forest_args);
%profile viewer;
%%
for oo = 0:20:340
    row = 128; col = 128; halfwidth = 4; contrast = 16;
    [image_out, label, label_centre] = create_ellipse_spicule_mb(halfwidth, contrast, oo, row, col);
    
    figure;
    subplot(1,3,1); imagesc(image_out); axis image; hold on; 
    subplot(1,3,2); imagesc(label); axis image; hold on;
    subplot(1,3,3); imagesc(label_centre); axis image; hold on;
end
%%
random_forest = mb_combine_rfs('rf_dir', '\\isbe-san1\mberks\dev\classification\rf\',...
    'tree_dir', 'C:\isbe\dev\classification\rf\rf_centre_ns\', 'delete_trees', 0,...
    'replace_tree_root', '\\isbe-san1\mberks\dev\',...
    'save_path', 'C:\isbe\dev\classification\rf\rf_centre_ns.mat');
%%
random_forest2 = mb_combine_rfs('rf_dir', 'C:\isbe\dev\classification\rf\',...
    'tree_dir', 'C:\isbe\dev\classification\rf\rf_centre_s\', 'delete_trees', 0,...
    'save_path', 'C:\isbe\dev\classification\rf\rf_centre_s.mat');
%%
%Generate test images
generate_bar_images('num_images', 100, 'image_dir', 'C:\isbe\dev\classification\data\bar_bg_128\', 'plot', 0);
%%
for ii = 11:100
    load(['C:\isbe\dev\classification\data\bar_bg_128\bar', zerostr(ii,3), '.mat']);
    [probability_image] = classify_image(bar_image, random_forest);
    save(['C:\isbe\dev\classification\data\bar_bg_128_results_1\bar', zerostr(ii,3), '_results.mat'], 'probability_image');
%     figure;
%     subplot(1,3,1); imagesc(bar_image); axis image;
%     subplot(1,3,2); imagesc(probability_image); axis image;
    
    [probability_image] = classify_image(bar_image, random_forest2);
    save(['C:\isbe\dev\classification\data\bar_bg_128_results_2\bar', zerostr(ii,3), '_results.mat'], 'probability_image');
%     subplot(1,3,3); imagesc(probability_image); axis image;
end
%%
%%
[roc_pts1,auc1,tp_counts1,fp_counts1] = compute_roc_image_set(...
    [mberksroot, 'classification/data/bar_bg_128'],...
    [mberksroot, 'classification/data/bar_bg_128_results_1']);
%
[roc_pts2,auc2,tp_counts2,fp_counts2] = compute_roc_image_set(...
    [mberksroot, 'classification/data/bar_bg_128'],...
    [mberksroot, 'classification/data/bar_bg_128_results_2']);
%%
figure; 
hold on;
plot(roc_pts1(:,1), roc_pts1(:,2));
plot(roc_pts2(:,1), roc_pts2(:,2), 'g');
%%
random_forest1 = u_load('C:\isbe\dev\classification\rf\rf_centre_ns.mat');
random_forest2 = u_load('C:\isbe\dev\classification\rf\rf_centre_s.mat');

bg_list = dir('C:\isbe\dev\background\images\normal512\*.bmp');
for ii = 1:10
    
    figure;
    
    mammo_bg = double(imread(['C:\isbe\dev\background\images\normal512\' bg_list(ii).name]));
    [probability_image] = classify_image(mammo_bg, random_forest1);
    save(['C:\isbe\dev\classification\data\bg_512_results_1\bg', zerostr(ii,3), '_results.mat'], 'probability_image');
    subplot(1,2,1); imagesc(probability_image); axis image;
    
    [probability_image] = classify_image(mammo_bg, random_forest2);
    save(['C:\isbe\dev\classification\data\bg_512_results_2\bg', zerostr(ii,3), '_results.mat'], 'probability_image');
    subplot(1,2,2); imagesc(probability_image); axis image;
end
%%
bg_list = dir('C:\isbe\dev\background\images\normal512\*.bmp');
for ii = 1:10
    
    figure;
    
    mammo_bg = double(imread(['C:\isbe\dev\background\images\normal512\' bg_list(ii).name]));
    write_im_from_colormap(mammo_bg, ['C:\isbe\dev\classification\figures\norm', zerostr(ii,3), '.bmp'], gray(256));
  
    load(['C:\isbe\dev\classification\data\bg_512_results_1\bg', zerostr(ii,3), '_results.mat'], 'probability_image');
    write_im_from_colormap(1-probability_image, ['C:\isbe\dev\classification\figures\norm', zerostr(ii,3), '_prob1.bmp'], gray(256));
    
    load(['C:\isbe\dev\classification\data\bg_512_results_2\bg', zerostr(ii,3), '_results.mat'], 'probability_image');
    write_im_from_colormap(1-probability_image, ['C:\isbe\dev\classification\figures\norm', zerostr(ii,3), '_prob2.bmp'], gray(256));
end
%%
mass_list = dir('C:\isbe\dev\masses1024x1024\*.mat');
for ii = 1:length(mass_list)
    load(['C:\isbe\dev\masses1024x1024\', mass_list(ii).name]);
    mass_roi = imresize(double(mass.mass_ROI), .5, 'bicubic');
    save(['C:\isbe\dev\image_data\masses512x512\mass', zerostr(ii,3), '.mat'], 'mass_roi');
end
%%
norm_list = dir('C:\isbe\dev\background\images\normal1024\*.bmp');
for ii = 1:length(norm_list)
    norm_roi = double(imread(['C:\isbe\dev\background\images\normal1024\', norm_list(ii).name]));
    norm_roi = imresize(norm_roi, .5, 'bicubic');
    save(['C:\isbe\dev\image_data\normal_512\norm', zerostr(ii,3), '.mat'], 'norm_roi');
end
%%
for ii = [28 47]%length(mass_list)
    load(['C:\isbe\dev\image_data\masses512x512\mass', zerostr(ii,3), '.mat'], 'mass_roi');
    figure; imagesc(mass_roi); axis image; colormap(gray(256));
end
%% 
mass_roi1 = u_load('C:\isbe\dev\image_data\masses512x512\mass028.mat');
mass_roi2 = u_load('C:\isbe\dev\image_data\masses512x512\mass047.mat');
norm_roi1 = u_load('M:\chen\data\normal_512\bg006.mat');
norm_roi2 = u_load('M:\chen\data\normal_512\bg016.mat');
random_forest = u_load('C:\isbe\dev\classification\rf\rf_centre_ns.mat');
[mass_roi1_prob] = classify_image(mass_roi1, random_forest);
[mass_roi2_prob] = classify_image(mass_roi2, random_forest);
[norm_roi1_prob] = classify_image(norm_roi1, random_forest);
[norm_roi2_prob] = classify_image(norm_roi2, random_forest);
save('C:\isbe\dev\image_data\mass_roi1_prob.mat', 'mass_roi1_prob');
save('C:\isbe\dev\image_data\mass_roi2_prob.mat', 'mass_roi2_prob');
save('C:\isbe\dev\image_data\norm_roi1_prob.mat', 'norm_roi1_prob');
save('C:\isbe\dev\image_data\norm_roi2_prob.mat', 'norm_roi2_prob');
%%
figure; 
subplot(1,2,1); imagesc(1-norm_roi1_prob); axis image; colormap(gray(256));
subplot(1,2,2); imagesc(norm_roi1); axis image; colormap(gray(256));

figure;
subplot(1,2,1); imagesc(1-norm_roi2_prob); axis image; colormap(gray(256));
subplot(1,2,2); imagesc(norm_roi2); axis image; colormap(gray(256));

figure;
subplot(1,2,1); imagesc(1-mass_roi1_prob); axis image; colormap(gray(256));
subplot(1,2,2); imagesc(mass_roi1); axis image; colormap(gray(256));

figure;
subplot(1,2,1); imagesc(1-mass_roi2_prob); axis image; colormap(gray(256));
subplot(1,2,2); imagesc(mass_roi2); axis image; colormap(gray(256));

%%
for ii = 1:55
    f1 = openfig(['temp_figs\fig' zerostr(ii,2) '.fig']);
    set(f1, 'windowstyle', 'normal', 'position', get(0,'ScreenSize'));
    frame1 = getframe(f1);
    gif1 = frame2im(frame1);
    [gif1a map] = rgb2ind(gif1, 2^8);
    
    if ii == 1
        imwrite(gif1a, map, 'temp_figs\dt_responses.gif', 'gif', 'WriteMode', 'overwrite', 'Loop', Inf, 'DelayTime', 0.5);
    else
        imwrite(gif1a, map, 'temp_figs\dt_responses.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
    end
    close(f1);
end