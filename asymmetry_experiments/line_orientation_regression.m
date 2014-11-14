%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Script for regressing orientation and applying non-maximal suppression
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Build a regression forest
sampling_method_args.num_samples = 2e4;
sampling_method_args.num_levels = 5;
sampling_method_args.bg_dir = [mberksroot, 'classification\data\normal_smooth512\'];

forest_args.sampling_method = 'sample_orientation_training_data';
forest_args.sampling_method_args = sampling_method_args;
forest_args.d = 20;
forest_args.n_trees = 200;
forest_args.mod = 180;
forest_args.tree_dir = [mberksroot, 'classification\rf\rf_reg_ori_', zerostr(2,2), '_trees\'];
forest_args.save_path = [mberksroot, 'classification\rf\rf_reg_ori_', zerostr(2,2), '.mat'];
random_forest = mb_random_forest_reg_train(forest_args);
%%
%--------------------------------------------------------------------------
% Test that the regression works

% 1. With a set of test data
load C:\isbe\dev\classification\rf\rf_reg_ori_01.mat
[test_data test_oris parameters] = sample_orientation_training_data('num_samples', 5e4, 'num_levels', 5);
test_oris_d = test_oris;

[test_oris_f] = mb_random_forest_reg_predict(random_forest, test_data);

max(abs(mb_mod(test_oris_f - test_oris_d, 180)))
mean(abs(mb_mod(test_oris_f - test_oris_d, 180)))
figure; hist(abs(mb_mod(test_oris_f - test_oris_d, 180)), 100);
%%
% 1. With a set of test data
load C:\isbe\dev\classification\rf\rf_reg_ori_02.mat
test_oris_d_2 = test_oris;

[test_oris_f_2] = mb_random_forest_reg_predict(random_forest, test_data);

max(abs(mb_mod(test_oris_f_2 - test_oris_d_2, 180)))
mean(abs(mb_mod(test_oris_f_2 - test_oris_d_2, 180)))
figure; hist(abs(mb_mod(test_oris_f_2 - test_oris_d_2, 180)), 100);
%%
% 2. By creating bars at a variety of angles
reg_forest = u_load('C:\isbe\dev\classification\rf\rf_reg_ori_02.mat');
class_forest = u_load('M:\chen\data\DTCWT_rf_fulltrees_W3L5\tree_combine\random_forest.mat');

bg = u_load('C:\isbe\dev\background\images\normal_smooth128\bg037.mat');
for ori = 0:36:175
    [bar label] = create_ellipse_bar(4, 5, ori, 128, 128, 64, 64);
    bar_ori = classify_image(...
        'image_in', bar+bg,...
        'forest', reg_forest,...
        'forest_type', 'regression',...
        'num_levels', 5,...
        'mask', [ones(128,64) zeros(128,64)]);
    
    bar_prob = classify_image(...
        'image_in', bar+bg,...
        'forest', class_forest,...
        'forest_type', 'isbe',...
        'num_levels', 5,...
        'mask', [ones(128,64) zeros(128,64)]);

    counts = hist(bar_ori(label), 0:179);
    [dummy max_ori] = max(counts);
    mean_ori = mean(bar_ori(label));
    display(['True angle = ', num2str(ori), ' mean angle = ', num2str(mean_ori), ' max angle = ', num2str(max_ori-1)]);
    
    figure; 
    subplot(1,2,1); imagesc(bar_ori); axis image; colormap(hsv(180)); colorbar;
    subplot(1,2,2); imagesc(bar_prob); axis image; colormap(hsv(180)); colorbar;
    title(['Bar angle = ', num2str(ori)]);
    
end
%%
%--------------------------------------------------------------------------
% Now try computing the orientations of a mass

%load the mass and line probability image
% mass_idx = 47;
% load C:\isbe\dev\classification\rf\rf_reg_ori_01.mat
% roi = u_load(['C:\isbe\dev\image_data\masses512x512\mass', zerostr(mass_idx,3), '.mat']);
% roi_prob = u_load(['C:\isbe\dev\image_data\predict_masses512x512_chen\probability_image', zerostr(mass_idx,3), '.mat']);
% roi_ori =
% u_load(['C:\isbe\dev\classification\data\masses512x512_orientation\mass_orientation', zerostr(mass_idx,3), '.mat']);
norm_idx = 16;
load C:\isbe\dev\classification\rf\rf_reg_ori_01.mat
roi = u_load(['M:\chen\data\line_detection_mammo\bg', zerostr(norm_idx,3), '.mat']);
roi_prob = u_load(['M:\chen\data\line_detection_mammo\bg', zerostr(norm_idx,3), '_line.mat']);
roi_prob = 1 - roi_prob;
roi_ori = u_load(['C:\isbe\dev\classification\data\normal_512_orientation\norm_orientation', zerostr(norm_idx,3), '.mat']);
%compute line orientations
% %[roi_ori] = classify_image(...
%     'image_in', roi,...
%     'forest', random_forest,...
%     'forest_type', 'regression',...
%     'num_levels', 5);


%combine the orientation and line strength into a single image
roi_ori = pi*roi_ori/180;
roi_combined = roi_prob .* exp(i*roi_ori);

%use the orientations to perform NMS on the line strength
roi_nms = mb_non_maximal_supp(roi_prob, roi_ori);

%perform hysterisis thresholding on the NMS image
[roi_nms_hyst] = hysterisis(roi_nms,[.3 .5]);
roi_final = roi_nms;
roi_final(~roi_nms_hyst) = 0;

%Combine final output with orientation into a single image
roi_combined_final = roi_final .* exp(i*roi_ori);
roi_combined_nms = roi_nms .* exp(i*roi_ori);

line_idx = find(roi_nms);
[line_y line_x] = ind2sub(size(roi), line_idx);

%Display the various images
figure; imagesc(roi); axis image; colormap(gray(180));
figure; imagesc(roi_prob); axis image; colormap(gray(180));
figure; imagesc(roi_ori); axis image; colormap(hsv(180));
figure; image(complex2rgb(roi_combined.^2)); axis image;
figure; image(complex2rgb(roi_combined.^2)); axis image; hold on;
quiver(1:8:512, (1:8:512)', cos(roi_ori(1:8:512,1:8:512)), -sin(roi_ori(1:8:512,1:8:512)), 'w');

figure; imagesc(roi_nms); axis image; colormap(gray(180));
figure; imagesc(roi_nms_hyst); axis image; colormap(gray(180));

figure; image(complex2rgb(roi_combined_nms .^2)); axis image;
figure; image(complex2rgb(roi_combined_nms .^2)); axis image; hold on;
quiver(line_x(1:2:end), line_y(1:2:end), real(roi_combined_nms(line_idx(1:2:end))), -imag(roi_combined_nms(line_idx(1:2:end))), 1.5, 'c');

figure; image(complex2rgb(roi_combined_final .^2)); axis image;
figure; image(complex2rgb(roi_combined_final .^2)); axis image; hold on;
quiver(line_x, line_y, real(roi_combined_final(line_idx)), -imag(roi_combined_final(line_idx)), 'w');


%
%Save image output
% write_im_from_colormap(roi, ['C:\isbe\dev\classification\figures\mass_', zerostr(mass_idx,3), '.bmp'], gray(256));
% imwrite(complex2rgb(roi_combined.^2), ['C:\isbe\dev\classification\figures\mass_', zerostr(mass_idx,3), '_ori.bmp']);
% write_im_from_colormap(roi_nms, ['C:\isbe\dev\classification\figures\mass_', zerostr(mass_idx,3), '_nms.bmp'], gray(256));
% write_im_from_colormap(roi_nms_hyst, ['C:\isbe\dev\classification\figures\mass_', zerostr(mass_idx,3), '_hyst.bmp'], gray(256));
% imwrite(complex2rgb(roi_combined_final.^2), ['C:\isbe\dev\classification\figures\mass_', zerostr(mass_idx,3), '_final.bmp']);
% imwrite(complex2rgb(roi_combined_nms.^2), ['C:\isbe\dev\classification\figures\mass_', zerostr(mass_idx,3), '_nms_ori.bmp']);

write_im_from_colormap(roi, ['C:\isbe\dev\classification\figures\norm_', zerostr(norm_idx,3), '.bmp'], gray(256));
imwrite(complex2rgb(roi_combined.^2), ['C:\isbe\dev\classification\figures\norm_', zerostr(norm_idx,3), '_ori.bmp']);
write_im_from_colormap(roi_nms, ['C:\isbe\dev\classification\figures\norm_', zerostr(norm_idx,3), '_nms.bmp'], gray(256));
write_im_from_colormap(roi_nms_hyst, ['C:\isbe\dev\classification\figures\norm_', zerostr(norm_idx,3), '_hyst.bmp'], gray(256));
imwrite(complex2rgb(roi_combined_final.^2), ['C:\isbe\dev\classification\figures\norm_', zerostr(norm_idx,3), '_final.bmp']);
imwrite(complex2rgb(roi_combined_nms.^2), ['C:\isbe\dev\classification\figures\norm_', zerostr(norm_idx,3), '_nms_ori.bmp']);

f = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', [50 50 512 512],...
    'PaperPositionMode','auto');
axes(...
    'Units', 'pixels',...
    'position', [0 0 512 512]); 

axis ij equal off; hold on;

image(complex2rgb(roi_combined.*roi_combined));
quiver(12:8:500, (12:8:500)', cos(roi_ori(12:8:500,12:8:500)), -sin(roi_ori(12:8:500,12:8:500)), 'c');
%print('-dtiff', '-noui', '-painters', f, '-r300', ['C:\isbe\dev\classification\figures\mass_', zerostr(mass_idx,3), '_ori_label.tif']);
print('-dtiff', '-noui', '-painters', f, '-r300', ['C:\isbe\dev\classification\figures\norm_', zerostr(norm_idx,3), '_ori_label.tif']);

f = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', [50 50 512 512],...
    'PaperPositionMode','auto');
axes(...
    'Units', 'pixels',...
    'position', [0 0 512 512]); 

axis ij equal off; hold on;

image(complex2rgb(roi_combined_nms .^2)); axis image; hold on;
quiver(line_x(1:4:end), line_y(1:4:end), real(roi_combined_nms(line_idx(1:4:end))), -imag(roi_combined_nms(line_idx(1:4:end))), 2, 'c');
%print('-dtiff', '-noui', '-painters', f, '-r300', ['C:\isbe\dev\classification\figures\mass_', zerostr(mass_idx,3), '_nms_ori_label.tif']);
print('-dtiff', '-noui', '-painters', f, '-r300', ['C:\isbe\dev\classification\figures\norm_', zerostr(norm_idx,3), '_nms_ori_label.tif']);

%%
%--------------------------------------------------------------------------
% Load in a forest we constructed on Hydra and compute orienation in all
% masses

load('C:\isbe\dev\classification\rf\rf_reg_ori_181424.mat') 

for ii = 1:179
    
    %load in mass
    mass_roi = u_load(['C:\isbe\dev\image_data\masses512x512\mass', zerostr(ii,3), '.mat']);

    %compute line orientations
    [mass_roi_orientation] = classify_image(...
        'image_in', mass_roi,...
        'forest', random_forest,...
        'forest_type', 'regression',...
        'num_levels', 5);
    
    save(['C:\isbe\dev\classification\data\masses512x512_orientation\mass_orientation', zerostr(ii,3), '.mat'], 'mass_roi_orientation');
    clear mass*
end
%%
%--------------------------------------------------------------------------
% Use the line orientation in combination to the line probability images to
% non-maximally suppress the lines
mkdir C:\isbe\dev\classification\data\masses512x512_nms\

for ii = 107:111
    
    %load in mass line strength and orientation
    mass_roi_prob = u_load(['C:\isbe\dev\image_data\predict_masses512x512_chen\probability_image', zerostr(ii,3), '.mat']);
    mass_roi_prob = 1 - mass_roi_prob;
    mass_roi_ori = u_load(['C:\isbe\dev\classification\data\masses512x512_orientation\mass_orientation', zerostr(ii,3), '.mat']);
    
    %Convert orientation to radians
    mass_roi_ori = pi*mass_roi_ori/180;

    %use the orientations to perform NMS on the line strength
    mass_roi_nms = mb_non_maximal_supp(mass_roi_prob, mass_roi_ori);
    
    %save the NMS line strength
    save(['C:\isbe\dev\classification\data\masses512x512_nms\mass_nms', zerostr(ii,3), '.mat'], 'mass_roi_nms');
    
end
    
    