tic;
[training_data training_oris parameters] = sample_orientation_training_data('num_samples', 2e4);
toc;
%%
[test_data test_oris parameters] = sample_orientation_training_data('num_samples', 5e4, 'num_levels', 5, 'idx_ratio', 0.05);
%test_oris_d = 180*angle(test_oris) / pi;
test_oris_d = test_oris;
%%
[tree] = mb_tree_reg_train(training_data, training_oris);
[y_fit, nodes] = mb_tree_predict(tree, test_data);
%%
sampling_method_args.num_samples = 2e4;
sampling_method_args.num_levels = 5;
sampling_method_args.bg_dir = [mberksroot, 'classification\data\normal_smooth512\'];

forest_args.sampling_method = 'sample_orientation_training_data';
forest_args.sampling_method_args = sampling_method_args;
forest_args.d = 20;
forest_args.n_trees = 200;
forest_args.tree_dir = [mberksroot, 'classification\rf\rf_reg_ori_', zerostr(1,2), '_trees\'];
forest_args.save_path = [mberksroot, 'classification\rf\rf_reg_ori_', zerostr(1,2), '.mat'];
random_forest = mb_random_forest_reg_train(forest_args);
%%
load C:\isbe\dev\classification\rf\rf_reg_ori_03.mat
[y_fit] = mb_random_forest_reg_predict(random_forest, test_data);
%test_oris_f = 180*angle(y_fit) / pi;
test_oris_f = y_fit;

max(abs(mb_mod(test_oris_f - test_oris_d, 180)))
mean(abs(mb_mod(test_oris_f - test_oris_d, 180)))
figure; hist(abs(mb_mod(test_oris_f - test_oris_d, 180)), 100);
%%
bg = u_load('C:\isbe\dev\background\images\normal_smooth128\bg037.mat');
for ori = 0:7:175
    [bar label] = create_ellipse_bar_mb(4, 5, ori, 128, 128);
    bar_ori = classify_image(bar+bg, random_forest, 'regression', 5);
    counts = hist(bar_ori(label), 0:179);
    [dummy max_ori] = max(counts);
    mean_ori = mean(bar_ori(label));
    display(['True angle = ', num2str(ori), ' mean angle = ', num2str(mean_ori), ' max angle = ', num2str(max_ori-1)]);
    
    figure; imagesc(bar_ori); axis image; colormap(hsv(180)); colorbar;
    title(['Bar angle = ', num2str(ori)]);
    
end
%%
load C:\isbe\dev\classification\rf\rf_reg_ori_01.mat
roi2 = u_load('C:\isbe\dev\image_data\masses512x512\mass028.mat');
roi2_prob = u_load('C:\isbe\dev\image_data\predict_masses512x512_chen\probability_image028.mat');
roi2_prob = 1 - roi2_prob;

%
[roi2_ori] = classify_image(...
    'image_in', roi2,...
    'forest', random_forest,...
    'forest_type', 'regression',...
    'num_levels', 5);
roi2_ori = pi*roi2_ori/180;
roi2_combined = roi2_prob .* exp(i*roi2_ori);

figure; 
imagesc(roi2_ori); axis image; colormap(hsv(180));

figure; image(complex2rgb(roi2_combined.*roi2_combined)); axis image; hold on;
quiver(1:8:512, (1:8:512)', cos(roi2_ori(1:8:512,1:8:512)), -sin(roi2_ori(1:8:512,1:8:512)), 'w');
%%
f2 = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', [50 50 512 512],...
    'PaperPositionMode','auto');
axes(...
    'Units', 'pixels',...
    'position', [0 0 512 512]); 

axis ij equal off; hold on;

image(complex2rgb(roi2_combined.*roi2_combined));
%print('-dtiff', '-noui', '-painters', f2, '-r300', 'C:\isbe\dev\classification\figures\mass_028_ori.tif');

quiver(12:8:500, (12:8:500)', cos(roi2_ori(12:8:500,12:8:500)), -sin(roi2_ori(12:8:500,12:8:500)), 'c');
print('-dtiff', '-noui', '-painters', f2, '-r300', 'C:\isbe\dev\classification\figures\mass_028_ori_label.tif');

%%
mass_roi2_nms = mb_non_maximal_supp(mass_roi2_prob, mass_roi2_ori);
figure; imagesc(mass_roi2_nms); axis image; colormap(gray(256));
%%
mass_roi_nms_thresh = mass_roi2_nms;
mass_roi_nms_thresh(mass_roi_nms_thresh < .3) = 0;
figure; imagesc(mass_roi_nms_thresh); axis image; colormap(gray(256));
%%
[mass_roi_nms_hyst] = hysterisis(mass_roi2_nms,[.3 .5]);
figure; imagesc(mass_roi_nms_hyst); axis image; colormap(gray(256));
%%
figure; imagesc(mass_roi2); axis image; colormap(gray(256)); hold on;
[line_pts_y line_pts_x] = find(mass_roi_nms_hyst);
plot(line_pts_x, line_pts_y, 'r.', 'markersize', 2);
%%
num_angles = 24; % number of orientations for detecting lines
min_scale = 5; % smallest scale for detecting linear segments;
max_scale = 15; % largest scale for detecting linear segments;
[mass_roi2_linop, mass_roi2_ori_linop] = line_operator_conv(mass_roi2, num_angles,min_scale, max_scale, 'degrees');
figure; imagesc(mass_roi2_linop); axis image; colormap(gray(256)); hold on;
quiver(1:8:512, (1:8:512)', cos(mass_roi2_ori_linop(1:8:512,1:8:512)), -sin(mass_roi2_ori_linop(1:8:512,1:8:512)), 'r');
figure; imagesc(mass_roi2_ori_linop); axis image; colormap(hsv(180));
mass_roi2_linop_nms = mb_non_maximal_supp(mass_roi2_linop, mass_roi2_ori_linop, 1);
figure; imagesc(mass_roi2_linop_nms); axis image; colormap(gray(256));
%%
random_forest = mb_combine_rfs('rf_dir', '\\isbe-san1\mberks\dev\classification\rf_reg_ori_181424\',...
    'tree_dir', 'C:\isbe\dev\classification\rf\rf_reg_ori_181424\', 'delete_trees', 0,...
    'replace_tree_root', '\\isbe-san1\mberks\dev\',...
    'save_path', 'C:\isbe\dev\classification\rf\rf_reg_ori_181424.mat');
%%
for ii = 1:30
    norm = u_load(['C:\isbe\dev\image_data\normal_512\bg', zerostr(ii,3), '.mat']);
    figure; imagesc(norm); axis image; colormap(gray(256));
end
%%
a = randn(1e4,1)/4 + .5;
a(a>1 | a<0) = [];
pdf = hist(a,100) / length(a);

cdf = cumsum(pdf);
theta = linspace(0, 1, length(cdf));
new_sample = interp1(cdf, theta, rand(1e4,1));

selection_probs = interp1(theta, pdf, mass_nms);

%
new_sample2 = zeros(1e4,1);
ii = 1;
while ii <= 1e4
    valid_idx = find(selection_probs >= rand);
    if ~isempty(valid_idx)
        new_sample2(ii) = mass_nms(valid_idx(ceil(rand*length(valid_idx))));
        ii = ii+1;
    end
end
%%
figure; hold on;
plot(1:100, hist(a,100), 'r');
plot(1:100, hist(new_sample,100), 'g');
plot(1:100, hist(new_sample2,100), 'b');