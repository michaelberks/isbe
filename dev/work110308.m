
mkdir C:\isbe\asymmetry_project\data\synthetic_data\real512_dt\test\
for ii = 1:100
    test_image = u_load(['C:\isbe\asymmetry_project\data\synthetic_lines\real512\image' zerostr(ii,3) '.mat']);
    test_label = load(['C:\isbe\asymmetry_project\data\synthetic_lines\real512\labels\label' zerostr(ii,3) '.mat']);
    
    [rows cols] = find(test_label.label_centre & test_label.label < 2);
    idx = find(test_label.label_centre & test_label.label < 2);
    dt = compute_dual_tree(test_image, 5, 0);
            
    %Sample DT coefficients from specified rows and cols according to
    %sampling arguments
    X = sample_dt_data(dt, rows, cols, ...
        'feature_shape', 'rect',...
        'feature_type', 'complex',...
        'levels', 1:5,...
        'do_max', 0,...
        'rotate', 0,...
        'win_size', 3,...
        'pca', []);
    X = reshape(X, length(rows), 9, 6, 5);
    
    y = complex(cosd(test_label.label_orientation(idx)), sind(test_label.label_orientation(idx)));
    
    save(['C:\isbe\asymmetry_project\data\synthetic_data\real512_dt\test\X_' zerostr(ii,3) '.mat'], 'X');
    save(['C:\isbe\asymmetry_project\data\synthetic_data\real512_dt\test\y_' zerostr(ii,3) '.mat'], 'y');
end
%%
mask_list = dir('C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\*.mat');
for ii = 1:length(mask_list)
    mask = u_load(['C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\' mask_list(ii).name]);
    n_pts(ii) = sum(mask(:)) / 4;
end
%%
[errors_1] = compute_image_orientation_errors(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512',...
    'Z:\data\synthetic_lines\real512\results\286715\');

[errors_2] = compute_image_orientation_errors(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512',...
    'Z:\data\synthetic_lines\real512\results\306079\200\');

[errors_3] = compute_image_orientation_errors(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512',...
    'Z:\data\synthetic_lines\real512\results\305894\') ;
%
ori_errors_1 = sort(abs(errors_1(:,1)));
ori_errors_2 = sort(abs(errors_2(:,1)));
ori_errors_3 = sort(abs(errors_3(:,1)));

mag_errors_1 = sort(abs(errors_1(:,2)));
mag_errors_2 = sort(abs(errors_2(:,2)));
mag_errors_3 = sort(abs(errors_3(:,2)));

com_errors_1 = sort(prod(abs(errors_1),2) / mean(mag_errors_1));
com_errors_2 = sort(prod(abs(errors_2),2) / mean(mag_errors_2));
com_errors_3 = sort(prod(abs(errors_3),2) / mean(mag_errors_3));

ori_cdf_1 = zeros(101,1);
ori_cdf_2 = zeros(101,1);
ori_cdf_3 = zeros(101,1);

mag_cdf_1 = zeros(101,1);
mag_cdf_2 = zeros(101,1);
mag_cdf_3 = zeros(101,1);

com_cdf_1 = zeros(101,1);
com_cdf_2 = zeros(101,1);
com_cdf_3 = zeros(101,1);

for ii = 1:100
    x = ceil(ii*size(ori_errors_3, 1)/100);
    ori_cdf_1(ii+1) = ori_errors_1(x,1);
    ori_cdf_2(ii+1) = ori_errors_2(x,1);
    ori_cdf_3(ii+1) = ori_errors_3(x,1);

    mag_cdf_1(ii+1) = mag_errors_1(x,1);
    mag_cdf_2(ii+1) = mag_errors_2(x,1);
    mag_cdf_3(ii+1) = mag_errors_3(x,1);

    com_cdf_1(ii+1) = com_errors_1(x,1);
    com_cdf_2(ii+1) = com_errors_2(x,1);
    com_cdf_3(ii+1) = com_errors_3(x,1);
end
%
figure; hold on;
plot(ori_cdf_1, 0:100, 'g');
plot(ori_cdf_2, 0:100, 'r');
plot(ori_cdf_3, 0:100, 'b');
legend({...
    ['1x1 windows old: mean = ' num2str(mean(ori_errors_1)) ', median = ' num2str(median(ori_errors_1))],...
    ['1x1 windows: mean = ' num2str(mean(ori_errors_2)) ', median = ' num2str(median(ori_errors_2))],...
    ['3x3 windows: mean = ' num2str(mean(ori_errors_3)) ', median = ' num2str(median(ori_errors_3))]},...
    'location', 'southeast');

figure; hold on;
plot(mag_cdf_1, 0:100, 'g');
plot(mag_cdf_2, 0:100, 'r');
plot(mag_cdf_3, 0:100, 'b');
legend({...
    ['1x1 windows old: mean = ' num2str(mean(mag_errors_1)) ', median = ' num2str(median(mag_errors_1))],...
    ['1x1 windows: mean = ' num2str(mean(mag_errors_2)) ', median = ' num2str(median(mag_errors_2))],...
    ['3x3 windows: mean = ' num2str(mean(mag_errors_3)) ', median = ' num2str(median(mag_errors_3))]},...
    'location', 'southeast');

figure; hold on;
plot(com_cdf_1, 0:100, 'g');
plot(com_cdf_2, 0:100, 'r');
plot(com_cdf_3, 0:100, 'b');
legend({...
    ['1x1 windows old: mean = ' num2str(mean(com_errors_1)) ', median = ' num2str(median(com_errors_1))],...
    ['1x1 windows: mean = ' num2str(mean(com_errors_2)) ', median = ' num2str(median(com_errors_2))],...
    ['3x3 windows: mean = ' num2str(mean(com_errors_3)) ', median = ' num2str(median(com_errors_3))]},...
    'location', 'southeast');
%%
%forest_list = {'243311', '286715', '305891', '305894', '305897', '299400', '299402', '299407' '306079'};%,'299430','299442',
%forest_list = {'306060','306076','306079'};
line_type = 'centre_line';

forest_list = {'286713', '286714', '286715'};
ori_leg = cell(length(forest_list),1);
mag_leg = cell(length(forest_list),1);
com_leg = cell(length(forest_list),1);
pct_leg = cell(length(forest_list),1);
r_idx = ceil(437523*rand(1e4,1));
wstyle = 'docked';
figure('windowstyle', wstyle); a1 = axes; hold all; 
title('CDF of orientation error'); 
xlabel('Orientation error');
figure('windowstyle', wstyle); a2 = axes; hold all; 
title('CDF of dispersion magnitudes');
xlabel('Dispersion magnitudes');
figure('windowstyle', wstyle); a3 = axes; hold all; 
title('CDF of orientation errors weighted by dispersion magnitudes');
xlabel('Weighted orientation errors');
figure('windowstyle', wstyle); a4 = axes; hold all; 
title('Mean orientation error for the Nth percentile of dispersion magnitudes');
ylabel('Mean orientation error')
xlabel('Percentile of samples sorted by dispersion magnitude');
%%
for ii = 1:length(forest_list)
    try
        [errors] = compute_image_orientation_errors(...
            'C:\isbe\asymmetry_project\data\synthetic_lines\real512',...
            ['Z:\data\synthetic_lines\real512\results\' forest_list{ii} '\'], line_type);
%         errors = u_load(['Z:\data\synthetic_lines\real512\results\' forest_list{ii} '\errors\ori_errors.mat']);
    catch
        [errors] = compute_image_orientation_errors(...
            'C:\isbe\asymmetry_project\data\synthetic_lines\real512',...
            ['Z:\data\synthetic_lines\real512\results\' forest_list{ii} '\200\'], line_type);
%         errors = u_load(['Z:\data\synthetic_lines\real512\results\' forest_list{ii} '\200\errors\ori_errors.mat']);
    end
        
    ori_errors = sort(abs(errors(:,1)));
    mag_errors = sort(errors(:,2));
    com_errors = sort(prod(abs(errors),2) / mean(mag_errors));
    errors = sortrows(errors, -2);
    
    ori_cdf = zeros(101,1);
    mag_cdf = zeros(101,1);
    com_cdf = zeros(101,1);
    mean_pct = zeros(100,1);
    for jj = 1:100
        x = ceil(jj*size(ori_errors, 1)/100);
        ori_cdf(jj+1) = ori_errors(x,1);
        mag_cdf(jj+1) = mag_errors(x,1);
        com_cdf(jj+1) = com_errors(x,1);
        mean_pct(jj) = mean(abs(errors(1:x,1)));
    end
    
    plot(a1, ori_cdf, (0:100)/100, 'linewidth', 2);
    plot(a2, mag_cdf, (0:100)/100, 'linewidth', 2);
    plot(a3, com_cdf, (0:100)/100, 'linewidth', 2);
    plot(a4, 1:100, mean_pct, 'linewidth', 2);
    
    ori_leg{ii+3} = ...
        [forest_list{ii} ': (mean, median) = (' num2str(mean(ori_errors),4) ', ' num2str(median(ori_errors),4) ')'];
    mag_leg{ii+3} = ...
        [forest_list{ii} ': (mean, median) = (' num2str(mean(mag_errors),4) ', ' num2str(median(mag_errors),4) ')'];
    com_leg{ii+3} = ...
        [forest_list{ii} ': (mean, median) = (' num2str(mean(com_errors),4) ', ' num2str(median(com_errors),4) ')'];
    pct_leg{ii+3} = [forest_list{ii} ': mean at 50th pcntile = ' num2str(mean_pct(50),4)];
%     figure;
%     plot(errors(r_idx,1), errors(r_idx,2), 'rx');
%     title(['Orientations errors vs dispersion magnitudes. Forest ' forest_list{ii}]);
%     xlabel('Orientation errors');
%     ylabel('Dispersion magnitudes');
end
legend(a1, ori_leg, 'location', 'southeast');
legend(a2, mag_leg, 'location', 'southeast');
legend(a3, com_leg, 'location', 'southeast');
legend(a4, pct_leg, 'location', 'southeast');
%%
random_forest = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\305894\random_forest.mat');
random_forest.tree_root = 'C:\isbe\asymmetry_project\data\line_orientation_rfs\';
save('C:\isbe\asymmetry_project\data\line_orientation_rfs\305894\random_forest.mat', 'random_forest');

random_forest = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\306060\random_forest.mat');
random_forest.tree_root = 'C:\isbe\asymmetry_project\data\line_orientation_rfs\';
save('C:\isbe\asymmetry_project\data\line_orientation_rfs\306060\random_forest.mat', 'random_forest');

random_forest = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\306079\random_forest.mat');
random_forest.tree_root = 'C:\isbe\asymmetry_project\data\line_orientation_rfs\';
save('C:\isbe\asymmetry_project\data\line_orientation_rfs\306079\random_forest.mat', 'random_forest');
%%
roi = u_load('C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\mass_roi\024RCC_roi.mat');

rf_3_4 = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\243311\random_forest.mat');

args.image_in = roi;
args.forest_type = 'orientation';

args.sampling_args.do_max = 0;
args.sampling_args.rotate = 0;
args.sampling_args.use_nag = 0;
args.sampling_args.feature_shape = 'rect';
args.sampling_args.feature_type = 'conj';

%
args.forest = rf_3_4;
args.sampling_args.num_levels = 4;
args.sampling_args.win_size = 3;
ori_map_3_4 = classify_image(args);