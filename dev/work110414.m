ret_dir = [asymmetryroot 'data/retinograms/DRIVE/'];
mkdir([ret_dir 'training\images_extended\']);
mkdir([ret_dir 'test\images_extended\']);
for ii = 1:40
    display(['processing image ' num2str(ii)]);
    if ii < 21
        data = 'test';
    else
        data = 'training';
    end
    
    %Load in ground truth
    gt = logical(imread([ret_dir data '\1st_manual\' zerostr(ii,2) '_manual1.gif']));
    
    %Load in retinogram
    ret = imread([ret_dir data '\images\' zerostr(ii,2) '_' data '.tif']);
    
    %load in mask of edge of image
    mask = logical(imread([ret_dir data '\mask\' zerostr(ii,2) '_' data '_mask.gif']));
    mask_inner = imerode(mask, strel('disk', 10));
    mask_inner2 = imerode(mask, strel('disk', 20));
    mask_ring = mask_inner & ~mask_inner2;
    
    for ch = 1:3
        channel = ret(:,:,ch);
        
        %Process the retinogram to get rid of the hard edge
        [y_int x_int] = find(~mask_inner);
        [y x] = find(mask_ring);
        z = double(channel(mask_ring));
        z_int = griddata(x, y, z, x_int, y_int, 'nearest');
        channel(~mask_inner) = z_int;
        ret(:,:,ch) = channel;
    end
    save([ret_dir data '\images_extended\' zerostr(ii,2) '_' data '_ext.mat'], 'ret');
    
%     %Load in the orientation mask
%     ori_map = u_load([ret_dir data '\orientations\' zerostr(ii,2) '_ori1.mat']);
%     
%     %Get row and columns of pixels in ground truth
%     [rows cols] = find(gt);
%     
%     %Compute dual-tree of retinogram
%     dt = compute_dual_tree(ret, 6, 0);
%             
%     %Sample DT coefficients from specified rows and cols according to
%     %sampling arguments
%     X = sample_dt_data(dt, rows, cols, sampling_args);
%     X = reshape(X, size(X,1), 9, 6, 6);
%     
%     %Sample orientations
%     y = ori_map(gt);
%     
%     %Save X and y
%     save([asymmetryroot 'data/precomputed_data/drive/dt/' data '/X_' zerostr(ii,2) '.mat'], 'X');
%     save([asymmetryroot 'data/precomputed_data/drive/dt/' data '/y_' zerostr(ii,2) '.mat'], 'y');
    
end
%%
[training_data training_labels] = sample_vessel_dt_data(...
    'num_samples', 2e3,...
    'saved_data_dir', 'C:\isbe\asymmetry_project\data\synthetic_data\drive\dt\training',...
    'channel', 'all',...
    'win_size', 1,...
    'num_levels', 4,...
    'feature_type', 'conj',...
    'feature_shape', 'rect',...
    'rotate', 0,...
    'do_max', 0,...
    'data_list', [], ...
    'label_list', [], ...
    'data_name', 'X',...
    'label_name', 'y',...
    'save_path', []);
%%
for ii = 1:40
    display(['processing image ' num2str(ii)]);
    if ii < 21
        data = 'test';
    else
        data = 'training';
    end
    movefile(...
        [asymmetryroot 'data/synthetic_data/drive/dt/' data '/X_' zerostr(ii,2) '.mat'],...
        [asymmetryroot 'data/synthetic_data/drive/dt/' data '/X_rgb_' zerostr(ii,2) '.mat']);
    movefile(...
        ['Z:\data/synthetic_data/drive/dt/' data '/X_' zerostr(ii,2) '.mat'],...
        ['Z:\data/synthetic_data/drive/dt/' data '/X_rgb_' zerostr(ii,2) '.mat']);
end
%%
ret_dir = [asymmetryroot 'data/retinograms/DRIVE/'];

[X_rgb y] = generate_vessel_classification_data(... % the user's input
    'num_samples', 2e3,...
    'image_dir', [ret_dir 'training\images_extended\'],...
    'foveal_mask_dir', [ret_dir 'training\foveal_masks\'],...
    'vessel_mask_dir', [ret_dir 'training\vessel_masks\'],... % the mandatory arguments
    'selected_images', [], ...
    'rgb_channel', 'rgb',...
    'win_size', 1,...
    'num_levels', 4,...
    'feature_type', 'conj',...
    'feature_shape', 'rect',...
    'rotate', 0,...
    'do_max', 0,...
    'decomp_type', 'dt',...
    'bg_ratio', 1);
%%
ret_dir = [asymmetryroot 'data/retinograms/DRIVE/'];
mkdir([ret_dir 'training\foveal_masks\']);
mkdir([ret_dir 'training\vessel_masks\']);
mkdir([ret_dir 'test\foveal_masks\']);
mkdir([ret_dir 'test\vessel_masks\']);
for ii = 1:40
    display(['processing image ' num2str(ii)]);
    if ii < 21
        data = 'test';
    else
        data = 'training';
    end
    
    %Load in vessel and foeval mask, convert to logical and resave
    vessel_mask = logical(imread([ret_dir data '\1st_manual\' zerostr(ii,2) '_manual1.gif']));
    foveal_mask = logical(imread([ret_dir data '\mask\' zerostr(ii,2) '_' data '_mask.gif']));
    
    save([ret_dir data '\vessel_masks\' zerostr(ii,2) '_' data '_v_mask.mat'], 'vessel_mask');
    save([ret_dir data '\foveal_masks\' zerostr(ii,2) '_' data '_f_mask.mat'], 'foveal_mask');
end
%%
ret_dir = [asymmetryroot 'data/retinograms/DRIVE/'];
% v_rf_g_1 = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\313394\random_forest.mat');
% v_rf_g_1.tree_root = 'C:\isbe\asymmetry_project\data\line_orientation_rfs\';
% v_rf_g_3 = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\313395\random_forest.mat');
% v_rf_g_3.tree_root = 'C:\isbe\asymmetry_project\data\line_orientation_rfs\';
% v_rf_rgb = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\313393\random_forest.mat');
% v_rf_rgb.tree_root = 'C:\isbe\asymmetry_project\data\line_orientation_rfs\';
v_rf_rgb_3 = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\313534\random_forest.mat');
v_rf_rgb_3.tree_root = 'C:\isbe\asymmetry_project\data\line_orientation_rfs\';
v_rf_all = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\313392\random_forest.mat');
v_rf_all.tree_root = 'C:\isbe\asymmetry_project\data\line_orientation_rfs\';

args.forest_type = 'orientation';
args.sampling_args.do_max = 0;
args.sampling_args.rotate = 0;
args.sampling_args.use_nag = 0;
args.sampling_args.feature_shape = 'rect';
args.sampling_args.feature_type = 'conj';
args.sampling_args.num_levels = 4;

for ii = 18:20
    display(['processing image ' num2str(ii)]);
    
    %Load in vessel and foeval mask, convert to logical and resave
    %Load in retinogram
    ret = u_load([ret_dir 'test\images_extended\' zerostr(ii,2) '_test_ext.mat']);
    mask = u_load([ret_dir 'test\foveal_masks\' zerostr(ii,2) '_test_f_mask.mat']);
    
%     args.forest = v_rf_g_1;
%     args.sampling_args.win_size = 1;
%     args.image_in = ret(:,:,2);
%     ori_map = classify_image(args);
%     save([ret_dir 'test\orientations\rf\g_1\' zerostr(ii,2) '_ori.mat'], 'ori_map');
%     
%     args.forest = v_rf_g_3;
%     args.sampling_args.win_size = 3;
%     args.image_in = ret(:,:,2);
%     ori_map = classify_image(args);
%     save([ret_dir 'test\orientations\rf\g_3\' zerostr(ii,2) '_ori.mat'], 'ori_map');
%     
%     args.forest = v_rf_rgb;
%     args.sampling_args.win_size = 1;
%     args.image_in = rgb2gray(ret);
%     ori_map = classify_image(args);
%     save([ret_dir 'test\orientations\rf\rgb_1\' zerostr(ii,2) '_ori.mat'], 'ori_map');

    args.forest = v_rf_rgb_3;
    args.sampling_args.win_size = 3;
    args.image_in = rgb2gray(ret);
    ori_map = classify_image(args);
    save([ret_dir 'test\orientations\rf\rgb_3\' zerostr(ii,2) '_ori.mat'], 'ori_map');
    
    args.forest = v_rf_all;
    args.sampling_args.win_size = 1;
    args.image_in = ret;
    ori_map = classify_image_rgb(args);
    save([ret_dir 'test\orientations\rf\all_1\' zerostr(ii,2) '_ori.mat'], 'ori_map');
end
%%
%make labels
ret_dir = [asymmetryroot 'data/retinograms/DRIVE/'];
mkdir([ret_dir 'training\labels\']);
mkdir([ret_dir 'test\labels\']);
for ii = 1:40
    display(['processing image ' num2str(ii)]);
    if ii < 21
        data = 'test';
    else
        data = 'training';
    end
    label = u_load([ret_dir data '\vessel_masks\' zerostr(ii,2) '_' data '_v_mask.mat']);
    label_centre = bwmorph(label, 'skel', 'inf');
    label_orientation = u_load([ret_dir data '\orientations\' zerostr(ii,2) '_ori1.mat']);
    label_orientation = mod(90*angle(label_orientation)/pi, 180);
    save([ret_dir data '\labels\' zerostr(ii,2) '_label.mat'],...
        'label', 'label_centre', 'label_orientation');
end
%%
ret_dir = [asymmetryroot 'data/retinograms/DRIVE/'];
[ori_er_c_g1 mie_c_g1] =...
    compute_image_orientation_errors([ret_dir 'test\'], [ret_dir 'test\orientations\rf\g_1\'], 'centre_line', 0, 1);
[ori_er_c_g3 mie_c_g3] =...
    compute_image_orientation_errors([ret_dir 'test\'], [ret_dir 'test\orientations\rf\g_3\'], 'centre_line', 0, 1);
[ori_er_c_rgb mie_c_rgb] =...
    compute_image_orientation_errors([ret_dir 'test\'], [ret_dir 'test\orientations\rf\rgb_1\'], 'centre_line', 0, 1);
[ori_er_a_g1 mie_a_g1] =...
    compute_image_orientation_errors([ret_dir 'test\'], [ret_dir 'test\orientations\rf\g_1\'], 'all_line', 0, 1);
[ori_er_a_g3 mie_a_g3] =...
    compute_image_orientation_errors([ret_dir 'test\'], [ret_dir 'test\orientations\rf\g_3\'], 'all_line', 0, 1);
[ori_er_a_rgb mie_a_rgb] =...
    compute_image_orientation_errors([ret_dir 'test\'], [ret_dir 'test\orientations\rf\rgb_1\'], 'all_line', 0, 1);
%%
%--------------------------------------------------------------------------
%Now do the same for orientations
warning('off', 'load_uint8:missing_variables');

wstyle = 'docked';
label_type = 'all_line';

f1 = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', [50 50 800 800],...
    'PaperPositionMode','auto');
a1 = axes; hold all; 
title('\fontsize{18} \bf  CDF of orientation estimation errors'); 
xlabel('\fontsize{16} \bf Orientation error (degrees)');
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

ori_dir = {'g_1', 'g_3', 'rgb_1', 'rgb_3', 'all_1'};
test_label = {'G, 1x1', 'G, 3x3', 'RGB, 1x1', 'RGB, 3x3', 'All, 1x1'};

ori_leg = cell(length(ori_dir),1);
mag_leg = cell(length(ori_dir),1);
com_leg = cell(length(ori_dir),1);
pct_leg = cell(length(ori_dir),1);

ret_dir = [asymmetryroot 'data/retinograms/DRIVE/test/'];
%
for ii = 1:length(ori_dir)    
    %Save results
    errors = compute_image_orientation_errors(...
        ret_dir, [ret_dir 'orientations\rf\' ori_dir{ii} '\'], label_type);
    
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
        mean_pct(jj) = naNmean(abs(errors(1:x,1)));
    end
    
    plot(a1, ori_cdf, (0:100)/100, 'linewidth', 2);
    plot(a2, mag_cdf, (0:100)/100, 'linewidth', 2);
    plot(a3, com_cdf, (0:100)/100, 'linewidth', 2);
    plot(a4, 1:100, mean_pct, 'linewidth', 2);
    
    ori_leg{ii} = ...
        ['\fontsize{16} ' test_label{ii} ': (mean, median) = ('...
        num2str(naNmean(ori_errors),4) ', ' num2str(naNmedian(ori_errors),4) ')'];
    mag_leg{ii} = ...
        ['\fontsize{16} ' test_label{ii} ': (mean, median) = ('...
        num2str(naNmean(mag_errors),4) ', ' num2str(naNmedian(mag_errors),4) ')'];
    com_leg{ii} = ...
        ['\fontsize{16} ' test_label{ii} ': (mean, median) = ('...
        num2str(naNmean(com_errors),4) ', ' num2str(naNmedian(com_errors),4) ')'];
    pct_leg{ii} = ['\fontsize{16} ' test_label{ii} ': mean at 50th pcntile = ' num2str(mean_pct(50),4)];
end

legend(a1, ori_leg, 'location', 'southeast');
legend(a2, mag_leg, 'location', 'southeast');
legend(a3, com_leg, 'location', 'southeast');
legend(a4, pct_leg, 'location', 'southeast');
%%
ret_dir = [asymmetryroot 'data/retinograms/DRIVE/'];
mkdir([ret_dir 'training\labels\']);
mkdir([ret_dir 'test\labels\']);
data = 'test';
for ii = 1:20
    display(['processing image ' num2str(ii)]);

    load([ret_dir 'test\orientations\rf\g_1\' zerostr(ii,2) '_ori.mat'], 'ori_map');
    imwrite(complex2rgb(ori_map.^2),...
        [ret_dir 'test\display_images\orientation_maps\rf\' zerostr(ii,2) '_ori.bmp']);
    
end
%%
for ii = 19:20
    classify_image_set(...
        '313534',...
        'retinograms/DRIVE/test/images_extended', ... % non-strict mode
        'task_id', ii, ...
        'num_jobs', 20, ...
        'forest_dir', 'line_orientation_rfs', ...
        'mask_dir', 'retinograms/DRIVE/test/foveal_masks',...
        'use_nag', 0);
end
for ii = 18:20
    classify_image_set(...
        '313392',...
        'retinograms/DRIVE/test/images_extended', ... % non-strict mode
        'task_id', ii, ...
        'num_jobs', 20, ...
        'forest_dir', 'line_orientation_rfs', ...
        'mask_dir', 'retinograms/DRIVE/test/foveal_masks',...
        'use_nag', 0);
end
%%
mkdir Z:\data\retinograms\DRIVE\test\vessel_maps
for ii = 1:20
    copyfile(...
        ['Z:\data\retinograms\DRIVE\test\images_extended\results\313538\' zerostr(ii,2) '_test_ext_class.mat'],...
        ['Z:\data\retinograms\DRIVE\test\vessel_maps\' zerostr(ii,2) '_test_vessels.mat']);
end
%%
ret_dir = [asymmetryroot 'data/retinograms/DRIVE/test/'];
ori_errors = compute_image_orientation_errors(...
    ret_dir, [ret_dir 'orientations/rf/rgb_3/'], 'centre_line');
vessel_class = get_image_classifications(...
    ret_dir, [ret_dir 'vessel_maps/rf/rgb_3/'], 'centre_line');

figure('windowstyle', 'normal'); 
plot(ori_errors(r_idx(1:1e4),2), vessel_class(r_idx(1:1e4)), 'r.'); 
axis([0 1 0 1]); axis equal;
title('Scatter plot of classification probability vs dispersion magnitude at vessel centres');
xlabel('Dispersion magnitude');
ylabel('Classification probability');