%%
%Compute Matlab responses
im = imread('C:\isbe\nailfold\test_im.png');
g2dd = compute_gaussian_2nd_derivatives_d(im, 1, 3);
h2dd = compute_hilbert_2nd_derivatives_d(im, 1, 3);
g2di = steer_gaussian_2nd_derivatives(g2dd, [], 6);
h2di = steer_hilbert_2nd_derivatives(h2dd, [], 6);

[rr cc] = meshgrid(1:size(im,1), 1:size(im,2));
[g2ds] = interpolate_filter_responses(g2di, cc(:), rr(:), 'win_size', 1, 'interp_method', 'cubic');%bilinear
[h2ds] = interpolate_filter_responses(h2di, cc(:), rr(:), 'win_size', 1, 'interp_method', 'cubic');

switch_cols = [1:3:16 2:3:17 3:3:18];
g2ds = g2ds(:, switch_cols);
h2ds = h2ds(:, switch_cols);

mag_m = sqrt(g2ds.^2 + h2ds.^2);
%%
%Compute C responses
bob = load('C:/isbe/nailfold/im_output.txt');
mag_c = [bob(:,1:2:end) bob(:,2:2:end)];
%%
diff1 = mag_m(:,1:6) - mag_c(:,1:6);
diff2 = mag_m(:,7:12) - mag_c(:,7:12);
diff3 = mag_m(:,13:18) - mag_c(:,13:18);

mask1 = mag_m(:,1:6) == 0;
mask2 = mag_m(:,7:12) == 0;
mask3 = mag_m(:,13:18) == 0;

figure; hist(diff1(~mask1), -100:2:100);
figure; hist(diff2(~mask2), -100:2:100);
figure; hist(diff3(~mask3), -100:2:100);

display(['Mean absolute difference in level 1: ' num2str(mean(abs(diff1(~mask1))))]);
display(['Mean absolute difference in level 2: ' num2str(mean(abs(diff2(~mask2))))]);
display(['Mean absolute difference in level 3: ' num2str(mean(abs(diff3(~mask3))))]);

display(['Max absolute difference in level 1: ' num2str(max(abs(diff1(~mask1))))]);
display(['Max absolute difference in level 2: ' num2str(max(abs(diff2(~mask2))))]);
display(['Max absolute difference in level 3: ' num2str(max(abs(diff3(~mask3))))]);
%%
for i_col = 7:18
    im_m = reshape(mag_m(:,i_col), 256, 256);
    im_c = reshape(mag_c(:,i_col), 256, 256);
    
    figure; 
    subplot(1,3,1); imgray(im_m); colorbar('location', 'NorthOutside');
    subplot(1,3,2); imgray(im_c); colorbar('location', 'NorthOutside');
    subplot(1,3,3); imgray(im_m - im_c); colorbar('location', 'NorthOutside');
end
%%
lev = 2;
for i_col = 7:18
    im_m = reshape(mag_m(:,i_col), 256, 256);
    im_c = reshape(mag_c(:,i_col), 256, 256);   
    
    band = i_col - 6*(lev-1);
    
    pixel_wdth = 2^(lev-1);
    row = 64 / pixel_wdth;
    offset = (pixel_wdth - 1) / (2^lev);

    mag1 = sqrt(g2di{lev}(row,:,band).^2 + h2di{lev}(row,:,band).^2);
    
    figure; hold all;
    plot(offset + (1:pixel_wdth:256), mag1, 'r', 'linewidth', 2);
    
    plot(1:256, im_m(64,:), 'b--', 'linewidth', 2);
    plot(1:256, im_c(64,:), 'g-.', 'linewidth', 2);
    
    if rem(i_col, 6) == 0
        lev = lev + 1;
    end
end
%%
figure; hold all;
for i_col = [1 7 13]
    im_m = reshape(mag_m(:,i_col), 256, 256);
    im_c = reshape(mag_c(:,i_col), 256, 256);
    
    plot(im_m(50, 100:130));
    plot(im_c(50, 100:130));
end

lc = {'M1', 'C1', 'M2', 'C2', 'M3', 'C3'};
legend(lc);
%%
im1m = imresize(double(im), 0.5, 'lanczos2');
im2m = imresize(im1m, 0.5, 'lanczos2');

im1c = load('C:\isbe\nailfold\im_l1.txt');
im2c = load('C:\isbe\nailfold\im_l2.txt');
%%
figure; hold all;
plot(64:1:128, im(64, 64:128));
plot(64:2:128, im1m(32, 32:64));
plot(64:2:128, im1c(32, 32:64));
plot(64:4:128, im2m(16, 16:32));
plot(64:4:128, im2c(16, 16:32));
legend({'I0', 'M1', 'C1', 'M2', 'C2'});
figure;
subplot(1,2,1); imgray(im1m); colorbar('location', 'NorthOutside');
subplot(1,2,2); imgray(im1m - im1c); colorbar('location', 'NorthOutside');

figure;
subplot(1,2,1); imgray(im2m); colorbar('location', 'NorthOutside');
subplot(1,2,2); imgray(im2m - im2c); colorbar('location', 'NorthOutside');
%%
g2dc = cell(3,1);
for i_lev = 1:3
    g2dc{i_lev}(:,:,1) = load(['C:\isbe\nailfold\im_Ixy_' num2str(i_lev-1) '.txt']);
    g2dc{i_lev}(:,:,2) = load(['C:\isbe\nailfold\im_Ixx_' num2str(i_lev-1) '.txt']);
    g2dc{i_lev}(:,:,3) = load(['C:\isbe\nailfold\im_Iyy_' num2str(i_lev-1) '.txt']);
end
%
h2dc = cell(3,1);
for i_lev = 1:3
    h2dc{i_lev}(:,:,1) = load(['C:\isbe\nailfold\im_Ia_' num2str(i_lev-1) '.txt']);
    h2dc{i_lev}(:,:,2) = load(['C:\isbe\nailfold\im_Ib_' num2str(i_lev-1) '.txt']);
    h2dc{i_lev}(:,:,3) = load(['C:\isbe\nailfold\im_Ic_' num2str(i_lev-1) '.txt']);
    h2dc{i_lev}(:,:,4) = load(['C:\isbe\nailfold\im_Id_' num2str(i_lev-1) '.txt']);
end
%%
for i_lev = 1:3
    for i_band = 1:3
        figure;
        subplot(1,2,1); imgray(g2dc{i_lev}(:,:,i_band)); colorbar('location', 'NorthOutside');
        subplot(1,2,2); imgray(g2dc{i_lev}(:,:,i_band) - g2dd{i_lev}(:,:,i_band)); colorbar('location', 'NorthOutside');
    end
end
%
for i_lev = 1:3
    for i_band = 1:4
        figure;
        subplot(1,2,1); imgray(h2dc{i_lev}(:,:,i_band)); colorbar('location', 'NorthOutside');
        subplot(1,2,2); imgray(h2dc{i_lev}(:,:,i_band) - h2dd{i_lev}(:,:,i_band)); colorbar('location', 'NorthOutside');
    end
end
%%
vessel_patch = double(rgb2gray(imread('C:\isbe\nailfold\images\BII_comp\10598c_small.png')));
vessel_patch = imresize(vessel_patch, 0.5, 'lanczos2');

rf_ori = u_load('C:\isbe\nailfold\models\vessel\orientation\rf_regression\296621\predictor.mat');
rf_ori.tree_root = 'C:\isbe\nailfold\models\vessel\orientation\rf_regression/';
ori_args = u_load('C:\isbe\nailfold\models\vessel\orientation\rf_regression\296621\job_args.mat');

rf_det = u_load('C:\isbe\nailfold\models\vessel\detection\rf_classification\296655\predictor.mat');
rf_det.tree_root = 'C:\isbe\nailfold\models\vessel\detection\rf_classification/';
det_args = u_load('C:\isbe\nailfold\models\vessel\detection\rf_classification\296655\job_args.mat');

rf_wid = u_load('C:\isbe\nailfold\models\vessel\width\rf_regression\297037\predictor.mat');
rf_wid.tree_root = 'C:\isbe\nailfold\models\vessel\width\rf_regression/';
wid_args = u_load('C:\isbe\nailfold\models\vessel\width\rf_regression\297037\job_args.mat');
%%
tic;
[patch_ori] = predict_image(...
    'image_in', vessel_patch,...
    'decomposition_args', ori_args.decomposition_args,...
    'predictor', rf_ori, ...
    'prediction_type', 'rf_regression',...
    'output_type', 'orientation',...
    'use_probs', 0,...
    'mask', [],...
    'tree_mask', [], ...
    'num_trees', [], ...
    'max_size', 128,...
    'incremental_results', 0);
toc;
tic;
[patch_ori50] = predict_image(...
    'image_in', vessel_patch,...
    'decomposition_args', ori_args.decomposition_args,...
    'predictor', rf_ori, ...
    'prediction_type', 'rf_regression',...
    'output_type', 'orientation',...
    'use_probs', 0,...
    'mask', [],...
    'tree_mask', [], ...
    'num_trees', 50, ...
    'max_size', 128,...
    'incremental_results', 0);
toc;
%%
tic;
[patch_det] = predict_image(...
    'image_in', vessel_patch,...
    'decomposition_args', det_args.decomposition_args,...
    'predictor', rf_det, ...
    'prediction_type', 'rf_classification',...
    'output_type', 'detection',...
    'use_probs', 0,...
    'mask', [],...
    'tree_mask', [], ...
    'num_trees', [], ...
    'max_size', 128,...
    'incremental_results', 0);
toc;
%%
tic;
[patch_det50] = predict_image(...
    'image_in', vessel_patch,...
    'decomposition_args', det_args.decomposition_args,...
    'predictor', rf_det, ...
    'prediction_type', 'rf_classification',...
    'output_type', 'detection',...
    'use_probs', 0,...
    'mask', [],...
    'tree_mask', [], ...
    'num_trees', 50, ...
    'max_size', 128,...
    'incremental_results', 0);
toc;
tic;
[patch_det50p] = predict_image(...
    'image_in', vessel_patch,...
    'decomposition_args', det_args.decomposition_args,...
    'predictor', rf_det, ...
    'prediction_type', 'rf_classification',...
    'output_type', 'detection',...
    'use_probs', 1,...
    'mask', [],...
    'tree_mask', [], ...
    'num_trees', 50, ...
    'max_size', 128,...
    'incremental_results', 0);
toc;
%%
[patch_wid] = predict_image(...
    'image_in', vessel_patch,...
    'decomposition_args', wid_args.decomposition_args,...
    'predictor', rf_wid, ...
    'prediction_type', 'rf_regression',...
    'output_type', 'width',...
    'use_probs', 0,...
    'mask', [],...
    'tree_mask', [], ...
    'num_trees', [], ...
    'max_size', 128,...
    'incremental_results', 0);
%%
figure; imgray(complex2rgb(patch_ori));
figure; imgray(patch_det);
figure; imgray(two_map_colour(patch_det, patch_wid, [], [15 50]));
%%
[patch_det2] = predict_image(...
    'image_in', vessel_patch,...
    'decomposition_args', det_args.decomposition_args,...
    'predictor', rf_det, ...
    'prediction_type', 'rf_classification',...
    'output_type', 'detection',...
    'use_probs', 0,...
    'mask', [],...
    'tree_mask', [], ...
    'num_trees', [], ...
    'max_size', 128,...
    'incremental_results', 0);
%%
im = double(rgb2gray(imread('C:\isbe\nailfold\images\BII_comp\10598c_small.png')));

tic;
g2dd = compute_gaussian_2nd_derivatives_d(im, 4, 3);
h2dd = compute_hilbert_2nd_derivatives_d(im, 4, 3);
g2di = steer_gaussian_2nd_derivatives(g2dd, [], 6);
h2di = steer_hilbert_2nd_derivatives(h2dd, [], 6);

[rr cc] = meshgrid(1:size(im,1), 1:size(im,2));
[g2ds] = interpolate_filter_responses(g2di, rr(:), cc(:), 'win_size', 1, 'interp_method', 'cubic');%bilinear
[h2ds] = interpolate_filter_responses(h2di, rr(:), cc(:), 'win_size', 1, 'interp_method', 'cubic');

mag_m = sqrt(g2ds.^2 + h2ds.^2);
phase_m = atan2(h2ds, g2ds);
toc;
%%
[vessel_centre] = ...
    extract_vessel_centres(patch_det, patch_ori, patch_wid,...
    'prob_sigma',           2,...
    'ori_sigma',            2,...
    'width_sigma',            2,...
    'strong_vessel_thresh', 0.25,...
    'weak_vessel_thresh',   0);

[vessel_centre50] = ...
    extract_vessel_centres(patch_det50, patch_ori50, patch_wid,...
    'prob_sigma',           2,...
    'ori_sigma',            2,...
    'width_sigma',            2,...
    'strong_vessel_thresh', 0.25,...
    'weak_vessel_thresh',   0);
%%
base_dir = 'C:\isbe\nailfold\data\rsa_study\set12g_half\';
training_ims = dir([base_dir 'images\*.mat']);
training_masks = dir([base_dir 'vessel_masks\*.mat']);

max_mag_vals_v = zeros(450, 6, 5);
max_mag_vals_b = zeros(450, 6, 5);

for i_im = 451:900
    im = u_load([base_dir 'images\' training_ims(i_im).name]);
    mask = u_load([base_dir 'vessel_masks\' training_masks(i_im).name]);
    
    g2dd = compute_gaussian_2nd_derivatives_d(im, 1, 5);
    h2dd = compute_hilbert_2nd_derivatives_d(im, 1, 5);
    g2di = steer_gaussian_2nd_derivatives(g2dd, [], 6);
    h2di = steer_hilbert_2nd_derivatives(h2dd, [], 6);
    
    for i_level = 1:5
        if i_level > 1
            mask = imresize(mask, size(g2di{i_level}(:,:,1)));
        end
        
        for i_band = 1:6
            band = sqrt(g2di{i_level}(:,:,i_band).^2 + h2di{i_level}(:,:,i_band).^2);
            
            if any(mask(:))
                max_mag_vals_v(i_im-450, i_band, i_level) = max(band(mask));
            end
            if any(~mask(:))
                max_mag_vals_b(i_im-450, i_band, i_level) = max(band(~mask));
            end
        end
        
    end
end
%%
for i_level = 1:5
    figure;
        
    for i_band = 1:6
        subplot(2,3,i_band); hold all;
        
        plot(sort(max_mag_vals_b(:, i_band, i_level)), (1:450)/450, 'b', 'linewidth', 2);
        plot(sort(max_mag_vals_v(:, i_band, i_level)), (1:450)/450, 'r', 'linewidth', 2);
    end
end
%%
base_dir = 'C:\isbe\nailfold\data\rsa_study\set12g_half\';
training_ims = dir([base_dir 'images\*.mat']);
training_masks = dir([base_dir 'vessel_masks\*.mat']);
training_widths = dir([base_dir 'width_maps\*.mat']);
    
all_vessel_widths = [];

for i_im = 451:900
    width_map = u_load([base_dir 'width_maps\' training_widths(i_im).name]);
    mask = u_load([base_dir 'vessel_masks\' training_masks(i_im).name]);
    all_vessel_widths = [all_vessel_widths; width_map(mask)]; %#ok
end
%%
figure; 
subplot(1,2,1); hist(all_vessel_widths);
subplot(1,2,2); hist(log(all_vessel_widths));
%%
log_vessel_widths = unique(log(all_vessel_widths));
max_width = log_vessel_widths(end);
idx = 1:length(log_vessel_widths);

num_samples = 1e4;

for i_rp = 1:5
    samples = sort((max_width-0.5)*rand(1e4,1)+0.5);
    curr_idx = 1;
    sample_idx = zeros(num_samples,1);
    for i_s = 1:num_samples

        while(log_vessel_widths(curr_idx) < samples(i_s))
            curr_idx = curr_idx + 1;
        end
        sample_idx(i_s) = curr_idx;
    end
    figure; hist(log_vessel_widths(sample_idx));
end
%%
DECOMP_TYPE="{'g2dia','h2dia'}" NUM_ANGLES=6 WIN_SIZE=1 SIGMA_RANGE="[2 5]" OUTPUT_TYPE="orientation" PREDICTION_TYPE="rf_regression" NUM_TREES=10 IMAGE_TYPE="real" DATA_ROOT="scratch/nailfold/" IMAGE_ROOT="data/rsa_study/set12g_half" NUM_SAMPLES=50000 MODEL_ROOT="models/vessel" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="{'g2dia','h2dia'}" NUM_ANGLES=6 WIN_SIZE=1 SIGMA_RANGE="[2 5]" OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" NUM_TREES=10 IMAGE_TYPE="real" DATA_ROOT="scratch/nailfold/" IMAGE_ROOT="data/rsa_study/set12g_half" NUM_SAMPLES=50000 MODEL_ROOT="models/vessel" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="{'g2dia','h2dia'}" NUM_ANGLES=6 WIN_SIZE=1 SIGMA_RANGE="[2 5]" OUTPUT_TYPE="width" PREDICTION_TYPE="rf_regression" NUM_TREES=10 IMAGE_TYPE="real" DATA_ROOT="scratch/nailfold/" IMAGE_ROOT="data/rsa_study/set12g_half" NUM_SAMPLES=50000 MODEL_ROOT="models/vessel" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh

DATA_ROOT="scratch/nailfold/models/" MODEL_ROOT="vessel/orientation/rf_regression" MODEL_PATH="675752" MAKE_SAMPLED_MAPS=1 qsub -V -l twoday -hold_jid 675752 matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
DATA_ROOT="scratch/nailfold/models/" MODEL_ROOT="vessel/detection/rf_classification" MODEL_PATH="675753" MAKE_SAMPLED_MAPS=1 qsub -V -l twoday -hold_jid 675753 matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
DATA_ROOT="scratch/nailfold/models/" MODEL_ROOT="vessel/width/rf_regression" MODEL_PATH="675754" MAKE_SAMPLED_MAPS=1 qsub -V -l twoday -hold_jid 675754 matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh

DATA_ROOT="scratch/nailfold/" MODEL_ROOT="models/vessel/orientation/rf_regression" MODEL_PATH="675752" NUM_JOBS=50 IMAGE_ROOT="data/rsa_study/master_set" USE_SAMPLED_MAPS=0 OVERWRITE=0 qsub -l twoday -V -t 1-46 matlab_code/trunk/hydra/cuc/predict_image_set.sh
DATA_ROOT="scratch/nailfold/" MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="675753" NUM_JOBS=50 IMAGE_ROOT="data/rsa_study/master_set" USE_SAMPLED_MAPS=0 OVERWRITE=0 qsub -l twoday -V -t 1-46 matlab_code/trunk/hydra/cuc/predict_image_set.sh
DATA_ROOT="scratch/nailfold/" MODEL_ROOT="models/vessel/width/rf_regression" MODEL_PATH="675754" NUM_JOBS=50 IMAGE_ROOT="data/rsa_study/master_set" USE_SAMPLED_MAPS=0 OVERWRITE=0 qsub -l twoday -V -t 1-46 matlab_code/trunk/hydra/cuc/predict_image_set.sh

DATA_ROOT="scratch/nailfold/" MODEL_ROOT="models/apex" MODEL_PATH="set12g_half_296655" NUM_JOBS=455 IMAGE_ROOT="data/rsa_study/master_set" IMAGE_DIR="predictions/detection/rf_classification/675753" MASK_DIR="fov_masks" FG_MASK_DIR="vessel_centres/w1" OVERWRITE=0 THRESH=0.5 MODEL_NAME="rf" MAX_SIZE=1000 qsub -V -t 1 -l twoday matlab_code/trunk/hydra/cuc/predict_apex_offsets_set.sh

