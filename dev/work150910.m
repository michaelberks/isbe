vessel_pred_c = load('C:\isbe\nailfold\data\rsa_study\master_set\detected_capillaries_cxx\10598c_vessels_v_pred.txt');
vessel_pred_m = u_load('C:\isbe\nailfold\data\rsa_study\master_set\predictions\detection\rf_classification\296655\10598c_pred.mat');

ori_pred_c = read_complex_txt('C:\isbe\nailfold\data\rsa_study\master_set\detected_capillaries_cxx\10598c_vessels_o_pred.txt');
ori_pred_m = u_load('C:\isbe\nailfold\data\rsa_study\master_set\predictions\orientation\rf_regression\296621\10598c_pred.mat');

width_pred_c = load('C:\isbe\nailfold\data\rsa_study\master_set\detected_capillaries_cxx\10598c_vessels_w_pred.txt');
width_pred_m = u_load('C:\isbe\nailfold\data\rsa_study\master_set\predictions\width\rf_regression\297037\10598c_pred.mat');
%%
figure;
a1 = subplot(2,1,1); imgray(vessel_pred_m);
a2 = subplot(2,1,2); imgray(vessel_pred_c);
linkaxes([a1 a2]);

figure;
a1 = subplot(2,1,1); imgray(complex2rgb(ori_pred_m));
a2 = subplot(2,1,2); imgray(complex2rgb(ori_pred_c));
linkaxes([a1 a2]);

figure;
a1 = subplot(2,1,1); imgray(width_pred_m);
a2 = subplot(2,1,2); imgray(width_pred_c);
linkaxes([a1 a2]);

figure;
a1 = subplot(2,1,1); imgray(vessel_pred_m);
a2 = subplot(2,1,2); imgray(vessel_pred_m - vessel_pred_c);
linkaxes([a1 a2]);

figure;
a1 = subplot(2,1,1); imgray(complex2rgb(ori_pred_m));
a2 = subplot(2,1,2); imgray(complex2rgb(ori_pred_m .* conj(ori_pred_c)));
linkaxes([a1 a2]);

figure;
a1 = subplot(2,1,1); imgray(width_pred_m);
a2 = subplot(2,1,2); imgray(width_pred_m - width_pred_c);
linkaxes([a1 a2]);
%%
vc_m = u_load('C:\isbe\nailfold\data\rsa_study\master_set\vessel_centres\full_centres\10598c_vc.mat');
vc_c = load('C:\isbe\nailfold\data\rsa_study\master_set\detected_capillaries_cxx1\10598c_vessel_centres.txt');
vc_cc.x = vc_c(:,1)+1;
vc_cc.y = vc_c(:,2)+1;
vc_cc.prob = vc_c(:,3);
vc_cc.ori = complex(vc_c(:,4), vc_c(:,5));
vc_cc.width = vc_c(:,6);


vc_in_m_not_c = setdiff([vc_m.x vc_m.y], vc_c(:,1:2)+1, 'rows');
vc_in_c_not_m = setdiff(vc_c(:,1:2)+1, [vc_m.x vc_m.y], 'rows');

figure;
a1 = subplot(2,1,1); imgray(vessel_pred_m);
plot(vc_in_m_not_c(:,1), vc_in_m_not_c(:,2), 'r.');

a2 = subplot(2,1,2); imgray(vessel_pred_c);
plot(vc_in_c_not_m(:,1), vc_in_c_not_m(:,2), 'r.');

linkaxes([a1 a2]);
%%
aoi_c = load('C:\isbe\nailfold\data\rsa_study\master_set\detected_capillaries_cxx1\10598c_apex_o_image.txt');
aoc_c = load('C:\isbe\nailfold\data\rsa_study\master_set\detected_capillaries_cxx1\10598c_apex_c_prediction.txt')';
aoo_c = load('C:\isbe\nailfold\data\rsa_study\master_set\detected_capillaries_cxx1\10598c_apex_o_prediction.txt');

aoi_m = load('C:\isbe\nailfold\data\rsa_study\master_set\apex_maps\set12g_half_296655(new)\10598c_pred.mat');

figure;
a1 = subplot(2,1,1); imgray(aoi_m.apex_offset_map); caxis([0 0.5]);
a2 = subplot(2,1,2); imgray(aoi_c); caxis([0 0.5]);
linkaxes([a1 a2]);
%%
[vc_both idx_m idx_c] = intersect([vc_m.x vc_m.y], vc_c(:,1:2)+1, 'rows');

idx_m2 = false(size(vc_m.x));
idx_m2(idx_m) = 1;

idx_c2 = false(size(vc_c,1),1);
idx_c2(idx_c) = 1;

figure; plot((aoi_m.apex_class_pred(idx_m)+aoc_c(idx_c))/2, aoi_m.apex_class_pred(idx_m) - aoc_c(idx_c), 'r.');
hold on;
plot([0 1], [0 0], 'k');

figure; plot((vc_m.prob(idx_m)+vc_c(idx_c,3))/2, vc_m.prob(idx_m) - vc_c(idx_c,3), 'r.');
figure; plot(vc_m.prob(idx_m), vc_c(idx_c,3), 'r.'); axis equal;

figure; plot((vc_m.width(idx_m)+vc_c(idx_c,6))/2, vc_m.width(idx_m) - vc_c(idx_c,6), 'r.');
figure; plot(vc_m.width(idx_m), vc_c(idx_c,6), 'r.'); axis equal;
%%
aox_c = aoo_c(idx_c,1);
aox_m = aoi_m.apex_offset_x_pred(idx_m);
aoy_c = aoo_c(idx_c,2);
aoy_m = aoi_m.apex_offset_y_pred(idx_m);
aoc_m = aoi_m.apex_class_pred(idx_m);
aoc_ci = aoc_c(idx_c);

figure; plot(aoc_m + rand(size(aoc_m))/100, aoc_ci+ rand(size(aoc_ci))/100, 'r.'); axis equal;
figure; plot(aox_m(aoc_m > 0.5), aox_c(aoc_m > 0.5), 'r.'); axis equal;
figure; plot(aoy_m(aoc_m > 0.5), aoy_c(aoc_m > 0.5), 'r.'); axis equal;
%%
rf_c = u_load('C:\isbe\nailfold\models\apex\classification\frog\rf.mat');
rf_x = u_load('C:\isbe\nailfold\models\apex\offset_x\set12g_half_296655\rf.mat');
rf_y = u_load('C:\isbe\nailfold\models\apex\offset_y\set12g_half_296655\rf.mat');

% rf_c = u_load('C:\isbe\nailfold\models\apex\classification\frog\rf.mat');
% rf_x = u_load('C:\isbe\nailfold\models\apex\offset_x\frog\rf.mat');
% rf_y = u_load('C:\isbe\nailfold\models\apex\offset_y\frog\rf.mat');

hog_args.cell_sz = [8 8];
hog_args.block_sz = [2 2];
hog_args.num_ori_bins = 9;
hog_args.norm_method = 'none';
hog_args.block_spacing = 8;
hog_args.gradient_operator = [-1 0 1];
hog_args.spatial_sigma = 0;
hog_args.angle_wrap = 1;

[apex_offset_map_m apex_class_pred_m apex_offset_x_pred_m apex_offset_y_pred_m] = ...
    predict_apex_offsets(...
    'apex_class_rf', rf_c,...
    'apex_offset_x_rf', rf_x,...
    'apex_offset_y_rf', rf_y,...
    'vessel_feature_im', vessel_pred_m,...
    'vessel_centre', vc_m,...
    'smoothing_sigma', 1,...
    'num_cells', 8,...
    'hog_args', hog_args,...
    'xy', [],...
    'apex_class_thresh', 0.5,...
    'separate_trees', 0,...
    'max_size', 1000,...
    'base_width', 20,...
    'include_pts', []);

[apex_offset_map_c apex_class_pred_c apex_offset_x_pred_c apex_offset_y_pred_c] = ...
    predict_apex_offsets(...
    'apex_class_rf', rf_c,...
    'apex_offset_x_rf', rf_x,...
    'apex_offset_y_rf', rf_y,...
    'vessel_feature_im', vessel_pred_c,...
    'vessel_centre', vc_cc,...
    'smoothing_sigma', 1,...
    'num_cells', 8,...
    'hog_args', hog_args,...
    'xy', [],...
    'apex_class_thresh', 0.5,...
    'separate_trees', 0,...
    'max_size', 1000,...
    'base_width', 20,...
    'include_pts', []);
%%
figure; plot(aoi_m.apex_class_pred, apex_class_pred_m, 'r.');
figure; plot(aoi_m.apex_offset_x_pred(aoi_m.apex_class_pred>0.5),...
    apex_offset_x_pred_m(aoi_m.apex_class_pred>0.5), 'r.');
figure; plot(aoi_m.apex_offset_y_pred(aoi_m.apex_class_pred>0.5),...
    apex_offset_y_pred_m(aoi_m.apex_class_pred>0.5), 'r.');

figure; plot(aoc_c, apex_class_pred_c, 'r.'); axis equal;
figure; plot(aoo_c(apex_class_pred_c>0.5,1), apex_offset_x_pred_c(apex_class_pred_c>0.5), 'r.'); axis equal;
figure; plot(aoo_c(apex_class_pred_c>0.5,2), apex_offset_y_pred_c(apex_class_pred_c>0.5), 'r.'); axis equal;
%%
v_hog = u_load('C:\isbe\nailfold\data\rsa_study\set12g_half\vessel_hogs\296655\enlargedapex0140_vessel_hog_part_0001.mat');
v_centre = u_load('C:\isbe\nailfold\data\rsa_study\set12g_half\vessel_centres\full_centres\enlargedapex0140_vessel_vc.mat');
v_pred = u_load('C:\isbe\nailfold\data\rsa_study\set12g_half\predictions\detection\rf_classification\296655\enlargedapex0140_vessel_pred.mat');
v_ori = u_load('C:\isbe\nailfold\data\rsa_study\set12g_half\predictions\orientation\rf_regression\296621\enlargedapex0140_vessel_pred.mat');
v_wid = u_load('C:\isbe\nailfold\data\rsa_study\set12g_half\predictions\width\rf_regression\297037\enlargedapex0140_vessel_pred.mat');
%%
compute_vessel_centre_hogs( ... % non-strict mode
    'task_id',              10, ...
    'num_jobs',             450, ...
    'data_dir',             [nailfoldroot 'data/rsa_study/set12g_half/'],...
    'feature_im_dir',       'predictions/detection/rf_classification/296655/',...
    'centre_dir',           'vessel_centres/full_centres/',...
    'hog_dir',              'bob/',...
    'max_size',             inf,...
    'smoothing_sigma',      1,...
    'num_cells',            8,...
    'cell_sz',              8,... %Size of HoG cells in blocks
    'block_sz',             [2 2],...%Size of blocks in cells
    'num_ori_bins',         9,... %Number of bins in orientation histograms
    'norm_method',          'none',... %Method for local normalisation
    'block_spacing',        8,... %Separation between blocks in pixels - controls level of overlap
    'gradient_operator',    [-1 0 1],...
    'spatial_sigma',        0, ...
    'angle_wrap',           1,...
    'base_width',           20, ...
    'make_parts_folders',   0);
%%
rf_c = u_load('C:\isbe\nailfold\models\apex\classification\frog\rf.mat');
rf_x = u_load('C:\isbe\nailfold\models\apex\offset_x\frog\rf.mat');
rf_y = u_load('C:\isbe\nailfold\models\apex\offset_y\frog\rf.mat');

% rf_c = u_load('C:\isbe\nailfold\models\apex\classification\frog\rf.mat');
% rf_x = u_load('C:\isbe\nailfold\models\apex\offset_x\frog\rf.mat');
% rf_y = u_load('C:\isbe\nailfold\models\apex\offset_y\frog\rf.mat');

hog_args.cell_sz = [8 8];
hog_args.block_sz = [2 2];
hog_args.num_ori_bins = 9;
hog_args.norm_method = 'none';
hog_args.block_spacing = 8;
hog_args.gradient_operator = [-1 0 1];
hog_args.spatial_sigma = 0;
hog_args.angle_wrap = 1;

vessel_pred_m = u_load('C:\isbe\nailfold\data\rsa_study\master_set\predictions\detection\rf_classification\296655\11053c_pred.mat');
vc_m = u_load('C:\isbe\nailfold\data\rsa_study\master_set\vessel_centres\full_centres\11053c_vc.mat');
aoi_m = load('C:\isbe\nailfold\data\rsa_study\master_set\apex_maps\frog\full_centres\11053c_pred.mat');

[apex_offset_map apex_class_pred apex_offset_x_pred apex_offset_y_pred] = ...
    predict_apex_offsets(...
    'apex_class_rf', rf_c,...
    'apex_offset_x_rf', rf_x,...
    'apex_offset_y_rf', rf_y,...
    'vessel_feature_im', vessel_pred_m,...
    'vessel_centre', vc_m,...
    'smoothing_sigma', 1,...
    'num_cells', 8,...
    'hog_args', hog_args,...
    'xy', [],...
    'apex_class_thresh', 0.5,...
    'separate_trees', 0,...
    'max_size', 1000,...
    'base_width', 20,...
    'include_pts', []);
figure; plot(aoi_m.apex_class_pred, apex_class_pred, 'r.'); axis equal;
figure; plot(aoi_m.apex_offset_x_pred, apex_offset_x_pred, 'r.'); axis equal;
figure; plot(aoi_m.apex_offset_y_pred, apex_offset_y_pred, 'r.'); axis equal;
%%
[apex_offset_map_cm] = ...
    transform_apex_offset_preds(aoc_c, aoo_c(:,1), aoo_c(:,2),...
        vc_cc, size(vessel_pred_c,1), size(vessel_pred_c,2), 20, aoc_c>0.5, 0);

%%    
%Compute local maxima and save
im_name = '22652c';
root_dir = 'C:\isbe\nailfold\data\rsa_study\master_set\';
nailfold_image = u_load([root_dir 'images\' im_name '.mat']);
nailfold_mask = u_load([root_dir 'fov_masks\' im_name '_f_mask.mat']);

vc_c = load([root_dir 'detected_capillaries_cxx\' im_name '_vessel_centres.txt']);
vessel_centre.x = vc_c(:,1)+1;
vessel_centre.y = vc_c(:,2)+1;
vessel_centre.prob = vc_c(:,3);
vessel_centre.ori = complex(vc_c(:,4), vc_c(:,5));
vessel_centre.width = vc_c(:,6);
clear vc_c;

vessel_pred_c = load([root_dir 'detected_capillaries_cxx\' im_name '_vessels_v_pred.txt']);
ori_pred_c = read_complex_txt([root_dir 'detected_capillaries_cxx\' im_name '_vessels_o_pred.txt']);
width_pred_c = load([root_dir 'detected_capillaries_cxx\' im_name '_vessels_w_pred.txt']);

vessel_predictions = cat(3, vessel_pred_c, ori_pred_c, width_pred_c);
clear vessel_pred_c ori_pred_c width_pred_c;
%
apex_offset_map = load([root_dir 'detected_capillaries_cxx\' im_name '_apex_o_image.txt']);
apex_class_pred = load([root_dir 'detected_capillaries_cxx\' im_name '_apex_c_prediction.txt'])';
aoo_c = load([root_dir 'detected_capillaries_cxx\' im_name '_apex_o_prediction.txt']);
apex_offset_x_pred = aoo_c(:,1);
apex_offset_y_pred = aoo_c(:,2);
clear aoo_c;

[candidate_xy candidate_scores] = ...
    local_image_maxima(apex_offset_map, 20, nailfold_mask, 0, 0);
    
save_path = [root_dir 'detected_capillaries_cxx\' im_name '_dc.mat'];

save(save_path, 'vessel_predictions', 'vessel_centre', 'apex_offset_map',...
            'apex_class_pred', 'apex_offset_x_pred', 'apex_offset_y_pred',...
            'candidate_xy', 'candidate_scores');
%%    
detect_capillaries(nailfold_image, 4:6, ...
    'nailfold_mask',	nailfold_mask,...
    'save_path',        save_path);
%%
load(save_path, 'candidates_hogs', 'candidate_oris', 'candidate_widths',...
            'candidate_rescores', 'candidate_displacements',...
            'selected_distal', 'selected_non_distal', 'candidate_class_probs',...
            'candidate_class', 'fell_at_the_last');
        
cxx_candidates = load([root_dir 'detected_capillaries_cxx\' im_name '_apex_candidates.txt']);
candidate_xy_c = cxx_candidates(:,1:2)+1;
candidate_scores_c = cxx_candidates(:,6);
candidate_oris_c = complex(cxx_candidates(:,4), cxx_candidates(:,5));
candidate_widths_c = cxx_candidates(:,3);
candidate_rescores_c = cxx_candidates(:,7);
candidate_displacements_c = cxx_candidates(:,8);
selected_distal_c = cxx_candidates(:,11);
selected_non_distal_c = cxx_candidates(:,12);
candidate_class_probs_c = cxx_candidates(:,10);
candidate_class_c = cxx_candidates(:,9);
merged_with_c = cxx_candidates(:,13);
fell_at_the_last_c = cxx_candidates(:,14);
%%
[candidates_in_m_not_c in_m_not_c_idx] = setdiff(candidate_xy, candidate_xy_c, 'rows');
[candidates_in_c_not_m in_c_not_m_idx] = setdiff(candidate_xy_c, candidate_xy, 'rows');
%%
v_hog = load('C:\isbe\nailfold\hog_responses.txt');


vessel_feature_patch_c = load('C:\isbe\nailfold\hog_patch.txt');
hog_args.cell_sz = [8 8];
hog_args.block_sz = [2 2];
hog_args.num_ori_bins = 12;
hog_args.norm_method = 'none';%'l1-sqrt';
hog_args.block_spacing = 8;
hog_args.gradient_operator = [-1 0 1];
hog_args.spatial_sigma = 0;
hog_args.angle_wrap = 1;
[v_hog_c] = compute_HoG(vessel_feature_patch_c, hog_args);
[v_hog_m] = compute_HoG(vessel_feature_patch, hog_args);
%%
rf_r = u_load('C:\isbe\nailfold\models\apex\rescoring\miccai_all\rf.mat');
[~,votes] = random_forest_class_predict(rf_r, v_hog);
%%
im_name = '22652c';
root_dir = 'C:\isbe\nailfold\data\rsa_study\master_set\';
nailfold_image = u_load([root_dir 'images\' im_name '.mat']);
nailfold_mask = u_load([root_dir 'fov_masks\' im_name '_f_mask.mat']);

vc_c = load([root_dir 'detected_capillaries_cxx\' im_name '_vessel_centres.txt']);
vessel_centre.x = vc_c(:,1)+1;
vessel_centre.y = vc_c(:,2)+1;
vessel_centre.prob = vc_c(:,3);
vessel_centre.ori = complex(vc_c(:,4), vc_c(:,5));
vessel_centre.width = vc_c(:,6);
clear vc_c;

vessel_pred_c = load([root_dir 'detected_capillaries_cxx\' im_name '_vessels_v_pred.txt']);
ori_pred_c = read_complex_txt([root_dir 'detected_capillaries_cxx\' im_name '_vessels_o_pred.txt']);
width_pred_c = load([root_dir 'detected_capillaries_cxx\' im_name '_vessels_w_pred.txt']);

vessel_predictions = cat(3, vessel_pred_c, ori_pred_c, width_pred_c);
clear vessel_pred_c ori_pred_c width_pred_c;
%
apex_offset_map = load([root_dir 'detected_capillaries_cxx\' im_name '_apex_o_image.txt']);
apex_class_pred = load([root_dir 'detected_capillaries_cxx\' im_name '_apex_c_prediction.txt'])';
aoo_c = load([root_dir 'detected_capillaries_cxx\' im_name '_apex_o_prediction.txt']);
apex_offset_x_pred = aoo_c(:,1);
apex_offset_y_pred = aoo_c(:,2);
clear aoo_c;

cxx_candidates = load([root_dir 'detected_capillaries_cxx\' im_name '_apex_candidates.txt']);
candidate_xy = cxx_candidates(:,1:2)+1;
candidate_scores = cxx_candidates(:,6);
candidate_oris = complex(cxx_candidates(:,4), cxx_candidates(:,5));
candidate_widths = cxx_candidates(:,3);
candidate_rescores = cxx_candidates(:,7);
candidate_displacements = cxx_candidates(:,8);

selected_distal_c = cxx_candidates(:,11)>0;
selected_non_distal_c = cxx_candidates(:,12)>0;
candidate_class_probs_c = cxx_candidates(:,10);
candidate_class_c = cxx_candidates(:,9);
merged_with_c = cxx_candidates(:,13);
fell_at_the_last_c = cxx_candidates(:,14)>0;

save_path = [root_dir 'detected_capillaries_cxx\' im_name '_dc.mat'];

save(save_path, 'vessel_predictions', 'vessel_centre', 'apex_offset_map',...
            'apex_class_pred', 'apex_offset_x_pred', 'apex_offset_y_pred',...
            'candidate_xy', 'candidate_scores', 'candidate_oris', 'candidate_widths',...
            'candidate_rescores', 'candidate_displacements');
%%    

detect_capillaries(nailfold_image, 6, ...
    'nailfold_mask',	nailfold_mask,...
    'save_path',        save_path);
load(save_path, 'selected_distal', 'selected_non_distal', 'candidate_class_probs',...
            'candidate_class', 'fell_at_the_last', 'merged_with');
        






