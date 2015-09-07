sample_data = load('C:\isbe\nailfold\models\vessel\sample_data.txt');
mask_dir = 'C:\isbe\nailfold\data\rsa_study\cxx2\vessel_masks\';
mask_list = [dir([mask_dir 'r_*.png']); dir([mask_dir 'e*.png']); dir([mask_dir 'g*.png']); dir([mask_dir 'n*.png']);];

d = size(sample_data,2) - 4;
%%
for i_im = 0 + (1:10)
    
    mask = imread([mask_dir mask_list(i_im).name]);
    im_pts = sample_data(:,d+1) == i_im - 1;
    x = sample_data(im_pts,d+2) + 1;
    y = sample_data(im_pts,d+3) + 1;
    c_s = sample_data(im_pts,d+4) > 0;
    
    idx = sub2ind([256 256], y, x);
    c_i = mask(idx) > 0;
    
    figure; imgray(mask);
    plot(x(c_s), y(c_s), 'rx');
    plot(x(~c_s), y(~c_s), 'gx');
    
    if any(c_i ~= c_s)
        
        display([num2str(i_im) ' = bugger!']);
    end
end
%%
sample_data = load('C:\isbe\nailfold\models\apex\sample_data.txt');
mask_dir = 'C:\isbe\nailfold\data\rsa_study\cxx2\vessel_masks\';
mask_list = [dir([mask_dir 'e*.png']); dir([mask_dir 'g*.png']); dir([mask_dir 'n*.png']);];

d = size(sample_data,2) - 4;
%
for i_im = 0+(1:10)
    
    mask = imread([mask_dir mask_list(i_im).name]);
    im_pts = sample_data(:,d+1) == i_im - 1;
    x = sample_data(im_pts,d+2) + 1;
    y = sample_data(im_pts,d+3) + 1;
    c_s = sample_data(im_pts,d+4) > 0;
    
    idx = sub2ind([256 256], y, x);
    c_i = mask(idx) > 0;
    
    figure; imgray(mask);
    plot(x(c_s), y(c_s), 'rx');
    plot(x(~c_s), y(~c_s), 'gx');
    
    %if any(c_i ~= c_s)
    %    
    %    display([num2str(i_im) ' = bugger!']);
    %end
end
%%
sample_data = load('C:\isbe\nailfold\models\apex\sample_data.txt');
mask_dir = 'C:\isbe\nailfold\data\rsa_study\cxx2\vessel_masks\';
mask_list = [dir([mask_dir 'e*.png']); dir([mask_dir 'g*.png']); dir([mask_dir 'n*.png']);];

d = size(sample_data,2) - 5;
%
for i_im = (1:10)
    
    mask = imread([mask_dir mask_list(i_im).name]);
    im_pts = sample_data(:,d+1) == i_im - 1;
    x = sample_data(im_pts,d+2) + 1;
    y = sample_data(im_pts,d+3) + 1;
    offset = sample_data(im_pts,d+4:5);
    
    
    figure; imgray(mask);
    plot(x, y, 'rx');
    
    %if any(c_i ~= c_s)
    %    
    %    display([num2str(i_im) ' = bugger!']);
    %end
end
%%
prob_dir = 'C:\isbe\nailfold\data\rsa_study\cxx2\predictions\detection\rf_classification\296655\';
prob_list = [dir([prob_dir 'e*.txt']); dir([prob_dir 'g*.txt']); dir([prob_dir 'n*.txt']);];

for i_im = 1:10
    prob_im = load([prob_dir prob_list(i_im).name]);
    
    im_pts = sample_data(:,d+1) == i_im - 1;
    x = sample_data(im_pts,d+2) + 1;
    y = sample_data(im_pts,d+3) + 1;
end
%%
prob_dir = 'C:\isbe\nailfold\data\rsa_study\cxx2\predictions\detection\rf_classification\296655\';
vessel_prob = load([prob_dir 'enlargedapex0140_vessel_pred.txt']);
vessel_class = load('C:\isbe\nailfold\data\rsa_study\cxx2\vessel_hogs\enlargedapex0140_vessel_apex_c_prediction.txt');
vessel_centres = load('C:\isbe\nailfold\data\rsa_study\cxx2\vessel_hogs\enlargedapex0140_vessel_apex_v_centres.txt');
vessel_class_im = zeros(size(vessel_prob));
idx = sub2ind(size(vessel_prob), vessel_centres(:,2)+1, vessel_centres(:,1)+1);
vessel_class_im(idx) = vessel_class;
figure;
subplot(1,2,1); imgray(vessel_prob);
plot(vessel_centres(:,1), vessel_centres(:,2), 'r.');
subplot(1,2,2); imgray(vessel_class_im);
%%
predict_apex_offsets_set(...
    'task_id',              10, ...
    'num_jobs',             450, ...
    'model_id',             'set12g_half_296655\',...
    'model_root',           [nailfoldroot,'models/apex'], ...
    'model_name',           'rf',...
    'data_dir',             'C:\isbe\nailfold\data\rsa_study\cxx2\',...
    'feature_im_dir',       'predictions\detection\rf_classification\296655\',...
    'feature_im_ext',       '*.txt',...
    'fov_mask_dir',         [],...
    'centre_dir',           'vessel_centres',...
    'apex_map_dir',         'apex_maps',...
    'max_size',             1000,...
    'separate_trees',       0,...
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
    'base_width',           20,...
    'apex_class_thresh',    0,...
    'overwrite',            1);
%%
vessel_prob = load('C:\isbe\nailfold\data\rsa_study\cxx2\predictions\detection\rf_classification\296655\enlargedapex0140_vessel_pred.txt');
load('C:\isbe\nailfold\data\rsa_study\cxx2\apex_maps\set12g_half_296655\enlargedapex0140_vessel_pred.mat')
load('C:\isbe\nailfold\data\rsa_study\cxx2\vessel_centres\enlargedapex0140_vessel_vc.mat');
vessel_class_im = zeros(size(vessel_prob));
idx = sub2ind(size(vessel_prob), vessel_centre.y, vessel_centre.x);
vessel_class_im(idx) = apex_class_pred;
figure;
imgray(vessel_class_im);
figure;
imgray(vessel_prob);
    
%%
args.vessel_feature_im = vessel_prob;
args.vessel_centre.x = vessel_centres(:,1);
args.vessel_centre.y = vessel_centres(:,2);
args.vessel_centre.prob = vessel_centres(:,3);
args.vessel_centre.ori = complex(vessel_centres(:,4), vessel_centres(:,5));
args.vessel_centre.width = vessel_centres(:,6);

args.smoothing_sigma = 2;
args.num_cells = 8;
args.hog_args.cell_sz = [8 8];
args.hog_args.num_ori_bins = 9;
args.apex_class_thresh = 0.5;
args.base_width = 20;

patch_sz = args.num_cells*args.hog_args.cell_sz(1);
patch_sz = patch_sz + 2; %Account for padding
patch_sz2 = (patch_sz - 1)/2;

%Set up x,y coordinates for template patch
x = repmat(-patch_sz2:patch_sz2, patch_sz, 1);
y = x';
args.xy = [x(:) y(:)];

hog_sz = args.hog_args.num_ori_bins * args.num_cells^2;
num_pts = length(args.vessel_centre.x);

num_preds = 1;
include_pts = true(num_pts,1);

apex_class_pred = zeros(num_pts,1);
apex_offset_x_pred = zeros(num_pts,num_preds);
apex_offset_y_pred = zeros(num_pts,num_preds);

%Create smoothing kernel for feature image
if args.smoothing_sigma 
    if length(args.smoothing_sigma) == 1
        g = gaussian_filters_1d(args.smoothing_sigma);
        g = g / sum(g);
    else
        g = args.smoothing_sigma;
    end
    
    %Smooth the vessel probs
    args.vessel_feature_im = conv2(g', g, args.vessel_feature_im, 'same');
end

vessel_hog = zeros(num_pts, hog_sz);
for i_pt = 1:num_pts

    %Get predicted scale and orientation at this point
    vxc = args.vessel_centre.x(i_pt);
    vyc = args.vessel_centre.y(i_pt);
    ori_c = angle(args.vessel_centre.ori(i_pt))/2;
    width_c = args.vessel_centre.width(i_pt);

    %Get scale relative to base width a make rotation matrix
    rot = [cos(ori_c) -sin(ori_c); sin(ori_c) cos(ori_c)];
    scale = width_c / args.base_width;

    %Transform points given scale and angle and translate to
    %candidate position
    xya = args.xy * rot * scale;
    xa = reshape(xya(:,1) + vxc, [patch_sz patch_sz]);
    ya = reshape(xya(:,2) + vyc, [patch_sz patch_sz]);

    %Sample vessel prob patch
    vessel_feature_patch = interp2(args.vessel_feature_im, xa, ya, '*linear', 0);
    [hog] = compute_HoG(vessel_feature_patch, args.hog_args);       
    vessel_hog(i_pt,:) = hog(:)';
end

% 
%         %Classify whether points in part can point to an apex
%         [~,votes] = random_forest_class_predict(args.apex_class_rf, vessel_hog);
%         apex_class_pred_i = votes(:,2) / length(args.apex_class_rf.trees); clear votes;
% 
%         %Threshold these class scores
%         include_pts_i = apex_class_pred_i > args.apex_class_thresh;
% 
%         %Save these predictions in the main containers
%         apex_class_pred(part_idx,:) = apex_class_pred_i; 
%         include_pts(part_idx,:) = include_pts(part_idx,:) & include_pts_i; 
% 
%         %Now predict the apex offsets for points above the threshold
%             apex_offset_x_pred(part_idx(include_pts_i),:) = ...
%                 random_forest_reg_predict(args.apex_offset_x_rf, vessel_hog(include_pts_i,:));
%             apex_offset_y_pred(part_idx(include_pts_i),:) = ...
%                 random_forest_reg_predict(args.apex_offset_y_rf, vessel_hog(include_pts_i,:));
%%

%hog_patch = double(x.^2 + y.^2 < 25^2);
[hog_mat] = compute_HoG(hog_patch, 'num_ori_bins', 9, 'spatial_sigma', 2);
vessel_hog = hog_mat(:)';

angles_c = linspace(-pi/2, pi/2, 10);
angles_m = linspace(0, pi, 10);
nx = 8; ny = 8; na = 9;

plot_n = 1;
figure;
for i_r = 0:7
    for i_c = 0:7
        
        if (plot_n == 65)
            figure; plot_n = 1;
        end
        
        subplot(8,8,plot_n); hold all; axis off; 
        plot_n = plot_n + 1;
        for i_o = 0:8
            
            col_c = na*(i_r*nx + i_c) + i_o + 1;
            col_m = ny*(i_o*nx + i_c) + i_r + 1;
            
            u_c = test_data(81,col_m)*[-cos(angles_c(i_o+1)) cos(angles_c(i_o+1))];
            v_c = test_data(81,col_m)*[sin(angles_c(i_o+1)) -sin(angles_c(i_o+1))];
            
            u_m = vessel_hog(1,col_m)*[-cos(angles_m(i_o+1)) cos(angles_m(i_o+1))];
            v_m = vessel_hog(1,col_m)*[sin(angles_m(i_o+1)) -sin(angles_m(i_o+1))];
            
            plot(u_c, v_c, 'r--');
            plot(u_m, v_m, 'g-');
        end
    end
end
%%
x = test_data(:,577) + 1;
y = test_data(:,578) + 1;
vessel_class = load('C:\isbe\nailfold\data\rsa_study\cxx2\vessel_hogs\enlargedapex0140_vessel_apex_c_prediction.txt');

idx = sub2ind([256 256], y, x);
vessel_class_im = zeros(256);
vessel_class_im(idx) = 1 - vessel_class;
figure; imgray(vessel_class_im);
%%
na = 9;
colours = jet(na);
theta = linspace(-pi/2, pi/2);

sa = 1/pi;
%sa = 1/(2*pi);
bin_idx = mod(floor(na*(theta*sa + 1.5)),na);


figure; axis equal; hold on;
for i_a = 0:na-1
    idx = bin_idx == i_a;
    x = cos(theta(idx));
    y = sin(theta(idx));
    plot(x, y, '.', 'markeredgecolor', colours(i_a+1,:));
end
%%
apex_class_pred = random_forest_reg_predict(apex_class_rf, test_data(:,1:576));
apex_offset_x_pred = random_forest_reg_predict(apex_offset_x_rf, test_data(:,1:576));
apex_offset_y_pred = random_forest_reg_predict(apex_offset_y_rf, test_data(:,1:576));

apex_class_predc = load('C:\isbe\nailfold\data\rsa_study\cxx2\vessel_hogs\enlargedapex0140_vessel_apex_c_prediction.txt');
apex_class_predo = load('C:\isbe\nailfold\data\rsa_study\cxx2\vessel_hogs\enlargedapex0140_vessel_apex_o_prediction.txt');

figure; 
subplot(1,3,1); plot(apex_class_pred, 1-apex_class_predc, 'r.'); axis equal;
subplot(1,3,2); plot(apex_offset_y_pred, apex_class_predo(:,2), 'r.'); axis equal;
subplot(1,3,3); plot(apex_offset_x_pred, apex_class_predo(:,1), 'r.'); axis equal;

x = test_data(:,577) + 1;
y = test_data(:,578) + 1;

idx = sub2ind([256 256], y, x);
vessel_class_im = zeros(256);
vessel_class_im(idx) = apex_class_pred;
figure; imgray(vessel_class_im);
%%
write_forest_txt(apex_class_rf, [apex_class_rf.tree_dir 'txt']);
write_forest_txt(apex_offset_x_rf, [apex_offset_x_rf.tree_dir 'txt']);
write_forest_txt(apex_offset_y_rf, [apex_offset_y_rf.tree_dir 'txt']);
%%
extract_apex_map_maxima(... % the user's input
    'image_names',       {'enlargedapex0140_vessel'},...
    'data_dir',             [nailfoldroot 'data/rsa_study/cxx2/'],...
    'fov_mask_dir',         '',...
    'centre_dir',           'vessel_centres',...
    'apex_map_dir',         'apex_maps/set12g_half_296655',...
    'candidates_dir',       'apex_maps/set12g_half_296655/version1',...
    'exclusion_zone',       20,...
    'transform_scores',     0,...
    'separate_trees',  0,...
    'apex_class_thresh',    0.5,...
    'discard_pts',          0,...
    'base_width',           20,...
    'use_island_max',       1,...
    'island_thresh',        0,...
    'overwrite',            1);
%%
run_ncm_batch('rsa_study/cxx2', 2-7, ...
    'im_names',             {'enlargedapex0140_vessel'},...
    'im_idx',               [],...
    'prob_dir',             'rf_classification/296655/',... 296655
    'ori_dir',              'rf_regression/296621/',...
    'width_dir',            'rf_regression/297037/',...
    'fov_mask_dir',         '',...
    'centre_dir',           'vessel_centres/',...
    'apex_map_dir',         'apex_maps/set12g_half_296655',...
    'candidates_dir',       'version0',...
    'hog_dir',              'hogs',...
    'rescore_dir',          'rescores',...
    'displacement_dir',     'displacements',...
    'selected_dir',         'selected_apexes',...
    'metrics_dir',          'apex_metrics',...
    'aam_dir',              'aam',...
    'feature_im_dir',       'images',...
    'save_dir',             'results',....
    'model_dir',            [nailfoldroot 'data/rsa_study/models/apex_templates/'],...
    'apex_class_rf',        [nailfoldroot 'models/apex/rescoring/miccai_all/rf.mat'],...        
    'class_map',            [nailfoldroot,'models/apex/final_MAP/miccai_class_MAP/class_map.mat'],...
    'prob_sigma',           1,...
    'ori_sigma',            0,...
    'width_sigma',          1,...
    'exclusion_zone',       20,...
    'transform_scores',     0,...
    'apex_class_thresh',    0.5,...
    'discard_pts',          0,...
    'base_width',           20,...
    'use_island_max',       0,...
    'island_thresh',        0,...
    'feature_sigma',        0,...
    'num_cells',            8,...
    'cell_sz',              8,... %Size of HoG cells in blocks
    'block_sz',             [2 2],...%Size of blocks in cells
    'num_ori_bins',         12,... %Number of bins in orientation histograms
    'norm_method',          'l1-sqrt',... %Method for local normalisation
    'block_spacing',        8,... %Separation between blocks in pixels - controls level of overlap
    'gradient_operator',    [-1 0 1],...
    'spatial_sigma',        0, ...
    'angle_wrap',           1,...
    'dist_thresh',          24^2,...
    'min_candidates',       3,...
    'initial_thresh',       0.3,...
    'grid_spacing',         8,...
    'strong_vessel_thresh', 0.8,...
    'angle_discard_thresh', 75*pi/180,...
    'do_final_cull',        1,...
    'do_post_merge',        1,...
    'merge_dist_thresh',    60,...
    'merge_connect_thresh', 0.25,...
    'merge_n_pts',          20,...
    'overwrite',            0,...
    'do_aam',               0,...
    'width_predictor',      [],...[nailfoldroot 'models/apex/width/rf.mat'],...
    'do_distal',            1,...
    'do_nondistal',         1,...
    'selected_features',    [],... %Do them all for now
    'plot_distal',          1,...
    'plot_nondistal',       1, ...
    'um_per_pix',           1.25,...
    'aam_thresh',           -2e4,...
    'plot_rejected',        1,...
    'plot_r',               1,...
    'plot_c',               1);
