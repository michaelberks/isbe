
args.num_jobs = 450;
args.model_id = '20131029T134052';
%args.model_id = 'frog';
args.model_root = [nailfoldroot 'models/apex'];
args.model_name = 'rf';
args.data_dir= [nailfoldroot 'data/rsa_study/set12g/'];
args.apex_class_thresh = 0;
%     'feature_im_dir',       'predictions/detection/rf_classification/257273/',...
%     'centre_dir',           'vessel_centres/',...
%     'hog_dir',              'vessel_hogs/',...
%     'max_size',             unixenv('MAX_SIZE', 1000),...
%     'smoothing_sigma',      2,...
%     'num_cells',            8,...
%     'cell_sz',              8,... %Size of HoG cells in blocks
%     'block_sz',             [2 2],...%Size of blocks in cells
%     'num_ori_bins',         9,... %Number of bins in orientation histograms
%     'norm_method',          'l1-sqrt',... %Method for local normalisation
%     'block_spacing',        8,... %Separation between blocks in pixels - controls level of overlap
%     'gradient_operator',    [-1 0 1],...
%     'spatial_sigma',        0, ...
%     'angle_wrap',           1,...
%     'base_width',           20);

for ii = 1:10
    args.task_id = ii;
    predict_apex_offsets_set(args);
end
%%
im_list = dir(['C:\isbe\nailfold\data\rsa_study\set12g\apex_maps\' args.model_id '\*.mat']);
for ii = 1:10
    load(['C:\isbe\nailfold\data\rsa_study\set12g\apex_maps\frog\' im_list(ii).name]);
    figure; imgray(apex_offset_map);
end
%%
args.num_jobs = 602;
args.model_id = '20131029T134052';
%args.model_id = 'frog';
args.model_root = [nailfoldroot 'models/apex'];
args.model_name = 'rf';
args.data_dir= [nailfoldroot 'data/rsa_study/test/'];
args.task_id = 3;
profile on; predict_apex_offsets_set(args); profile viewer;
profsave(profile('info'),'profile_results4')

%%

%Load in the classification and offset models
rf_class_dir = [args.model_root '/classification/' args.model_id '/'];
rf_offset_x_dir = [args.model_root '/offset_x/' args.model_id '/'];
rf_offset_y_dir = [args.model_root '/offset_y/' args.model_id '/'];
apex_class_rf = u_load([rf_class_dir args.model_name '.mat']);
apex_offset_x_rf = u_load([rf_offset_y_dir args.model_name '.mat']);
apex_offset_y_rf = u_load([rf_offset_x_dir args.model_name '.mat']);

vessel_hog = u_load('C:\isbe\nailfold\data\rsa_study\set12g\vessel_hogs\enlargedapex0004_vessel_hog.mat');
%%
tic;
apex_offset_x_pred = random_forest_reg_predict(apex_offset_x_rf, vessel_hog);
toc;
num_pts = size(vessel_hog,1);
tic;
for i_pt = 1:num_pts
    apex_offset_x_pred_i = random_forest_reg_predict(apex_offset_x_rf, vessel_hog(i_pt,:));
end
toc;
%%
part_num = 1;
hog_dir = 'C:\isbe\nailfold\data\rsa_study\set12g\vessel_hogs\';
im_list = dir('C:\isbe\nailfold\data\rsa_study\set12g\vessel_centres\*.mat');
for i_hog = 1:length(im_list)
    im_name = im_list(i_hog).name(1:end-7);
    movefile(...
        [hog_dir im_name '_hog/hog_part_' zerostr(part_num,4) '.mat'],...
        [hog_dir im_name '_hog_part_' zerostr(part_num,4) '.mat']);
    rmdir([hog_dir im_name '_hog']);
            
end
%%
DATA_ROOT="scratch/nailfold/" MODEL_ROOT="models/apex" MODEL_PATH="frog" NUM_JOBS=602 IMAGE_ROOT="data/rsa_study/test" IMADE_DIR="predictions/detection/rf_classification/257273" OVERWRITE=0 qsub -l short -V -t 3 matlab_code/trunk/hydra/cuc/predict_apex_offsets_set.sh
%%
im_list = dir('C:\isbe\nailfold\data\rsa_study\test\images\*.mat');

apex_map_dir = 'C:\isbe\nailfold\data\rsa_study\test\apex_maps\frog\';
apex_cluster_dir = 'C:\isbe\nailfold\data\rsa_study\test\apex_clusters\';
apex_gt_dir = 'C:\isbe\nailfold\data\rsa_study\test\apex_gt\';
vessel_centres_dir = 'C:\isbe\nailfold\data\rsa_study\test\vessel_centres\';
fov_mask_dir = 'C:\isbe\nailfold\data\rsa_study\test\fov_masks\';

apex_class_thresh = 0.5;
base_width = 20;
for i_im = 1:10
    im_name = im_list(i_im).name(1:6);
    
    load([apex_map_dir im_name '_pred.mat']);
    load([apex_gt_dir im_name '_gt.mat']);
    %load([apex_cluster_dir im_name '_apex_clusters.mat']);
    load([vessel_centres_dir im_name '_vc.mat']);
    f_mask = u_load([fov_mask_dir im_name '_f_mask.mat']);
    
    [discard_pts] = discard_edge_preds(vessel_centre, f_mask);
    include_pts = ~discard_pts & (apex_class_pred > apex_class_thresh);
    
    [apex_offset_map] = ...
        transform_apex_offset_preds(apex_class_pred, apex_offset_x_pred, apex_offset_y_pred,...
            vessel_centre, nrows, ncols, base_width, include_pts);
        
    [candidates_xy candidate_scores] = ...
        local_image_maxima(apex_offset_map, 40, f_mask, 0);
    
    fig_h = figure; imgray(apex_offset_map);
    %plot_apex_clusters(vessels, fig_h);
    plot(candidates_xy(:,1), candidates_xy(:,2), 'rx');
    plot(apex_xy(is_distal,1), apex_xy(is_distal,2), 'go');
end
%%
image_dir = 'C:\isbe\nailfold\data\rsa_study\test\images\';
ori_dir = 'C:\isbe\nailfold\data\rsa_study\test\predictions\orientation\rf_regression\259076\';
width_dir = 'C:\isbe\nailfold\data\rsa_study\test\predictions\width\rf_regression\257847\';
candidates_dir = 'C:\isbe\nailfold\data\rsa_study\test\apex_maps\frog\local_maxima\';
model_dir = 'C:\isbe\nailfold\data\rsa_study\models\apex_templates\';
    
im_num = '15373c';

load([model_dir 'mean_shape.mat'], 'mean_shape');
load([candidates_dir im_num '_candidates.mat'], 'candidate_xy'); 
nailfold = u_load([image_dir im_num '.mat']);
vessel_ori = u_load([ori_dir im_num '_pred.mat']);
vessel_width = u_load([width_dir im_num '_pred.mat']);

g = gaussian_filters_1d(2);
g = g / sum(g);
vessel_width = conv2(g', g, vessel_width, 'same');

initialise_aam_candidates([image_dir im_num '.mat'], vessel_width, vessel_ori, candidate_xy, 'debug', 1, 'mean_shape', mean_shape);
%%
figure; imgray(complex2rgb(ori_patch));
plot(vessel_xy(:,1), vessel_xy(:,2));

[nxy] = compute_spline_normals(vessel_xy);

plot([vessel_xy(:,1)+10*nxy(:,1) vessel_xy(:,1)-10*nxy(:,1)]', [vessel_xy(:,2)+10*nxy(:,2) vessel_xy(:,2)-10*nxy(:,2)]');
%%
clear
load('temp.mat', 'width_patch', 'ori_patch', 'image_patch', 'mean_shape', 'vessel_xy');
cx = vessel_xy(16,1);
cy = vessel_xy(16,2);
vessel_xy_c = [vessel_xy(:,1)-cx vessel_xy(:,2)-cy];

[nxy] = compute_spline_normals(vessel_xy);
vessel_xy_ori = exp(-2i*atan(nxy(:,1)./nxy(:,2)));

for i_theta = linspace(-pi/4, pi/4, 20)
    figure; imgray(image_patch);
    rot = [cos(i_theta) -sin(i_theta); sin(i_theta) cos(i_theta)];  
    vessel_xy_i = bsxfun(@plus, vessel_xy_c*rot*1.25, [cx cy]);

    vessel_xy_ori_i = vessel_xy_ori * exp(-2i*i_theta);

    patch_ori_i = interp2(ori_patch, vessel_xy_i(:,1), vessel_xy_i(:,2));
    ori_diff = abs(angle(patch_ori_i .* vessel_xy_ori_i)/2);
    
    fitting_error = sum(abs(patch_ori_i));% - sum(ori_diff);
    %fitting_error = -sum(ori_diff);
    
    plot(vessel_xy_i(:,1), vessel_xy_i(:,2));
    
    title(num2str(sum(fitting_error)));
    
end
%%
clear
load('temp.mat', 'width_patch', 'ori_patch', 'image_patch', 'mean_shape', 'vessel_xy');

cx = vessel_xy(16,1);
cy = vessel_xy(16,2);
vessel_xy_c = [vessel_xy(:,1)-cx vessel_xy(:,2)-cy];

[nxy] = compute_spline_normals(vessel_xy);
vessel_xy_ori = exp(-2i*atan(nxy(:,1)./nxy(:,2)));

num_thetas = 10;
num_scales = 5;
thetas = linspace(-pi/8, pi/8, num_thetas);
scales = linspace(0.75, 1.25, num_scales);

best_fit = -inf;
for i_theta = 1:num_thetas
    theta = thetas(i_theta);
    rot = [cos(theta) -sin(theta); sin(theta) cos(theta)];  
    vessel_xy_ori_i = vessel_xy_ori * exp(-2i*theta);
    
    for i_scale = 1:num_scales
        scale = scales(i_scale);
        
        vessel_xy_i = bsxfun(@plus, vessel_xy_c*rot*scale, [cx cy]);

        patch_ori_i = interp2(ori_patch, vessel_xy_i(:,1), vessel_xy_i(:,2));
        ori_diff = abs(angle(patch_ori_i .* vessel_xy_ori_i)/2);
        
        fit_score = sum(abs(patch_ori_i)) - sum(ori_diff);
        %fit_score = sum(abs(patch_ori_i));
        %fit_score = -sum(ori_diff);
        
        if fit_score > best_fit
            best_fit = fit_score;
            best_xy = vessel_xy_i;
            best_theta = theta;
            best_scale = scale;
        end
    end
end
figure; imgray(image_patch);
plot(best_xy(:,1), best_xy(:,2));  
title(['Best \theta = ' num2str(round(180*best_theta/pi)) ', best scale =' num2str(best_scale), ', score =' num2str(best_fit)]);
%%
