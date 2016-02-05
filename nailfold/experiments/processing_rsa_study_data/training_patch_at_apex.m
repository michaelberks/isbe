function [all_rf_widths all_gt_widths profile_features model_scores g2d_features] =...
    training_patch_at_apex(start_i, end_i, do_plot)
%Looking at remeasuring apex width

if ~exist('start_i', 'var') || isempty(start_i)
    start_i = 1;
end
if ~exist('end_i', 'var') || isempty(end_i)
    end_i = start_i;
end

args.data_dir = [nailfoldroot 'data/rsa_study/set12g_half/'];
args.image_dir = 'images';
args.prob_dir = 'rf_classification/296655';
args.ori_dir = 'rf_regression/296621';
args.width_dir = 'rf_regression/297037';
args.contour_dir = 'vessel_contours';
args.model_dir = [nailfoldroot 'data/rsa_study/models/apex_templates/'];
args.aam_dir = 'aam/';
args.aam_name = 'aam/orig/2/vessel_apex_orig.smd';
args.prob_sigma = 2;
args.ori_sigma = 0;
args.width_sigma = 2;
args.base_width = 20;

image_dir = [args.data_dir args.image_dir '/'];

%ori_dir = [args.data_dir 'predictions/orientation/' args.ori_dir '/'];
width_dir = [args.data_dir 'predictions/width/' args.width_dir '/'];
width_gt_dir = [args.data_dir 'width_maps/'];

%contour_dir = [args.data_dir args.contour_dir '/'];
aam_dir = [args.data_dir args.aam_dir '/'];

pred_list = dir([width_dir '*.mat']);
% num_images = length(pred_list);

% %1. Workout number of images in job
% job_size = ceil(num_images / args.num_jobs);
% 
% %2. Workout start and end indices for job
% start_i	= (args.task_id-1)*job_size + 1;
% end_i	= min(args.task_id*job_size, num_images);
display(['Extracting vessel centres from images ' num2str(start_i) ' to ' num2str(end_i)]);

load([args.model_dir 'mean_shape.mat'], 'mean_shape'); 

predictor.trees = cell(100,1);
for i_tree = 1:100
    predictor.trees{i_tree} = ['rf_tree' zerostr(i_tree,4) '.mat'];
end
predictor.D = 30;
predictor.d = 5;
predictor.tree_dir = 'C:\isbe\nailfold\models\apex\width\trees\';
predictor.tree_root = [];
predictor.regression_method = 'rf_regression';
predictor.sampled_data_dir = [];

num_images = end_i - start_i + 1;
curr_im = 1;

all_rf_widths = zeros(num_images, 31);
all_gt_widths = zeros(num_images, 31);
model_scores = zeros(num_images, 1);
profile_features = zeros(31, 300, num_images);

decomposition_args.decomp_type = {'g2dia', 'h2dia'}; %Downsampling/Interpolating form of Gaussian derivatives
decomposition_args.win_size = 3;            %Window size about pixel to sample features
decomposition_args.sigma_range = [2 5];     %Finest scale Gaussian sigma and number of levels
decomposition_args.num_angles = 6;
decomposition_args.do_max = 0;              %Keep all 6 sub-bands (0) or only keep the maximum (1)
decomposition_args.rotate = 0;              %Try and make rotation invariant (best left switched off)
decomposition_args.normalise = 0;
feature_size = 3*3*6*2*5;
g2d_features = zeros(31, feature_size, num_images);

%Loop though each image
for i_im = start_i:end_i
    im_name = pred_list(i_im).name(1:end-9);
    
    display(['Processing image ' num2str(i_im) ', ' datestr(now)]);  
    
    %Load in images and vessel markup
    vessel_im = u_load([image_dir im_name '.mat']);
%     vessel_ori = u_load([ori_dir pred_list(i_im).name]);
    vessel_width = u_load([width_dir pred_list(i_im).name]);
    vessel_width_gt = u_load([width_gt_dir im_name '_width.mat']);
%     s = load([contour_dir im_name '_vc.mat']);   
%     
%     %Smooth the vessel probs
%     if args.ori_sigma
%         g_ori = gaussian_filters_1d(args.ori_sigma);
%         g_ori = g_ori / sum(g_ori);
%         vessel_ori = conv2(g_ori', g_ori, vessel_ori, 'same');
%     end 
%     if args.width_sigma
%         g_width = gaussian_filters_1d(args.width_sigma);
%         g_width = g_width / sum(g_width);
%         vessel_width_sm = conv2(g_width', g_width, vessel_width, 'same');
%     end   

%     candidate_xy = s.vessel_centre(s.apex_idx(1),:);
%     initialise_aam_candidates([image_dir im_name '.mat'], vessel_width_sm, vessel_ori, candidate_xy, ...
%             'aam_dir', [aam_dir im_name '/'],...
%             'mean_shape', mean_shape,...
%             'base_width', 15,...
%             'width_sigma', 2,...
%             'ori_sigma', 0,...
%             'max_num_candidates', inf,...
%             'thetas', linspace(-pi/8, pi/8, 20),...
%             'scales', [0.8 0.9 1.0 1.1 1.2],...
%             'debug', 0);
%         
%     fit_aam_to_candidates(...
%         'aam_dir', [aam_dir im_name '/'],...
%         'aam_exe', 'ncm_sandpit_mb',...
%         'aam_path', [args.model_dir args.aam_name],...
%         'delete_candidate_patches', 0);
    
    load([aam_dir im_name '\apex' zerostr(1, 4) '_aam.mat'], 'apex_candidate');
    
    %apex_x = apex_candidate.vessel_xy(:,1) + apex_candidate.sc;
    %apex_y = apex_candidate.vessel_xy(:,2) + apex_candidate.sr;
    fapex_x = apex_candidate.fitted_vessel_xy(:,1) + apex_candidate.sc;
    fapex_y = apex_candidate.fitted_vessel_xy(:,2) + apex_candidate.sr;
    model_scores(i_im) = apex_candidate.model_score;
    
    fapex_i = sub2ind(size(vessel_width), round(fapex_y), round(fapex_x));
    apex_widths = vessel_width(fapex_i);
    ap_gt_widths = vessel_width_gt(fapex_i);
    
    %apex_widths = interp2(vessel_width, fapex_x, fapex_y, '*linear');
%     [predicted_width_ratios profile_features(:,:,i_im)] = predict_apex_width_image(...
%        vessel_im, [fapex_x, fapex_y], apex_widths, predictor, 'num_profile_pts', 50, 'scale_offsets', [0.5 1 2]);
   
    [responses] = compute_filter_responses(vessel_im, decomposition_args);
    g2d_features(:,:,i_im) = sample_image_features(responses, round(fapex_y), round(fapex_x), decomposition_args);

%     v_gt_widths = sqrt(sum((s.outer_edge - s.inner_edge).^2,2));
%     ap_gt_widths = zeros(31,1);
    
    if do_plot
        figure;
        imgray(vessel_im);
        %plot(apex_x, apex_y, 'r', 'linewidth', 2);
        plot(fapex_x, fapex_y, 'g--');
        title(['Model score: ' num2str(apex_candidate.model_score,3)]);
        
%         normal_xy = compute_spline_normals([fapex_x, fapex_y]);
    end    
    
    for i_pt = 1:31
%         if do_plot
%             pt_width = apex_widths(i_pt) * predicted_width_ratios(i_pt);
%             nx = fapex_x(i_pt) + [-.5 .5]*normal_xy(i_pt,1)*pt_width;
%             ny = fapex_y(i_pt) + [-.5 .5]*normal_xy(i_pt,2)*pt_width;
%             plot(nx, ny, '.');
%         end
%         d = sum(bsxfun(@minus, s.vessel_centre, [fapex_x(i_pt) fapex_y(i_pt)]).^2,2);
%         [~,min_i] = min(d);
%         ap_gt_widths(i_pt) = v_gt_widths(min_i);
    end
    
    all_rf_widths(curr_im,:) = apex_widths';
    all_gt_widths(curr_im,:) = ap_gt_widths';
    
    curr_im = curr_im + 1;
    
end