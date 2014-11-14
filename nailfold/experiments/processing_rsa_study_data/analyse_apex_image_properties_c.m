function [v_stats] = analyse_apex_image_properties_c(start_idx, end_idx)
%ANALYSE_APEX_IMAGE_PROPERTIES *Insert a one line summary here*
%   [] = analyse_apex_image_properties()
%
% Inputs:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 18-Jun-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
if isdeployed    
    task_idx = unixenv('SGE_TASK_ID', 1);
    start_idx = 20*(task_idx-1) + 1;
    end_idx = 20*task_idx;
else
    if ~exist('start_idx', 'var'); start_idx = 1; end
    if ~exist('end_idx', 'var'); end_idx = 10; end
end

%Args (hardcoded for now)
apex_radius = 5;
max_apex_guess = 500;
vessel_prob_smoothing_sigma = 2;
curvature_smoothing_sigma = 2;
strong_vessel_thresh = 0.25;
curv_max = 0.5;

%Set up image lists
if ispc
    rsa_dir = 'rsa_study/';
else
    rsa_dir = [];
end

model_dir = [nailfoldroot 'data/' rsa_dir 'models/apex_templates/'];
image_dir = [nailfoldroot 'data/' rsa_dir 'test/images/'];
fov_mask_dir = [nailfoldroot 'data/' rsa_dir 'test/fov_masks/'];
prob_dir = [nailfoldroot 'data/' rsa_dir 'test/predictions/detection/rf_classification/222836/'];
ori_dir = [nailfoldroot 'data/' rsa_dir 'test/predictions/orientation/rf_regression/222835/'];
%width_dir = [nailfoldroot 'data/' rsa_dir 'test/predictions/width/rf_regression/182367/'];
markup_dir = [nailfoldroot 'data/' rsa_dir 'markup/tmoore/'];
centre_dir = [nailfoldroot 'data/' rsa_dir 'test/vessel_centres/'];
corr_scores_dir = [nailfoldroot 'data/' rsa_dir 'test/corr_scores/'];
shape_scores_dir = [nailfoldroot 'data/' rsa_dir 'test/shape_scores/'];

create_folder(centre_dir);
create_folder(corr_scores_dir);
create_folder(shape_scores_dir);

pred_list = dir([prob_dir '*.mat']);

num_images = length(pred_list);

%Pre-allocate outputs
v_stats.vessel_prob.apex_mean = zeros(max_apex_guess,1);
v_stats.vessel_prob.apex_max = zeros(max_apex_guess,1);
v_stats.vessel_prob.apex_val = zeros(max_apex_guess,1);
v_stats.vessel_prob.background = zeros(num_images,100);
%
v_stats.vessel_curv.apex_mean = zeros(max_apex_guess,1);
v_stats.vessel_curv.apex_max = zeros(max_apex_guess,1);
v_stats.vessel_curv.apex_val = zeros(max_apex_guess,1);
v_stats.vessel_curv.background = zeros(num_images,100);
%
v_stats.vessel_corr.apex_mean = zeros(max_apex_guess,1);
v_stats.vessel_corr.apex_max = zeros(max_apex_guess,1);
v_stats.vessel_corr.apex_val = zeros(max_apex_guess,1);
v_stats.vessel_corr.background = zeros(num_images,100);
%
% v_stats.vessel_width.apex_mean = zeros(max_apex_guess,1);
% v_stats.vessel_width.apex_val = zeros(max_apex_guess,1);
v_stats.vessel_width.apex_truth = zeros(max_apex_guess,1);
%
v_stats.vessel_centres.apex_d2c = zeros(max_apex_guess,1);
%
v_stats.vessel_properties.shape = cell(max_apex_guess,1);
v_stats.vessel_properties.size = cell(max_apex_guess,1);
v_stats.vessel_properties.im_num = zeros(max_apex_guess,1);
v_stats.vessel_properties.apex_pos = zeros(max_apex_guess,2);
%
v_stats.vessel_prob_corr.background = zeros(num_images,100);
v_stats.vessel_prob_curv.background = zeros(num_images,100);
%
apex_count = 1;
%
vessel_prob_bins = linspace(0, 1, 100);
vessel_corr_bins = linspace(-1, 1, 100);
vessel_curv_bins = linspace(-curv_max, curv_max, 100);

%Load in model templates
load([model_dir 'apex_templates.mat']);

%Create circular mask for each apex
patch_sz = 2*apex_radius+1;
circ_x = repmat(-apex_radius:apex_radius, patch_sz, 1);
circ_y = circ_x';
circ_mask = (circ_x.^2 + circ_y.^2) < apex_radius^2;
ncirc_mask = ~circ_mask;

%Create smoothing kernel for vessel probs
g = gaussian_filters_1d(vessel_prob_smoothing_sigma);
g = g / sum(g);

%Set up priors for apex streaming
load([model_dir 'apex_curv_stats.mat']);
ori_hist_r_smoothed = imfilter(curv_stats.r.pred_ori_by_step.hist, conv(g,g), 'circular');
ori_hist_r_smoothed = imfilter(ori_hist_r_smoothed, conv(g,g)', 'replicate');

ori_hist_l_smoothed = imfilter(curv_stats.l.pred_ori_by_step.hist, conv(g,g), 'circular');
ori_hist_l_smoothed = imfilter(ori_hist_l_smoothed, conv(g,g)', 'replicate');

apex_prior = cat(3, ori_hist_l_smoothed, ori_hist_r_smoothed);
num_streams = 100;

if ~exist('start_idx', 'var'); start_idx = 1; end;
if ~exist('end_idx', 'var'); end_idx = num_images; end;

%Loop though each image
for i_im = start_idx:end_idx
    im_name = pred_list(i_im).name(1:end-9);
    
    %if exist([centre_dir im_name '_vc.mat'], 'file');
    %    continue;
    %end
    
    display(['Processing image ' num2str(i_im) ', ' datestr(now)]);
    
    
    
    %Load in images and vessel markup
    try
        vessel_prob = u_load([prob_dir pred_list(i_im).name]);
        vessel_ori = u_load([ori_dir pred_list(i_im).name]);
%         vessel_width = u_load([width_dir pred_list(i_im).name]);
        vessel_im = u_load([image_dir im_name '.mat']);
        fov_mask = u_load([fov_mask_dir im_name '_f_mask.mat']);
    catch last_err
        display(last_err.message);
        continue;
    end
    
    %Smooth the vessel probs
    vessel_prob = conv2(g', g, vessel_prob, 'same');
    vessel_ori = conv2(g', g, vessel_ori, 'same');
    
    %Equalize image intensities
    vessel_im = equalise_nailfold_intensity(vessel_im);
    
    %Compute curvature
    vessel_curv = -(complex_gaussian_curvature(vessel_ori ./ (abs(vessel_ori)+1e-6), curvature_smoothing_sigma));
    
    %Compute correlation with g2d template
    [~, ~, vessel_corr] = template_match_apexes(vessel_im, ...
        {'g2d', apex_template.g2}, 'threshold', 1);
        
    %Compute NMS centrelines
    vessel_nms = mb_non_maximal_supp(vessel_prob, angle(vessel_ori)/2);
    strong_vessels = vessel_nms > strong_vessel_thresh;
    if any(strong_vessels(:))
        [rstrong cstrong] = find(strong_vessels);
        vessel_centre_mask = bwselect(vessel_nms > 0, cstrong, rstrong, 8);
    else
        vessel_centre_mask = strong_vessels;
    end
    [vessel_centre_y vessel_centre_x] = find(vessel_centre_mask);
 
    %Now see if we have markup we can process
    vessel_markup_list = dir([markup_dir '*' im_name '*.txt']);    
    
    if ~isempty(vessel_markup_list)
        
        %Read in markup
        vessel_markup = read_markup_from([markup_dir vessel_markup_list(end).name]);
        
        %Initialise the background mask
        background_mask = vessel_centre_mask;
        background_mask = padarray(background_mask, [apex_radius apex_radius]);

        %Loop through each marked vessel apex
        num_vessels = length(vessel_markup.vessels);
        for i_v = 1:num_vessels

            %Check this is a valid vessel
            anchor_xy = vessel_markup.vessels(i_v).anchor;

            if isempty(anchor_xy); continue; end

            %Check if distal
            is_distal = vessel_markup.vessels(i_v).ncm_vessel_properties.is_distal;

            if is_distal
                num_apices = length(vessel_markup.vessels(i_v).apices);
                for i_a = 1:num_apices
                    if isempty(vessel_markup.vessels(i_v).apices(i_a).inner_point)                    
                        %Use the anchor
                        axy = round(anchor_xy);

                    else
                        %Compute the centre of the apex
                        apex_xy =  ...
                            [ vessel_markup.vessels(i_v).apices(i_a).outer_point;...
                              vessel_markup.vessels(i_v).apices(i_a).inner_point];
                        apex_width = sqrt(sum(diff(apex_xy).^2));

                        v_stats.vessel_properties.apex_pos(apex_count,:) = ...
                            mean(apex_xy);

                        apex_xy = mean(apex_xy);

                        %Find the nearest centreline pixel to the apex
                        c_dists = sqrt((vessel_centre_x-apex_xy(1)).^2 + (vessel_centre_y-apex_xy(2)).^2);
                        [min_dist, min_i] = min(c_dists);
                        axy = [vessel_centre_x(min_i) vessel_centre_y(min_i)];


                        %Sample patches from curvature, correlation, vessel
                        %prob etc
                        curv_patch = sample_window(vessel_curv, patch_sz, axy(2), axy(1));
                        corr_patch = sample_window(vessel_corr, patch_sz, axy(2), axy(1));
                        prob_patch = sample_window(vessel_prob, patch_sz, axy(2), axy(1));
    %                     width_patch = sample_window(vessel_width, patch_sz, axy(2), axy(1));
                        centre_patch = sample_window(vessel_centre_mask, patch_sz, axy(2), axy(1));

                        %Compute outputs for apex
                        v_stats.vessel_prob.apex_mean(apex_count,:) = ...
                            mean(prob_patch(centre_patch & circ_mask));
                        v_stats.vessel_prob.apex_max(apex_count,:) = ...
                            max(prob_patch(centre_patch & circ_mask));
                        v_stats.vessel_prob.apex_val(apex_count,:) = ...
                            vessel_prob(axy(2), axy(1));

                        %
                        v_stats.vessel_curv.apex_mean(apex_count,:) = ...
                            mean(curv_patch(centre_patch & circ_mask));
                        v_stats.vessel_curv.apex_max(apex_count,:) = ...
                            max(curv_patch(centre_patch & circ_mask));
                        v_stats.vessel_curv.apex_val(apex_count,:) = ...
                            vessel_curv(axy(2), axy(1));

                        %
                        v_stats.vessel_corr.apex_mean(apex_count,:) = ...
                            mean(corr_patch(centre_patch & circ_mask));
                        v_stats.vessel_corr.apex_max(apex_count,:) = ...
                            max(corr_patch(centre_patch & circ_mask));
                        v_stats.vessel_corr.apex_val(apex_count,:) = ...
                            vessel_corr(axy(2), axy(1));

                        %
    %                     v_stats.vessel_width.apex_mean(apex_count,:) = ...
    %                         mean(width_patch(centre_patch & circ_mask));
                        v_stats.vessel_width.apex_truth(apex_count,:) = ...
                            apex_width;
    %                     v_stats.vessel_width.apex_val(apex_count,:) = ...
    %                         vessel_width(axy(2), axy(1));

                        %
                        v_stats.vessel_centres.apex_d2c(apex_count,:) = ...
                            min_dist;

                        %
                        v_stats.vessel_properties.shape{apex_count,:} = ...
                            vessel_markup.vessels(i_v).ncm_vessel_properties.shape;
                        v_stats.vessel_properties.size{apex_count,:} = ...
                            vessel_markup.vessels(i_v).ncm_vessel_properties.size;
                        v_stats.vessel_properties.im_num(apex_count,:) = ...
                            i_im;                                       

                        %Increment the vessel counts
                        apex_count = apex_count + 1;

                    end
                    %Mask out the region around the apex
                    axy = axy + apex_radius;
                    background_mask(axy(2)+circ_y(:,1), axy(1)+circ_x(1,:)) = ...
                        background_mask(axy(2)+circ_y(:,1), axy(1)+circ_x(1,:)) & ncirc_mask;
                end
            else
                %Mask out the region around the anchor
                axy = round(anchor_xy) + apex_radius;
                background_mask(axy(2)+circ_y(:,1), axy(1)+circ_x(1,:)) = ...
                    background_mask(axy(2)+circ_y(:,1), axy(1)+circ_x(1,:)) & ncirc_mask;
            end
        end

        %remove backgorund padding
        background_mask([1:apex_radius end-apex_radius+1:end],:) = [];
        background_mask(:,[1:apex_radius end-apex_radius+1:end]) = [];

        %Compute outputs for background
        v_stats.vessel_prob.background(i_im,:) = ...
            hist(vessel_prob(background_mask), vessel_prob_bins);
        v_stats.vessel_corr.background(i_im,:) = ...
            hist(vessel_corr(background_mask), vessel_corr_bins);
        v_stats.vessel_curv.background(i_im,:) = ...
            hist(vessel_curv(background_mask), vessel_curv_bins);

        %Combination of vessel prob and curv/corr
        v_stats.vessel_prob_corr.background(i_im,:) = ...
            hist(vessel_corr(background_mask).* vessel_prob(background_mask), vessel_corr_bins);
        v_stats.vessel_prob_curv.background(i_im,:) = ...
            hist(vessel_curv(background_mask).* vessel_prob(background_mask), vessel_curv_bins);
    end
    
    %----------------------------------------------------------------------
    % through the centre pts and mark which would be discarded depending at
    % each stage of the analysis workflow
    num_pts = length(vessel_centre_x);    
    vessel_centre_discards = zeros(num_pts, 1);
    remaining_pts = true(num_pts,1);
    apex_corr_max_points = false(num_pts,1);
    
    [nrows ncols] = size(vessel_prob); %#ok
    vessel_centre_prob = vessel_prob(vessel_centre_mask);
    vessel_centre_ori = vessel_ori(vessel_centre_mask);
    vessel_centre_curv = vessel_curv(vessel_centre_mask);
    vessel_centre_corr = vessel_corr(vessel_centre_mask);
    
    centre_labels = bwlabel(vessel_centre_mask, 8);
    centre_labels = centre_labels(vessel_centre_mask);
    
    % #1, #2, #3
    %Discard centre lines near the edge of the image, and any long straight
    %horizontal lines
    border_mask = fov_mask & ~imerode(fov_mask, strel('disk', 32));
    border_centre = border_mask(vessel_centre_mask);
    [~, border_ori] = gaussian_1st_derivative_gradient(fov_mask, 16);   
    border_ori = exp(2i*border_ori(vessel_centre_mask));
    
    for i_lab = 1:max(centre_labels(:))
        label_mask = centre_labels == i_lab;
        
        label_ori = vessel_centre_ori(label_mask);        
        label_border = border_centre(label_mask);
        label_border_ori = border_ori(label_mask);
        label_corr = vessel_centre_corr(label_mask);
        mean_ori = mean(label_ori);
        
        %Pixels near the border, with matching orientation - probably picke
        %dup edge sof the image frames
        if any(label_border) && ...
            mean(abs(ori_error(label_border_ori(label_border), label_ori(label_border)))) < (pi*10/180)
            remaining_pts(label_mask) = 0;
            vessel_centre_discards(label_mask) = 1;
        
        %Straight horizontal lines through the image
        %elseif sum(label_mask(:)) > 20 && ...
        %        (abs(angle(mean_ori)/2) < pi/6) && (abs(mean_ori) > 0.25) 
        %            remaining_pts(label_mask) = 0;  
        %            vessel_centre_discards(label_mask) = 2;
        
        %Lines connected to pixels with strong apex correlation - the pixel
        %with max correlation is kept and assumed to be an apex, the rest
        %discarded
        elseif any(label_corr > 0.4)

            label_x = vessel_centre_x(label_mask);
            label_y = vessel_centre_x(label_mask);
            
            label_apex_corr_max_points = false(size(label_corr));
            [~,max_idx] = max(label_corr);
            
            discard_mask = ((label_x - label_x(max_idx)).^2 + ...
                (label_y - label_y(max_idx)).^2) < 20^2;
            discard_mask(max_idx) = 1;
            
            label_apex_corr_max_points(max_idx) = 1;
            apex_corr_max_points(label_mask) = label_apex_corr_max_points;  
            remaining_pts(discard_mask) = 0;
            vessel_centre_discards(discard_mask) = 3;
        end
        
    end

    %Take a copy of the pixels with max apex correlation to use in the next
    %stage of analysis
    max_corr_apex_x = vessel_centre_x(apex_corr_max_points);
    max_corr_apex_y = vessel_centre_y(apex_corr_max_points);
    max_corr_vals = vessel_centre_corr(apex_corr_max_points); 
    
    [maxima_corr_xy, maxima_corr_vals] = ...
        apply_local_exclusion([max_corr_apex_x max_corr_apex_y],...
        max_corr_vals, 50);
    
    mid_y = median(maxima_corr_xy(:,2));
    
    % #4
    %Discard pixels, more than a set distance away from the median height of the
    %selected apex pixels
    out_of_bounds_pts = (vessel_centre_y < mid_y-250) | ...
        (vessel_centre_y > mid_y+200);
    
    vessel_centre_discards(remaining_pts & out_of_bounds_pts) = 4;
    remaining_pts = remaining_pts & ~out_of_bounds_pts;
    
    % #5
    %Discard pixels with low probability
    low_prob_pts = vessel_centre_prob < 0.4;
    vessel_centre_discards(remaining_pts & low_prob_pts) = 5;
    remaining_pts = remaining_pts & ~low_prob_pts;
    
    % #6
    %Discard pixels with low curvature
    %low_curv_pts = abs(vessel_centre_curv) < 0.01;
    %vessel_centre_discards(remaining_pts & low_curv_pts) = 6;
    %remaining_pts = remaining_pts & ~low_curv_pts;

    
    %Take a copy of the remaining pts and apply probabilistic stream to
    %compute shape prior scores
    apex_x = vessel_centre_x(remaining_pts);
    apex_y = vessel_centre_y(remaining_pts);
    
    [apex_prior_scores] = apex_prior_prob_track_mult(...
                    vessel_prob, vessel_ori,...
                    apex_y, apex_x,...
                    'num_streams', num_streams, 'stopping_prob', 0.1,...
                    'apex_prior', apex_prior);
    apex_prior_scores_sum = squeeze(sum(sum(apex_prior_scores,2)));
   
    %Compute the highest shape prior for each connected region
    shape_prior_mask = vessel_centre_mask;
    shape_prior_mask(vessel_centre_mask) = remaining_pts;

    centre_labels = bwlabel(shape_prior_mask, 8);
    apex_shape_max_points = false(num_pts,1);
    
    num_labels = max(centre_labels(:));
    max_shape_vals = zeros(num_labels,1);
    max_shape_apex_x = zeros(num_labels,1);
    max_shape_apex_y = zeros(num_labels,1);

    for i_lab = 1:num_labels
        label_mask = centre_labels == i_lab;      
        label_scores = apex_prior_scores_sum(label_mask(shape_prior_mask));
        [label_y label_x] = find(label_mask);

        label_apex_shape_max_points = false(size(label_scores));
        [max_shape_vals(i_lab), max_idx] = max(label_scores);
        max_shape_apex_x(i_lab) = label_x(max_idx);
        max_shape_apex_y(i_lab) = label_y(max_idx);
        label_apex_shape_max_points(max_idx) = 1;
        apex_shape_max_points(label_mask(vessel_centre_mask)) = label_apex_shape_max_points;  

        remaining_pts(label_mask(vessel_centre_mask)) = 0;
        vessel_centre_discards(label_mask(vessel_centre_mask)) = 7;

    end
    
    [maxima_shape_xy, maxima_shape_vals] = ...
        apply_local_exclusion([max_shape_apex_x max_shape_apex_y],...
        max_shape_vals, 50);
       
    save([centre_dir im_name '_vc.mat'], 'vessel_centre_*', 'nrows', 'ncols');
    save([corr_scores_dir im_name '_cs.mat'], 'maxima_corr_xy', 'maxima_corr_vals', 'nrows', 'ncols');
    save([shape_scores_dir im_name '_ss.mat'], 'apex_x', 'apex_y', 'apex_prior_scores', 'maxima_shape_xy', 'maxima_shape_vals', 'nrows', 'ncols');
        
end

%Discard any additional pre-allocated space in the outputs
if apex_count <= max_apex_guess
    v_stats.vessel_prob.apex_mean(apex_count:end,:) = [];
    v_stats.vessel_prob.apex_max(apex_count:end,:) = [];
    v_stats.vessel_prob.apex_val(apex_count:end,:) = [];
    %
    v_stats.vessel_curv.apex_mean(apex_count:end,:) = [];
    v_stats.vessel_curv.apex_max(apex_count:end,:) = [];
    v_stats.vessel_curv.apex_val(apex_count:end,:) = [];
    %
    v_stats.vessel_corr.apex_mean(apex_count:end,:) = [];
    v_stats.vessel_corr.apex_max(apex_count:end,:) = [];
    v_stats.vessel_corr.apex_val(apex_count:end,:) = [];
    %
%     v_stats.vessel_width.apex_mean(apex_count:end,:) = [];
    v_stats.vessel_width.apex_truth(apex_count:end,:) = [];
%     v_stats.vessel_width.apex_val(apex_count:end,:) = [];
    %
    v_stats.vessel_centres.apex_d2c(apex_count:end,:) = [];
    %
    v_stats.vessel_properties.shape(apex_count:end,:) = [];
    v_stats.vessel_properties.size(apex_count:end,:) = [];
    v_stats.vessel_properties.im_num(apex_count:end,:) = [];
    v_stats.vessel_properties.apex_pos(apex_count:end,:) = [];
end
if isdeployed
    save([nailfoldroot 'data/' rsa_dir 'models/apex_image_statistics_c_' zerostr(task_idx, 2) '.mat'],...
        'v_stats');
end    
    
    
