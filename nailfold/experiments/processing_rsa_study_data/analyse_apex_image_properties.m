function [image_stats] = analyse_apex_image_properties(start_idx, end_idx)
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

%Args (hardcoded for now)
apex_radius = 5;
max_apex_guess = 8000;
vessel_prob_smoothing_sigma = 2;
curvature_smoothing_sigma = 2;
strong_vessel_thresh = 0.25;
curv_max = 0.2;
d2c_max = 100;

%Set up image lists
if ispc
    rsa_dir = 'rsa_study/';
else
    rsa_dir = [];
end

model_dir = [nailfoldroot 'data/' rsa_dir 'models/apex_templates/'];
image_dir = [nailfoldroot 'data/' rsa_dir 'test/images/'];
prob_dir = [nailfoldroot 'data/' rsa_dir 'test/predictions/detection/rf_classification/222836/'];
ori_dir = [nailfoldroot 'data/' rsa_dir 'test/predictions/orientation/rf_regression/222835/'];
%width_dir = [nailfoldroot 'data/' rsa_dir 'test/predictions/width/rf_regression/182367/'];
markup_dir = [nailfoldroot 'data/' rsa_dir 'markup/tmoore/'];

pred_list = dir([prob_dir '*.mat']);

num_images = length(pred_list);

%Pre-allocate outputs
image_stats.vessel_prob.apex_mean = zeros(max_apex_guess,1);
image_stats.vessel_prob.apex_max = zeros(max_apex_guess,1);
image_stats.vessel_prob.background = zeros(num_images,100);
%
image_stats.vessel_curv.apex_mean = zeros(max_apex_guess,1);
image_stats.vessel_curv.apex_max = zeros(max_apex_guess,1);
image_stats.vessel_curv.background = zeros(num_images,100);
%
image_stats.vessel_corr.apex_mean = zeros(max_apex_guess,1);
image_stats.vessel_corr.apex_max = zeros(max_apex_guess,1);
image_stats.vessel_corr.background = zeros(num_images,100);
%
% image_stats.vessel_width.apex_mean = zeros(max_apex_guess,1);
% image_stats.vessel_width.apex_truth = zeros(max_apex_guess,1);
%
image_stats.vessel_centres.apex_d2c = zeros(max_apex_guess,1);
image_stats.vessel_centres.background = zeros(num_images,100);
%
image_stats.vessel_properties.shape = cell(max_apex_guess,1);
image_stats.vessel_properties.size = cell(max_apex_guess,1);
image_stats.vessel_properties.im_num = zeros(max_apex_guess,1);
%
image_stats.vessel_prob_corr.background = zeros(num_images,100);
image_stats.vessel_prob_curv.background = zeros(num_images,100);
%
apex_count = 1;
%
vessel_prob_bins = linspace(0, 1, 100);
vessel_corr_bins = linspace(-1, 1, 100);
vessel_curv_bins = linspace(0, curv_max, 100);
vessel_d2c_bins = linspace(0, d2c_max, 100);

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

if ~exist('start_idx', 'var'); start_idx = 1; end;
if ~exist('end_idx', 'var'); end_idx = num_images; end;

%Loop though each image
for i_im = start_idx:end_idx
    display(['Processing image ' num2str(i_im) ', ' datestr(now)]);
    
    im_name = pred_list(i_im).name(1:end-9);
    
    %Load in images and vessel markup
    try
        vessel_prob = u_load([prob_dir pred_list(i_im).name]);
        vessel_ori = u_load([ori_dir pred_list(i_im).name]);
        %vessel_width = u_load([width_dir pred_list(i_im).name]);
        vessel_im = u_load([image_dir im_name '.mat']);

        vessel_markup_list = dir([markup_dir '*' im_name '*.txt']);    
        vessel_markup = read_markup_from([markup_dir vessel_markup_list(end).name]);
    catch last_err
        display(last_err.message);
        continue;
    end
    
    %Smooth the vessel probs and vessel oris
    vessel_prob = conv2(g', g, vessel_prob, 'same');
    vessel_ori = conv2(g', g, vessel_ori, 'same');
    
    %Equalize image intensities
    vessel_im = equalise_nailfold_intensity(vessel_im);
    
    %Compute curvature
    vessel_curv = abs(complex_gaussian_curvature(vessel_ori, curvature_smoothing_sigma));
    
    %Compute correlation with g2d template
    [~, ~, vessel_corr] = template_match_apexes(vessel_im, ...
        {'g2d', apex_template.g2}, 'threshold', 1);
        
    %Compute NMS centrelines
    vessel_nms = mb_non_maximal_supp(vessel_prob, angle(vessel_ori)/2);
    strong_vessels = vessel_nms > strong_vessel_thresh;
    if any(strong_vessels(:))
        [rstrong cstrong] = find(strong_vessels);
        vessel_centre = bwselect(vessel_nms > 0, cstrong, rstrong, 8);
    else
        vessel_centre = strong_vessels;
    end
    
    %Compute distance transform to centrelines
    vessel_d2c = bwdist(vessel_centre);
    
    %Initialise the background mask
    background_mask = vessel_centre;
    
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
                if isempty(vessel_markup.vessels(i_v).apices.inner_point)                    
                    %Use the anchor
                    axy = round(anchor_xy);
                    
                else
                    %Compute the centre of the apex
                    apex_xy =  ...
                        [ vessel_markup.vessels(i_v).apices(i_a).outer_point;...
                          vessel_markup.vessels(i_v).apices(i_a).inner_point];                     
                    axy = round(mean(apex_xy));
%                     apex_width = sqrt(sum(diff(apex_xy).^2));
                      
                    %Sample patches from curvature, correlation, vessel
                    %prob etc
                    curv_patch = sample_window(vessel_curv, patch_sz, axy(2), axy(1));
                    corr_patch = sample_window(vessel_corr, patch_sz, axy(2), axy(1));
                    prob_patch = sample_window(vessel_prob, patch_sz, axy(2), axy(1));
                    %width_patch = sample_window(vessel_width, patch_sz, axy(2), axy(1));
                    
                    %Compute outputs for apex
                    image_stats.vessel_prob.apex_mean(apex_count,:) = ...
                        mean(prob_patch(circ_mask));
                    image_stats.vessel_prob.apex_max(apex_count,:) = ...
                        max(prob_patch(circ_mask));
                    
                    %
                    image_stats.vessel_curv.apex_mean(apex_count,:) = ...
                        mean(curv_patch(circ_mask));
                    image_stats.vessel_curv.apex_max(apex_count,:) = ...
                        max(curv_patch(circ_mask));
                    
                    %
                    image_stats.vessel_corr.apex_mean(apex_count,:) = ...
                        mean(corr_patch(circ_mask));
                    image_stats.vessel_corr.apex_max(apex_count,:) = ...
                        max(corr_patch(circ_mask));
                    
                    %
%                     image_stats.vessel_width.apex_mean(apex_count,:) = ...
%                         mean(width_patch(circ_mask));
%                     image_stats.vessel_width.apex_truth(apex_count,:) = ...
%                         apex_width;
                    
                    %
                    image_stats.vessel_centres.apex_d2c(apex_count,:) = ...
                        vessel_d2c(axy(2), axy(1));
                    
                    %
                    image_stats.vessel_properties.shape{apex_count,:} = ...
                        vessel_markup.vessels(i_v).ncm_vessel_properties.shape;
                    image_stats.vessel_properties.size{apex_count,:} = ...
                        vessel_markup.vessels(i_v).ncm_vessel_properties.size;
                    image_stats.vessel_properties.im_num(apex_count,:) = ...
                        i_im;                    
                    
                    
                    %Increment the vessel counts
                    apex_count = apex_count + 1;
                      
                end
                %Mask out the region around the apex
                background_mask(axy(2)+circ_y(:,1), axy(1)+circ_x(1,:)) = ...
                    background_mask(axy(2)+circ_y(:,1), axy(1)+circ_x(1,:)) & ncirc_mask;
            end
        else
            %Mask out the region around the anchor
            axy = round(anchor_xy);
            background_mask(axy(2)+circ_y(:,1), axy(1)+circ_x(1,:)) = ...
                background_mask(axy(2)+circ_y(:,1), axy(1)+circ_x(1,:)) & ncirc_mask;
        end
    end
    
    %Compute outputs for background
    image_stats.vessel_prob.background(i_im,:) = ...
        hist(vessel_prob(background_mask), vessel_prob_bins);
    image_stats.vessel_corr.background(i_im,:) = ...
        hist(vessel_corr(background_mask), vessel_corr_bins);
    image_stats.vessel_curv.background(i_im,:) = ...
        hist(vessel_curv(background_mask), vessel_curv_bins);
    image_stats.vessel_centres.background(i_im,:) = ...
        hist(vessel_d2c(background_mask), vessel_d2c_bins);
    
    %Combination of vessel prob and curv/corr
    image_stats.vessel_prob_corr.background(i_im,:) = ...
        hist(vessel_corr(background_mask).* vessel_prob(background_mask), vessel_corr_bins);
    image_stats.vessel_prob_curv.background(i_im,:) = ...
        hist(vessel_curv(background_mask).* vessel_prob(background_mask), vessel_curv_bins);
    
    save([nailfoldroot 'data/' rsa_dir 'models/apex_image_statistics.mat'], 'image_stats');
        
end

%Discard any additional pre-allocated space in the outputs
if apex_count <= max_apex_guess
    image_stats.vessel_prob.apex_mean(apex_count:end,:) = [];
    image_stats.vessel_prob.apex_max(apex_count:end,:) = [];
    %
    image_stats.vessel_curv.apex_mean(apex_count:end,:) = [];
    image_stats.vessel_curv.apex_max(apex_count:end,:) = [];
    %
    image_stats.vessel_corr.apex_mean(apex_count:end,:) = [];
    image_stats.vessel_corr.apex_max(apex_count:end,:) = [];
    %
%     image_stats.vessel_width.apex_mean(apex_count:end,:) = [];
%     image_stats.vessel_width.apex_truth(apex_count:end,:) = [];
    %
    image_stats.vessel_centres.apex_d2c(apex_count:end,:) = [];
    %
    image_stats.vessel_properties.shape(apex_count:end,:) = [];
    image_stats.vessel_properties.size(apex_count:end,:) = [];
    image_stats.vessel_properties.im_num(apex_count:end,:) = [];
end
save([nailfoldroot 'data/' rsa_dir 'models/apex_image_statistics.mat'], 'image_stats');    
    
    
