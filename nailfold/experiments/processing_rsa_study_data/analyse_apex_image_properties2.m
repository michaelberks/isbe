function [] = analyse_apex_image_properties2(start_idx, end_idx)
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
vessel_prob_smoothing_sigma = 2;

%Set up image lists
if ispc
    rsa_dir = 'rsa_study/';
else
    rsa_dir = [];
end

model_dir = [nailfoldroot 'data/' rsa_dir 'models/apex_templates/'];
prob_dir = [nailfoldroot 'data/' rsa_dir 'test/predictions/detection/rf_classification/222836/'];
ori_dir = [nailfoldroot 'data/' rsa_dir 'test/predictions/orientation/rf_regression/222835/'];
centre_dir = [nailfoldroot 'data/' rsa_dir 'test/vessel_centres/'];
% shape_context_dir = [nailfoldroot 'data/' rsa_dir 'test/shape_contexts/'];
% create_folder(shape_context_dir);

apex_prior_dir = [nailfoldroot 'data/' rsa_dir 'test/apex_priors/'];
create_folder(apex_prior_dir);

pred_list = dir([prob_dir '*.mat']);

%Create smoothing kernel for vessel probs
g = gaussian_filters_1d(vessel_prob_smoothing_sigma);
g = g / sum(g);

load([model_dir 'apex_curv_stats.mat']);
ori_hist_r_smoothed = imfilter(curv_stats.r.pred_ori_by_step.hist, conv(g,g), 'circular');
ori_hist_r_smoothed = imfilter(ori_hist_r_smoothed, conv(g,g)', 'replicate');

ori_hist_l_smoothed = imfilter(curv_stats.l.pred_ori_by_step.hist, conv(g,g), 'circular');
ori_hist_l_smoothed = imfilter(ori_hist_l_smoothed, conv(g,g)', 'replicate');

apex_prior = cat(3, ori_hist_l_smoothed, ori_hist_r_smoothed);

%Loop though each image
for i_im = start_idx:end_idx
    display(['Processing image ' num2str(i_im) ', ' datestr(now)]);
    
    im_name = pred_list(i_im).name(1:6);
    
    %Load in images and vessel markup
    try
        vessel_prob = u_load([prob_dir pred_list(i_im).name]);
        vessel_ori = u_load([ori_dir pred_list(i_im).name]);
        load([centre_dir im_name '_vc.mat'], 'vessel_centre*');
    catch last_err
        display(last_err.message);
        continue;
    end
    
    %Smooth the vessel probs
    vessel_prob = conv2(g', g, vessel_prob, 'same');
    %vessel_ori = conv2(g', g, vessel_ori, 'same');
    
    [nrows ncols] = size(vessel_prob);
    vessel_centre_mask = false(nrows, ncols);
    vessel_centre_idx = sub2ind([nrows ncols], vessel_centre_y, vessel_centre_x);
    vessel_centre_mask(vessel_centre_idx) = 1;
    
    vessel_curv = zeros(nrows, ncols);
    vessel_curv(vessel_centre_mask) = abs(vessel_centre_curv);
    
    %Update the centres to include only high curvatur epoints
    apex_mask = vessel_centre_mask & (vessel_curv > 0.025) & (vessel_prob > 0.8);
        
    %Now take only the highest curvature from each connected section
    centre_labels = bwlabel(apex_mask, 8);
    for i_lab = 1:max(centre_labels(:))
        label_mask = centre_labels == i_lab;
        label_curv = vessel_curv(label_mask);
        label_bw = false(size(label_curv));
        [~, max_idx] = max(label_curv);
        label_bw(max_idx) = 1;
        apex_mask(label_mask) = label_bw;
    end
    include_pts = apex_mask(vessel_centre_mask);
    apex_x = vessel_centre_x(include_pts);
    apex_y = vessel_centre_y(include_pts);
    
%     [sc_pi pc_pi] = shape_context_prob_track_mult(...
%          vessel_ori.*vessel_prob, vessel_prob, vessel_ori,...
%          vessel_centre_y(include_pts), vessel_centre_x(include_pts),...
%          'num_streams', 1e3, 'stopping_prob', 0.2); %#ok
%      
%      save([shape_context_dir im_name '_sc.mat'], 'sc_pi', 'pc_pi');
     
     [apex_prior_scores] = apex_prior_prob_track_mult(...
            vessel_prob, vessel_ori,...
            apex_y, apex_x,...
            'num_streams', 1e3, 'stopping_prob', 0.1,...
            'apex_prior', apex_prior); %#ok
     
     save([apex_prior_dir im_name '_ap.mat'], 'apex_prior_scores', 'apex_x', 'apex_y');
            
end  
    
    
