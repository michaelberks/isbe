function [vessel_contour_probs bg_contour_probs] = compute_matched_vessel_contour_probs()
%COMPUTE_MATCHED_VESSEL_CONTOUR_PROBS *Insert a one line summary here*
%   [] = compute_matched_vessel_contour_probs()
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
% Created: 28-Jun-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

%Define constants
vessel_prob_smoothing_sigma = 2;
g = gaussian_filters_1d(vessel_prob_smoothing_sigma);
g = g / sum(g);
do_plot = 1;
rsa_dir = 'rsa_study';
max_length = 100;

%Define directories
vessel_size = 'normal';
vessel_dir = [nailfoldroot 'data/' rsa_dir '/apexes/' vessel_size '/'];
prob_dir = [nailfoldroot 'data/' rsa_dir '/test/predictions/detection/rf_classification/182321/'];
ori_dir = [nailfoldroot 'data/' rsa_dir '/test/predictions/orientation/rf_regression/182263/'];

%Get list of vessel names
v_list = dir([vessel_dir, '*vessel.mat']);

curr_nailfold_name = 'frog';

vessel_contour_probs = zeros(101, 101);
bg_contour_probs = zeros(101, 101);

v_count = 1;

%Loop through each vessel
for i_ve = 1:length(v_list)

    %load vessel
    vessel_name = [vessel_dir, v_list(i_ve).name];
    vessel_struc = u_load(vessel_name);
    
    %check it is processable
    if vessel_struc.quality < 2
        continue;
    end
    
    %check whether we've already loaded the main image maps for this vessel
    nailfold_name = vessel_struc.vessel_properties.nailfold_name;
    [~, nailfold_name] = fileparts(nailfold_name);
    
    %If not load them
    if ~strcmpi(curr_nailfold_name, nailfold_name)
        curr_nailfold_name = nailfold_name;
               
        if exist([prob_dir nailfold_name '_pred.mat'], 'file');
            display(['Loading ' nailfold_name]);
            
            vessel_prob = u_load([prob_dir nailfold_name '_pred.mat']);
            vessel_ori = u_load([ori_dir nailfold_name '_pred.mat']);
            vessel_prob = conv2(g', g, vessel_prob, 'same');
        end
    end
    
    %If there definietly are matching image maps
    if exist([prob_dir nailfold_name '_pred.mat'], 'file');
        display(['Processing ' vessel_name ' from ' nailfold_name]);
        
        %Workout which bit of the main ori and prob maps to extract
        sr = vessel_struc.vessel_properties.sr;
        sc = vessel_struc.vessel_properties.sc;
        er = sr + size(vessel_struc.vessel_patch, 1) - 1;
        ec = sc + size(vessel_struc.vessel_patch, 2) - 1;
        
        prob_patch = vessel_prob(sr:er, sc:ec);
        ori_patch = vessel_ori(sr:er, sc:ec);
        
        %Load in the apex
        apex_struc = load([vessel_dir v_list(i_ve).name(1:8) '.mat']);
        apex_xy(:,1) = apex_struc.apex_xy(:,1) +...
            apex_struc.apex_properties.sc - vessel_struc.vessel_properties.sc;
        apex_xy(:,2) = apex_struc.apex_xy(:,2) +...
            apex_struc.apex_properties.sr - vessel_struc.vessel_properties.sr;
        
        %Resample the vessel pts to be equalled spaced, 1 pixel apart
        v_pts = spline_contour(vessel_struc.v_pts, [], 1);
        v_pts = correct_pts(prob_patch, v_pts);
        
        %Find the vessel pt nearets the apex        
        dists = sum(bsxfun(@minus, v_pts, mean(apex_xy)).^2, 2);
        [~, apex_idx] = min(dists);
        
        %Set max length?
        
        %Now separate out the left and right halves of the vessel contour
        v_pts_l = flipud(v_pts(1:apex_idx,:));
        v_pts_r = v_pts(apex_idx:end,:);
        
        v_pts_l(max_length:end,:) = [];
        v_pts_r(max_length:end,:) = [];
        
        
        %Compute the contour probs for the left and right halves
        [contour_probs_l] = compute_contour_probability(v_pts_l, prob_patch, ori_patch, 101, linspace(0,1,101));
        [contour_probs_r] = compute_contour_probability(v_pts_r, prob_patch, ori_patch, 101, linspace(0,1,101));
        
        %Add to the main sum
        vessel_contour_probs = vessel_contour_probs +...
            contour_probs_l + contour_probs_r;
        
        %Now compute matched probs for the backgorund contours...?
        bg_pts_l = v_pts_l; bg_pts_l(:,2) = bg_pts_l(:,2) + 40;
        bg_pts_r = v_pts_l; bg_pts_r(:,2) = bg_pts_r(:,2) + 40;
        
        [bg_probs_l] = compute_contour_probability(bg_pts_l, prob_patch, ori_patch, 101, linspace(0,1,101));
        [bg_probs_r] = compute_contour_probability(bg_pts_r, prob_patch, ori_patch, 101, linspace(0,1,101));
        bg_contour_probs = bg_contour_probs +...
            bg_probs_l + bg_probs_r;
        
        if do_plot && v_count <= 10
            figure;
            subplot(1,4,1); imgray(vessel_struc.vessel_patch);
            plot(v_pts_l(:,1), v_pts_l(:,2), 'g');
            plot(v_pts_r(:,1), v_pts_r(:,2), 'y');
            plot(apex_xy(:,1), apex_xy(:,2), 'r-x');
            
            subplot(1,4,2); imgray(prob_patch);
            plot(v_pts_l(:,1), v_pts_l(:,2), 'g');
            plot(v_pts_r(:,1), v_pts_r(:,2), 'y');
            plot(apex_xy(:,1), apex_xy(:,2), 'r-x');

            subplot(1,4,3); imgray(complex2rgb(ori_patch));
            plot(vessel_struc.v_pts(:,1), vessel_struc.v_pts(:,2), 'w');
            plot(apex_xy(:,1), apex_xy(:,2), 'k-x');
            
            [~,~, path_map] = shape_context_prob_track_mult(...
                ori_patch.*prob_patch, prob_patch, ori_patch,...
                round(v_pts(apex_idx,2)), round(v_pts(apex_idx,1)),...
                'num_streams', 1e4, 'stopping_prob', 0.03);
            
            subplot(1,4,4); imgray(path_map);
            plot(v_pts_l(:,1), v_pts_l(:,2), 'g');
            plot(v_pts_r(:,1), v_pts_r(:,2), 'y');
            plot(apex_xy(:,1), apex_xy(:,2), 'r-x');
            
            figure; 
            subplot(2,2,1); plot(contour_probs_l);
            subplot(2,2,2); plot(contour_probs_r);
            subplot(2,2,3); plot(bg_probs_l);
            subplot(2,2,4); plot(bg_probs_r);
        end
        
        v_count = v_count + 1;
        
        
        
    end
end
    %v_pts = vessel_struc.v_pts;
    
function v_pts = correct_pts(vessel_prob, v_pts)

    v_samp = spline_contour(v_pts, [], 5);
    v_normal = compute_spline_normals(v_samp);
    for i_pt = 1:size(v_samp,1)
        [nx ny np] = improfile(vessel_prob, ...
            [-5 5]*v_normal(i_pt,1) + v_samp(i_pt,1),...
            [-5 5]*v_normal(i_pt,2) + v_samp(i_pt,2));
        [~,max_i] = max(np);
        v_samp(i_pt,:) = [nx(max_i) ny(max_i)];
    end
    
    v_pts = spline_contour(v_samp, [], 1);
        
        
