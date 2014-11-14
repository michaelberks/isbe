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
max_length = 1000;

%Define directories
prob_dir = [nailfoldroot 'data/' rsa_dir '/set12g_half/predictions/detection/rf_classification/304666/'];
ori_dir = [nailfoldroot 'data/' rsa_dir '/set12g_half/predictions/orientation/rf_regression/296621/'];
contour_dir = [nailfoldroot 'data/' rsa_dir '/set12g_half/vessel_contours/'];
centre_dir = [nailfoldroot 'data/' rsa_dir '/set12g_half/vessel_centres/'];

vessel_contour_probs = zeros(101, 101);
bg_contour_probs = zeros(101, 101);

%Loop through each vessel
for i_sz = {'normal', 'enlarged', 'giant'}
    
    vessel_size = i_sz{1};
    apex_dir = [nailfoldroot 'data/' rsa_dir '/apexes/' vessel_size '/'];
    

    %Get list of vessel names
    v_list = dir([prob_dir, vessel_size '*.mat']);
    
    for i_ve = 1:10%length(v_list)%1:20% 

        
        %load vessel
        apex_name = v_list(i_ve).name(length(vessel_size)+(1:8));
        vessel_struc = u_load([apex_dir apex_name '_vessel.mat']);
        contour_struc = load([contour_dir vessel_size apex_name '_vessel_vc.mat'], 'apex_idx', 'vessel_centre');
        centre_struc = load([centre_dir vessel_size apex_name '_vessel_vc.mat'], 'vessel_centre');
        
        %Load in the prob and ori maps
        prob_patch = u_load([prob_dir v_list(i_ve).name]);
        prob_patch = conv2(g', g, prob_patch, 'same');
        ori_patch = u_load([ori_dir v_list(i_ve).name]);
        ori_patch = conv2(g', g, ori_patch, 'same');
        
        display(['Processing ' apex_name ' from ' vessel_size]);  

         
        vcx = centre_struc.vessel_centre.x;
        vcy = centre_struc.vessel_centre.y; clear centre_struc;

        %Resample the vessel pts to be equalled spaced, 1 pixel apart
        v_pts = spline_contour(contour_struc.vessel_centre, [], 1);
        apex_idx = contour_struc.apex_idx;
        apex_xy = [...
            contour_struc.vessel_centre(apex_idx,1)...
            contour_struc.vessel_centre(apex_idx,2)];
        
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
        bg_pts_l = v_pts_l; %bg_pts_l(:,2) = bg_pts_l(:,2) + 40;
        bg_pts_r = v_pts_l; %bg_pts_r(:,2) = bg_pts_r(:,2) + 40;

        [bg_probs_l] = compute_contour_probability(bg_pts_l, rot90(prob_patch,2), rot90(ori_patch,2), 101, linspace(0,1,101));
        [bg_probs_r] = compute_contour_probability(bg_pts_r, rot90(prob_patch,2), rot90(ori_patch,2), 101, linspace(0,1,101));
        bg_contour_probs = bg_contour_probs +...
            bg_probs_l + bg_probs_r;

        if do_plot && i_ve <= 10
            figure;
            a1 = subplot(1,4,1); imgray(imresize(vessel_struc.vessel_patch,0.5));
            plot(v_pts_l(:,1), v_pts_l(:,2), 'g');
            plot(v_pts_r(:,1), v_pts_r(:,2), 'y');
            plot(apex_xy(:,1), apex_xy(:,2), 'r-x');

            a2 = subplot(1,4,2); imgray(prob_patch);
            plot(v_pts_l(:,1), v_pts_l(:,2), 'g');
            plot(v_pts_r(:,1), v_pts_r(:,2), 'y');
            plot(apex_xy(:,1), apex_xy(:,2), 'r-x');
            plot(vcx, vcy, 'b.', 'markersize', 2);
            
            a3 = subplot(1,4,3); imgray(complex2rgb(ori_patch));
            plot(v_pts_l(:,1), v_pts_l(:,2), 'w');
            plot(v_pts_r(:,1), v_pts_r(:,2), 'w');
            plot(apex_xy(:,1), apex_xy(:,2), 'k-x');

            [~,~, path_map] = shape_context_prob_track_mult(...
                ori_patch.*prob_patch, prob_patch, ori_patch,...
                round(v_pts(apex_idx,2)), round(v_pts(apex_idx,1)),...
                'num_streams', 1e4, 'stopping_prob', 0.1);

            a4 = subplot(1,4,4); imgray(path_map);
            plot(v_pts_l(:,1), v_pts_l(:,2), 'g');
            plot(v_pts_r(:,1), v_pts_r(:,2), 'y');
            plot(apex_xy(:,1), apex_xy(:,2), 'r-x');
            plot(vcx, vcy, 'm.', 'markersize', 4);

            linkaxes([a1 a2 a3 a4]);
%             figure; 
%             subplot(2,2,1); plot(contour_probs_l);
%             subplot(2,2,2); plot(contour_probs_r);
%             subplot(2,2,3); plot(bg_probs_l);
%             subplot(2,2,4); plot(bg_probs_r);
        end
    end
end
    
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