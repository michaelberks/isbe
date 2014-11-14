function [] = compute_vessel_corrected_path(angle_with_dist_lc, angle_with_dist_rc)
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
rsa_dir = 'rsa_study';

%Define directories
ori_dir = [nailfoldroot 'data/' rsa_dir '/set12/predictions/orientation/rf_regression/222835/'];

%Loop through each vessel
for i_sz = {'normal', 'enlarged'}
    
    vessel_size = i_sz{1};
    apex_dir = [nailfoldroot 'data/' rsa_dir '/apexes/' vessel_size '/'];
    contour_dir = [nailfoldroot 'data/' rsa_dir '/vessel_contours/' vessel_size '/'];

    %Get list of vessel names
    v_list = dir([ori_dir, vessel_size '*.mat']);
    
    
    for i_ve = 1:length(v_list)%1:20% 

        
        %load vessel
        apex_name = v_list(i_ve).name(length(vessel_size)+(1:8));
        contour_struc = load([contour_dir apex_name '_vessel_contour.mat']);

        %Load in the prob and ori maps
        ori_patch = u_load([ori_dir v_list(i_ve).name]);
        ori_patch = conv2(g', g, ori_patch, 'same');
        
        display(['Processing ' apex_name ' from ' vessel_size]);

        %Resample the vessel pts to be equalled spaced, 1 pixel apart
        v_pts = spline_contour(contour_struc.vessel_centre, [], 1);

        apex_struc = load([apex_dir apex_name '.mat']);
        vessel_struc = u_load([apex_dir apex_name '_vessel.mat']);

        %Correct apex coordinates frame
        apex_xy(:,1) = apex_struc.apex_xy(:,1) +...
            apex_struc.apex_properties.sc - vessel_struc.vessel_properties.sc;
        apex_xy(:,2) = apex_struc.apex_xy(:,2) +...
            apex_struc.apex_properties.sr - vessel_struc.vessel_properties.sr;

        %Find the vessel pt nearets the apex        
        dists = sum(bsxfun(@minus, v_pts, mean(apex_xy)).^2, 2);
        [~, apex_idx] = min(dists);
        
        %Check orientation of apex
        apex_vec = v_pts(apex_idx-1,:) - v_pts(apex_idx+1,:);
        apex_angle = atan2(-apex_vec(2), apex_vec(1));
        if (apex_angle > -pi/2) && (apex_angle < pi/2);
            v_pts = flipud(v_pts);
            apex_idx = size(v_pts,1) - apex_idx + 1;
        end

        %Get predicted orientations
        predicted_ori = interp2(ori_patch, v_pts(2:end-1,1), v_pts(2:end-1,2), 'nearest');
         
        theta_orig = angle(predicted_ori(apex_idx-1))/2;
        path_length = 50;
        

        path_template = zeros(path_length,2);
        path_template(1,:) = v_pts(apex_idx,:);
        
        projected_path_l = ...
            project_path(ori_patch, -cos(theta_orig), sin(theta_orig), path_template, path_length, []);
        
        projected_path_lc = ...
            project_path(ori_patch, -cos(theta_orig), sin(theta_orig), path_template, path_length, angle(angle_with_dist_lc)/2);
        
        projected_path_r = ...
            project_path(ori_patch, cos(theta_orig), -sin(theta_orig), path_template, path_length, []);
        
        projected_path_rc = ...
            project_path(ori_patch, cos(theta_orig), -sin(theta_orig), path_template, path_length, angle(angle_with_dist_rc)/2);
                                      
        
        plot_num = rem(i_ve - 1, 12) + 1;
        if plot_num == 1
            figure;
        end

        subplot(3, 4, plot_num); axis equal ij; hold on;
        plot(v_pts(:,1), v_pts(:,2), 'k');
        plot(v_pts(apex_idx,1), v_pts(apex_idx,2), 'rx');

        plot(projected_path_l(:,1), projected_path_l(:,2), 'r', 'linewidth', 2);
        plot(projected_path_r(:,1), projected_path_r(:,2), 'g', 'linewidth', 2);
        plot(projected_path_lc(:,1), projected_path_lc(:,2), 'm', 'linewidth', 2);
        plot(projected_path_rc(:,1), projected_path_rc(:,2), 'b', 'linewidth', 2);

    end
end

function projected_path = ...
    project_path(ori_map, xi, yi, projected_path, path_length, correction_table)

for i_step = 2:path_length
    theta = angle(ori_map(...
        round(projected_path(i_step-1,2)),...
        round(projected_path(i_step-1,1)))) / 2;
    
    if ~isempty(correction_table)
        theta = theta + correction_table(i_step-1);
    end
    
    xj = cos(theta);
    yj = -sin(theta);

    if (xi*xj) + (yi*yj) < 0
        xi = -xj;
        yi = -yj;
    else
        xi = xj;
        yi = yj;
    end
    projected_path(i_step,1) = projected_path(i_step-1,1) + xi;
    projected_path(i_step,2) = projected_path(i_step-1,2) + yi;
end