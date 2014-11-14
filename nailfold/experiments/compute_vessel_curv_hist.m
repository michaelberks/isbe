function [curv_stats] = compute_vessel_curv_hist()
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
debug = 0;
max_length = 100;
vessel_prob_smoothing_sigma = 2;
g = gaussian_filters_1d(vessel_prob_smoothing_sigma);
g = g / sum(g);
rsa_dir = 'rsa_study';

%Define directories
ori_dir = [nailfoldroot 'data/' rsa_dir '/set12/predictions/orientation/rf_regression/222835/'];

n_ori = 180;
n_mags = 50;

%Set up containers
for field = {'l', 'r'}
    
    curv_stats.(field{1}).pred_error_by_step.hist = zeros(max_length, n_ori);
    curv_stats.(field{1}).pred_error_by_step.dist = zeros(max_length, 1);

    curv_stats.(field{1}).marked_ori_by_step.hist = zeros(max_length, n_ori);
    curv_stats.(field{1}).marked_ori_by_step.dist = zeros(max_length, 1);

    curv_stats.(field{1}).pred_ori_by_step.hist = zeros(max_length, n_ori);
    curv_stats.(field{1}).pred_ori_by_step.dist = zeros(max_length, 1);

    curv_stats.(field{1}).pred_mag_by_step.hist = zeros(max_length, n_mags);

    curv_stats.(field{1}).marked_curv_by_step.hist = zeros(max_length, n_ori);
    curv_stats.(field{1}).marked_curv_by_step.dist = zeros(max_length, 1);

    curv_stats.(field{1}).pred_curv_by_step.hist = zeros(max_length, n_ori);
    curv_stats.(field{1}).pred_curv_by_step.dist = zeros(max_length, 1);

    curv_stats.(field{1}).map_curv_by_step.hist = zeros(max_length, n_ori);
    curv_stats.(field{1}).map_curv_by_step.dist = zeros(max_length, 1);
end

%Loop through each vessel
for i_sz = {'normal', 'enlarged'}
    
    vessel_size = i_sz{1};
    apex_dir = [nailfoldroot 'data/' rsa_dir '/apexes/' vessel_size '/'];
    contour_dir = [nailfoldroot 'data/' rsa_dir '/vessel_contours/' vessel_size '/'];

    %Get list of vessel names
    v_list = dir([ori_dir, vessel_size '*.mat']);
    
    for i_ve = 1:length(v_list)

        
        %load vessel
        apex_name = v_list(i_ve).name(length(vessel_size)+(1:8));
        contour_struc = load([contour_dir apex_name '_vessel_contour.mat']);

        %Load in the prob and ori maps
        ori_patch = u_load([ori_dir v_list(i_ve).name]);
        ori_patch = conv2(g', g, ori_patch, 'same');
        
        curv_patch = complex_gaussian_curvature(ori_patch ./ (abs(ori_patch) + 1e-6), 1);
        
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
            
        v_pts_l = flipud(v_pts(1:apex_idx,:));
        v_pts_r = v_pts(apex_idx:end,:);
        
        curv_stats.l = ...
            compute_stats(ori_patch, curv_patch, curv_stats.l, v_pts_l, max_length, n_ori, n_mags, debug && (i_ve <= 10));
        
        curv_stats.r = ...
            compute_stats(ori_patch, curv_patch, curv_stats.r, v_pts_r, max_length, n_ori, n_mags, debug && (i_ve <= 10));
        
        
        if debug && i_ve <= 40
            theta_orig = angle(ori_patch(round(v_pts(apex_idx,2)), round(v_pts(apex_idx,2))))/2;
            path_length = 50;


            path_template = zeros(path_length,2);
            path_template(1,:) = v_pts(apex_idx,:);

            projected_path_l = ...
                project_path(ori_patch, -cos(theta_orig), sin(theta_orig), path_template, path_length, []);

            projected_path_lc = ...
                project_path(ori_patch, -cos(theta_orig), sin(theta_orig), path_template, path_length, angle(curv_stats.pred_error_by_step_lci)/2);

            projected_path_r = ...
                project_path(ori_patch, cos(theta_orig), -sin(theta_orig), path_template, path_length, []);

            projected_path_rc = ...
                project_path(ori_patch, cos(theta_orig), -sin(theta_orig), path_template, path_length, angle(curv_stats.pred_error_by_step_rci)/2);


            plot_num = rem(i_ve - 1, 12) + 1;
            if plot_num == 1
                figure;
            end

            subplot(3, 4, plot_num); axis equal ij; hold on;
            plot(v_pts(:,1), v_pts(:,2), 'k');
            plot(v_pts(apex_idx,1), v_pts(apex_idx,2), 'rx');

            plot(projected_path_l(:,1), projected_path_l(:,2), 'r', 'linewidth', 2);
            plot(projected_path_r(:,1), projected_path_r(:,2), 'g', 'linewidth', 2);
            plot(projected_path_lc(:,1), projected_path_lc(:,2), 'm--', 'linewidth', 2);
            plot(projected_path_rc(:,1), projected_path_rc(:,2), 'b--', 'linewidth', 2);
        end
    end
end

%Normalise distribution by hist counts

for field = {'l', 'r'}
    
    curv_stats.(field{1}).pred_error_by_step.dist = ...
        curv_stats.(field{1}).pred_error_by_step.dist ./ ...
        sum(curv_stats.(field{1}).pred_error_by_step.hist,2);

    curv_stats.(field{1}).marked_ori_by_step.dist = ...
        curv_stats.(field{1}).marked_ori_by_step.dist ./ ...
        sum(curv_stats.(field{1}).marked_ori_by_step.hist,2);

    curv_stats.(field{1}).pred_ori_by_step.dist = ...
        curv_stats.(field{1}).pred_ori_by_step.dist ./ ...
        sum(curv_stats.(field{1}).pred_ori_by_step.hist,2);

    curv_stats.(field{1}).marked_curv_by_step.dist = ...
        curv_stats.(field{1}).marked_curv_by_step.dist ./ ...
        sum(curv_stats.(field{1}).marked_curv_by_step.hist,2);

    curv_stats.(field{1}).pred_curv_by_step.dist = ...
        curv_stats.(field{1}).pred_curv_by_step.dist ./ ...
        sum(curv_stats.(field{1}).pred_curv_by_step.hist,2);

    curv_stats.(field{1}).map_curv_by_step.dist= ...
        curv_stats.(field{1}).map_curv_by_step.dist ./ ...
        sum(curv_stats.(field{1}).map_curv_by_step.hist,2);
end

%--------------------------------------------------------------------------
function curv_stats_i = compute_stats(ori_map, curv_map, curv_stats_i, v_pts, max_length, n_ori, n_mags, debug)

%Workout orientation and change of orientation of each step in the
%marked vessel contour, normalise to orientation of 1st step
xy_steps = diff(v_pts);
marked_ori = exp(2i*atan(-xy_steps(:,2) ./ xy_steps(:,1)));
normalising_ori = conj(marked_ori(1));

marked_ori = marked_ori * normalising_ori;
marked_ori_bin_idx = ceil(.5*n_ori*(angle(marked_ori) + pi)/pi);

marked_curv = marked_ori(2:end) .* conj(marked_ori(1:end-1));
marked_curv_bin_idx = ceil(.5*n_ori*(angle(marked_curv) + pi)/pi);

%Get predicted orientations - normalise to orientation of 1st step
pred_ori = interp2(ori_map, v_pts(1:end-1,1), v_pts(1:end-1,2), 'nearest');
pred_ori = pred_ori * normalising_ori;
pred_ori_bin_idx = ceil(.5*n_ori*(angle(pred_ori) + pi)/pi);

pred_mag = abs(pred_ori);
pred_mag_bin_idx = ceil(n_mags*pred_mag);

%GEt predicted curvature (both from the map, and computed from the sampled
%orientations)
map_curv = interp2(curv_map, v_pts(1:end-1,1), v_pts(1:end-1,2), 'nearest');
map_curv_bin_idx = ceil(n_ori*(map_curv + pi/2)/pi);

pred_curv = pred_ori(2:end) .* conj(pred_ori(1:end-1));
pred_curv_bin_idx = ceil(.5*n_ori*(angle(pred_curv) + pi)/pi);

%Compute difference between predicted and actual orientations
pred_error = marked_ori .* exp(-1i*angle(pred_ori));
pred_error_bin_idx = ceil(.5*n_ori*(angle(pred_error) + pi)/pi);

for i_pt = 1:min(size(v_pts,1)-2, max_length)

    %Increment histogram counts
    curv_stats_i.pred_error_by_step.hist(i_pt, pred_error_bin_idx(i_pt)) =...
        curv_stats_i.pred_error_by_step.hist(i_pt, pred_error_bin_idx(i_pt)) + 1;
    
    curv_stats_i.marked_ori_by_step.hist(i_pt, marked_ori_bin_idx(i_pt)) =...
        curv_stats_i.marked_ori_by_step.hist(i_pt, marked_ori_bin_idx(i_pt)) + 1;    
    
    curv_stats_i.pred_ori_by_step.hist(i_pt, pred_ori_bin_idx(i_pt)) =...
        curv_stats_i.pred_ori_by_step.hist(i_pt, pred_ori_bin_idx(i_pt)) + 1;
    
    curv_stats_i.pred_mag_by_step.hist(i_pt, pred_mag_bin_idx(i_pt)) =...
        curv_stats_i.pred_mag_by_step.hist(i_pt, pred_mag_bin_idx(i_pt)) + 1;

    curv_stats_i.marked_curv_by_step.hist(i_pt, marked_curv_bin_idx(i_pt)) =...
        curv_stats_i.marked_curv_by_step.hist(i_pt, marked_curv_bin_idx(i_pt)) + 1;
    
    curv_stats_i.pred_curv_by_step.hist(i_pt, pred_curv_bin_idx(i_pt)) =...
        curv_stats_i.pred_curv_by_step.hist(i_pt, pred_curv_bin_idx(i_pt)) + 1;
    
    curv_stats_i.map_curv_by_step.hist(i_pt, map_curv_bin_idx(i_pt)) =...
        curv_stats_i.map_curv_by_step.hist(i_pt, map_curv_bin_idx(i_pt)) + 1;
    
    
    %Add to distributions
    curv_stats_i.pred_error_by_step.dist(i_pt) = ...
        curv_stats_i.pred_error_by_step.dist(i_pt) + pred_error(i_pt);
   
    curv_stats_i.marked_ori_by_step.dist(i_pt) = ...
        curv_stats_i.marked_ori_by_step.dist(i_pt) + marked_ori(i_pt);

    curv_stats_i.pred_ori_by_step.dist(i_pt) = ...
        curv_stats_i.pred_ori_by_step.dist(i_pt) + pred_ori(i_pt);

    curv_stats_i.marked_curv_by_step.dist(i_pt) = ...
        curv_stats_i.marked_curv_by_step.dist(i_pt) + marked_curv(i_pt);

    curv_stats_i.pred_curv_by_step.dist(i_pt) = ...
        curv_stats_i.pred_curv_by_step.dist(i_pt) + pred_curv(i_pt);

    curv_stats_i.map_curv_by_step.dist(i_pt) = ...
        curv_stats_i.map_curv_by_step.dist(i_pt) + map_curv(i_pt);

end

%Plot some stuff if we're debugging
if debug
         
    figure;
    a1 = subplot(1,2,1); hold off; plot(v_pts(:,1), v_pts(:,2)); axis ij equal; hold on;
    a2 = subplot(1,2,2); hold off; plot(v_pts(:,1), v_pts(:,2)); axis ij equal; hold on;     

    for i_pt = 1:size(v_pts,1)-2

        if angle(pred_error(i_pt)) > .5
            plot(a2, ...
                v_pts(i_pt+1,1)+[-1 1]*cos(angle(pred_ori(i_pt))/2),... %cos(marked_ori_change(i_pt)/2) + 
                v_pts(i_pt+1,2)-[-1 1]*sin(angle(pred_ori(i_pt))/2), 'r');
        elseif angle(ori_diff(i_pt)) > -.5
            plot(a2, ...
                v_pts(i_pt+1,1)+[-1 1]*cos(angle(pred_ori(i_pt))/2),... %cos(marked_ori_change(i_pt)/2) + 
                v_pts(i_pt+1,2)-[-1 1]*sin(angle(pred_ori(i_pt))/2), 'y');
        else
            plot(a2, ...
                v_pts(i_pt+1,1)+[-1 1]*cos(angle(pred_ori(i_pt))/2),... %cos(marked_ori_change(i_pt)/2) + 
                v_pts(i_pt+1,2)-[-1 1]*sin(angle(pred_ori(i_pt))/2), 'g');
        end
        text(v_pts(i_pt+1,1), v_pts(i_pt+1,2), num2str(pred_error_bin_index(i_pt)));

        if angle(marked_curv(i_pt))/2 > 0.2
            plot(a1, ...
                v_pts(i_pt+1,1)+[0 xy_steps(i_pt,1)],... %cos(marked_ori_change(i_pt)/2) + 
                v_pts(i_pt+1,2)+[0 xy_steps(i_pt,2)], 'r');  %sin(marked_ori_change(i_pt)/2) + 
        elseif angle(marked_curv(i_pt))/2 > -0.2
            plot(a1, ...
                v_pts(i_pt+1,1)+[0 xy_steps(i_pt,1)],... %cos(marked_ori_change(i_pt)/2) + 
                v_pts(i_pt+1,2)+[0 xy_steps(i_pt,2)], 'y');  %sin(marked_ori_change(i_pt)/2) +
        else
            plot(a1, ...
                v_pts(i_pt+1,1)+[0 xy_steps(i_pt,1)],... %cos(marked_ori_change(i_pt)/2) + 
                v_pts(i_pt+1,2)+[0 xy_steps(i_pt,2)], 'g'); %sin(marked_ori_change(i_pt)/2) + 
        end
    end
    linkaxes([a1 a2]);

end

%--------------------------------------------------------------------------
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