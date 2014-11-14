function [outer_edge_xy, inner_edge_xy, vessel_mask] = detect_vessel_edge(vessel_patch, vessel_pts, apex_xy, varargin)
%DETECT_VESSEL_EDGE *Insert a one line summary here*
%   [outer_edge, inner_edge, vessel_mask] = detect_vessel_edge(vessel_patch, vessel_pts)
%
% Inputs:
%      vessel_patch - *Insert description of input variable here*
%
%      vessel_pts - *Insert description of input variable here*
%
%
% Outputs:
%      outer_edge - *Insert description of input variable here*
%
%      inner_edge - *Insert description of input variable here*
%
%      vessel_mask - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 22-Mar-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

args = u_packargs(varargin,... % the user's input
    '0', ... % strict mode
    'sigma0', 2,...
    'sigma1', 2,...
    'alpha', 0.5,...
    'beta', 1,...
    'gamma', 3,...
    'resolution', 1,...
    'initial_spacing', 2,...
    'search_width_fraction', 0.1,...
    'iterate', true,...
    'plot', 1,...
    'quiet', 0);
 clear varargin;

%Constants
sigma0_w = round(5*args.sigma0);
num_pts = size(vessel_pts,1);

%Smooth vessel path
g = gaussian_filters_1d(args.sigma0, sigma0_w);
g = g / sum(g);

vessel_pts = [...
    repmat(vessel_pts(1,:), sigma0_w, 1);...
    vessel_pts;...
    repmat(vessel_pts(end,:), sigma0_w, 1)];

v_pts_sm = [conv(vessel_pts(:,1), g, 'same') conv(vessel_pts(:,2), g, 'same')];
v_pts_sm = v_pts_sm(sigma0_w+(1:num_pts),:);

%Fit spline to vessel path to get evenly spaced points
[v_pts_sm] = spline_contour(v_pts_sm, [], args.initial_spacing);
num_pts = size(v_pts_sm,1);

%Project normals out from path, and use profiles to detect outer edge
normal_xy = compute_spline_normals(v_pts_sm);

%Find the apex pt and check the y component for
%this points up (is negative). Otherwise reverse the normals
[~,~,apex_pt] = line_cross(v_pts_sm, apex_xy);
if normal_xy(apex_pt,2) > 0
    normal_xy = -normal_xy;
end

apex_width = sqrt(sum(diff(apex_xy).^2));
apex_halfwidth = floor(apex_width/2);
centre_prof_width = [-round(apex_width) round(apex_width)];
[g1d_responses] = compute_gaussian_1st_derivatives(vessel_patch, apex_width/2);
[g2d_responses] = compute_gaussian_2nd_derivatives(vessel_patch, apex_width/2);
[h2d_responses] = compute_hilbert_2nd_derivatives_sep(vessel_patch, apex_width/2);

%Make an estimated mask of the vessel based on the initial annotation and
%the appex width
initial_vessel_mask = false(size(vessel_patch));
vessel_idx = sub2ind(size(vessel_patch), round(vessel_pts(:,2)), round(vessel_pts(:,1)));
initial_vessel_mask(vessel_idx) = 1;
initial_vessel_mask = imdilate(initial_vessel_mask, strel('disk', floor(0.8*apex_width/2)));

% Pre-allocate containers for the normal profiles and the associated
% sampling points
norm_width = diff(centre_prof_width)+1;
normal_p = zeros(num_pts, norm_width);
normal_x = zeros(num_pts, norm_width);
normal_y = zeros(num_pts, norm_width);
normal_mask = false(num_pts, norm_width);
normal_g2d = zeros(num_pts, norm_width);
normal_h2d = zeros(num_pts, norm_width);
normal_g1d = zeros(num_pts, norm_width);

for i_n = 1:num_pts  
    %Get normal profile
    n_x = v_pts_sm(i_n,1)+normal_xy(i_n,1)*centre_prof_width;
    n_y = v_pts_sm(i_n,2)+normal_xy(i_n,2)*centre_prof_width;
    
    theta = atan2(-normal_xy(i_n,2), normal_xy(i_n,1));
    
    g1d_responses_p = zeros(1,norm_width,1,2);
    for i_g = 1:2
        g1d_responses_p(:,:,:,i_g) = improfile(g1d_responses(:,:,1,i_g), n_x, n_y, norm_width, 'bicubic')';
    end
    g2d_responses_p = zeros(1,norm_width,1,3);
    for i_g = 1:3
        g2d_responses_p(:,:,:,i_g) = improfile(g2d_responses(:,:,1,i_g), n_x, n_y, norm_width, 'bicubic')';
    end
    h2d_responses_p = zeros(1,norm_width,1,4);
    for i_h = 1:4
        h2d_responses_p(:,:,:,i_h) = improfile(h2d_responses(:,:,1,i_h), n_x, n_y, norm_width, 'bicubic')';
    end
        
    [cx, cy, cp] = improfile(vessel_patch, n_x, n_y, norm_width, 'bilinear');
    normal_p(i_n, :) = cp';
    normal_x(i_n, :) = cx';
    normal_y(i_n, :) = cy';
    
    mask_p = improfile(initial_vessel_mask, n_x, n_y, norm_width, 'nearest');
    mask_p(isnan(mask_p)) = 0;
    normal_mask(i_n, :) = mask_p;
    
    normal_g1d(i_n, :) = steer_gaussian_1st_derivatives(g1d_responses_p, theta);
    normal_g2d(i_n, :) = steer_gaussian_2nd_derivatives(g2d_responses_p, theta);
    normal_h2d(i_n, :) = steer_hilbert_2nd_derivatives(h2d_responses_p, theta);
end

centre_col = ceil(norm_width/2);
edge_mask = false(num_pts, norm_width);
centre_mask = false(num_pts, norm_width);

centre_xy = zeros(num_pts,2);
centre_idx = NaN(num_pts,1);
outer_edge_xy = NaN(num_pts,2);
inner_edge_xy = NaN(num_pts,2);

if args.plot
    figure; 
    axes1 = subplot(2,3,1); imgray(initial_vessel_mask); %#ok
    axes2 = subplot(2,3,2); imgray(normal_mask);
    axes3 = subplot(2,3,3); imgray(normal_g2d);
end
for i_n = 1:num_pts
    profile_phase = atan(normal_g2d(i_n,:) ./ normal_h2d(i_n,:));
    poss_edges = abs(profile_phase) < 0.2;
    
    left_edge_i = find(~normal_mask(i_n,centre_col:-1:1), 1, 'first');
    if isempty(left_edge_i)
        left_edge_i = 1;
    else
        left_edge_i = centre_col - left_edge_i + 2;
    end
    initial_left_edge = false(1,norm_width);
    initial_left_edge(max(1,left_edge_i-apex_halfwidth):left_edge_i+apex_halfwidth)=1;
    
    right_edge_i = find(~normal_mask(i_n,centre_col:end), 1, 'first');
    if isempty(right_edge_i)
        right_edge_i = norm_width;
    else
        right_edge_i = centre_col + right_edge_i - 2;
    end
    initial_right_edge = false(1,norm_width);
    initial_right_edge(right_edge_i-apex_halfwidth:min(norm_width,right_edge_i+apex_halfwidth))=1;   
    
    left_edge = find(poss_edges(1:centre_col-1)&initial_left_edge(1:centre_col-1), 1, 'last');
    if ~isempty(left_edge)
        edge_mask(i_n,1:left_edge) = 1;
    end
    right_edge = find(poss_edges(centre_col+1:end)&initial_right_edge(centre_col+1:end), 1, 'first');
    if ~isempty(right_edge)
        right_edge = right_edge + centre_col;
        
        plot(axes2, right_edge, i_n, 'rx');
        plot(axes3, right_edge, i_n, 'rx');
        
        edge_mask(i_n,right_edge:end) = 1;
        
        outer_edge_xy(i_n,1) = normal_x(i_n, right_edge);
        outer_edge_xy(i_n,2) = normal_y(i_n, right_edge);

        centre_mask(i_n,:) = false(1,norm_width);
        centre_mask(i_n,centre_col+(-apex_halfwidth:apex_halfwidth)) = 1;
        centre_mask(i_n,:) = centre_mask(i_n,:) & (profile_phase > (pi/2 - 0.1));

        poss_centre = find(diff(normal_g2d(i_n,right_edge:-1:round(right_edge-apex_width))) <= 0, 1, 'first');
        if ~isempty(poss_centre) && (poss_centre > 1)
            centre_idx(i_n) = right_edge-poss_centre+1;
            centre_xy(i_n,1) = normal_x(i_n, centre_idx(i_n));
            centre_xy(i_n,2) = normal_y(i_n, centre_idx(i_n));
            
            plot(axes2, centre_idx(i_n), i_n, 'yx');
            plot(axes3, centre_idx(i_n), i_n, 'yx');
            
            left_edge = 2*centre_idx(i_n)-right_edge;
            if (left_edge > 1) && initial_left_edge(left_edge)
                plot(axes2, 2*centre_idx(i_n)-right_edge, i_n, 'gx');
                plot(axes3, 2*centre_idx(i_n)-right_edge, i_n, 'gx');

                inner_edge_xy(i_n,:) = 2*centre_xy(i_n,:) - outer_edge_xy(i_n,:);
            end
        end                   
    end
end

%Get list of points that have been found
outer_idx_i = ~any(isnan(outer_edge_xy),2);
inner_idx_i = ~any(isnan(inner_edge_xy),2);

% %Median filter the points to get rid of rogue outliers
% outer_edge_xy(outer_idx_i,:) = medfilt1(outer_edge_xy(outer_idx_i,:), 5);
% inner_edge_xy(inner_idx_i,:) = medfilt1(inner_edge_xy(inner_idx_i,:), 5);

%Remove outer edge points that cross
[outer_idx_j] = remove_overlaps(centre_xy(outer_idx_i,:), outer_edge_xy(outer_idx_i,:), 10);
outer_idx = false(num_pts,1);
outer_idx(outer_idx_i) = outer_idx_j;
outer_edge_xy(~outer_idx,:) = NaN;

%Remove inner points that cross
[inner_idx_j] = remove_overlaps(centre_xy(inner_idx_i,:), inner_edge_xy(inner_idx_i,:), 10);
inner_idx = false(num_pts,1);
inner_idx(inner_idx_i) = inner_idx_j;
inner_edge_xy(~inner_idx,:) = NaN;
%
dists_i = cumsum([0; sum(diff(outer_edge_xy(outer_idx,:)).^2,2)]);
dists = zeros(num_pts,1);
dists(outer_idx,:) = dists_i;
evenly_spaced_dists = linspace(0, dists_i(end), round(dists_i(end)/2)+1);

%First lets pair up any inner point that is missing an outer point
inner_pair_xy = ...
    interp1(dists(inner_idx&outer_idx,:), inner_edge_xy(inner_idx&outer_idx,:), dists(outer_idx,:));

%Now median filter the widths to remove any rogue outliers


outer_pair_xy = ...
    interp1(dists(outer_idx,:), outer_edge_xy(outer_idx,:), evenly_spaced_dists);
inner_pair_xy = ...
    interp1(dists(inner_idx&outer_idx,:), inner_edge_xy(inner_idx&outer_idx,:), evenly_spaced_dists);

discard_idx = any(isnan(outer_pair_xy),2) | any(isnan(inner_pair_xy),2);
outer_pair_xy(discard_idx,:) = [];
inner_pair_xy(discard_idx,:) = [];

centre_pair_xy = (outer_pair_xy + inner_pair_xy)/2;
[~,~,paired_apex_pt] = line_cross(centre_pair_xy, apex_xy);

%Get widths at each pair
profile_vector = centre_pair_xy - outer_pair_xy;
profile_widths = sqrt(sum(profile_vector.^2,2));

%Now build a vessel mask by placing a disk of appropriate width at each
%centre point
vessel_mask = false(size(vessel_patch));
for i_p = 1:size(centre_pair_xy,1)
    disk_r = floor(profile_widths(i_p));
    vessel_mask_i = false(size(vessel_patch));
    vessel_mask_i(round(centre_pair_xy(i_p,2)), round(centre_pair_xy(i_p,1))) = 1;
    vessel_mask = vessel_mask | imdilate(vessel_mask_i, strel('disk', disk_r));
end
    


% %Smooth the widths of the outer points
% profile_vector = centre_xy(outer_idx,:) - outer_edge_xy(outer_idx);
% profile_widths = sqrt(sum(profile_vector.^2,2));
% profile_vector_u = bsxfun(@rdivide, profile_vector, profile_widths);
% 
% profile_widths_sm = medfilt1(profile_widths, 5);
% profile_vector_sm = bsxfun(@times, profile_vector_u, profile_widths_sm);
% inner_pair_xy = outer_pair_xy - profile_vector_sm;
% 
% %Smooth the widths of the inner points



% paired_apex_pt = sum(paired_idx(1:apex_pt));
% outer_apex_pt = sum(outer_idx(1:apex_pt));
% inner_apex_pt = sum(inner_idx(1:apex_pt));
% 
% inner_pair_xy = inner_edge_xy(paired_idx,:);
% outer_pair_xy = outer_edge_xy(paired_idx,:);
% centre_pair_xy = centre_xy(paired_idx,:);

%Median filter widths of normal profiles to remove rogue profiles and
%smooth the vessel edge
% profile_vector = outer_pair_xy - inner_pair_xy;
% profile_widths = sqrt(sum(profile_vector.^2,2));
% profile_vector_u = bsxfun(@rdivide, profile_vector, profile_widths);
% 
% profile_widths_sm = medfilt1(profile_widths, 5);
% profile_vector_sm = bsxfun(@times, profile_vector_u, profile_widths_sm);
% inner_pair_xy = outer_pair_xy - profile_vector_sm;

% %Make a BW mask of the vessel using each half of the vessel independently
% vessel_xy1 = [flipud(inner_pair_xy(1:paired_apex_pt,:)); outer_pair_xy(1:paired_apex_pt,:)];
% vessel_xy2 = [flipud(inner_pair_xy(paired_apex_pt:end,:)); outer_pair_xy(paired_apex_pt:end,:)];
% vessel_mask = ...
%     poly2mask(vessel_xy1(:,1), vessel_xy1(:,2),...
%         size(vessel_patch,1), size(vessel_patch,2)) | ...
%     poly2mask(vessel_xy2(:,1), vessel_xy2(:,2),...
%         size(vessel_patch,1), size(vessel_patch,2));

if args.plot 
    subplot(2,3,4); imgray(vessel_patch);
    plot(vessel_pts(:,1), vessel_pts(:,2), 'y-.');   
    
    subplot(2,3,5); imgray(vessel_patch);
    plot(vessel_pts(:,1), vessel_pts(:,2), 'y-.');
    plot([inner_pair_xy(:,1) outer_pair_xy(:,1)]',...
        [inner_pair_xy(:,2) outer_pair_xy(:,2)]');
    plot(centre_pair_xy(:,1), centre_pair_xy(:,2), 'g.');

    plot(outer_edge_xy(:,1), outer_edge_xy(:,2), 'rx');
    plot(inner_edge_xy(:,1), inner_edge_xy(:,2), 'bx');
    
    subplot(2,3,6); imgray(vessel_mask);
end

%--------------------------------------------------------------------------
function [keep_idx] = remove_overlaps(contour1, contour2, neighbourhood)

num_pts = size(contour1,1);
normals_cross = false(num_pts, num_pts);
for i_p = 1: num_pts - neighbourhood
    
    vec1 = [contour1(i_p,:); contour2(i_p,:)];
    for j_p = i_p + (1:neighbourhood)
        
        vec2 = [contour1(j_p,:); contour2(j_p,:)];
        
        normals_cross(i_p, j_p) = line_cross(vec1, vec2);
        normals_cross(j_p, i_p) = normals_cross(i_p, j_p);
    end
end

keep_pts = (1:num_pts)';
while any(normals_cross(:))
    [~, remove_i] = max(sum(normals_cross));
    %display(remove_i);
    
    normals_cross(remove_i,:) = [];
    normals_cross(:,remove_i) = [];
    keep_pts(remove_i) = [];
end

keep_idx = false(num_pts,1);
keep_idx(keep_pts) = 1;
