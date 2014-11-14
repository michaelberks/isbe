function [outer_edge_xy, inner_edge_xy, vessel_mask] = detect_vessel_edge_snake(vessel_patch, vessel_pts, apex_xy, varargin)
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
    'search_width_fraction', 0.1,...
    'iterate', true,...
    'plot', 1,...
    'quiet', 0);
 clear varargin;
 
%Constants
sigma0_w = round(5*args.sigma0);
sigma1_w = round(5*args.sigma1);

initial_spacing = 5;
outer_prof_widths = [-20 50];
%centre_prof_width = [-40 40];

neighbourhood = 10;

num_pts = size(vessel_pts,1);

%Display the patch
if args.plot
    figure; imgray(vessel_patch);
    plot(vessel_pts(:,1), vessel_pts(:,2), 'r-');
    plot(vessel_pts(:,1), vessel_pts(:,2), 'r.');
end

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
[v_pts_sm] = spline_contour(v_pts_sm, [], initial_spacing);
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
centre_prof_width = [-round(apex_width) round(apex_width)];
%sym_filt = [ones(1,round(apex_width/2)) -ones(1,round(apex_width/2))] / (2*round(apex_width/2));
[g1d_responses] = compute_gaussian_1st_derivatives(vessel_patch, apex_width/2);
[g2d_responses] = compute_gaussian_2nd_derivatives(vessel_patch, apex_width/2);
[h2d_responses] = compute_hilbert_2nd_derivatives_sep(vessel_patch, apex_width/2);
[g2d_mag, g2d_ori] = gaussian_2nd_derivative_line(g2d_responses);


% Pre-allocate containers for the normal profiles and the associated
% sampling points
norm_width = diff(centre_prof_width)+1;
normal_p = zeros(num_pts, norm_width);
%normal_sym = zeros(num_pts, norm_width);
normal_g2d = zeros(num_pts, norm_width);
normal_h2d = zeros(num_pts, norm_width);
normal_g1d = zeros(num_pts, norm_width);
normal_g2m = zeros(num_pts, norm_width);

normal_x = zeros(num_pts, norm_width);
normal_y = zeros(num_pts, norm_width);

for i_n = 1:num_pts  
    %Get normal profile
    n_x = v_pts_sm(i_n,1)+normal_xy(i_n,1)*centre_prof_width;
    n_y = v_pts_sm(i_n,2)+normal_xy(i_n,2)*centre_prof_width;
    
    theta = atan(-normal_xy(i_n,2)/normal_xy(i_n,1));
    
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
    normal_g1d(i_n, :) = steer_gaussian_1st_derivatives(g1d_responses_p, theta);
    normal_g2d(i_n, :) = steer_gaussian_2nd_derivatives(g2d_responses_p, theta);
    normal_h2d(i_n, :) = steer_hilbert_2nd_derivatives(h2d_responses_p, theta);
    %normal_sym(i_n, :) = conv(cp', sym_filt, 'same');
    normal_g2m(i_n, :) = improfile(g2d_mag, n_x, n_y, norm_width, 'bicubic')';
end

centre_col = ceil(norm_width/2);
edge_mask = false(num_pts, norm_width);
centre_mask = false(num_pts, norm_width);

centre_xy = zeros(num_pts,2);
centre_idx = NaN(num_pts,1);
outer_edge_xy = NaN(num_pts,2);
inner_edge_xy = NaN(num_pts,2);
for i_n = 1:num_pts
    profile_phase = atan(normal_g2d(i_n,:) ./ normal_h2d(i_n,:));
    poss_edges = abs(profile_phase) < 0.2;
    
    
    left_edge = find(poss_edges(1:centre_col-1), 1, 'last');
    if ~isempty(left_edge)
        edge_mask(i_n,1:left_edge) = 1;
    end
    right_edge = find(poss_edges(centre_col+1:end), 1, 'first');
    if ~isempty(right_edge)
        right_edge = right_edge + centre_col;
        edge_mask(i_n,right_edge:end) = 1;
        
        outer_edge_xy(i_n,1) = normal_x(i_n, right_edge);
        outer_edge_xy(i_n,2) = normal_y(i_n, right_edge);

        centre_mask(i_n,:) = (profile_phase > (pi/2 - 0.1));
%         poss_centres = find(...
%             centre_mask(i_n,:) & ~edge_mask(i_n,:));
%         
%         if ~isempty(poss_centres)
%             [~,min_idx] = min(abs(poss_centres - centre_col));
%             centre_idx(i_n) = poss_centres(min_idx);
%             centre_xy(i_n,1) = normal_x(i_n, centre_idx(i_n));
%             centre_xy(i_n,2) = normal_y(i_n, centre_idx(i_n));
%         else
%             centre_xy(i_n,:) = NaN;
%         end
        poss_centre = find(diff(normal_g2d(i_n,right_edge:-1:1)) <= 0, 1, 'first');
        if ~isempty(poss_centre)
            centre_idx(i_n) = right_edge-poss_centre+1;
            centre_xy(i_n,1) = normal_x(i_n, centre_idx(i_n));
            centre_xy(i_n,2) = normal_y(i_n, centre_idx(i_n));
            
            inner_edge_xy(i_n,:) = 2*centre_xy(i_n,:) + outer_edge_xy(i_n,:);
        end                   
    end
end

normal_score = normal_g2d;
normal_score(edge_mask) = -inf;

%
initial_edge = [ones(num_pts,1)*centre_col (1:num_pts)'];
search_width = norm_width;

%Run first iteration of the snake
[snake_edge,snake_energy] = mb_snake_normal(initial_edge, ...
    args.alpha, args.beta, search_width, args.resolution, normal_score, normal_x, normal_y);

centre_xy = zeros(num_pts,2);
for i_n = 1:num_pts
    centre_xy(i_n,1) = normal_x(i_n, snake_edge(i_n,1));
    centre_xy(i_n,2) = normal_y(i_n, snake_edge(i_n,1));
end

return;
norm_width = diff(outer_prof_widths)+1;
%idx_outer = (-outer_prof_widths(1)+1):norm_width;
idx_outer = 1:norm_width;
outer_width = length(idx_outer);

normal_p = zeros(num_pts, norm_width);
normal_sym = zeros(num_pts, outer_width);
normal_g2d = zeros(num_pts, outer_width);

normal_x = zeros(num_pts, outer_width);
normal_y = zeros(num_pts, outer_width);


apex_width = sqrt(sum(diff(apex_xy).^2));
sym_filt = [ones(1,round(apex_width/2)) -ones(1,round(apex_width/2))] / (2*round(apex_width/2));
[g2d_responses] = compute_gaussian_2nd_derivatives(vessel_patch, apex_width/2);

for i_n = 1:num_pts  
    %Get normal profile
    n_x = v_pts_sm(i_n,1)+normal_xy(i_n,1)*outer_prof_widths;
    n_y = v_pts_sm(i_n,2)+normal_xy(i_n,2)*outer_prof_widths;
    
    theta = atan(-normal_xy(i_n,2)/normal_xy(i_n,1));
    cc = cos(theta).^2;
    ss = sin(theta).^2;
    s2 = sin(2*theta);
    
    Ixy = improfile(g2d_responses(:,:,1,1), n_x, n_y, norm_width, 'bicubic');
    Ixx = improfile(g2d_responses(:,:,1,2), n_x, n_y, norm_width, 'bicubic');
    Iyy = improfile(g2d_responses(:,:,1,3), n_x, n_y, norm_width, 'bicubic');
        
    [cx, cy, cp] = improfile(vessel_patch, n_x, n_y, norm_width, 'bilinear');
    normal_p(i_n, :) = cp';
    normal_x(i_n, :) = cx(idx_outer)';
    normal_y(i_n, :) = cy(idx_outer)';
    normal_g2d(i_n, :) = (Ixx.*cc + Iyy.*ss + Ixy.*s2)';
    normal_sym(i_n, :) = conv(cp', sym_filt, 'same');
end

apex_xy = sortrows(apex_xy,2);
wide_apex_xy = bsxfun(@plus, 2*subtract_mean(apex_xy), mean(apex_xy));
apex_profile = improfile(vessel_patch, wide_apex_xy(:,1), wide_apex_xy(:,2), 'bilinear');
half_width = floor(length(apex_profile)/2);

norm_edge = normxcorr2(apex_profile', normal_p);
norm_edge = norm_edge(:,half_width+idx_outer);

%Take derivatives of profile to look for edge the profiles
%[~, dg] = gaussian_filters_1d(args.sigma1, sigma1_w);
%norm_edge = imfilter(normal_p, -dg, 'replicate');
%norm_edge = norm_edge(:,idx_outer);

%Compute the maximal edge strength on each row and use as the intial points
%for the snake
initial_edge = [zeros(num_pts,1) (1:num_pts)'];
[snake_edge,snake_energy] = mb_snake_normal(initial_edge, ...
    args.alpha, args.beta, search_width, args.resolution, norm_edge, normal_x, normal_y);

% 
% for i_n = 1:num_pts
%     %Get initial edge as first maxima
%     max_outer = find(diff(norm_edge(i_n,:)) <= 0, 1, 'first');       
%     initial_edge(i_n,1) = max_outer;
% end
[~, initial_edge(:,1)] = max(norm_edge,[],2);

search_width = round(args.search_width_fraction*norm_width);

%Run first iteration of the snake
[snake_edge,snake_energy] = mb_snake_normal(initial_edge, ...
    args.alpha, args.beta, search_width, args.resolution, norm_edge, normal_x, normal_y);

%Repeat iterations of the snake until the overall energy increases again
%(note this may find a local minima)
if args.iterate
    go_on = true;
    i_iter = 2;
    while go_on
        [snake_edge,snake_energy(i_iter)] = mb_snake_normal(snake_edge, ...
            args.alpha, args.beta, search_width, args.resolution, norm_edge, normal_x, normal_y);
        go_on = snake_energy(i_iter) < snake_energy(i_iter-1);
        if ~args.quiet
            display(['Iteration: ', num2str(i_iter)]);
        end
        i_iter = i_iter+1;
    end
end

outer_edge_xy = zeros(num_pts,2);
for i_n = 1:num_pts
    outer_edge_xy(i_n,1) = normal_x(i_n, snake_edge(i_n,1));
    outer_edge_xy(i_n,2) = normal_y(i_n, snake_edge(i_n,1));
end
if args.plot
    plot(outer_edge_xy(:,1), outer_edge_xy(:,2), 'g-');
    plot(outer_edge_xy(:,1), outer_edge_xy(:,2), 'g.');
end

return;

%Remove any overlapping points
keep_idx = remove_overlaps(outer_edge_xy, v_pts_sm, neighbourhood);
outer_edge_xy = outer_edge_xy(keep_idx,:);
num_pts = size(outer_edge_xy,1);

%Smooth outer edge
outer_edge_xy = [...
    repmat(outer_edge_xy(1,:), sigma0_w, 1);...
    outer_edge_xy;...
    repmat(outer_edge_xy(end,:), sigma0_w, 1)];
outer_edge_xy = medfilt1(outer_edge_xy, neighbourhood, [], 1);
outer_edge_xy = [conv(outer_edge_xy(:,1), g, 'same') conv(outer_edge_xy(:,2), g, 'same')];
outer_edge_xy = outer_edge_xy(sigma0_w+(1:num_pts),:);
if args.plot
    plot(outer_edge_xy(:,1), outer_edge_xy(:,2), 'b-');
end
    
%Compute normals from outer edge - check again if we need to reverse
normal_xy = compute_spline_normals(outer_edge_xy);
if normal_xy(apex_pt,2) > 0
    normal_xy = -normal_xy;
end

%Find inner edge points from normal profiles
inner_edge_xy = zeros(num_pts,2);
idx_inner = (inner_prof_width+2):(2*inner_prof_width+1);
for i_n = 1:num_pts
    
    %Get normal profile
    [norm_xp norm_yp norm_profile] = improfile(vessel_patch,...
        v_pts_sm(i_n,1)+normal_xy(i_n,1)*inner_prof_width*[1 -1],...
        v_pts_sm(i_n,2)+normal_xy(i_n,2)*inner_prof_width*[1 -1],...
        2*inner_prof_width+1, 'bilinear');
    
    %Compute derivative of profile
    [~, dg] = gaussian_filters_1d(args.sigma1, sigma1_w);
    normal_g1 = conv(norm_profile, dg, 'same');

    %Get edge as first maxima
    max_inner = find(diff(normal_g1(idx_inner)) <= 0, 1, 'first');
       
    x_edge_inner = norm_xp(idx_inner(max_inner));
    y_edge_inner = norm_yp(idx_inner(max_inner));
    inner_edge_xy(i_n,:) = [x_edge_inner y_edge_inner];
    
    if args.plot
        plot(x_edge_inner, y_edge_inner, 'yx');
    end
end

%Check no points overlap
keep_idx = remove_overlaps(outer_edge_xy, inner_edge_xy, neighbourhood);
inner_edge_xy = inner_edge_xy(keep_idx,:);
num_pts = size(inner_edge_xy,1);

%Smooth outer edge
inner_edge_xy = [...
    repmat(inner_edge_xy(1,:), sigma0_w, 1);...
    inner_edge_xy;...
    repmat(inner_edge_xy(end,:), sigma0_w, 1)];
inner_edge_xy = medfilt1(inner_edge_xy, neighbourhood, [], 1);
inner_edge_xy = [conv(inner_edge_xy(:,1), g, 'same') conv(inner_edge_xy(:,2), g, 'same')];
inner_edge_xy = inner_edge_xy(sigma0_w+(1:num_pts),:);

if args.plot
    plot(inner_edge_xy(:,1), inner_edge_xy(:,2), 'm-');
end

%Finally fit a spline curve to both inner an outer points at equal
%distances
[outer_edge_xy] = spline_contour(outer_edge_xy);
[inner_edge_xy] = spline_contour(inner_edge_xy, size(outer_edge_xy,1));

total_vessel_xy = [outer_edge_xy; flipud(inner_edge_xy)];
vessel_mask = poly2mask(total_vessel_xy(:,1), total_vessel_xy(:,2),...
    size(vessel_patch,1), size(vessel_patch,2));



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

keep_idx = (1:num_pts)';
while any(normals_cross(:))
    [~, remove_i] = max(sum(normals_cross));
    %display(remove_i);
    
    normals_cross(remove_i,:) = [];
    normals_cross(:,remove_i) = [];
    keep_idx(remove_i) = [];
end
