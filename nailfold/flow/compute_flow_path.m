function [vessel_centre particles_x particles_y particles_p particles_d hit_sink smooth_ves] =...
    compute_flow_path(varargin)
%COMPUTE_FLOW_PATH *Insert a one line summary here*
%   [] = compute_flow_path(varargin)
%
% COMPUTE_FLOW_PATH uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Optional Arguments:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 26-Feb-2016
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, '0', ...
    'vessel_name', [],...
    'flow_data_dir',            'N:\Nailfold Capillaroscopy\wellcome\flow_data\',...
    'flow_results_dir',         'N:\Nailfold Capillaroscopy\wellcome\flow_results\',...
    'model_root',               'C:\isbe\nailfold\models\',...
    'rf_ori',                   'vessel\orientation\rf_regression\296621\predictor.mat',...
    'rf_vessels',               'vessel\detection\rf_classification\296655\predictor.mat',...
    'rf_width',                 'vessel\width\rf_regression\297037\predictor.mat',...
    'rf_vessels_args',          'vessel\detection\rf_classification\296655\job_args.mat',...
    'vessel_mask_thresh',       0.25,...
    'vessel_flow',              [],...
    'vessel_predictions',       [],...
    'apex_xy',                  [],...
    'patch_contrast_scale',     70,...
    'vessels_sigma',            2,...
    'refind_apex',              1,...
    'hog_class_forest_path',    'apex\classification\set12g_half_296655\rf.mat',...
    'hog_off_x_forest_path',    'apex\offset_x\set12g_half_296655\rf.mat',...
    'hog_off_y_forest_path',    'apex\offset_y\set12g_half_296655\rf.mat',...
    'separate_trees',           0,...
    'apex_class_thresh',        0.1,...
    'base_width',               20,...
    'num_cells',                8,...
    'cell_sz',                  8,... %Size of HoG cells in blocks
    'block_sz',                 [2 2],...%Size of blocks in cells
    'num_ori_bins',             9,... %Number of bins in orientation histograms
    'norm_method',              'none',... %Method for local normalisation 'none'
    'block_spacing',            8,... %Separation between blocks in pixels - controls level of overlap
    'gradient_operator',        [-1 0 1],...
    'spatial_sigma',            0, ...
    'angle_wrap',               1,...
    'd', 1,...
    'M', 200,...
    'normal_max_step_length', 4,...
    'junction_max_step_length', 25,...
    'junction_angle_spread',    2*pi/3,...
    'N', 400, ...
    'double_angle', 0,...
    'use_angle_weights', 1,...
    'plot', 0);
clear varargin;

%Load frames and flow predictions
if isempty(args.vessel_flow)
    load([args.flow_results_dir args.vessel_name], 'flow_results');
    patch_flow = flow_results.flowPyramidEst{1}; 
else
    patch_flow = args.vessel_flow;
end
[rows cols] = size(patch_flow);
[xx yy] = meshgrid(1:cols, 1:rows);

if isempty(args.apex_xy)
    load([args.flow_results_dir args.vessel_name], 'apex_xy');
    apex_xy = round(apex_xy); %#ok
else
    apex_xy = args.apex_xy;
end


if isempty(args.vessel_predictions)
    
    %Re-predict vessel probs
    load([args.flow_data_dir args.vessel_name],...
        'cropped_frames');
    
    %Load forests
    rfs = cell(3,1);
    rfs{1} = u_load([args.model_root args.rf_vessels]);
    rfs{2} = u_load([args.model_root args.rf_ori]);
    rfs{3} = u_load([args.model_root args.rf_width]);
    rf_args = u_load([args.model_root args.rf_vessels_args]);

    patch = args.patch_contrast_scale*mean(cropped_frames,3);
    vessel_predictions = predict_image(... % non-strict mode
        'image_in', patch,...
        'decomposition_args', rf_args.decomposition_args,...
        'predictor', rfs, ...
        'prediction_type', {'rf_classification', 'rf_regression', 'rf_regression'},...
        'output_type', {'detection', 'orientation', 'width'},...
        'use_probs', 0,...
        'mask', [],...
        'tree_mask', [], ...
        'num_trees', [], ...
        'max_size', 1024,...
        'incremental_results', 0);   
    
    patch_ves = vessel_predictions(:,:,1);
    patch_ori = vessel_predictions(:,:,2);
    patch_wid = vessel_predictions(:,:,3);
else
    patch_ves = args.vessel_predictions(:,:,1);
    patch_ori = args.vessel_predictions(:,:,2);
    patch_wid = args.vessel_predictions(:,:,3);
end   

patch_ves = mb_pad(patch_ves, [rows cols] - size(patch_ves), 'replicate', 'post');
patch_ori = mb_pad(patch_ori, [rows cols] - size(patch_ori), 'replicate', 'post');

%Smooth to get a continous ridge at vessel centre
g_prob = gaussian_filters_1d(args.vessels_sigma);
g_prob = g_prob / sum(g_prob);
smooth_ves = conv2(g_prob', g_prob, patch_ves, 'same');
smooth_flow = conv2(g_prob', g_prob, patch_flow, 'same');
smooth_ori = conv2(g_prob', g_prob, patch_ori, 'same');
smooth_wid = conv2(g_prob', g_prob, patch_wid, 'same');

%make a vessel mask connected to the apex
vessel_mask = smooth_ves > args.vessel_mask_thresh;
if vessel_mask(apex_xy(2), apex_xy(1))
    vessel_mask = bwselect(vessel_mask, apex_xy(1), apex_xy(2));
end

%Get the connected centre mask by thinning to infinity, then thicken by 1
%dilation
thin_mask = bwmorph(vessel_mask, 'thin', 'inf');
medium_mask = imdilate(thin_mask, strel('disk', 2));

%Get coordinates for thin and medium mask
[cy cx] = find(thin_mask);
[vy vx] = find(vessel_mask);
[my mx] = find(medium_mask);

%Re-predict apex location
if args.refind_apex
    apex_class_rf = u_load([args.model_root args.hog_class_forest_path]);
    apex_offset_x_rf = u_load([args.model_root args.hog_off_x_forest_path]);
    apex_offset_y_rf = u_load([args.model_root args.hog_off_y_forest_path]);

    % Get HoG args from main args
    hog_args.cell_sz = [args.cell_sz args.cell_sz];
    hog_args.block_sz = args.block_sz;
    hog_args.num_ori_bins = args.num_ori_bins;
    hog_args.norm_method = args.norm_method;
    hog_args.block_spacing = args.block_spacing;
    hog_args.gradient_operator = args.gradient_operator;
    hog_args.spatial_sigma = args.spatial_sigma;
    hog_args.angle_wrap = args.angle_wrap;

    % Get patch size and form template x,y coordinates for the patch
    patch_sz = args.num_cells*args.cell_sz;
    patch_sz = patch_sz + 2; %Account for padding
    patch_sz2 = (patch_sz - 1)/2;

    % Get hog size from the hog_args struct
    % Set up x,y coordinates for template patch
    hx = repmat(-patch_sz2:patch_sz2, patch_sz, 1);
    hy = hx';
    hxy = [hx(:) hy(:)];
    
    vessel_centre.x = cx;
    vessel_centre.y = cy;
    vessel_centre.ori = smooth_ori(thin_mask);
    vessel_centre.width = smooth_wid(thin_mask);

    [apex_offset_map] = ...
        predict_apex_offsets(...
            'apex_class_rf', apex_class_rf,...
            'apex_offset_x_rf', apex_offset_x_rf,...
            'apex_offset_y_rf', apex_offset_y_rf,...
            'vessel_feature_im', smooth_ves, ...
            'vessel_centre', vessel_centre, ...
            'separate_trees', args.separate_trees,...
            'smoothing_sigma', 0,...
            'num_cells', args.num_cells,...
            'hog_args', hog_args,...
            'xy', hxy,...
            'apex_class_thresh', args.apex_class_thresh,...
            'base_width', args.base_width,...
            'include_pts', []);

    [~, max_idx] = max(apex_offset_map(:));
    [maxr maxc] = ind2sub([rows cols], max_idx);
    apex_xy = [maxc maxr];
    
    %Update the vessel mask given the new apex
    vessel_mask_new = bwselect(vessel_mask, apex_xy(1), apex_xy(2));
    if any(vessel_mask_new(:) ~= vessel_mask(:))
        vessel_mask = vessel_mask_new;
        thin_mask = thin_mask & vessel_mask;
        medium_mask = medium_mask & vessel_mask;

        [cy cx] = find(thin_mask);
        [vy vx] = find(vessel_mask);
        [my mx] = find(medium_mask);
    end    
end

%Correct the predicted orientation given flow direction
patch_theta = angle(smooth_ori)/2;
patch_theta_c = complex(cos(patch_theta), sin(patch_theta));
angle_diff = smooth_flow .* patch_theta_c;
swap_mask = abs(angle(angle_diff)) > pi/2;
smooth_flow = abs(smooth_ori) .* conj(patch_theta_c);
smooth_flow(swap_mask) = smooth_flow(swap_mask) .* exp(1i*pi);

%Find point on thin maks closest to apex and set this region by 0 distance
dists = (cx - apex_xy(1)).^2 + (cy - apex_xy(:,2)).^2;
[~, min_i] = min(dists);
x1 = cx(min_i);
y1 = cy(min_i);
apex_mask = false(rows, cols);
apex_mask(y1, x1) = 1;
apex_mask = imdilate(apex_mask, strel('disk', 2));

%Make a distance map for each point on the medium mask to the vessel apex
max_dist = sum(thin_mask(:));
dist_mask_0 = false(rows, cols);
dist_mask_0(apex_mask) = 1;
dist_map = inf(rows, cols);
dist_map(apex_mask) = 0;
for i_d = 1:max_dist
    dist_mask_1 = imdilate(dist_mask_0, strel('disk', 1)) & medium_mask;
    dist_map(dist_mask_1 & ~dist_mask_0) = i_d;
    if all(dist_mask_0(:) == dist_mask_1(:))
        break;
    else
        dist_mask_0 = dist_mask_1;
    end
end

%Now extrapolate the distance from the medium mask to all vessel pixels
%and smooth
dist_map(vessel_mask) = griddata(cx, cy, dist_map(thin_mask),...
    vx, vy, 'nearest');
dist_map(isinf(dist_map)) = max(dist_map(~isinf(dist_map))) + 1;
dist_map = conv2(g_prob', g_prob, dist_map, 'same');
figure; imgray(dist_map); colorbar;

%Next work out whether moving in the direction of flow increases or
%decreases on the distance map - points moving toward the apex should
%decrease, points away from the apex should increase - more robust to do
%this only for the medium mask points again, tidy this map to remove
%outliers, then extrapolate to full vessel
abs_f = abs(smooth_flow(medium_mask));
myf = my + 2*imag(smooth_flow(medium_mask))./abs_f;
mxf = mx + 2*real(smooth_flow(medium_mask))./abs_f;

m_dist0 = dist_map(medium_mask);
m_dist1 = interp2(dist_map, mxf, myf, '*linear*');
valid = ~isinf(m_dist1);

m_directions = sign(m_dist1(valid) - m_dist0(valid));

%Tidy up small regions in the sign map
sign_map = zeros(rows,cols);
sign_map(medium_mask) = m_directions;
s = bwconncomp(sign_map == -1);
for i_c = 1:length(s.PixelIdxList)
    if length(s.PixelIdxList{i_c}) < 100
        sign_map(s.PixelIdxList{i_c}) = 1;
    end
end
s = bwconncomp(sign_map == 1);
for i_c = 1:length(s.PixelIdxList)
    if length(s.PixelIdxList{i_c}) < 100
        sign_map(s.PixelIdxList{i_c}) = -1;
    end
end

%Extrapolate to full vessel
m_directions = sign_map(medium_mask);
v_directions = griddata(mx(valid), my(valid), m_directions(valid),...
    vx, vy, 'nearest');

%Subtract the mmin value (i.e. the largest negative distance) so that the
%vessel origin (rather than the apex) has zero distance, set all non-vessel
%pixels to inf
dist_map(vessel_mask) = dist_map(vessel_mask) .* v_directions;
%dist_map = dist_map - min(dist_map(medium_mask)); - don't bother with this
dist_map(~vessel_mask) = inf;
vessel_mask = ~isinf(dist_map);
medium_mask = medium_mask & vessel_mask;
thin_mask = thin_mask & vessel_mask;
figure; imgray(dist_map); colorbar;

%
% Locate any junctions and make junction mask - this can be conservatively
% large as the vessel path should still track finewhere the junction mask
% overlaps normal unidirectional vessels
[junction_mask] =...
    get_junction_mask(vessel_mask, thin_mask,'junction_detection');
    
figure; imgray(junction_mask);

for direction = 1:2
    
    %Start at the apex, and define a sink at the end of the limb - for the
    %forward direction this is the point with highest distance
    if direction == 1
        [ysink xsink] = find(dist_map > (max(dist_map(medium_mask))-10) & medium_mask);
        
    else
        %For the reverse direction this is the point with smallest distance
        [ysink xsink] = find(dist_map < (min(dist_map(medium_mask))+10) & medium_mask);
        
        %Also need to reverse distances and flow map
        dist_map = -dist_map;
        smooth_flow = -smooth_flow;
    end

    vessel_sink = false(rows, cols);
    for i_pt = 1:length(xsink)
        vessel_sink = vessel_sink | ...
            (xx-xsink(i_pt)).^2 + (yy-ysink(i_pt)).^2 < 25;
    end
    vessel_sink = vessel_sink & vessel_mask;
    figure; imgray(vessel_sink);
    plot(cx, cy, 'rx');
    plot(x1, y1, 'go');

    %
    % Apply flow jetstream algorithm to get candidate path
    [particles_xi particles_yi particles_pi particles_di hit_sink] =...
        jetstream_flow(smooth_ves, smooth_flow, dist_map, vessel_sink, [x1 y1], ...
        'd', 1,...
        'M', args.M,...
        'normal_max_step_length', args.normal_max_step_length,...
        'junction_max_step_length', args.junction_max_step_length,...
        'junction_angle_spread', args.junction_angle_spread,...
        'use_angle_weights', args.use_angle_weights,...
        'I_jun', junction_mask,...
        'N', args.N, ...
        'double_angle', 0,...
        'plot', 0);

    if any(hit_sink)
        particles_xi(~hit_sink,end) = mean(particles_xi(hit_sink,end));
        particles_yi(~hit_sink,end) = mean(particles_yi(hit_sink,end));
    end
    
    if direction == 1
        particles_x = particles_xi;
        particles_y = particles_yi;
        particles_p = particles_pi;
        particles_d = particles_di;
    else
        particles_x = [fliplr(particles_xi) particles_x]; %#ok
        particles_y = [fliplr(particles_yi) particles_y]; %#ok
        particles_p = [fliplr(particles_pi) particles_p]; %#ok
        particles_d = [fliplr(particles_di) particles_d]; %#ok
    end
        
end

%Take robust average of particles as intial path
particle_path = medfilt1([median(particles_x)' median(particles_y)']);
particle_path(1,:) = [mean(particles_x(:,1)) mean(particles_y(:,1))];
particle_path(end,:) = [mean(particles_x(:,end)) mean(particles_y(:,end))];

figure; imgray(smooth_ves);
plot(particles_x, particles_y, '.');
plot(particle_path(:,1), particle_path(:,2), 'r', 'linewidth', 2);

%Now pick the best centre line given the particles

%Discard any duplicates from the path then upsample coordinates to be pixel
%spaced, and downsample to have a coarser (but evenly spaced) path
duplicates = [false; ~any(diff(particle_path),2)];
particle_path(duplicates,:) = [];
[particle_path_hi] = spline_contour(particle_path, [], 1);
[particle_path_lo] = spline_contour(particle_path_hi(1:5:end,:), [], 5);

%Move points perpindicular to the vessel path to get a smooth contour that
%locks onto the vessel centre line ridge
num_pts = size(particle_path_lo,1);
normal_xy = compute_spline_normals(particle_path_lo);

outer_prof_widths = [-5 5];
norm_width = diff(outer_prof_widths)+1;

normal_p = zeros(num_pts, norm_width);
normal_x = zeros(num_pts, norm_width);
normal_y = zeros(num_pts, norm_width);

for i_n = 1:num_pts  
    %Get normal profile
    n_x = particle_path_lo(i_n,1)+normal_xy(i_n,1)*outer_prof_widths;
    n_y = particle_path_lo(i_n,2)+normal_xy(i_n,2)*outer_prof_widths;
     
    [cx, cy, cp] = improfile(smooth_ves, n_x, n_y, norm_width, 'bilinear');
    normal_p(i_n, :) = cp;
    normal_x(i_n, :) = cx';
    normal_y(i_n, :) = cy';

end

initial_edge = [zeros(num_pts,1) (1:num_pts)'];
[snake_edge] = mb_snake_normal(initial_edge, ...
    0.01, 0.01, norm_width, 1, normal_p, normal_x, normal_y);

vessel_centre = zeros(num_pts,2);
for i_n = 1:num_pts
    vessel_centre(i_n,1) = normal_x(i_n, snake_edge(i_n,1));
    vessel_centre(i_n,2) = normal_y(i_n, snake_edge(i_n,1));
end

figure; imgray(smooth_ves);
plot(vessel_centre(:,1), vessel_centre(:,2));
