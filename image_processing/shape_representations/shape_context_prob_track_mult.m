function [shape_contexts path_contexts path_map] =...
    shape_context_prob_track_mult(feature_map, prob_map, orientation_map, r_idx, c_idx, varargin)
%SHAPE_CONTEXT *Insert a one line summary here*
%   [] = shape_context()
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
% Created: 07-Jun-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

args = u_packargs(varargin, 0,...
    'weighting_map', [],...
    'prior_left', [],...
    'prior_right', [],...
    'num_rdist', 5,...
    'num_theta', 12, ...
    'log_r', 1.5,...
    'min_rdist', 10,...
    'num_streams', 1e3,...
    'step_length', 1,...
    'stopping_prob', 0.2);
clear varargin;

%Apply weighting to features if set
if ~isempty(args.weighting_map)    
    feature_map = feature_map.*args.weighting_map;
end

num_pts = length(r_idx);
shape_contexts = zeros(args.num_rdist, args.num_theta, num_pts);
path_contexts = zeros(args.num_rdist, args.num_theta, num_pts);

%Make radial masks
radial_masks = make_log_radial_masks(args.num_rdist, args.min_rdist, args.log_r, 0);
mask_sz = size(radial_masks,1);
mask_sz2 = floor(mask_sz/2);
patch_rows = 1:mask_sz;
patch_cols = 1:mask_sz;
image_rows = -mask_sz2:mask_sz2;
image_cols = -mask_sz2:mask_sz2;
[rows cols] = size(feature_map);

%Make map of radial indices
radial_idx = zeros(mask_sz, mask_sz);
for i_rdist = 1:args.num_rdist
    radial_idx(radial_masks(:,:,i_rdist)) = i_rdist;
end

%Make map of angles in mask
max_rdist = floor(mask_sz/2);
xy = repmat(-max_rdist:max_rdist, mask_sz, 1);
theta_map = atan2(-xy', xy);

ori_theta = angle(orientation_map) / 2;
ori_D = abs(orientation_map);
%ori_D = prob_map;
ori_sigma = sqrt(-2*log(ori_D));

if nargout > 2
    path_map = zeros(rows, cols);
end

for i_pt = 1:num_pts
    r_i = r_idx(i_pt);
    c_i = c_idx(i_pt);
    
    %%Pre-allocate shape context for this point as a sparse matrix
    %shape_context_i = sparse(args.num_rdist, args.num_theta);
    
    %Make shape context masks oriented to point's direction
    theta0 = ori_theta(r_i, c_i);
    theta_idx = ceil( .5*args.num_theta*mod(theta_map-theta0, 2*pi)/pi);
    
    %Get feature and selection patch about local point
    image_rows_i = r_i+image_rows;
    image_cols_i = c_i+image_cols;
    patch_rows_i = patch_rows;
    patch_cols_i = patch_cols;
    
    discard_rows = (image_rows_i < 1) | (image_rows_i > rows);
    discard_cols = (image_cols_i < 1) | (image_cols_i > cols);
    image_rows_i(discard_rows) = [];
    image_cols_i(discard_cols) = [];
    patch_rows_i(discard_rows) = [];
    patch_cols_i(discard_cols) = [];
    
    feature_patch = zeros(mask_sz);
    feature_patch(patch_rows_i, patch_cols_i) = feature_map(image_rows_i, image_cols_i);
    
    ori_theta_patch = zeros(mask_sz);
    ori_theta_patch(patch_rows_i, patch_cols_i) = ori_theta(image_rows_i, image_cols_i);
    
    ori_sigma_patch = zeros(mask_sz);
    ori_sigma_patch(patch_rows_i, patch_cols_i) = ori_sigma(image_rows_i, image_cols_i);
    
    prob_patch = zeros(mask_sz);
    prob_patch(patch_rows_i, patch_cols_i) = prob_map(image_rows_i, image_cols_i);
    
    %Now track in this patch
    [shape_context_i path_context_i path_map_i] = ...  %
        shape_context_prob_track(...
            prob_patch, ori_theta_patch, ori_sigma_patch, feature_patch, max_rdist+1, max_rdist+1, ...
            args.step_length, args.num_streams, ...
            mask_sz, radial_idx, theta_idx, args.num_rdist, args.num_theta, args.stopping_prob, args.prior_left, args.prior_right);
        
    if exist('path_map', 'var')
        path_map(image_rows_i, image_cols_i) = path_map(image_rows_i, image_cols_i) +...
            path_map_i(patch_rows_i, patch_cols_i);
    end
    
    %Save the sparse mat back into the main container and derotate all the
    %angles by theta0
    shape_contexts(:,:,i_pt) = shape_context_i .* exp(-2i*theta0);
    path_contexts(:,:,i_pt) = path_context_i;
end