function [apex_prior_scores path_map] =...
    apex_prior_prob_track_mult(prob_map, orientation_map, r_idx, c_idx, varargin)
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
    'apex_prior', [],...
    'max_steps', 50,...
    'num_streams', 1e3,...
    'step_length', 1,...
    'stopping_prob', 0.2);
clear varargin;

num_pts = length(r_idx);
apex_prior_scores = zeros(args.max_steps, 2, num_pts);

%Set up coords for sampling
mask_sz = 2*args.max_steps + 1;
patch_rows = 1:mask_sz;
patch_cols = 1:mask_sz;
image_rows = -args.max_steps:args.max_steps;
image_cols = -args.max_steps:args.max_steps;
[rows cols] = size(prob_map);

ori_theta = angle(orientation_map) / 2;
ori_D = abs(orientation_map);
%ori_D = prob_map;
ori_sigma = sqrt(-2*log(ori_D));
ori_sigma(~ori_D) = 1e3;

if nargout > 1
    path_map = zeros(rows, cols);
end

for i_pt = 1:num_pts
    r_i = r_idx(i_pt);
    c_i = c_idx(i_pt);
    
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
    
    ori_theta_patch = zeros(mask_sz);
    ori_theta_patch(patch_rows_i, patch_cols_i) = ori_theta(image_rows_i, image_cols_i);
    
    ori_sigma_patch = zeros(mask_sz);
    ori_sigma_patch(patch_rows_i, patch_cols_i) = ori_sigma(image_rows_i, image_cols_i);
    
    prob_patch = zeros(mask_sz);
    prob_patch(patch_rows_i, patch_cols_i) = prob_map(image_rows_i, image_cols_i);
    
    %Now track in this patch
        
    [apex_prior_score path_map_i] = ...
        apex_prior_prob_track(...
            prob_patch, ori_theta_patch, ori_sigma_patch, args.max_steps + 1, args.max_steps + 1,...
            args.step_length, args.num_streams, ...
            args.max_steps, args.stopping_prob, args.apex_prior);
        
    if exist('path_map', 'var')
        path_map(image_rows_i, image_cols_i) = path_map(image_rows_i, image_cols_i) +...
            path_map_i(patch_rows_i, patch_cols_i);
    end
    
    %Save the sparse mat back into the main container and derotate all the
    %angles by theta0
    apex_prior_scores(:,:,i_pt) = apex_prior_score;
end
apex_prior_scores = apex_prior_scores / args.num_streams;