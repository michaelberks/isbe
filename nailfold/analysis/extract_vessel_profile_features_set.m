function [profile_features displacements width_ratios] = extract_vessel_profile_features_set(contour_dir, vessel_dir, varargin)
%EXTRACT_VESSEL_PROFILE_FEATURES *Insert a one line summary here*
%   [profile_features] = extract_vessel_profile_features(outer_edge, inner_edge, vessel_centre)
%
% Inputs:
%      outer_edge, inner_edge - *Insert description of input variable here*
%
%      vessel_centre - *Insert description of input variable here*
%
%
% Outputs:
%      profile_features - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 17-Apr-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

args = u_packargs(varargin,... % the user's input
    '0', ... % strict mode
    'selected_idx', [],...
    'num_profile_pts', 5,...
    'theta_offsets', pi*[-15 0 15]/180,...
    'scale_offsets', [0.75 1 1.25],...
    'num_displacements', 10,...
    'max_displacement', [],...
    'plot', 0,...
    'quiet', 0);
 clear varargin;
 
v_files = dir([contour_dir '*vessel_contour.mat']);
if isempty(args.selected_idx)
    selected_idx = 1:length(v_files);
else
    selected_idx = args.selected_idx;
end
args = rmfield(args, 'selected_idx');

%Loop through data once to add up total number of points
num_pts = 0;
for i_v = selected_idx
    vessel_struc = load([contour_dir v_files(i_v).name]);
    num_pts = num_pts + size(vessel_struc.outer_edge,1);
end
num_samples = num_pts*args.num_displacements;

num_angles = length(args.theta_offsets);
num_scales = length(args.scale_offsets);
feature_size = 2*num_scales*args.num_profile_pts;

%Pre-allocate outputs
profile_features = zeros(num_samples, feature_size);
displacements = zeros(num_samples, 1);
width_ratios = zeros(num_samples, 1);

%Now loop through each image to extract the features
sample_idx = 0;
for i_v = selected_idx
    
    display(['Extracting features from vessel ' num2str(i_v)]);
    
    %Load the data for this vessel
    vessel_struc = load([contour_dir v_files(i_v).name]);
    apex_struc = u_load([vessel_dir v_files(i_v).name(1:end-12) '.mat']);
     
    %Compute the equalised vessel patch
    vessel_patch = double(apex_struc.vessel_patch);
    g = gaussian_filters_1d(16, 48);
    g = g / sum(g);
    im_edges = conv2(g', g, ones(size(vessel_patch)), 'same');
    vessel_patch_smoothed = conv2(g', g, vessel_patch, 'same') ./ im_edges;
    vessel_patch_equalised = vessel_patch - vessel_patch_smoothed;

    %Sample the features from this image
    num_pts_i = size(vessel_struc.outer_edge,1);
    sample_idx = sample_idx(end)+(1:num_pts_i*args.num_displacements*num_angles);
    
    [profile_features(sample_idx,:) displacements(sample_idx) width_ratios(sample_idx)] = ...
        sample_vessel_profile_features_image(...
        vessel_patch_equalised,...
        vessel_struc.vessel_centre,...
        vessel_struc.outer_edge, vessel_struc.inner_edge, args);
end


