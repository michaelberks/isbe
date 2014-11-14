function [profile_features displacements width_ratios] = sample_vessel_profile_features_image(vessel_patch, vessel_centre, outer_edge, inner_edge, varargin)
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
    'num_profile_pts', 5,...
    'theta_offsets', pi*[-15 0 15]/180,...
    'scale_offsets', [0.75 1 1.25],...
    'num_displacements', 10,...
    'max_displacement', [],...
    'min_wr', 0.5,...
    'max_wr', 1.5,...
    'plot', 1,...
    'quiet', 0);
 clear varargin;

num_pts = size(vessel_centre,1);
num_angles = length(args.theta_offsets);
num_scales = length(args.scale_offsets);

num_samples = num_pts*args.num_displacements*num_angles;
feature_size = 2*num_scales*args.num_profile_pts;

%Pre-allocate outputs
profile_features = zeros(num_samples, feature_size);
displacements = zeros(num_samples, 1);
width_ratios = zeros(num_samples, 1);

%Compute vessel width_ratios from inner and outer edges
vessel_widths = sqrt(sum((inner_edge-outer_edge).^2,2));

%Compute edge normals from inner and outer edges
normal_xy = bsxfun(@rdivide, inner_edge - outer_edge, vessel_widths);

%Compute Gaussian (and its hilbert pair) responses
%Use the half the mean apex width as the central scale
mean_width = mean(vessel_widths);
[g2d_responses] = compute_gaussian_2nd_derivatives(vessel_patch, mean_width*args.scale_offsets/2);
[h2d_responses] = compute_hilbert_2nd_derivatives_sep(vessel_patch, mean_width*args.scale_offsets/2);

%Also use the mean vessel width to set the max displacement if not supplied
%by user
if isempty(args.max_displacement)
    args.max_displacement = mean_width/4;
end

%--------------------------------------------------------------------------
%Leave out the mask stuff for now
%Make an estimated mask of the vessel based on the initial annotation and
%the appex width
% initial_vessel_mask = false(size(vessel_patch));
% vessel_idx = sub2ind(size(vessel_patch), round(vessel_centre(:,2)), round(vessel_centre(:,1)));
% initial_vessel_mask(vessel_idx) = 1;
% initial_vessel_mask = imdilate(initial_vessel_mask, strel('disk', floor(0.8*apex_width2)));
% normal_mask = false(num_pts, args.num_displacements, num_scales, num_angles, args.num_profile_pts);
%--------------------------------------------------------------------------
    
if args.plot
    figure; imgray(vessel_patch);
end

sample_idx = 1;
for i_pt = 1:num_pts
        
    %Set up profile sample positions
    profile_pos = linspace(-vessel_widths(i_pt)/2, vessel_widths(i_pt)/2, args.num_profile_pts);
    
    %Compute angle from normal vector
    theta = atan2(-normal_xy(i_pt,2), normal_xy(i_pt,1));
    
    for i_dis = 1:args.num_displacements
               
        if i_dis == 1
            displacement = 0;
            %width_ratio = 1;
        else
            displacement = 2*args.max_displacement*rand - args.max_displacement;
            
        end
        width_ratio = args.min_wr + (args.max_wr - args.min_wr)*rand;
        profile_pos = profile_pos*width_ratio;
        
        for theta_offset = args.theta_offsets
            %Store the responses and their associated displacement in the
            %containers
            profile_features(sample_idx,:) = extract_vessel_profile_features(...
                g2d_responses, h2d_responses, vessel_centre(i_pt,:), profile_pos,...
                theta+theta_offset, ones(size(args.scale_offsets)), displacement);
            displacements(sample_idx,:) = displacement;
            width_ratios(sample_idx,:) = width_ratio;
            
            %Increment the sample index
            sample_idx = sample_idx + 1;
        end        
    end           
end


