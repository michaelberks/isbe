function [predicted_width_ratios profile_features] = ...
    predict_apex_width_image(vessel_patch, vessel_centre, apex_width, rf_predictor, varargin)
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
    'theta', [],...
    'num_profile_pts', 5,...
    'scale_offsets', [0.75 1 1.25],...
    'sigma_vote', 2,...
    'plot', 0,...
    'debug', 0,...
    'quiet', 0);
 clear varargin;

if args.plot
    figure; imgray(vessel_patch);
end
    
num_pts = size(vessel_centre,1);
num_scales = length(args.scale_offsets);
feature_size = 2*num_scales*args.num_profile_pts;

%Compute edge normals from inner and outer edges
%Project normals out from path, and use profiles to detect outer edge
if isempty(args.theta)
    if size(vessel_centre,3)
        normal_xy = compute_spline_normals(vessel_centre);
        theta = atan2(-normal_xy(:,2), normal_xy(:,1));
    else
        error('If a normal direction is supplied, vessel centre must contain at least 3 points to compute a normal direction');
    end
else
    theta = args.theta;
end

%Compute Gaussian (and its hilbert pair) responses
%Use the half the mean apex width as the central scale
sigmas = mean(apex_width)*args.scale_offsets/2;
[g2d_responses] = compute_gaussian_2nd_derivatives(vessel_patch, sigmas);
[h2d_responses] = compute_hilbert_2nd_derivatives_sep(vessel_patch, sigmas);
  
profile_features = zeros(num_pts, feature_size);
    
for i_pt = 1:num_pts
        
    %Set up profile sample positions
    profile_pos = linspace(-apex_width(i_pt)/1.5, apex_width(i_pt)/1.5, args.num_profile_pts);
    
    profile_features(i_pt,:) = extract_vessel_profile_features(...
        g2d_responses, h2d_responses, vessel_centre(i_pt,:), profile_pos,...
        theta(i_pt), ones(size(args.scale_offsets)), 0);

end
predicted_width_ratios = random_forest_reg_predict(rf_predictor, profile_features, 0);    
    
if args.plot
    %plot(predicted_centres(:,1), predicted_centres(:,2), 'x');
end
    


