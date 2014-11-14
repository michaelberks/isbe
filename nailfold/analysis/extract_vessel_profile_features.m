function [profile_features] = extract_vessel_profile_features(...
    g2d_responses, h2d_responses, profile_centre, profile_pos, theta, scale_offsets, displacement)
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

num_scales = length(scale_offsets);
num_profile_pts = length(profile_pos);

%Rotate normal vectors by angle offset and compute final
%sampling positions
normal_xy_i = [cos(theta) -sin(theta)];

% Pre-allocate containers for the normal profile features
normal_g2d = zeros(num_scales, num_profile_pts);
normal_h2d = zeros(num_scales, num_profile_pts);
            
for i_sc = 1:num_scales

    %Apply scale to sampling positions
    scale_offset = scale_offsets(i_sc);
    profile_pos_i = profile_pos*scale_offset + displacement;       

    n_x = profile_centre(1)+normal_xy_i(1)*profile_pos_i;
    n_y = profile_centre(2)+normal_xy_i(2)*profile_pos_i;

    %Sample profiles from each filter response map
    g2d_responses_p = zeros(1,num_profile_pts,1,3);
    for i_g = 1:3
        g2d_responses_p(:,:,:,i_g) = interp2(g2d_responses(:,:,i_sc,i_g), n_x, n_y, '*linear')';
    end
    h2d_responses_p = zeros(1,num_profile_pts,1,4);
    for i_h = 1:4
        h2d_responses_p(:,:,:,i_h) = interp2(h2d_responses(:,:,i_sc,i_h), n_x, n_y, '*linear')';
    end
    normal_g2d(i_sc, :) = steer_gaussian_2nd_derivatives(g2d_responses_p, theta);
    normal_h2d(i_sc, :) = steer_hilbert_2nd_derivatives(h2d_responses_p, theta);

    %%Sample profiles from mask?    
    %mask_p = interp2(initial_vessel_mask, n_x, n_y, 'nearest');
    %mask_p(isnan(mask_p)) = 0;
    %normal_mask(i_pt, i_dis, i_sc, i_th, :) = mask_p;

end
profile_features = [normal_g2d(:)' normal_h2d(:)'];


