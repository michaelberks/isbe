function [discard_map] = remove_circles(image_in, sigma_range, num_angles,...
    mag_thresh, phase_thresh, diff_thresh)
%REMOVE_CIRCLES *Insert a one line summary here*
%   [discard_map] = remove_circles(image_in, sigma_range, mag_thresh, phase_thresh)
%
% Inputs:
%      image_in - *Insert description of input variable here*
%
%      sigma_range - *Insert description of input variable here*
%
%      mag_thresh - *Insert description of input variable here*
%
%      phase_thresh - *Insert description of input variable here*
%
%
% Outputs:
%      discard_map - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 07-Nov-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

if ~exist('sigma_range', 'var') || isempty(sigma_range)
    sigma_range = [2 4];
end
if ~exist('num_angles', 'var') || isempty(num_angles)
    num_angles = 6;
end
if ~exist('mag_thresh', 'var') || isempty(mag_thresh)
    mag_thresh = 1;
end
if ~exist('phase_thresh', 'var') || isempty(phase_thresh)
    phase_thresh = pi/6;
end
if ~exist('diff_thresh', 'var') || isempty(diff_thresh)
    diff_thresh = 0.8;
end

separable_responses = compute_gaussian_2nd_derivatives(image_in, sigma_range);
r_g2d = steer_gaussian_2nd_derivatives(separable_responses, [], num_angles);

separable_responses = compute_hilbert_2nd_derivatives_sep(image_in, sigma_range);
r_h2d = steer_hilbert_2nd_derivatives(separable_responses, [], num_angles);

clear separable_responses;
 
r_mag = sqrt(r_g2d.^2 + r_h2d.^2);
r_phase = atan(r_h2d ./ r_g2d);
clear r_g2d r_h2d;

discard_map = false(size(image_in));

for i_sc = 1:length(sigma_range)
    discard_mag = all(r_mag(:,:,i_sc,:) > mag_thresh, 4);
    discard_phase = all(abs(r_phase(:,:,i_sc,:)) < phase_thresh, 4);
    discard_diff = (min(r_mag(:,:,i_sc,:),[],4) ./ max(r_mag(:,:,i_sc,:),[],4))  > diff_thresh;
    
    %discard_map_i = imdilate(discard_mag & discard_phase, strel('disk', sigma_range(i_sc)));
    discard_map_i = imdilate(discard_mag & discard_diff & discard_phase,...
        strel('disk', sigma_range(i_sc)));
    discard_map = discard_map | discard_map_i;
end


