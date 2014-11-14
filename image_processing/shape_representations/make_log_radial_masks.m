function [radial_masks radial_idx] = make_log_radial_masks(num_r, min_r, log_r, do_plot)
%MAKE_LOG_RADIAL_MASKS *Insert a one line summary here*
%   [masks] = make_shape_context_mask(num_r, min_r, num_theta)
%
% Inputs:
%      num_r - *Insert description of input variable here*
%
%      min_r - *Insert description of input variable here*
%
%      log_r - *Insert description of input variable here*
%
%
% Outputs:
%      masks - *Insert description of input variable here*
%
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
if ~exist('do_plot', 'var')
    do_plot = 0;
end

%Workout sizes of each scale bin
powers = 0:num_r-1;
scales = log_r.^powers;
rdists = [0 cumsum(scales*min_r)];
max_rdist = ceil(rdists(end));
mask_sz = 2*max_rdist + 1;

%Make dist map
xy = repmat(-max_rdist:max_rdist, mask_sz, 1);
dist_map = sqrt(xy.^2 + xy'.^2);

if do_plot; figure; num_plots = min(3, num_r); end

%Now make masks
radial_masks = false(mask_sz, mask_sz, num_r);
radial_idx = zeros(mask_sz, mask_sz);
for i_r = 1:num_r
        
    radial_masks(:,:,i_r) = ...
        (dist_map > rdists(i_r)) &...
        (dist_map <= rdists(i_r+1));
    
    radial_idx(radial_masks(:,:,i_r)) = i_r;

    if do_plot && (i_r <= num_plots);
        subplot(1, num_plots, p_num);
        imgray(radial_masks(:,:,i_r,i_theta));
        p_num = p_num+1;
    end
end
                    


