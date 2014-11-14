function [sc_masks] = make_shape_context_mask(num_r, min_r, num_theta, theta0, log_r, do_plot)
%MAKE_SHAPE_CONTEXT_MASK *Insert a one line summary here*
%   [masks] = make_shape_context_mask(num_r, min_r, num_theta)
%
% Inputs:
%      num_r - *Insert description of input variable here*
%
%      min_r - *Insert description of input variable here*
%
%      num_theta - *Insert description of input variable here*
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

%Work out limits of angular bins
thetas = linspace(0, 2*pi, num_theta+1);

%Make dist and theta map
xy = repmat(-max_rdist:max_rdist, mask_sz, 1);
dist_map = sqrt(xy.^2 + xy'.^2);
theta_map = mod(atan2(-xy', xy)-theta0,2*pi);

if do_plot; figure; p_num = 1; end
%Now make maps
sc_masks = false(mask_sz, mask_sz, num_r, num_theta);
for i_r = 1:num_r
    for i_theta = 1:num_theta
        
        sc_masks(:,:,i_r,i_theta) = ...
            (dist_map > rdists(i_r)) &...
            (dist_map <= rdists(i_r+1)) &...
            (theta_map >= thetas(i_theta)) &...
            (theta_map < thetas(i_theta+1));
        
        if do_plot;
            subplot(num_r, num_theta, p_num);
            imgray(sc_masks(:,:,i_r,i_theta));
            p_num = p_num+1;
        end
    end
end
                    


