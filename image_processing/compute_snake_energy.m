function [E_i, E_a, E_b] = compute_snake_energy(xy_pts, feat_img)
%COMPUTE_SNAKE_ENERGY *Insert a one line summary here*
%   [E_i, E_a, E_b] = compute_snake_energy(xy_pts)
%
% Inputs:
%      xy_pts - *Insert description of input variable here*
%
%
% Outputs:
%      E_i - *Insert description of input variable here*
%
%      E_a - *Insert description of input variable here*
%
%      E_b - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 23-Apr-2014
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
E_i = sum(interp2(feat_img, xy_pts(:,1), xy_pts(:,2), 'linear'));

xy_d1 = diff(xy_pts);
xy_d2 = diff(xy_d1);

%Sum of norm
%E_a = sum(sqrt(sum(xy_d1.^2,2)));
%E_b = sum(sqrt(sum(xy_d2.^2,2)));
E_a = sum(xy_d1(:).^2);
E_b = sum(xy_d2(:).^2);
