function [smoothed_vol, orientations] = ...
    tangential_smoothing_vol(vol_in)
%TANGENTIAL_SMOOTHING_VOL *Insert a one line summary here*
%   [smoothed_image] = tangential_smoothing_vol(image_in)
%
% Inputs:
%      image_in - *Insert description of input variable here*
%
%      orientations - *Insert description of input variable here*
%
%
% Outputs:
%      smoothed_image - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 27-Jul-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

%Apply to tangential smoothing to each slice
[n_y, n_x, n_z] = size(vol_in);

smoothed_vol = zeros(n_y, n_x, n_z);
compute_orientations = nargout > 1;

if compute_orientations   
    orientations = zeros(n_y, n_x, n_z);
end

for i_slice = 1:n_z
    if compute_orientations
        [smoothed_slice, slice_orientations] = ...
            tangential_smoothing(vol_in(:,:,i_slice), [], 0);
        orientations(:,:,i_slice) = slice_orientations;
    else
        [smoothed_slice] = ...
            tangential_smoothing(vol_in(:,:,i_slice), [], 0);
    end
    smoothed_vol(:,:,i_slice) = smoothed_slice;
end
smoothed_vol(isnan(smoothed_vol)) = 0;
if compute_orientations   
    orientations(isnan(smoothed_vol)) = 0;
end
