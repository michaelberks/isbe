function [cdf_data] = cdf_transform3d(data_in, roi_mask, apply_by_slice)
%CDF_TRANSFORM transform data array of any size containing arbitrary
%distribution of values to same size array P(x) where P is the CDF of
%observed values x
%   [cdf_data] = cdf_transform(data_in)
%
% Inputs:
%      data_in - numerical array, any size
%
%
% Outputs:
%      cdf_data - array, same size as data in where cdf_data(i) =
%      P(data_in(:)), for CDf P(x)
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 16-Oct-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
if nargin < 2 || isempty(roi_mask)
    roi_mask = true(size(data_in));
end
    
if nargin < 3
    apply_by_slice = true;
end

if apply_by_slice
    [n_y, n_x, n_z] = size(roi_mask);
    cdf_data = zeros(n_y, n_x, n_z);
    for i_z = 1:n_z
        cdf_data(:,:,i_z) = cdf_transform(data_in(:,:,i_z), roi_mask(:,:,i_z));
    end
else
    cdf_data = cdf_transform(data_in, roi_mask);
end
