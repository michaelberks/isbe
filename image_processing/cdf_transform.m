function [cdf_data] = cdf_transform(data_in, roi_mask)
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
    
[sorted_vals, sorted_idx] = sort(data_in(roi_mask));
n_vals = length(sorted_vals);
cdf_vals = (0:n_vals-1) / n_vals;
[~, cdf_idx] = sort(sorted_idx);
cdf_vals = cdf_vals(cdf_idx);
cdf_data = zeros(size(data_in));
cdf_data(roi_mask) = cdf_vals;
