function [index_map] = map_to_nearest_val(feature_map, feature_vals)
%MAP_TO_NEAREST_VAL *Insert a one line summary here*
%   [index_map] = map_to_nearest_val(feature_map, feature_val)
%
% Inputs:
%      feature_map - *Insert description of input variable here*
%
%      feature_vals - *Insert description of input variable here*
%
%
% Outputs:
%      index_map - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 01-May-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
[rows cols] = size(feature_map);

[~,index_map] = min( abs(bsxfun(@minus, feature_map(:), feature_vals(:)')), [], 2);
index_map = reshape(index_map, rows, cols);
