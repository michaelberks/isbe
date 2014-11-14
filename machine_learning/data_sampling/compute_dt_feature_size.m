function [n_pixels] = compute_dt_feature_size(sampling_args)
%COMPUTE_DT_FEATURE_SIZE *Insert a one line summary here*
%   [feature_size] = compute_dt_feature_size(sampling_args)
%
% Inputs:
%      sampling_args - *Insert description of input variable here*
%
%
% Outputs:
%      feature_size - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 14-Sep-2010
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
switch sampling_args.feature_shape
    case 'rect'
        n_pixels = sampling_args.win_size*sampling_args.win_size;

    case 'clock'
        n_pixels = (sampling_args.win_size+1);

    otherwise
        warning(['Feature shape: ', sampling_args.feature_shape, ' not recognised, using square windows (feature_type = ''rect'')']); %#ok
        n_pixels = sampling_args.win_size*sampling_args.win_size;
end