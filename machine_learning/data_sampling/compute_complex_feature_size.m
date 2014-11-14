function [num_features] = compute_complex_feature_size(sampling_args)
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
switch sampling_args.feature_type
    case {'all', 'conj', 'ilp', 'icp', 'real_imag', 'real_abs_imag'}
        num_features = 2;

    case {'real', 'mag', 'phase', 'complex', 'imag'}    
        num_features = 1;

    otherwise
        warning(['Feature type: ', sampling_args.feature_type, ' not recognised, using phase and magnitude (feature_type = ''all'')']); %#ok
        num_features = 2;
end