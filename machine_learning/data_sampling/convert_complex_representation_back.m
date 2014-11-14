function [complex_samples band_reduction] = convert_complex_representation_back(complex_samples, feature_type)
%CONVERT_COMPLEX_REPRESENTATION *Insert a one line summary here*
%   [] = convert_dt_representation(varargin)
%
% CONVERT_COMPLEX_REPRESENTATION uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Optional Arguments:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 22-Feb-2011
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
n_cols = size(complex_samples, 2);
switch feature_type
    case 'all'
        %Convert from mag and phase
        complex_samples = complex_samples(:,1:n_cols/2) .* exp(1i * complex_samples(:,1+n_cols/2:end) );
        band_reduction = 2;
    case 'real_imag'
        complex_samples = complex(complex_samples(:,1:n_cols/2), complex_samples(:,1+n_cols/2:end) );
        band_reduction = 2;
     case 'complex'
        %We're ok - do nothing
        band_reduction = 1;
    case {'mag', 'phase', 'real', 'imag', 'conj', 'ilp', 'icp', 'real_abs_imag'}
        error(['Cannot convert back to a full complex representation from ' feature_type ' data.']);

    otherwise
        error(['Complex feature type ', feature_type, ' not recognized']);
end 
