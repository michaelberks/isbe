function [out_var, scaling] = save_uint8(filepath, in_var)
%SAVE_UINT8 Save a double array as a unit8, keeping a record of the scaling
%required to produce values in range [0, 255]
%   [out_var, scaling] = save_uint8(filepath, in_var)
%
% Inputs:
%      filepath - to save variable
%
%      in_var - double array to be saved
%
%
% Outputs:
%      out_var - input array now in unit8 format, with the range of in_var
%      maximally spread over [0, 255]
%
%      scaling - struc containing fields 'm', 'c' such that: c =
%      min(in_var), m = max(in_var) - min(in_var), and thus 
%       out_var = unit8( 255 * (in_var - c) / m) )
%
%
% Example:
%
% Notes:
%
% See also: LOAD_UINT8
%
% Created: 27-Oct-2010
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

if isnumeric(in_var)
    
    if ~isreal(in_var)
        scaling.complex = true;
        in_var = [real(in_var); imag(in_var)];
    end
    zero_mask = abs(in_var) > 0;
    c = min(in_var(zero_mask));
    m = max(in_var(zero_mask)) - c;

    out_var = uint8(255 * (in_var(zero_mask) - c) / m );
    scaling.m = m;
    scaling.c = c;

    save(filepath, 'out_var', 'scaling', 'zero_mask');
else
    warning('save_uint8:incorrect_type', 'Variable in_var is non-numeric and will be saved using SAVE');
    save(filepath, 'in_var');
end
        
        


