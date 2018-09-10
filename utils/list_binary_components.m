function [bin_comps] = list_binary_components(dec_in)
%LIST_BINARY_COMPONENTS *Insert a one line summary here*
%   [bin_comps] = list_binary_components(dec_in)
%
% Inputs:
%      dec_in - *Insert description of input variable here*
%
%
% Outputs:
%      bin_comps - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 12-Feb-2018
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
max_factor = floor(log2(dec_in));
bin_comps = [];
for pow = 0:max_factor
    bin_comp = 2^pow;
    if rem(dec_in, 2*bin_comp)
        bin_comps(1,end+1) = bin_comp; %#ok
        dec_in = dec_in - bin_comp;
    end
end
        

