function [g_pyr] = mb_g_pyramid(image_in, num_levels, filt)
%MB_G_PYRAMID *Insert a one line summary here*
%   [g_pyr] = mb_g_pyramid(image_in,num_levels,filt)
%
% Inputs:
%      image_in- *Insert description of input variable here*
%
%      num_levels- *Insert description of input variable here*
%
%      filt- *Insert description of input variable here*
%
%
% Outputs:
%      g_pyr- *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 06-Apr-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

if nargin < 3
    filt = [1 4 6 4 1] / 16;
end

half_size = (length(filt)-1) / 2;

g_pyr = cell(num_levels, 1);

g_pyr{1} = image_in; clear image_in;
lo = g_pyr{1};
for lev = 2:num_levels
    lo = conv2(padarray(lo, [half_size 0], 'symmetric'), filt', 'valid');
    lo(2:2:end,:) = [];

    lo = conv2(padarray(lo, [0 half_size], 'symmetric'), filt, 'valid');
    lo(:,2:2:end) = [];
    
    g_pyr{lev} = lo;
end
