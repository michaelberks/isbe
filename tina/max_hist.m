function [max_val, max_x] = max_hist(hist_struc)
%MAX_HIST *Insert a one line summary here*
%   [max_val, max_x] = max_hist(hist_struc)
%
% Inputs:
%      hist_struc - *Insert description of input variable here*
%
%
% Outputs:
%      max_val, max_x - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 29-Jun-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
[max_val max_i] = max(hist_struc.xcounts);

if max_i > 1 && max_i < hist_struc.nbins;
    h0  = max_val;
    hp1 = hist_struc.xcounts(max_i+1);
    hm1 = hist_struc.xcounts(max_i-1);
    % assuming parabola fit to 3 centre bins */
    centre = (hm1-hp1)/(2*(hp1+hm1)-4*h0);
    max_x = hist_struc.xmin+(max_i-0.5)*hist_struc.xincr;
    max_x = centre*(hist_struc.xincr) + max_x;
else  
    max_x = hist_struc.xbins(max_i);
end
