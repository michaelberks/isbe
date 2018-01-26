function [hist_val] = hist_val(hist_struc, x)
%HIST_VAL *Insert a one line summary here*
%   [hist_val] = hist_val(hist_struc, x)
%
% Inputs:
%      hist_struc, x - *Insert description of input variable here*
%
%
% Outputs:
%      hist_val - *Insert description of input variable here*
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
bin_idx = floor((x - hist_struc.xmin)/ hist_struc.xincr + 0.5) + 1;
if bin_idx >=1 && bin_idx <= hist_struc.nbins
    hist_val = hist_struc.xcounts(bin_idx);
else
    hist_val = 0;
end
