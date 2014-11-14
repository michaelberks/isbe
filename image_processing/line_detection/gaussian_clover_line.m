function [line_strength, orientation, scale] = gaussian_clover_line(im, scales, g_width)
%GAUSSIAN_2ND_DERIVATIVE_LINE *Insert a one line summary here*
%   [line_strength, orientation, scale] = gaussian_2nd_derivative_line(im, scales)
%
% Inputs:
%      im - *Insert description of input variable here*
%
%      scales - *Insert description of input variable here*
%
%
% Outputs:
%      line_strength - *Insert description of input variable here*
%
%      orientation - *Insert description of input variable here*
%
%      scale - *Insert description of input variable here*
%
%
% Example:
%
% Notes: This function is deprecated and now redirects to call GAUSSIAN_2ND_DERIVATIVE_LINE. Note that
% the 2 functions were always designed to produce identical results and
% only co-existed because of a development oversight. Redirecting as below
% ensures implementation between the functions remains consistent.
%
% See also: GAUSSIAN_2ND_DERIVATIVE_LINE
%
% Created: 01-Dec-2010
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
display('This function is deprecated, and now calls gaussian_2nd_derivative_line.');
if nargin < 5
    g_width = 5;
end
[line_strength, orientation, scale] = gaussian_2nd_derivative_line(im, scales, g_width);