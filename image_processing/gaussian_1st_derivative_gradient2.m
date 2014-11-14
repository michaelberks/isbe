function [line_strength, orientation, scale] = gaussian_1st_derivative_gradient2(im, scales)
%GAUSSIAN_1ST_DERIVATIVE_GRADIENT *Insert a one line summary here*
%   [grad_strength, orientation] = gaussian_1st_derivative_gradient(im, scale)
%
% Inputs:
%      im - *Insert description of input variable here*
%
%      scale - *Insert description of input variable here*
%
%
% Outputs:
%      grad_strength - *Insert description of input variable here*
%
%      grad_orientation - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 01-Dec-2010
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
display('This function is deprecated, and now calls gaussian_2nd_derivative_line.');
[line_strength, orientation, scale] = gaussian_1st_derivative_gradient(im, scales);