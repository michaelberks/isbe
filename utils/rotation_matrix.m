function [R] = rotation_matrix(theta, radians)
%ROTATION_MATRIX *Insert a one line summary here*
%   [R] = rotation_matrix(theta, radians)
%
% Inputs:
%      theta - *Insert description of input variable here*
%
%      radians - *Insert description of input variable here*
%
%
% Outputs:
%      R - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 17-Apr-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
if nargin < 2 || radians
    c = cos(theta);
    s = sin(theta);
else
    c = cosd(theta);
    s = sind(theta);
end
R = [c -s; s c];
