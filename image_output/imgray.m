function [h] = imgray(im_in, n)
%IMGRAY wrapper to function to display image with grayscale map
%   [fig] = imgray(im_in, n)
%
% Inputs:
%      im_in - image to display
%
%      n - number of grayscales to use
%
%
% Outputs:
%      fig - figure handle
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 17-Jan-2012
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
if nargin < 2
    n = 256;
end
h = imagesc(im_in); colormap(gray(n)); axis image; hold all;
