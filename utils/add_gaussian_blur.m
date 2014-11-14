function [im_r, im_c, gxy] = add_gaussian_blur(xc, yc, sigma, n_rows, n_cols, width)
%ADD_GAUSSIAN_BLUR *Insert a one line summary here*
%   [im_r, im_c, gxy] = add_gaussian_blur(xc, yc, sigma, n_rows, n_cols, width)
%
% Inputs:
%      xc - *Insert description of input variable here*
%
%      yc - *Insert description of input variable here*
%
%      sigma - *Insert description of input variable here*
%
%      n_rows - *Insert description of input variable here*
%
%      n_cols - *Insert description of input variable here*
%
%      width - *Insert description of input variable here*
%
%
% Outputs:
%      im_r - *Insert description of input variable here*
%
%      im_c - *Insert description of input variable here*
%
%      gxy - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 23-Oct-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
if ~exist('width', 'var') || isempty(width)
    width = ceil(3*sigma);
end

sigma2 = sigma^2;

xci = round(xc);
yci = round(yc);

dx = -width:width;

im_r = dx + yci;
im_c = dx + xci;

discard_rows = (im_r < 1) | (im_r > n_rows);
discard_cols = (im_c < 1) | (im_c > n_cols);

im_r(discard_rows) = [];
im_c(discard_cols) = [];

[x y] = meshgrid(dx(~discard_cols), dx(~discard_rows));

dxy2 = -(x.^2 + y.^2) / sigma2;
gxy = sigma*exp(dxy2) / (2*pi*sigma2);

