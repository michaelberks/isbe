function [centres_r, rot_mat, ncolsr, nrowsr, ncols, nrows, bad] = compute_mosiac_midline(fov_mask, varargin)
%COMPUTE_MOSIAC_MIDLINE *Insert a one line summary here*
%   [centres_r, rot_mat, ncolsr, nrowsr, ncols, nrows, bad] = compute_mosiac_midline(fov_mask, poly_n)
%
% Inputs:
%      fov_mask - *Insert description of input variable here*
%
%      poly_n - *Insert description of input variable here*
%
%
% Outputs:
%      centres_r - *Insert description of input variable here*
%
%      rot_mat - *Insert description of input variable here*
%
%      ncolsr - *Insert description of input variable here*
%
%      nrowsr - *Insert description of input variable here*
%
%      ncols - *Insert description of input variable here*
%
%      nrows - *Insert description of input variable here*
%
%      bad - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 04-Dec-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

args = u_packargs(varargin,... % the user's input
    0, ... % non-strict mode
    'ori_bins', -45:45,...
    'resize_factor', 8, ...
    'edge_ori_sigma', 8, ...
    'bad_mosaic_sz', 55,...
    'border_size', 2,...
    'min_depth', 470,...
    'poly_n', 5);
clear varargin;

%Get size of full size mask
[nrows ncols] = size(fov_mask);

%Reduce size of fov_mask, we don't need it so big
fov_mask4 = imresize(fov_mask, 1/args.resize_factor, 'nearest');
    
%See if the mask suggests a mis-registered mosaic
bad_mask = imopen(fov_mask4, strel('disk', args.bad_mosaic_sz));
bad = any(bad_mask(:));

%if bad we can bale out now
if bad
    centres_r = [];
    rot_mat = [];
    ncolsr = [];
    nrowsr = [];
    return;
end

%If not, get the anlge of the edges on the top and bottom of the mosaic
f_border = fov_mask & ~imerode(fov_mask, strel('disk', args.border_size));
[~, border_ori_map] = gaussian_1st_derivative_gradient(fov_mask, args.edge_ori_sigma);

border_ori = border_ori_map(f_border);
border_ori(border_ori < -pi/4) = [];
border_ori(border_ori > pi/4) = [];
border_ori = 180*border_ori/pi;

%Histogram these angles to work out what the rotation of the mosaic is
counts = hist(border_ori, args.ori_bins);
[~, max_idx] = max(counts);
max_ori = args.ori_bins(max_idx);
rot_mat = [cosd(max_ori) -sind(max_ori); sind(max_ori) cosd(max_ori)];

%Roate the mosaic so the frames are horizontal
fov_mask_r = imrotate(fov_mask, -max_ori, 'nearest', 'loose');
[nrowsr ncolsr] = size(fov_mask_r);

%Compute the height of the mosaic along it's width, and thus the midline
depth = zeros(ncolsr, 1);
centres_r = zeros(ncolsr, 2);
centres_r(:,1) = 1:ncolsr;
for i_col = 1:ncolsr
    depth(i_col) = sum(fov_mask_r(:,i_col));
    centres_r(i_col,2) = mean(find(fov_mask_r(:,i_col)));
end

%Discard the edges of the mosaic that aren't a full frame height 
start_col = find(depth > args.min_depth, 1);
end_col = find(depth > args.min_depth, 1, 'last');
centres_r = centres_r(start_col:end_col,:);

%Fit a smooth polynomial to midline
[pp, s, mu] = polyfit(centres_r(:,1), centres_r(:,2), args.poly_n);
centres_r(:,2) = polyval(pp, centres_r(:,1), s, mu);


