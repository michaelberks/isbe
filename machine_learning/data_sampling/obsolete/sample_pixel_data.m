function [im_samples] = sample_pixel_data(image_in, rows, cols, varargin)
%SAMPLE_PIXEL_DATA *Insert a one line summary here*
%   [dt_samples] = sample_dt_polar13(dt, rows, cols, rotate)
%
% Inputs:
%      image_in - image to decompose and sample
%
%      num_levels - num_levels to include
%
%      rows - row subscripts of sample locations
%
%      cols - column subscripts of sample locations
%
%      win_size - dimensions of window to sample
%
%      rotate - if true (default false) circular shift the columns of the
%      final representation so maximal orientation is in first row, thus
%      achieving approximate rotation invariance (+/- 15 degrees)
%
%
% Outputs:
%      linop_samples - *Insert description of input variable here*
%
%
% Example:
%
% Notes: Representation is based on that presented by Kingsbury in
% "Rotation-Invariant Local Feature Matching with Complex Wavelets"
%
% See also:
%
% Created: 23-Aug-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

% Unpack the arguments:
args = u_packargs(varargin, '0',...
    'win_size', 15,...
    'pca', []);

win_size = args.win_size; 

%Make copies of sample rows and cols at positions of local window patch
pad_w = floor(win_size/2);
win_idx = -pad_w:pad_w;
num_samples = length(rows);
rr = repmat(rows*ones(1,win_size) + ones(num_samples,1)*win_idx, 1, win_size);
cc = kron(cols*ones(1,win_size) + ones(num_samples,1)*win_idx, ones(1,win_size));

patch_idx = sub2ind(size(image_in), rr, cc);

im_samples = image_in(patch_idx);
im_samples = bsxfun(@minus, im_samples, mean(im_samples, 2));

if ~isempty(args.pca)
    %Transform sample using PCA modes
    im_samples = bsxfun(@minus, im_samples, args.pca.mean)*args.pca.modes;
end