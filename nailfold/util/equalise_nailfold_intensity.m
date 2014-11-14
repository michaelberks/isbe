function [nailfold_im_equalised] = ...
    equalise_nailfold_intensity(nailfold_im, smoothing_size)
%EQUALISE_NAILFOLD_INTENSITY *Insert a one line summary here*
%   [nailfold_equalised] = equalise_nailfold_intensity(nailfold_im, smoothing_size)
%
% Inputs:
%      nailfold_im - *Insert description of input variable here*
%
%      smoothing_size - *Insert description of input variable here*
%
%
% Outputs:
%      nailfold_equalised - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 01-May-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

if ~exist('smoothing_size', 'var'); smoothing_size = 32; end
if ~isa(nailfold_im, 'double')
    nailfold_im = double(nailfold_im);
end

%Make smoothing filter
g = gaussian_filters_1d(smoothing_size, 3*smoothing_size);
g = g / sum(g);

%Make mask to get region edge
nailfold_mask = make_nailfold_mosaic_mask(nailfold_im);
im_edges = conv2(g', g, double(nailfold_mask), 'same');

%Smooth the nailfold
nailfold_im(~nailfold_mask) = 0;
nailfold_im_smoothed = conv2(g', g, nailfold_im, 'same') ./ im_edges;
nailfold_im_equalised = nailfold_im - nailfold_im_smoothed;
nailfold_im_equalised(~nailfold_mask) = 0;
