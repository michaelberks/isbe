function [grad_strength, grad_orientation] = ...
    gaussian_1st_derivative_gradient_old(im, scales)
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

f_debug = (nargin==0 && nargout==0);
if f_debug
	datapath = 'A:\asym\data\mammograms\2004_screening_processed\mass_roi';
	load(fullfile(datapath,'045LML_roi.mat'));
	im = bg(401:700,151:450);
	
	scales = [1 2 4 8];
end

if ~isa(im, 'double')
    im = double(im);
end

min_grad_strength = inf*ones(size(im));
for ii = 1:length(scales)
	scale = scales(ii);
	
	ssq = scale ^ 2;
	width = round(6*scale);

	%make the 1st derivate filters - y direction is just transpose of x
	%direction
	[x,y] = meshgrid(-width:width, -width:width);
	dgau2D =-x.*exp(-(x.*x+y.*y)/(2*ssq)) * sqrt(2) / (ssq*sqrt(pi*ssq));
	%dgau2D = dgau2D / sum(dgau2D(:));

	% Convolve image with the directional derivatives
	ax = imfilter(im, dgau2D, 'conv','replicate');
	ay = imfilter(im, dgau2D', 'conv','replicate');

	%Get the magnitude of the filter response
	grad_strength = sqrt((ax.*ax) + (ay.*ay));

	%Now do the maximal suppression
	%compute orientations at each point
	grad_orientation = atan(-ay./ax);

	%We may have some NaNs where ax and ay are both zero, so set these to
	%zero
	grad_orientation(isnan(grad_orientation)) = 0;
	
	% get minimum gradient strength
	min_grad_strength = min(min_grad_strength,grad_strength);
end

if f_debug
	figure(1); colormap(gray(256));
		image(im); axis('image');
	figure(2); 
		imagesc(min_grad_strength); axis('image');
    clear;
end
