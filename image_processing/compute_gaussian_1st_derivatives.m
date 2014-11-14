function [responses] = compute_gaussian_1st_derivatives(im, scales)
% COMPUTE_GAUSSIAN_DERIVATIVES Compute responses to gaussian 1st derivative filters
%   [line_strength, orientation, scale] = compute_gaussian_derivatives(im, scales)
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
	im = rand(100,100);
	scales = [1 2 4 8];
	clc;
end

% pre-allocate output arguments
[r c] = size(im);
responses = zeros(r, c, length(scales), 2);
im = double(im);

for iscl = 1:length(scales)
    % Make 2nd order directional filters
    [g,dg,ddg] = gaussian_filters_1d(scales(iscl));

	% Filter the image with separable filters
	responses(:,:,iscl,1) = conv2(g',dg,im,'same'); % = Ix
	responses(:,:,iscl,2) = conv2(dg',g,im,'same'); % = Iy
end

if f_debug
	clear;
end

