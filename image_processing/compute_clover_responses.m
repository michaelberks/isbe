function [clover_responses] = compute_clover_responses(im, scales)
%COMPUTE_CLOVER_RESPOSES Compute responses to two 'cloverleaf' filters
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
clover_responses = zeros(r, c, length(scales), 3);
im = double(im);

if f_debug
	figure(1); clf; colormap(gray(256));
	figure(2); clf;
end

for iscl = 1:length(scales)
    % Make 2nd order directional filters
    [g,dg,ddg] = gaussian_filters_1d(scales(iscl));

	if f_debug
		figure(1);
			subplot(length(scales),2,(iscl-1)*2+1);
			imagesc(2*dg'*dg); axis('image');
			subplot(length(scales),2,(iscl-1)*2+2);
			imagesc(g'*ddg-ddg'*g); axis('image');
		f1 = dg'*dg;
		f2 = g'*ddg-ddg'*g;
		disp([sum(abs(g)) sum(abs(dg))]);
		disp([sum(abs(f1(:))) sum(abs(f2(:)))]);
	end
    
	% Filter the image with separable filters
	clover_responses(:,:,iscl,1) = conv2(g',ddg,im,'same'); % = Ixx
	clover_responses(:,:,iscl,2) = conv2(dg',dg,im,'same'); % = Ixy
	clover_responses(:,:,iscl,3) = conv2(ddg',g,im,'same'); % = Iyy
% 	g2d_responses(:,:,1) = imfilter(im, filters(:,:,jj), 'same', 'replicate');
end

if f_debug
	clear;
end
