function [image_in, label, label_centre, label_orientation, parameters] = ...
	generate_grain_image(varargin)
%GENERATE_GRAIN_IMAGE Generates an image with a 'grain' to it
%   [image_in, label, label_centre, label_orientation, parameters] = generate_grain_image(bg, args)
%
% Inputs:
%      
%      sigma - spread of the image in the direction chosen
% 
%      theta - hand-specified orientation
% 
%      image_size - width (and height) of the output image
% 
%      scale, offset - image out = (img*scale)+offset
%
% Outputs:
%      image_in - Resulting image
%
%      label - Whether a line exists (true for all pixels)
%
%      label_centre - 'Centre' of the 'line' (true for all pixels)
%
%      label_orientation - Orientation of line at that point
%
%      parameters - Parameters of grain (i.e. orientation)
%
%
% Example:
%
% Notes:
%
% See also:
%   GENERATE_LINE_IMAGE
%
% Created: 02-Feb-2011
% Author: Phil Tresadern 
% Email : philip.tresadern@manchester.ac.uk 
% Phone : +44 (0)161 275 5114 
% Copyright: (C) University of Manchester 

warning('off','ASYM:unexpectedArgument');
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
		'sigma',5,...
    'image_size', 128,...
    'theta', 180*rand,...
		'scale',255,...
		'offset',0);
warning('on','ASYM:unexpectedArgument');
		
% rip variables of interest
sigma = args.sigma;
image_size = args.image_size;
theta = args.theta;
	
% create 1D filter
x = -3*sigma:3*sigma;
if sigma>0
	g = exp(-0.5*(x.*x)./(sigma.*sigma)); g = g/sum(g);
else
	g = zeros(size(x)); g(3*sigma+1) = 1;
end

% rotate 1D filter to 2D
g2 = zeros(length(x),length(x));
g2(3*sigma+1,:) = g;
g2 = imrotate(g2,theta);

% create test image and filter it
N			= image_size+max(size(g2))-1;
grain	= conv2(rand(N,N),g2,'valid');

% return images and parameters
image_in = (args.scale*grain)+args.offset;
label = true(size(image_in));
label_centre = true(size(image_in));
label_orientation = theta*ones(size(image_in));
parameters.theta = theta;
parameters.sigma = sigma;

% if run as a script then show output
if nargout==0
	figure; imshow(uint8(image_in));
	clear image_in;
end
