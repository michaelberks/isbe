function [grain_im,theta] = create_grain_repeat(sig,theta,imsize)

if nargin<1,	sig = 5; end
if nargin<2,	theta = 180*rand; end
if nargin<3,	imsize = 128; end

% create 1D filter
x = -3*sig:3*sig;
if sig>0
	g = exp(-0.5*(x.*x)./(sig.*sig)); g = g/sum(g);
else
	g = zeros(size(x)); g(3*sig+1) = 1;
end

theta

% rotate 1D filter to 2D
g2 = zeros(length(x),length(x));
g2(3*sig+1,:) = g;
g2 = imrotate(g2,theta);

% create test image and filter it
a  = 255*rand(imsize,imsize);
a  = repmat(a,[3,3]);
a2 = conv2(a, g2, 'valid');
a2 = a2(1:imsize, 1:imsize);
a2 = repmat(a2, [3,3]);

% return fixed size image
grain_im = uint8(a2);

if nargout==0
	figure; imshow(grain_im);
	clear grain_im theta;
end