function [ddG,G,dG] = Wfilt(theta,sigma,halfwidth)

if nargin<1, theta = 0; end
if nargin<2, sigma = 3; end
if nargin<3, halfwidth = 3*sigma; end
	
% y is -x since the image axes are upside down ('ij' instead of 'xy')
x = -halfwidth:halfwidth; y = x;
[xx,yy] = meshgrid(x,y);

r = [cos(theta); sin(theta)];

k1 = 1/(2*pi*sigma*sigma);
k2 = -0.5/(sigma*sigma);
k3 = -1/(sigma*sigma);

% generate gaussian
G = zeros(length(y),length(x));
dG = zeros(length(y),length(x));
ddG = zeros(length(y),length(x));
for i = 1:length(x)
	for j = 1:length(y)
		G(j,i)		= k1*exp(k2*(xx(j,i)*xx(j,i) + (yy(j,i)*yy(j,i))));
		dG(j,i)		= k3*G(j,i)*(x(i)*r(1)+y(j)*r(2));
		ddG(j,i)	= k3*(dG(j,i)*(x(i)*r(1)+y(j)*r(2)) + G(j,i));
	end
end

if nargout==0,
	figure; imagesc(ddG); axis('image'); colormap(gray(256));
	clear ddG;
end


