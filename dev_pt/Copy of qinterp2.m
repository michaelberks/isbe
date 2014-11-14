function [img_out] = qinterp2(img_in,x,y,method)
% Faster version of interp2. Also returns the weights so that they can be
% reused.

% use bilinear interpolation by default
if nargin<4, method = 'linear'; end

% get integer limits and fractional bit in between
xf = floor(xx); xc = ceil(xx); xr = x-xf;
yf = floor(yy); yc = ceil(yy); yr = y-yf;

% convert to a grid
[xf,yf] = meshgrid(xf,yf);
[xc,yc] = meshgrid(xc,yc);
[xr,yr] = meshgrid(xr,yr);

switch method
	case 'linear',
		% need to define m here
		q1	= (1-xr(:)) .* (1-yr(:)) .* img_in(1+(xf(:)-1)*m+(yf(:)-1));
		q2	=    xr(:)  .* (1-yr(:)) .* img_in(1+(xc(:)-1)*m+(yf(:)-1));
		q3	= (1-xr(:)) .*    yr(:)  .* img_in(1+(xf(:)-1)*m+(yc(:)-1));
		q4	=    xr(:)  .*    yr(:)  .* img_in(1+(xc(:)-1)*m+(yc(:)-1));
		img_out = reshape(q1+q2+q3+q4,size(xf));
		
	case 'cubic',
		% do bicubic interpolation
end
