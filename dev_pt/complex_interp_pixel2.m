function z_out = complex_interp_pixel2(z_in, r, c, level, w, interpmethod)
%COMPLEX_INTERP_FULL_IMAGE interpolate a given band of a DT-CWT to the full pixel grid
%   [z_out] = complex_interp_pixel(z_in, r, c, level, w, interpmethod)
%
% Inputs:
%      z_in - 2D complex matrix to be interpolated
%
%      r - row of pixel location to interpolate to
%
%      c - column of pixel to interpolate to
%
%      level - the level of the DT-CWT the input band belongs to (so we
%      assume that M=m*2^level, N=n*2^level, where [m n] = size(z_in) )
%
%      w - the central frequencies of each sub-band (as set by Nick
%      Kingsbury)
%
%      interpmethod - string description of interpolation method passed to
%      Matlab's interp2 function
%
% Outputs:
%      z_out - interpolated complex array of size MxN
%
%
% Example:
%
% Notes: This code is based on the complex interpolation method kindly
% supplied by Nick Kingsbury (see below)
%
% See also: CPXINTERP2
%
% Created: 05-May-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

%Set default interpolation method
if nargin < 6 
    interpmethod = 'cubic'; 
end

% get size of complex input
[m n] = size(z_in);

% set up locations of full pixel grid
st = .5 - 2^-(level + 1);
xi = st + (2^-level)*c;
yi = st + (2^-level)*r;

%make central frequencies complex
jw = sqrt(-1) * w;

% compute pixels we'll need in order to interpolate (via mask creation)
x1 = floor(xi)*(m+2); y1 = floor(yi)+1;
x2 = ceil(xi)*(m+2);  y2 = ceil(yi)+1;
mask = zeros(m+2,n+2);
mask([x1+y1; x1+y2; x2+y1; x2+y2]) = 1;
mask = (conv2(mask,ones(3,3),'same')>0);
pixinds = find(mask);

[y,x] = ind2sub([m+2,n+2],pixinds);

exp_xjw = exp(-(x-1)*jw(2));
exp_yjw = exp(-(y-1)*jw(1));

% Pad the complex input by one/row column and unwrap the phase
z_unwrap = zeros(m+2,n+2);
z_unwrap(2:m+1,2:n+1) = z_in;
z_unwrap(pixinds) = z_unwrap(pixinds) .* exp_xjw .* exp_yjw;

% compute magnitude
z_mag = zeros(m+2,n+2);
z_mag(pixinds) = abs(z_unwrap(pixinds));

% Interpolate magnitude to the full image grid
[x,y] = meshgrid(0:n+1,0:m+1);
z_mag_int = interp2(x, y, z_mag, xi, yi, ['*',interpmethod]);

% Interpolate phase to the full image grid
z_tmp = zeros(m+2,n+2);
z_tmp(pixinds) = z_unwrap(pixinds) ./ (z_mag(pixinds)+1e-6);
z_phase_unwrap = angle(interp2(x, y, z_tmp, xi, yi, '*linear'));

% convert back to complex values
z_out = z_mag_int  .* exp(xi*jw(2)+yi*jw(1) + sqrt(-1)*z_phase_unwrap);

if 0
	% for debugging only
	z_out1 = complex_interp_pixel(z_in, r, c, level, w, interpmethod);
	err = z_out - z_out1;
	if norm(err,'fro')>1e-6
		display(err(1:5));
		keyboard
	end
end

return;
