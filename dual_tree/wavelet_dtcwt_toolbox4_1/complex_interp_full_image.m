function z_out = complex_interp_full_image(z_in, M, N, level, w, interpmethod)
%COMPLEX_INTERP_FULL_IMAGE interpolate a given band of a DT-CWT to the full pixel grid
%   [z_out] = complex_interp_full_image(z_in, M, N, level, w, interpmethod)
%
% Inputs:
%      z_in - 2D complex matrix to be interpolated
%
%      M - number of rows in full image
%
%      N - number of columns in full image
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

%get size of complex input
[m n] = size(z_in);

%check size is as expected if not give warning
if m*2^level ~= M || n*2^level ~= N
    warning('Matrix sizes do not match'); %#ok
end

%set up locations of sampling points at this location
x = repmat(0:n+1, m+2, 1);
y = repmat((0:m+1)', 1, n+2);

%set up locations of full pixel grid
st = .5 - 2^-(level + 1);
xi = repmat(st + (2^-level)*(1:N), M, 1);
yi = repmat(st + (2^-level)*(1:M)', 1, N);

%make central frequencies complex
jw = sqrt(-1) * w;

% Pad the complex input by one/row column
z_pad = zeros(m+2, n+2);
z_pad(2:m+1, 2:n+1) = z_in;

% Unwrap the phase
z_unwrap = z_pad .* exp(-x*jw(2)-y*jw(1));
%z_unwrap = z_pad;

% Interpolate to the full image grid
%z_out_unwrap = interp2(x, y, z_unwrap, xi, yi, interpmethod);
z_mag_unwrap = interp2(x, y, abs(z_unwrap), xi, yi, interpmethod);
z_phase_unwrap = interp2(x, y, z_unwrap ./ (abs(z_unwrap)+1e-6), xi, yi, 'linear');

%NAG code below is much faster but requires installing toolobox
% [px, py, lamda_r, mu_r, c_r] = e01da(x(1,:), y(:,1), real(z_unwrap(:)));
% [px, py, lamda_i, mu_i, c_i] = e01da(x(1,:), y(:,1), imag(z_unwrap(:)));
% z_out_unwrap = ...
%     reshape(complex(e02df(xi(1,:), yi(:,1), lamda_r, mu_r, c_r),...
%     e02df(xi(1,:), yi(:,1), lamda_i, mu_i, c_i)), size(xi));
    
    
% Rewrap the phases.
%z_out = z_out_unwrap .* exp(xi*jw(2)+yi*jw(1));
z_out = z_mag_unwrap  .* exp(xi*jw(2)+yi*jw(1) + i*angle(z_phase_unwrap));
%z_out = z_mag_unwrap .* exp(xi*jw(2)+yi*jw(1));
%z_out = z_out_unwrap;
return;
