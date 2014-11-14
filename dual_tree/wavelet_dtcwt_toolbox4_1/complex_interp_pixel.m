function z_out = complex_interp_pixel(z_in, r, c, level, varargin)
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

args = u_packargs(varargin, '0',...
    'interpmethod', 'cubic',...
    'band_frequency', [],...
    'interp_mag_phase', 1);

%Set default interpolation method
interpmethod = args.interpmethod; 
band_frequency = args.band_frequency;
interp_mag_phase = args.interp_mag_phase;
unwrap_phase = ~isempty(band_frequency);

%get size of complex input
[m n] = size(z_in);

%set up locations of sampling points at this location
x = repmat(0:n+1, m+2, 1);
y = repmat((0:m+1)', 1, n+2);

%set up locations of full pixel grid
st = .5 - 2^-(level + 1);
xi = st + (2^-level)*c;
yi = st + (2^-level)*r;

% Pad the complex input by one/row column
z_pad = zeros(m+2, n+2);
z_pad(2:m+1, 2:n+1) = z_in; clear z_in;

if unwrap_phase
    % Unwrap the phase
    z_pad = z_pad .* exp(-x*band_frequency(2) -y*band_frequency(1));
end    

%Now either interpolate the complex number (fast) or interpolate magnitude
%and phase separately
if interp_mag_phase
    z_mag = abs(z_pad);
    z_mag_int = interp2(x, y, z_mag, xi, yi, ['*',interpmethod]);
    z_phase_int = angle(interp2(x, y, z_pad ./ (z_mag+1e-6), xi, yi, '*linear'));
    z_out = z_mag_int  .* exp(1i*z_phase_int);        
else
    z_out = interp2(x, y, z_pad, xi, yi, ['*',interpmethod]);
end

if unwrap_phase
    % Rewrap the phases.
    z_out = z_out .* exp(xi*band_frequency(2)+yi*band_frequency(1));
end

return;
