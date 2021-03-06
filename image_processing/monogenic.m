% MONOGENIC - Michael Berks' test version of MONOFILT:
% Apply monogenic filters to an image to obtain 2D analytic signal
%
% See copyright below
%
% Usage: [local_amp, local_phase, local_ori] = ...
%       	monogenic(im, nscale, minWaveLength, mult, sigmaOnf)
%                             3         4           2     0.65
% Arguments:
% The convolutions are done via the FFT.  Many of the parameters relate 
% to the specification of the filters in the frequency plane.  
%
%   Variable       Suggested   Description
%   name           value
%  ----------------------------------------------------------
%    im                        Image to be convolved.
%    nscale          = 3;      Number of filter scales.
%    minWaveLength   = 4;      Wavelength of smallest scale filter.
%    mult            = 2;      Scaling factor between successive filters.
%    sigmaOnf        = 0.65;   Ratio of the standard deviation of the
%                              Gaussian describing the log Gabor filter's
%                              transfer function in the frequency domain
%                              to the filter center frequency. 
% 
% Returns:
%  
%        A  - cell array of monogenic energy responses.
%    theta  - cell array of phase orientation responses.
%      psi  - cell array of phase angle responses.
%    local_phase - cell array of local phase
%
%
% Experimentation with sigmaOnf can be useful depending on your application.
% I have found values as low as 0.2 (a filter with a *very* large bandwidth)
% to be useful on some occasions.
%
% See also: GABORCONVOLVE

% References:
% Michael Felsberg and Gerald Sommer. "A New Extension of Linear Signal
% Processing for Estimating Local Properties and Detecting Features"
% DAGM Symposium 2000, Kiel 
%
% Michael Felsberg and Gerald Sommer. "The monogenic signal" IEEE
% Transactions on Signal Processing, 49(12):3136-3144, December 2001

% Copyright (c) 2004-2005 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

%  October 2004 - Original version.
%  May     2005 - Orientation wrapping and code cleaned up.
%  August  2005 - Phase calculation improved.

function [local_amp local_phase local_ori] = ...
	monogenic(im, nscale, minWaveLength, mult, sigmaOnf, orientWrap)

    if nargin == 5
        orientWrap = 0;  % Default is no orientation wrapping
    end
    
    if ~isa(im,'double')
        im = double(im);
    end
    
    [rows,cols] = size(im);    
    IM = fft2(im);
    
    % Generate horizontal and vertical frequency grids that vary from
    % -0.5 to 0.5 
    [u1, u2] = meshgrid(([1:cols]-(fix(cols/2)+1))/(cols-mod(cols,2)), ...
			([1:rows]-(fix(rows/2)+1))/(rows-mod(rows,2))); %#ok

    u1 = ifftshift(u1);   % Quadrant shift to put 0 frequency at the corners
    u2 = ifftshift(u2);
    
    radius = sqrt(u1.^2 + u2.^2);    % Matrix values contain frequency
                                     % values as a radius from centre
                                     % (but quadrant shifted)
    
    % Get rid of the 0 radius value in the middle (at top left corner after
    % fftshifting) so that taking the log of the radius, or dividing by the
    % radius, will not cause trouble.
    radius(1,1) = 1;
    
    H1 = i*u1./radius;   % The two monogenic filters in the frequency domain
    H2 = i*u2./radius;
    
    f = im;
	h1f = real(ifft2(IM.*H1));
	h2f = real(ifft2(IM.*H2));
    
    local_amp = zeros(rows, cols, nscale+1);
	local_ori = zeros(rows, cols, nscale+1);
    local_phase = zeros(rows, cols, nscale+1);
    
	local_amp(:,:,1) = sqrt(f.^2 + h1f.^2 + h2f.^2);
    local_ori(:,:,1) = atan2(h1f, h2f);
    local_phase(:,:,1) = atan2(f, sqrt(h1f.^2 + h2f.^2));
    
    % The two monogenic filters H1 and H2 are oriented in frequency space
    % but are not selective in terms of the magnitudes of the
    % frequencies.  The code below generates bandpass log-Gabor filters
    % which are point-wise multiplied by H1 and H2 to produce different
    % bandpass versions of H1 and H2
    
    for s = 2:nscale+1
        wavelength = minWaveLength*mult^(s-1);
        fo = 1.0/wavelength;                  % Centre frequency of filter.
        logGabor = exp((-(log(radius/fo)).^2) / (2 * log(sigmaOnf)^2));  
        logGabor(1,1) = 0;                    % undo the radius fudge.	

        % Generate bandpass versions of H1 and H2 at this scale
        H1s = H1.*logGabor; 
        H2s = H2.*logGabor; 

        %  Apply filters to image in the frequency domain and get spatial
            %  results 
        f = real(ifft2(IM.*logGabor));    
        h1f = real(ifft2(IM.*H1s));
        h2f = real(ifft2(IM.*H2s));

        local_amp(:,:,s) = sqrt(f.^2 + h1f.^2 + h2f.^2);  % Magnitude of Energy.

        % Calculate the orientation and phase angles 
        local_ori(:,:,s) = atan2(h1f, h2f);              % Orientation.
        local_phase(:,:,s) = atan2(f, sqrt(h1f.^2 + h2f.^2));

    end
     
    % Here phase is measured relative to the h1f-h2f plane as an
    % 'elevation' angle that ranges over +- pi/2
    if orientWrap
        % Wrap orientation values back into the range 0-pi
        negind = find(local_ori < 0);
        local_ori(negind) = local_ori(negind) + pi;

        % Where orientation values have been wrapped we should
                % adjust phase accordingly **check**
%         local_phase(negind) = pi-local_phase(negind);
%         morethanpi = find(local_phase>pi);
%         local_phase(morethanpi) = local_phase(morethanpi)-2*pi;
    end