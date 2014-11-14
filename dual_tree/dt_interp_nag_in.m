function [knot_mag knot_im knot_re dt_dims] = dt_interp_nag_in(dual_tree, levels)
%DT_INTERP_NAG_IN *Insert a one line summary here*
%   [] = dt_interp_nag_in()
%
% Inputs:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 11-Aug-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester
if nargin < 2
    levels = 1:(size(dual_tree, 1) - 1);
end

%Get number of levels to put in full tree
num_levels = length(levels);

%Set central frequencies of bands for interp method
w = [-3 -1; -sqrt(5) -sqrt(5); -1 -3; 1 -3; sqrt(5) -sqrt(5); 3 -1]*pi/2.15;

%Pre-allocate space for interpolation knotpoints
knot_mag = cell(num_levels,1);
knot_im = cell(num_levels,1);
knot_re = cell(num_levels,1);
dt_dims = zeros(num_levels,2);
%for each sub-band, interpolate up the coefficients to the full size
for l = 1:num_levels
    
    %get size of complex input
    [m n] = size(dual_tree{levels(l)}(:,:,1));
    dt_dims(l,1) = m;
    dt_dims(l,2) = n;
    
    knot_mag{l} = zeros((m+2)*(n+2),1);
    knot_im{l} = zeros((m+2)*(n+2),1);
    knot_re{l} = zeros((m+2)*(n+2),1);

    %set up locations of sampling points at this location
    x = repmat(0:n+1, m+2, 1);
    y = repmat((0:m+1)', 1, n+2);
    
    for band = 1:6

        %make central frequencies complex
        jw = sqrt(-1) * w(band,:);

        % Pad the complex input by one/row column
        z_pad = zeros(m+2, n+2);
        z_pad(2:m+1, 2:n+1) = dual_tree{levels(l)}(:,:,band);

        % Unwrap the phase
        z_unwrap = z_pad .* exp(-x*jw(2)-y*jw(1)); clear z_pad;

        z_mag = abs(z_unwrap);
        % Compute augmented  knot points for interpolation using NAG library
        [dummy, dummy, dummy, dummy, knot_mag{l}(:,band)] = ...
            e01da(x(1,:), y(:,1), z_mag(:));
        [dummy, dummy, dummy, dummy, knot_re{l}(:,band)] = ...
            e01da(x(1,:), y(:,1), real(z_unwrap(:)) ./ (z_mag(:)+1e-6));
        [dummy, dummy, dummy, dummy, knot_im{l}(:,band)] = ...
            e01da(x(1,:), y(:,1), imag(z_unwrap(:)) ./ (z_mag(:)+1e-6));
        clear z_mag z_unwrap z_in;
    end
end
