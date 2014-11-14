function [line_response, orientation, scale] = gaussian_1st_derivative_gradient(im, scales)
%GAUSSIAN_1ST_DERIVATIVE_GRADIENT *Insert a one line summary here*
%   [grad_strength, orientation] = gaussian_1st_derivative_gradient(im, scale)
%
% Inputs:
%      im - *Insert description of input variable here*
%
%      scale - *Insert description of input variable here*
%
%
% Outputs:
%      grad_strength - *Insert description of input variable here*
%
%      grad_orientation - *Insert description of input variable here*
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

if nargin == 1
    %Get number of scales and image size from input responses
    responses = im; clear im;
    if iscell(responses)
        [r c] = size(responses{1}(:,:,1));
        num_scales = length(responses);
    else
        [r c num_scales] = size(responses(:,:,:,1));
    end
else
    %Set default filter width if none specified
    if nargin < 3
        g_width = 5;
    end
    
    %Get size of image then pad
    im = double(im);
    [r c] = size(im);
    padsize = round(g_width*max(scales));
    im = padarray(im, [padsize padsize], 'replicate');
    num_scales = length(scales);
end

%pre-allocate output arguments 
line_response = zeros(r, c);
orientation = zeros(r, c);
scale = zeros(r, c);

for i_scale = 1:num_scales
    
    if exist('responses', 'var')
        if iscell(responses)
            %Interpolate the responses back to the full pixel grid
            scaling = 2^(i_scale-1);
            offset = (i_scale-1) / (2^i_scale);
            rr = offset + (1:r)' / scaling;
            cc = offset + (1:c) / scaling;
            Ix = interp2(responses{i_scale}(:,:,1), cc, rr, 'cubic');
            Iy = interp2(responses{i_scale}(:,:,2), cc, rr, 'cubic');
        else
            Ix = responses(:,:,i_scale,1);
            Iy = -responses(:,:,i_scale,2);
        end
        sigma = 2^(i_scale-1);
    else
        sigma = scales(i_scale);
        %Make 1st order directional filters
        [g,dg] = gaussian_filters_1d(scales(i_scale), round(g_width*sigma));

        %Filter the image
        % Filter the image with separable filters and remove padding
        Ix = conv2(g',dg,im,'same'); Ix = Ix(padsize+(1:r), padsize+(1:c));
        Iy = conv2(dg',g,im,'same'); Iy = -Iy(padsize+(1:r), padsize+(1:c));
    end
    
    theta = atan2(Iy, Ix);
    %theta = atan(Iy ./ Ix);
    wo_theta = filter_output(Ix, Iy, theta);
    
    %if wo_theta at this scale so far is biggest, swap it into the main
    %line_strength, orientation and scale matrices - note the theta are
    %normal to the line direction
    swap_idx = abs(wo_theta) > abs(line_response);
    line_response(swap_idx) = wo_theta(swap_idx);
    orientation(swap_idx) = theta(swap_idx) + pi/2;
    scale(swap_idx) = sigma;
    
end

%wrap orientations back into [-pi pi] range
orientation(orientation > pi) = orientation(orientation > pi) - 2*pi;


function wo_theta = filter_output(Ix, Iy, theta)
%Compute the output of the filter for an arbitrary angle
wo_theta = Ix.*cos(theta) + Iy.*sin(theta);