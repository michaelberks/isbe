function [] = display_orientation(f, image_in, orientations, spacing, mask, cmap)
%DISPLAY_ORIENTATION *Insert a one line summary here*
%   [] = display_orientation(image_in, orientations, mask, spacing)
%
% Inputs:
%      image_in - *Insert description of input variable here*
%
%      orientations - *Insert description of input variable here*
%
%      mask - *Insert description of input variable here*
%
%      spacing - *Insert description of input variable here*
%
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 30-Sep-2010
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
if nargin < 6
    cmap = gray(256);
end
if nargin < 5
    mask = [];
end
if nargin < 4
    spacing = 8;
end
figure(f); imagesc(image_in); axis image; colormap(cmap); hold on;
[y x] = size(image_in);

is_complex = any(imag(orientations(:)));

if ~is_complex
    orientations = pi*orientations/180;
end
%quiver(1:spacing:x, 1:spacing:y, cos(orientations(1:spacing:y, 1:spacing:x)), -sin(orientations(1:spacing:y, 1:spacing:x)));

num_angles = 36;
ang_res = pi / num_angles;
colors = hsv(num_angles);

    
if isempty(mask)
    mask = true(y,x);
end

spacing_mask = false(size(mask));
spacing_mask(1:spacing:y, 1:spacing:x) = true;

if is_complex
    for ii = 1:num_angles
        theta = (ii - 0.5)*ang_res;
        %Get mask of pixels that have orientation within theta range
        angle_mask = mask & ...
                     spacing_mask & ...
                     (mod(angle(orientations),pi) > theta - 0.5*ang_res) &...
                     (mod(angle(orientations),pi) <= theta + 0.5*ang_res);
        [yy xx] = find(angle_mask);

        quiver(xx, yy, ...
            4*spacing*real(orientations(angle_mask)),...
           -4*spacing*imag(orientations(angle_mask)), 0, 'color', colors(ii,:));
    end
else
    mask = mask & spacing_mask;
    [yy xx] = find(mask);
    quiver(xx, yy, cos(orientations(mask)), -sin(orientations(mask)), 'r');    
end





