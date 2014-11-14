function display_weighted_orientation(line_strength, line_ori, spacing, mask, sigma)
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
if nargin < 5
    sigma = [];
end
if nargin < 4
    mask = [];
end
if nargin < 3
    spacing = 8;
end
[r c] = size(line_strength);

if ~any(imag(line_strength(:)));
    line_ori = pi*line_ori/180;
    line_strength = complex(line_strength .* cos(line_ori), line_strength .* sin(line_ori));
end

if ~isempty(sigma)
    line_strength = imfilter(line_strength, fspecial('gaussian', 6*sigma+1, sigma));
else
    sigma = 1;
end

if isempty(mask)
    %mask = abs(line_strength) > 0.25*max(abs(line_strength(:)));
    mask = mb_non_maximal_supp(abs(line_strength), angle(line_strength)) > 0;
end

figure;
a1 = axes;
%i1 = image(complex2rgb(line_strength .^2)); axis image; hold on;
i1 = imagesc(abs(line_strength .^2)); axis image; hold on;
set(i1, 'ButtonDownFcn', @display_hist);

spacing_mask = false(size(mask));
spacing_mask(1:spacing:r, 1:spacing:c) = true;

mask = mask & spacing_mask;

num_angles = 36;
ang_res = pi / num_angles;
colors = hsv(num_angles);
for ii = 0:num_angles
    theta = (ii - 0.5)*ang_res;
    %Get mask of pixels that have orientation within theta range
    angle_mask = mask & ...
                 spacing_mask & ...
                 (angle(line_strength) > theta - 0.5*ang_res) &...
                 (angle(line_strength) <= theta + 0.5*ang_res);
    [yy xx] = find(angle_mask);

    h = quiver(xx, yy, ...
        4*sigma*real(line_strength(angle_mask)),...
       -4*sigma*imag(line_strength(angle_mask)), 0, 'color', colors(mod(ii, num_angles)+1,:), 'linewidth', sigma);
   set(h, 'HitTest', 'off');
end

 
    function display_hist(hObject, eventdata) %#ok

        XYZ = get(a1, 'CurrentPoint');
        x = round(XYZ(1,1));
        y = round(XYZ(1,2));

        display([x y]);
        
        if (x > 1) && (x <= c) && (y > 1) && (y <= r)
            figure('name', ['Orientation histogram for x = ' num2str(x), ' y = ' num2str(y)]);
            for jj = 1:4
                g_width = 2^(2 + jj);
            
                %g_filt = fspecial('gaussian', 6*g_width+1, g_width);
                roi = sample_window(line_strength, 6*g_width + 1, y, x);
                
                subplot(2,2,jj);
                weighted_complex_rose(roi, 40);
                title(['\sigma = ' num2str(g_width)]);
            end  
        end

    end
end
%% Coloured arrows code
% num_angles = 36;
% ang_res = pi / num_angles;
% colors = hsv(num_angles);
% 
% for ii = 1:num_angles
%     theta = (ii - 0.5)*ang_res;
%     %Get mask of pixels that have orientation within theta range
%     angle_mask = mask & ...
%                  spacing_mask & ...
%                  (angle(orientations) > theta - 0.5*ang_res) &...
%                  (angle(orientations) <= theta + 0.5*ang_res);
%     [yy xx] = find(angle_mask);
% 
%     quiver(xx, yy, ...
%         4*spacing*real(orientations(angle_mask)),...
%        -4*spacing*imag(orientations(angle_mask)), 0, 'color', colors(ii,:));
% end



