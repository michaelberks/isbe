function [response, orientations] = line_operator_pyr(im, num_angles, ...
						  num_levels, angle_type, circle)
% LINE_OPERATOR_CONV - Find lines in 2D images
%   @[response, orientations] = line_operator_conv(im, num_angles,
%   min_scale, max_scale, angle_type)@
%
%   @im@ is the image, @num_angles@ is the number of angles at which
%   to apply the line (it is applied every 180/num_angles degrees),
%   @min_scale@ and @max_scale@ are the minimum and maximum scales to
%   use.  @angle_type@ determines how the orientations will be
%   reported, can be @'degrees'@ (the default) or @'radians'@.
%   Orientations are between 0 and 180 degrees or 0 and %$\pi$%
%   radians.

if nargin < 5
    circle = 0;
end

if nargin < 4
    angle_type = 'degrees';
end

% A few sanity checks.
if min(im(:)) < 0
    im = im + min(im(:));
end

%compute Gaussian pyramid of image
g_pyr = mb_g_pyramid(im, num_levels);

% Inititialise outputs
response = cell(num_levels, 1);
orientations = cell(num_levels, 1);  
  
for scale = 1:num_levels
    
    % Generate the lines and masks for this scale.
    [lines, neighbourhoods, angles] = make_masks(5, num_angles, angle_type);
    if circle
        neighbourhoods = make_neighbourhood2(5);
    end
    
    response{scale} = -inf(size(g_pyr{scale}));
    orientations{scale} = zeros(size(g_pyr{scale}));
    strength = zeros(size(g_pyr{scale}));
    
    % Calculate the line sums for each angle.
    for angle_index=1:length(lines),
      
        tmp_strength = imfilter(g_pyr{scale}, lines{angle_index}, 'symmetric', 'conv');
        inds = tmp_strength > strength;
        orientations{scale}(inds) = angles(angle_index);
        strength(inds) = tmp_strength(inds);
    end

    % Calculate the neighbourhood sums
    if ~circle
        for angle_index=1:length(lines),
            r = imfilter(g_pyr{scale}, neighbourhoods{angle_index}, 'symmetric', 'conv');
            inds = orientations{scale} == angles(angle_index);
            response{scale}(inds) = r(inds);
        end
    else
        response{scale} = imfilter(g_pyr{scale}, neighbourhoods ./ sum(neighbourhoods(:)), 'symmetric', 'conv');
    end
    
    response{scale} = strength - response{scale};
    response{scale} = ((response{scale} - min(response{scale}(:))) ...
		      ./ (max(response{scale}(:)) - min(response{scale}(:))));

    %Need to think about a comparison between scales here?
end

%   close(w);
  
  
function [lines, neighbourhoods, angles] = ...
    make_masks(scale, num_angles, angle_type)
% MAKE_MASKS - Makes line-op masks
%   The masks are quite sparse, but sparse matrices cannot be used
%   as input to @conv2@.
  
switch angle_type
    case 'degrees'
        angles = 0:180/num_angles:180;
    case 'radians'
        angles = 0:pi/num_angles:pi;
    otherwise
        error(['Unknown angle type ''' angle_type '''']);
end
  
%   line_ratio = 3;
%   hw = (scale - 1) / 2;
%   
%   [x_pts y_pts] = meshgrid(...
%       [-hw hw],...
%       [-hw ((-1e-3) - hw/line_ratio) -hw/line_ratio hw/line_ratio 1e-3+hw/line_ratio hw]);
%   
%   z_pts = [0 0 1 1  0 0; 0 0 1 1 0 0]';
%   [xi_pts yi_pts] = meshgrid(-hw:hw, -hw:hw);
%   
%     for theta=angles
%         xo_pts = reshape([cos(theta) sin(theta)] * [xi_pts(:)'; yi_pts(:)'], size(xi_pts)); 
%         yo_pts = reshape([-sin(theta) cos(theta)] * [xi_pts(:)'; yi_pts(:)'], size(xi_pts));
  
    line = zeros(scale);
    line((scale+1)/2,:) = 1;
    nhood = ones(scale);
    neighbourhood = double(nhood);
    i = 1;

    for theta=angles,
        l = imrotate(line, theta, 'bilinear', 'crop');
        lines{i} = l / sum(l(:)); %#ok
        n = imrotate(neighbourhood, theta, 'bilinear');
        neighbourhoods{i} = n ./ sum(n(:)); %#ok
        i = i + 1;
    end
  