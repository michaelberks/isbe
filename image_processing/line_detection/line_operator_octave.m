function [response, orientations] = line_operator_octave(im, num_angles, ...
						  num_scales, angle_type, if_plot)
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
    if_plot = 0;
end

if nargin < 4
    angle_type = 'degrees';
end

% A few sanity checks.
if min(im(:)) < 0
    im = im + min(im(:));
end

% Inititialise.
response = -inf;
response = response(ones(size(im)));
orientations = zeros(size(im));

% The neighbourhood response and its orientation at each scale.
response_scale = zeros(size(im));
orientations_scale = zeros(size(im));
  
for scale = 1:num_scales
    
    % Generate the lines and masks for this scale.
    [lines, neighbourhoods, angles] = make_masks(scale, num_angles, angle_type);
  
    strength = zeros(size(im));
    
    
    for angle_index=1:length(lines),
        % Calculate the line sums for each angle.
        tmp_strength = imfilter(im, lines{angle_index}, 'symmetric', 'conv');
        inds = tmp_strength > strength;
        orientations_scale(inds) = angles(angle_index);
        strength(inds) = tmp_strength(inds);

        % Calculate the neighbourhood sums.
        r = imfilter(im, neighbourhoods{angle_index}, 'symmetric', 'conv');
        inds = orientations_scale==angles(angle_index);
        response_scale(inds) = r(inds);
    end
    
    %Compute the final response as the difference between line strength and
    %neighbourhood response
    response_scale = strength - response_scale;
    
    %Normalise the response at this scale
    response_scale = ((response_scale - min(response_scale(:))) ...
		      ./ (max(response_scale(:)) - min(response_scale(:))));

    %Save the response for any location with maximal repsonse at this scale 
    inds = response < response_scale;
    response(inds) = response_scale(inds);
    orientations(inds) = orientations_scale(inds);
    
    if if_plot
        figure; imagesc(response); axis image; colormap(gray(256));
    end
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
    
    %Convert scale to filter length: 1->3, 2->5, 4->9 etc
    filter_length = (2^scale) + 1;
    
    centre = (filter_length + 1) / 2;
    xi = centre - (2^(scale-1) / 3);
    x1 = floor(xi);
    x2 = ceil(xi);
    
    line = zeros(filter_length);
    line(1:x1-1,:) = 0;
    line(x1,:) = (x2-xi)^2;
    line(x2,:) = 2*(x2-xi) - (x2-xi)^2 + 1;
    line(x2+1:centre,:) = 2;
    
    line(centre+1:end,:) = line(centre-1:-1:1,:);
    
    neighbourhood = ones(filter_length);
    
    for ii = 1:num_angles
        theta = angles(ii);
        
        l = imrotate(line, theta, 'bilinear', 'crop');
        lines{ii} = l / sum(l(:)); %#ok
        n = imrotate(neighbourhood, theta, 'bilinear');
        neighbourhoods{ii} = n ./ sum(n(:)); %#ok
        
    end
  