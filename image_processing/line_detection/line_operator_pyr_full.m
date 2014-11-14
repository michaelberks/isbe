function [response] = line_operator_pyr_full(im, num_angles, ...
						  num_levels, angle_type)
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
  
for scale = 1:num_levels
    
    [r c] = size(g_pyr{scale});
    
    % Generate the lines and masks for this scale.
    [lines, neighbourhoods] = make_masks(5, num_angles, angle_type);
    
    response{scale} = zeros(r, c, num_angles);
    
    % Calculate the neighbourhood sums
    for angle_index=1:length(lines)
        
            n_response = imfilter(g_pyr{scale}, neighbourhoods{angle_index}, 'symmetric', 'conv');
            l_response = imfilter(g_pyr{scale}, lines{angle_index}, 'symmetric', 'conv');
            response{scale}(:,:,angle_index) = l_response - n_response;
    end

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
  