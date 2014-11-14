function [response, orientations, scales] = line_operator_conv(im, num_angles, ...
						  min_scale, max_scale, ...
						  angle_type, mask, circle)
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

if nargin < 7
    circle = 0;
end
  
if nargin < 6
    mask = [];
end
  
if nargin < 5
    angle_type = 'degrees';
end

  % A few sanity checks.
  if min(im(:)) < 0
    im = im + min(im(:));
  end

  % We prefer odd numbered lengths for the masks.
  if mod (min_scale, 2) == 0
    min_scale = min_scale - 1;
  end
  
  if mod (max_scale, 2) ~= 0
    max_scale = max_scale + 1;
  end
  
  if min_scale > max_scale
    error('min_scale > max_scale');
  end

  % Avoid a divide by zero error.
  if min_scale == max_scale
    max_scale = min_scale + 1;
  end

  % Inititialise.
  response = -inf;
  response = response(ones(size(im)));
  orientations = zeros(size(im));
  scales = zeros(size(im));

  % The neighbourhood response and its orientation at each scale.
  response_scale = zeros(size(im));
  orientations_scale = zeros(size(im));
  
%   w = waitbar(0, 'Applying lin-op...');
  
  for scale=min_scale:2:max_scale,
    % Generate the lines and masks for this scale.
    [lines, neighbourhoods, angles] = make_masks(scale, num_angles, angle_type);
    if circle
      neighbourhoods = make_neighbourhood2(scale);
    end
  
    strength = zeros(size(im));
    % Calculate the line sums for each angle.
    for angle_index=1:length(lines),
      
      tmp_strength = imfilter(im, lines{angle_index}, 'symmetric', 'conv');
      inds = tmp_strength > strength;
      orientations_scale(inds) = angles(angle_index);
      strength(inds) = tmp_strength(inds);
    end

    % Calculate the neighbourhood sums.
    
    if ~circle
      for angle_index=1:length(lines),
        r = imfilter(im, neighbourhoods{angle_index}, 'symmetric', 'conv');
        inds = orientations_scale==angles(angle_index);
        response_scale(inds) = r(inds);
      end
    else
        response_scale = imfilter(im, neighbourhoods ./ sum(neighbourhoods(:)), 'symmetric', 'conv');
    end
    
    response_scale = strength - response_scale; 
    if ~isempty(mask)
        response_scale(~mask) = min(response_scale(mask));
    end
    response_scale = ((response_scale - min(response_scale(:))) ...
		      ./ (max(response_scale(:)) - min(response_scale(:))));

    inds = response < response_scale;
    response(inds) = response_scale(inds);
    orientations(inds) = orientations_scale(inds);
    scales(inds) = scale;
%     waitbar((scale - min_scale + 1) / (max_scale - min_scale), w);
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
  