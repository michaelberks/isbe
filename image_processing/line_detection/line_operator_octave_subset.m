function [response] = line_operator_octave_subset(im, num_angles, num_scales, rr, cc, if_plot)
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
if nargin < 6
    if_plot = 0;
end

angle_type = 'degrees';
sample_idx = sub2ind(size(im), rr, cc);

% Inititialise outputs
response = zeros([size(sample_idx), num_angles, num_scales]);
  
for scale = 1:num_scales
    
    % Generate the lines and masks for this scale.
    [lines, neighbourhoods] = make_masks(scale, num_angles, angle_type);   
    
    for angle_index=1:length(lines),
        % Calculate the line sums for each angle.
        l_strength = imfilter(im, lines{angle_index}, 'symmetric', 'conv');

        % Calculate the neighbourhood sums.
        n_strength = imfilter(im, neighbourhoods{angle_index}, 'symmetric', 'conv');
        
        response(:,:,angle_index,scale) = l_strength(sample_idx) - n_strength(sample_idx);
        
        if if_plot
            figure; imagesc(l_strength - n_strength); axis image; colormap(gray(256));
            hold on; plot(cc, rr, 'r.');
        end
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
  