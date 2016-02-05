function [output_sz, compound_transforms, dx, dy] = ...
    mosaic_limits(input_sz, compound_transforms)

% Return the corners of the mosaic, given an input image size and a list of 
% compound transforms (with respect to the first image).

x_min = 1;
x_max = input_sz(2);

y_min = 1;
y_max = input_sz(1);

n_frames = size(compound_transforms,3);
        
for i = 1:n_frames
    % Transform the corner points using the compounded transforms to update
    % the max size of the final mosaic
    corners = [ 1           1; 
                1           input_sz(1);
                input_sz(2) 1;
                input_sz(2) input_sz(1) ];     
            
    corners_t = (compound_transforms(:,:,i) * [corners ones(4,1)]')';

    x_min = min( x_min, min(corners_t(:,1)) );
    x_max = max( x_max, max(corners_t(:,1)) );
    
    y_min = min( y_min, min(corners_t(:,2)) );
    y_max = max( y_max, max(corners_t(:,2)) );
end

% Add a translation that will take the top left pixel to (1,1)
dx = 1 - x_min;
dy = 1 - y_min;

compound_transforms(1,3,:) = compound_transforms(1,3,:) + dx;
compound_transforms(2,3,:) = compound_transforms(2,3,:) + dy;

% Return size of output image corresponding to the compound transforms
% Note this is in (row,column) format.
output_sz = [ceil(y_max + dy), ceil(x_max + dx)];
