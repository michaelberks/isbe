function [img_label] = get_real_image_label(images, output_type)

% Load the appropriate label image
switch output_type
    case {'detection'}
        img_label = u_load(images.fg_mask);
        
    case {'centre_detection'}
        img_label = u_load(images.fg_mask);
        img_label = bwmorph(img_label, 'thin', 'inf');
        
    case {'class_label'}
        img_label = u_load(images.class_label);
        
    case {'orientation', 'centre_orientation'}
        img_label = u_load(images.ori);
        
    case {'mixed_orientation', 'mixed_centre_orientation'}
        % Could randomly sample one orientation per pixel and return this
        % as the ground truth label. This would, however, make it
        % impossible to use the same image location twice with different
        % labels.
        
        img_label = u_load(images.ori);

    case {'width'}
        % Load in width map
        img_label = u_load(images.width);
        
    otherwise
        error(['Output type ', output_type, ' not recognized']);
end
