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
        [nrows ncols] = size(images.fg_mask);
        nlayers = size(images.ori);
        img_label = images.ori(:,:,1);
        
        for i_r = 1:nrows
            for i_c = 1:ncols
                if images.fg_mask(i_r, i_c)
                    not_valid = true;
                    width = ceil(images.width(i_r, i_c) / 2);
                    width2 = 2*width + 1;
                    while not_valid
                        offset_r = i_r + ceil(width2*rand) - width - 1;
                        offset_c = i_c + ceil(width2*rand) - width - 1;
                        not_valid = ...
                            offset_r < 1 || offset_r > nrows || ...
                            offset_c < 1 || offset_c > ncols ...
                            || ~images.fg_mask(offset_r, offset_c);
                    end
                    offset_l = ceil(nlayers*rand);
                    img_label(i_r, i_c) = images.ori(offset_r,offset_c,offset_l);
                end
            end
        end
        %img_label = u_load(images.ori);

    case {'width'}
        % Load in width map
        img_label = u_load(images.width);
        
    otherwise
        error(['Output type ', output_type, ' not recognized']);
end
