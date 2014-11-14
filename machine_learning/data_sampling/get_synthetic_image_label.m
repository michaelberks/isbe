function img_label = get_synthetic_image_label(line_map, centre_map, ...
                                               orientation_map, width_map, ...
                                               output_type)

switch output_type
    case {'detection'}
        img_label = line_map;
        
    case {'centre_detection'}
        img_label = centre_map;
        
    case {'orientation', 'centre_orientation'}
        img_label = orientation_map;
        
    case {'width'}
        img_label = width_map;

    otherwise
        error(['Output type ', output_type, ' not recognized']);
end
