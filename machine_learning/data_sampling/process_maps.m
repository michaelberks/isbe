function [fg_map, bg_map] = process_maps(line_map, centre_map, fov_map, ...
                                         output_type, shrink_fov)

%If logical, assume line-map is a fg_mask we've loaded    
if islogical(line_map)    
    % Skeletonize the line map for certain prediction types if a centre map is
    % not given.
    if isempty(centre_map)
        switch output_type
            case {'centre_detection', ...
                  'junction_detection', 'junction_centre_detection', ...
                  'centre_orientation', 'mixed_centre_orientation',...
                  'centre_width'}
                % Create skeletonised vessel mask
                centre_map = bwmorph(line_map, 'thin', 'inf');
        end
    end
    
    %BG maps are just the negative of the line maps
    nonline_map = ~line_map;
    noncentre_map = ~centre_map;
    
else   
    %If not logical, assume it's a resampling map, which may have a second
    %frame giving bg sampling probs
    if size(line_map, 3) > 1
        nonline_map = line_map(:,:,2);
        line_map = line_map(:,:,1);
    else
        nonline_map = [];
    end
    %For centre methods, the line map IS the centre line sampling map
    centre_map = line_map;
    noncentre_map = nonline_map;
end

% If FoV mask defined then use it to mask line pixels
if ~isempty(fov_map)
    % Erode FoV mask for some methods (to remove pixels on the border?)
    if shrink_fov
        fov_map = fov_map & ~imerode(fov_map, strel('disk', 50));
    end
    
    line_map(~fov_map) = 0; %Remember line map may not be logical
    nonline_map(~fov_map) = 0;
    
    % If skeleton defined then mask this, too
    if ~isempty(centre_map)
        centre_map(~fov_map) = 0;
        noncentre_map(~fov_map) = 0;
    end
end 

% compute final foreground and background masks
switch output_type
    case {'detection', 'width'}
        fg_map = line_map;
        bg_map = nonline_map;

    case {'centre_detection'}
        fg_map = centre_map;
        bg_map = noncentre_map;
        
    case {'orientation', 'mixed_orientation'}
        fg_map = line_map;
        bg_map = nonline_map;
        
    case {'centre_orientation', 'mixed_centre_orientation', 'centre_width'}
        fg_map = centre_map;
        bg_map = [];

    case {'junction_detection', 'junction_centre_detection'}
        junction_mask = ...
            get_junction_mask(line_map, centre_map, output_type);

        % foreground pixels are those within the vessel and close to a junction
        % background are those within the vessel but not close to a junction
        fg_map = junction_mask;
        bg_map = (line_map & ~junction_mask);

    otherwise
        error(['Output type ', output_type, ' not recognized']);
end   
