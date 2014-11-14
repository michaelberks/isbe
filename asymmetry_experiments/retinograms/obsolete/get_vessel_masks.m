function [fg, bg] = get_vessel_masks(imlist, args)
% Compute foreground and background masks given a vessel image

% Load in vessel mask
vessel_mask = u_load([args.vessel_mask_dir imlist.v_mask]);

% Skeletonize the vessel mask for certain prediction types
switch args.prediction_type
    case {'centre_detection', ...
          'junction_detection', ...
          'junction_centre_detection', ...
          'centre_orientation', ...
          'mixed_centre_orientation'}
        %Create skeletonised vessel mask
        centre_mask = bwmorph(vessel_mask, 'thin', 'inf');
        
    otherwise
        centre_mask = [];
end

% Load foveal mask for some prediction types (why only the detection ones?)
switch args.prediction_type
    case {'detection', 'centre_detection', 'foveal_detection', ...
          'logistic_classification'}
        %Load in mask of foveal region
        foveal_mask = u_load([args.foveal_mask_dir imlist.f_mask]);
        
    otherwise
        foveal_mask = [];
end

% Erode (?) foveal mask for some methods
switch args.prediction_type
    case {'foveal_detection'}
        foveal_mask = foveal_mask & ~imerode(foveal_mask, strel('disk', 50));
end

% If foveal mask defined then use it to mask vessel pixels
if ~isempty(foveal_mask)
    vessel_mask = vessel_mask & foveal_mask;
    
    % If skeleton defined then erode this, too
    if ~isempty(centre_mask)
        centre_mask = centre_mask & foveal_mask;
    end
end

% compute final foreground and background masks
switch args.prediction_type
    case {'detection', 'foveal_detection', 'logistic_classification'}
        fg = (vessel_mask & foveal_mask);
        bg = (~vessel_mask & foveal_mask);

    case {'centre_detection'}
        fg = (centre_mask & foveal_mask);
        bg = (~centre_mask & foveal_mask);
        
    case {'orientation', 'mixed_orientation', ...
          'linear_regression', 'logistic_regression', 'boosted_regression', ...
          'width'}
        fg = vessel_mask;
        bg = [];
        
    case {'centre_orientation', 'mixed_centre_orientation'}
        fg = centre_mask;
        bg = [];
        
    case {'width'}
        % Load in width map
        width_map = u_load([args.width_dir imlist.width]);
        fg = width_map;
%         fg = f(width_map); % transform width_map to weighted map
        bg = [];

    case {'junction_detection', 'junction_centre_detection'}
        % Get xy coords of points
        [c_y c_x] = find(centre_mask);
        num_pts = length(c_x);

        % Distance transform of vessel mask - gives Euclidean distance of
        % every point to background segment
        v_thick = bwdist(~vessel_mask);
        
        % Create a mask of size ws x ws, with ones on the border and zeros
        % in the centre
        ws = 5;
        ws2 = floor(ws/2);
        offs = -ws2:ws2;
        count_mask = ones(ws); 
        count_mask(2:ws-1,2:ws-1) = 0;
        
        % Set up meshgrid and preallocate junction mask image
        [rows cols] = size(vessel_mask);
        [xx,yy] = meshgrid(1:cols, 1:rows);
        junction_mask = false(rows, cols);

        for ii = 1:num_pts
            patch = centre_mask(c_y(ii) + offs, c_x(ii) + offs);
            patch = bwselect(patch, ws2+1, ws2+1, 8) & count_mask;
            label = bwlabel(patch, 8);

            % if more than two vessel sprout from this point then add points
            % to the junction mask
            if max(label(:)) > 2
                if strcmp(args.prediction_type,'junction_detection')
                    % add a fat circle that lies just within the vessel
                    rad = v_thick(c_y(ii), c_x(ii));
                    circ = (xx - c_x(ii)).^2 + (yy - c_y(ii)).^2 < rad^2;
                    circ = circ & vessel_mask;
                    junction_mask = junction_mask | circ;
                else
                    % just add the centreline point
                    junction_mask(c_y(ii) + offs, c_x(ii) + offs) = true;
                end
            end
        end

        % foreground pixels are those within the vessel and close to a junction
        % background are those within the vessel but not close to a junction
        fg = junction_mask;
        bg = (~junction_mask & vessel_mask);

    otherwise
        if strcmp(get_username, 'ptresadern')
            error(['Unknown output type: ',args.prediction_type]);
        end
end   
