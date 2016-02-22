function [flowmap, mask, vessel_centre, widths] = ...
    create_flowmap_profile(vessel_centre, inner_edge, outer_edge, apex_idx, map_sz, cell_sz)

if ~exist('map_sz', 'var') || isempty(map_sz)
    min_x = min(min(outer_edge(:,1)), min(inner_edge(:,1)));
    min_y = min(min(outer_edge(:,2)), min(inner_edge(:,2)));
    
    if min_x > 32
        outer_edge(:,1) = outer_edge(:,1) - min_x + 32;
        inner_edge(:,1) = inner_edge(:,1) - min_x + 32;
        vessel_centre(:,1) = vessel_centre(:,1) - min_x + 32;
    end
    
    if min_y > 32
        outer_edge(:,2) = outer_edge(:,2) - min_x + 32;
        inner_edge(:,2) = inner_edge(:,2) - min_x + 32;
        vessel_centre(:,2) = vessel_centre(:,2) - min_x + 32;
    end
    
    max_x = max(max(outer_edge(:,1)), max(inner_edge(:,1)));
    max_y = max(max(outer_edge(:,2)), max(inner_edge(:,2)));
    
    map_sz = [max_y + 4 - rem(max_y,4), max_x + 32 - rem(max_x,4)];
end

if ~exist('cell_sz', 'var')
    cell_sz = 10;
end
    
%Check which is the top apex (i.e. the one with lowest y-coordinate)
num_apex = length(apex_idx);
if num_apex > 1
    [~,top_apex] = min(vessel_centre(apex_idx,2));
else
    top_apex = 1;
end

%Compute vessel widths from inner and outer edges
widths = sqrt(sum((outer_edge - inner_edge).^2,2));

%Check mean widths either side of apex, set the smaller to be the arterial
%side, and reverse edges and vessel path if necessary
w1 = mean(widths(1:apex_idx(top_apex)));
w2 = mean(widths(apex_idx(top_apex):end));

if w1 > w2
    widths = widths(end:-1:1,:);
    vessel_centre = vessel_centre(end:-1:1,:);
    %outer_edge = outer_edge(end:-1:1,:);
    %inner_edge = inner_edge(end:-1:1,:);
end

directions = zeros(size(vessel_centre));
directions(1,:) = vessel_centre(2,:) - vessel_centre(1,:);
directions(2:end-1,:) = vessel_centre(3:end,:) - vessel_centre(1:end-2,:);
directions(end,:) = vessel_centre(end,:) - vessel_centre(end-1,:);

mask = zeros(map_sz(1), map_sz(2), num_apex+1);
nans = nan(map_sz(1), map_sz(2), num_apex+1);
flowmap = complex(nans, nans);
n_points = size(vessel_centre, 1);

max_flow = 0;

%Now set ceiling for minimum cell width (which defines how wide the
%collision boundary is set)
wmin = min(min(widths), cell_sz);

apex_idx = [1; apex_idx(:); n_points];
for segment = 1:num_apex+1
    
    %Get vessel centre and direction for this layer
    pts = apex_idx(segment):apex_idx(segment+1);
    vessel_centre_i = vessel_centre(pts,:);
    widths_i = widths(pts);
    directions_i = directions(pts,:);
    
    for x = 1:map_sz(2)
        for y = 1:map_sz(1)
            diff = [vessel_centre_i(:,1)-x vessel_centre_i(:,2)-y];
            d = sum(diff.*diff, 2);

            [dist_squared, ind_nearest] = min(d);
            % Create a 'sink' around the endpoint where flow is zero.
            % (Cells are 'reborn' when they reach this part of the image.)
            if (segment==num_apex+1 && ind_nearest == length(pts))
                continue;
            end

            dist_nearest = sqrt(dist_squared);
            w = widths_i(ind_nearest);

            %PT code:
            % Assume the narrowest point of the arterial limb is the 
            % width of one blood cell. 
            % Therefore, anything more than (w-wmin)/2 from the centre 
            % must be colliding with a wall.
            %MB change - this only makes sense for normal width vessels.
            %For enlarged vessel better to assume cell size (even if this
            %means preselecting one) - see how wmin is set above
            is_inside_vessel = (dist_nearest < w/2);
            is_colliding = (dist_nearest >= (w-wmin)/2);

            if is_inside_vessel
                if is_colliding
                    % Add a component of flow toward the vessel centre
                    vec = vessel_centre_i(ind_nearest,:) - [x,y];
                    vec = vec / norm(vec);

                    wt = 1.0; % Arguable
                    flowdir = complex(directions_i(ind_nearest,1) + wt*vec(1), ...
                                      directions_i(ind_nearest,2) + wt*vec(2));
                else
                    % Set flow to act along the vessel centre
                    flowdir = complex(directions_i(ind_nearest,1), ...
                                      directions_i(ind_nearest,2));
                end

                % Normalize direction vector
                flowdir = flowdir / abs(flowdir);

                % Flow should be inversely proportional to w^2
                damping = 0.001;
                flowmag = 1 / (w*w) + damping;
                flowmap(y,x,segment) = flowmag * flowdir;

                if (flowmag > max_flow)
                    max_flow = flowmag;
                end

                % Update mask image
    %             mask(j,i) = 1; % Rectangle profile
                mask(y,x,segment) = cos(dist_nearest / (w/2) * pi/2); % Cosine profile

                % Make darkness proportional to thickness of the vessel
                % i.e. Enhance contrast of thicker limb of the vessel
    %             mask(j,i) = w * mask(j,i);
            end
        end
    end
end

flowmap = flowmap / max_flow;

% x = linspace(-3, 3, 9);
% g = exp(-0.5*x.*x);
% mask = conv2(g, g, mask, 'same');
mask = mask / max(mask(:));



        
        