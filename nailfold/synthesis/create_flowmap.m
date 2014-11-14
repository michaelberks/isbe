function [flowmap, mask] = ...
    create_flowmap(path_points, widths, dirs, map_sz)

if (nargin==0 && nargout==0), test_script(); return; end

mask = zeros(map_sz);
flowmap = complex(nan(map_sz), nan(map_sz));
n_points = size(path_points, 1);

max_flow = 0;
wmin = min(widths);

for x = 1:map_sz(2)
    for y = 1:map_sz(1)
        diff = [path_points(:,1)-x path_points(:,2)-y];
        d = sum(diff.*diff, 2);

        [dist_squared, ind_nearest] = min(d);
        % Create a 'sink' around the endpoint where flow is zero.
        % (Cells are 'reborn' when they reach this part of the image.)
        if (ind_nearest == n_points)
            continue;
        end
        
        dist_nearest = sqrt(dist_squared);
        w = widths(ind_nearest);
        
        % Assume the narrowest point of the arterial limb is the 
        % width of one blood cell. 
        % Therefore, anything more than (w-wmin)/2 from the centre 
        % must be colliding with a wall.
        is_inside_vessel = (dist_nearest < w/2);
        is_colliding = (dist_nearest >= (w-wmin)/2);
        
        if is_inside_vessel
            if is_colliding
                % Add a component of flow toward the vessel centre
                vec = path_points(ind_nearest,:) - [x,y];
                vec = vec / norm(vec);

                wt = 1.0; % Arguable
                flowdir = complex(dirs(ind_nearest,1) + wt*vec(1), ...
                                  dirs(ind_nearest,2) + wt*vec(2));
            else
                % Set flow to act along the vessel centre
                flowdir = complex(dirs(ind_nearest,1), ...
                                  dirs(ind_nearest,2));
            end
            
            % Normalize direction vector
            flowdir = flowdir / abs(flowdir);
            
            % Flow should be inversely proportional to w^2
            flowmag = 1 / (w*w);
            flowmap(y,x) = flowmag * flowdir;
            
            if (flowmag > max_flow)
                max_flow = flowmag;
            end
            
            % Update mask image
%             mask(j,i) = 1; % Rectangle profile
            mask(y,x) = cos(dist_nearest / (w/2) * pi/2); % Cosine profile
            
            % Make darkness proportional to thickness of the vessel
            % i.e. Enhance contrast of thicker limb of the vessel
%             mask(j,i) = w * mask(j,i);
        end
    end
end

flowmap = flowmap / max_flow;

% x = linspace(-3, 3, 9);
% g = exp(-0.5*x.*x);
% mask = conv2(g, g, mask, 'same');
mask = mask / max(mask(:));


function test_script()
clc;

s = [ 
   211418064
  2094606752
];
randn('state', s);

f_profile = true;
if (f_profile), profile clear; profile on; end
    pts = generate_path();

    [widths, dirs, inner, outer] = ...
        generate_edges(pts, [], 10);
    
    [pts, widths, inner, outer, imsz] = ...
        scale_vessel(pts, widths, inner, outer, 160);

    [flowmap, mask] = ...
        create_flowmap(pts, widths, dirs, imsz);
if (f_profile), profile off; profile report; end

figure(1); clf;    
subplot(1,3,1); hold on;
    plot(pts(:,1), pts(:,2), 'b.-');
    plot(pts(1,1), pts(1,2), 'go');
    plot(pts(end,1), pts(end,2), 'ro');
    plot(inner(:,1), inner(:,2), 'r.-');
    plot(outer(:,1), outer(:,2), 'g.-');
    axis('equal', 'ij', [1,imsz(2),2,imsz(1)]);
subplot(1,3,2);
    show_flow_as('rgb', flowmap);
subplot(1,3,3); colormap(gray(256));
    image(uint8(mask*255)); 
    axis('image');


        
        