function [cell_positions] = ...
    generate_cell_positions(flowmap, mask, ...
                            path_points, widths, ...
                            n_cells, n_frames, ...
                            sigma_f, samples_per_exposure)

if (nargin==0 && nargout==0), test_script(); return; end

if ~exist('sigma_f','var') || isempty(sigma_f), sigma_f = 0.4; end
if ~exist('samples_per_exposure','var') || isempty(samples_per_exposure), samples_per_exposure = 1; end

% Subsample in time for smoother motion
if samples_per_exposure ~= 1
    % More frames and a longer burn in
    n_frames = n_frames * samples_per_exposure;
    
    % Slower flow and smaller sigma
    flowmap = flowmap / samples_per_exposure;
    sigma_f = sigma_f / samples_per_exposure; % Probably wrong - need to check
end
    
imsz = size(flowmap);
[xx,yy] = meshgrid(1:imsz(2), 1:imsz(1));

cell_pos_t = zeros(n_cells, 3);
cell_positions = zeros(n_cells, 3, n_frames);

if ~isempty(path_points)
    path_points = complex(path_points(:,1), path_points(:,2));
end

if isempty(mask)
    % Generate initial states, spread evenly along the vessel centre
    for c = 1:n_cells
        ind = ceil(rand*size(path_points,1));
        p0 = path_points(ind); % create new cells along the vessel path
        cell_pos_t(c, 1) = p0 + sigma_f*complex(randn, randn);
        cell_pos_t(c, 2) = rand; % opacity
    end
else
    % Generate initial states, spread evenly within the vessel mask
    inds = find(mask);
    for c = 1:n_cells
        ind = inds(ceil(rand*length(inds)));
        cell_pos_t(c, 1) = complex(xx(ind), yy(ind));
        cell_pos_t(c, 2) = rand; % opacity
    end
end
cell_pos_t(:, 3) = 1; % 'incarnation' number

% Create list of cells that are 'live'
b_cells = ones(n_cells, 1);

% Create new cells close to vessel end from now on
if ~isempty(path_points)
    p0 = path_points(1); 
    rr_squared = (xx-real(p0)).^2 + (yy-imag(p0)).^2;
    r_mask = (rr_squared <= widths(1)^2) & mask;
else
    r_mask = ones(size(flowmap));
end

p0_inds = find(r_mask);

% Find the maximum displacement between two frames
max_flow = max(abs(flowmap(:)));

% Split the position update over n_its iterations so that the maximum
% displacement applied in any one iteration is limited (e.g. to one pixel 
% per iteration).
n_its = ceil(max_flow);

% Scale flowmap accordingly
flowmap = flowmap / n_its;
n_recorded = 0;
n_generated = 0;
burn_in = inf;

tb = timebar('title','Burning in...','limit',n_frames);
while n_recorded < n_frames
    for c = 1:n_cells
        % If a cell is not in use then create a new one
        if (b_cells(c) == 0)
            ind = p0_inds(ceil(rand*length(p0_inds)));
            cell_pos_t(c, 1) = complex(xx(ind), yy(ind));
            cell_pos_t(c, 2) = rand - 0; % opacity
            cell_pos_t(c, 3) = cell_pos_t(c, 3) + 1;

            b_cells(c) = 1;
        end

        % Compute the displacement over n_its iterations
        for it = 1:n_its
            % Quick approximation to interp2()
            dx = abs((1:imsz(2)) - real(cell_pos_t(c, 1)));
            xinds = (dx < 1);
            xwts = 1 - dx(xinds);

            dy = abs((1:imsz(1))' - imag(cell_pos_t(c, 1)));
            yinds = (dy < 1);
            ywts = 1 - dy(yinds);

            uv = (ywts' * flowmap(yinds,xinds) * xwts');

            if ~isnan(uv)
                % Displace by uv
                cell_pos_t(c, 1) = cell_pos_t(c, 1) + uv;
            else
                % Kill this cell and create a new one
                b_cells(c) = 0;
                
                % Skip to the next cell
                break; 
            end
        end % for it
    end % for c

    % Apply displacement to positions.
    cell_pos_t(:, 1) = cell_pos_t(:, 1) + ...
                       sigma_f * complex(randn(n_cells,1), randn(n_cells,1));
        
    % Only start recording cell positions once every cells has reached the
    % 'sink' in the vessel circuit.
    n_generated = n_generated + 1;
    if (all(cell_pos_t(c, 3) > 1)) || ...
       (n_generated >= burn_in)
        n_recorded = n_recorded + 1;
        cell_positions(:, :, n_recorded) = cell_pos_t;
        
        timebar(tb, 'title','Generating positions', ...
                    'advance');
    end
end
timebar(tb, 'close');
drawnow;

cell_positions = [real(cell_positions(:,1,:)) ...
                  imag(cell_positions(:,1,:)) ...
                  cell_positions(:,2,:) ...
                  cell_positions(:,3,:) ];
              

%% Test script
function test_script()
clc;

s = [ 
   211418064
  2094606752
]; randn('state', s);

s = randn('state'); disp(uint32(s));

f_profile = false;

timebar('closeall');

if (f_profile), profile clear; profile on; end
    pts = generate_path(500);

    [widths, dirs, inner, outer] = ...
        generate_edges(pts, [], 10, 1.6);

    [pts, widths, inner, outer, imsz] = ...
        scale_vessel(pts, widths, inner, outer, 160);

    [flowmap, mask] = ...
        create_flowmap(pts, widths, dirs, imsz);
    
    max_flow = 3.0;
    flowmap = flowmap * max_flow;

    n_cells = 100;
    n_frames = 200;
    samples_per_exposure = 10;
    cell_positions = generate_cell_positions(flowmap, mask, ...
                                             pts, widths, ...
                                             n_cells, n_frames+1, ...
                                             [], ...
                                             samples_per_exposure);
if (f_profile), profile off; profile report; end

figure(1); clf;
    animate_cell_positions(cell_positions(:,:,1:samples_per_exposure:end), imsz);

return

figure(1); clf; hold on;
    plot(inner(:,1), inner(:,2), 'r-');
    plot(outer(:,1), outer(:,2), 'r-');
    h = plot(cell_positions(:,1,1), cell_positions(:,2,1), 'k.', ...
             'markersize', 24);
    axis('ij', 'equal', [xx(1,1) xx(1,end) yy(1,1) yy(end,1)]);

set(h, 'xdatasource','cell_positions(:,1,i)',...
       'ydatasource','cell_positions(:,2,i)');
   
for i = 2:n_frames
    refreshdata(h, 'caller');
    pause(0.01);
end

