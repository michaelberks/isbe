function make_flow_movie()
clc; clear;

imgpath = 'U:\projects\nailfold\synthesis\showcase';
filename = fullfile(imgpath, 'flowmap.mat');
load(filename);

pts = -pts;
dirs = -dirs;
inner = -inner;
outer = -outer;

% flowmap = create_flowmap(pts, widths, dirs, [inner; outer], 160);
% return

create_flowmap_movie(pts, widths, dirs, [inner; outer], imgpath);
    
    
function create_flowmap_movie(path_points, widths, dirs, edge_points, imgroot)
        
map_width = 160;

% Compute border size
ymin = min(edge_points(:,2));
ymax = max(edge_points(:,2));
yrng = ymax - ymin;
border = 0.1 * yrng;

% Add border to limits
ymin = ymin - border;
ymax = ymax + border;
xmin = min(edge_points(:,1)) - border;
xmax = max(edge_points(:,1)) + border;

n_points = size(path_points,1);

% Now generate a pixelmap
gridstep = (xmax-xmin) / (map_width-1);
[xx,yy] = meshgrid(xmin:gridstep:xmax, ...
                   ymin:gridstep:ymax);

sz = size(xx) - mod(size(xx),2);
xx = xx(1:sz(1), 1:sz(2));
yy = yy(1:sz(1), 1:sz(2));

mask = zeros(size(xx));
flowmap = complex(nan(size(xx)), nan(size(xx)));

max_flow = 0;

% Get nearest path point
nearest_point = zeros(size(xx));
dist_mat = nan(size(xx));
for i = 1:size(xx,2)
    for j = 1:size(xx,1)
        diff = [path_points(:,1)-xx(j,i) path_points(:,2)-yy(j,i)];
        d = sum(diff.*diff, 2);

        [dist_squared, ind_nearest] = min(d);
        nearest_point(j,i) = ind_nearest;
        dist_mat(j,i) = dist_squared;
    end
end
dist_mat = sqrt(dist_mat);
width_mat = widths(nearest_point);

% Generate flow direction for all points
profile_mat = ones(size(dist_mat));
flowdir = nan(size(dist_mat));
flowmag = nan(size(dist_mat));

for i = 1:size(xx,2)
    for j = 1:size(xx,1)
        ind_nearest = nearest_point(j,i);
        dist_nearest = dist_mat(j,i);
        w = width_mat(j,i);
        
        is_inside_vessel = (dist_nearest < w);
        is_just_outside_vessel = (dist_nearest < 2*w);
        
        if is_inside_vessel
            profile_mat(j,i) = cos(dist_nearest / w); % Cosine profile
            
            % Make darkness proportional to thickness of the vessel
            % i.e. Enhance contrast of thicker limb of the vessel
%             mask(j,i) = w * mask(j,i);

            flowdir(j,i) = complex(dirs(ind_nearest,1), ...
                                   dirs(ind_nearest,2));
        elseif is_just_outside_vessel
            vec = path_points(ind_nearest,:) - [xx(j,i) yy(j,i)];
            vec = vec / norm(vec);
            
            flowdir(j,i) = complex(dirs(ind_nearest,1) + vec(1), ...
                                   dirs(ind_nearest,2) + vec(2));
            flowdir(j,i) = flowdir(j,i) / abs(flowdir(j,i)); % normalize vector
        else
            % Do nothing - no flow, no mask
            continue;
        end
        
        % Create a 'sink' around the endpoint where flow is zero.
        % Kill cells when they reach this part of the image.
        if (ind_nearest == n_points)
            flowdir(j,i) = complex(0,0);
            flowmag(j,i) = 0;
        else
            % Flow should be inversely proportional to w^2
            flowmag(j,i) = 1 / (w*w);
        end
    end
end

bgval = 1;
rgb = complex2rgb(flowdir, [],[],[], bgval);

nBlend = 100;

f = 1;
f_export = true;

imgpath = fullfile(imgroot, 'create_flowmap/interior');
if ~exist(imgpath, 'dir')
    mkdir(imgpath)
else
    delete(fullfile(imgpath,'*.png'));
end

interior_mask = (dist_mat < width_mat);
for p = unique([1, 1:4:n_points, n_points])
    mask = (interior_mask & nearest_point<=p);
    
    rgb2 = rgb;
    rgb2(repmat(~mask, [1,1,3])) = bgval;
    
    if f_export
        filename = sprintf('frame_%04d.png', f);
        imwrite(rgb2, fullfile(imgpath, filename));
        f = f + 1;
    else
        figure(1); clf;
            image(rgb2);
            axis('image','ij');
        drawnow;    
    end
end
rgb = rgb2;

f = 1;

imgpath = fullfile(imgroot, 'create_flowmap/feather');
if ~exist(imgpath, 'dir')
    mkdir(imgpath)
else
    delete(fullfile(imgpath,'*.png'));
end

% Blend in the vessel profile map
rgbw = complex2rgb(flowdir .* profile_mat, [],[],[], bgval);
for i = 1:nBlend
    alpha = i/nBlend;
    rgb2 = (1-alpha) * rgb + ...
              alpha  * rgbw;
    rgb2(repmat(~mask, [1,1,3])) = bgval;
    
    if f_export
        filename = sprintf('frame_%04d.png', f);
        imwrite(rgb2, fullfile(imgpath, filename));
        f = f + 1;
    else
        figure(1); clf;
            image(rgb2);
            axis('image','ij');
        drawnow;    
    end
end
rgb = rgbw;

f = 1;

imgpath = fullfile(imgroot, 'create_flowmap/exterior');
if ~exist(imgpath, 'dir')
    mkdir(imgpath)
else
    delete(fullfile(imgpath,'*.png'));
end

exterior_mask = (dist_mat < 2*width_mat);
for p = unique([1, 1:4:n_points, n_points])
    mask = (interior_mask | (exterior_mask & nearest_point<=p));
    
    rgb2 = rgb;
    rgb2(repmat(~mask, [1,1,3])) = bgval;
    
    if f_export
        filename = sprintf('frame_%04d.png', f);
        imwrite(rgb2, fullfile(imgpath, filename));
        f = f + 1;
    else
        figure(1); clf;
            image(rgb2);
            axis('image','ij');
        drawnow;    
    end
end
rgb = rgb2;

f = 1;

imgpath = fullfile(imgroot, 'create_flowmap/weight');
if ~exist(imgpath, 'dir')
    mkdir(imgpath)
else
    delete(fullfile(imgpath,'*.png'));
end

% Blend in the flow magnitude map
rgbw = complex2rgb(flowdir .* profile_mat .* flowmag, [],[],[], bgval);
for i = 1:nBlend
    alpha = i/nBlend;
    rgb2 = (1-alpha) * rgb + ...
              alpha  * rgbw;
    rgb2(repmat(~mask, [1,1,3])) = bgval;
    
    if f_export
        filename = sprintf('frame_%04d.png', f);
        imwrite(rgb2, fullfile(imgpath, filename));
        f = f + 1;
    else
        figure(1); clf;
            image(rgb2);
            axis('image','ij');
        drawnow;    
    end
end
