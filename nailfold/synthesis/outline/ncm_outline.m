function pts_out = ncm_outline(imgname, n_wanted)
if (nargin==0 && nargout==0), test(); return; end

if ~exist('n_wanted','var'), n_wanted = []; end

img = imread(imgname);
img = mean(img,3);

[r,c] = find(img == 0);

n_pts = length(r);
pts_in = [c, r];
pts_out = [];

pts_out(end+1,:) = pts_in(1,:);
pts_in(1,:) = [];

% Sort points into a sensible order
while ~isempty(pts_in)
    d_squared = sum((pts_in - ones(size(pts_in,1),1)*pts_out(end,:)).^2, 2);
    [shortest_d, nearest_ind] = min(d_squared);
    
    pts_out(end+1,:) = pts_in(nearest_ind,:);
    pts_in(nearest_ind,:) = [];
end

% Subsample
if ~isempty(n_wanted)
    inds = ceil(linspace(1, n_pts, n_wanted));
    pts_out = pts_out(inds,:);
end

return


function test()
clc;

imgpath = 'S:\projects\nailfold\synthesis';
imgname = fullfile(imgpath,'m_outline.png');

pp = ncm_outline(imgname);

[widths, dirs, inner, outer] = ...
    generate_edges(pp, 0.005, 20, 1.0);

figure(1); clf; hold on;
    plot(pp(:,1), pp(:,2), 'b.-');
    plot(pp(1,1), pp(1,2), 'go');
    plot(pp(end,1), pp(end,2), 'ro');
    plot(inner(:,1), inner(:,2), 'r.-');
    plot(outer(:,1), outer(:,2), 'g.-');
    axis('equal','ij');
    
[pp, widths, inner, outer, imsz] = ...
    scale_vessel(pp, widths, inner, outer, 400);

[flowmap, mask] = ...
    create_flowmap(pp, widths, dirs, imsz);

figure(2);
    show_flow_as('rgb', flowmap);

max_flow = 10.0;
flowmap = flowmap * max_flow;

% Choose cell size and number of cells
cell_sz = min(widths) * 0.8;
n_pixels_per_cell = cell_sz^2;
n_mask_pixels = sum(mask(:)>0);
n = 8;
n_cells = ceil(n * n_mask_pixels / n_pixels_per_cell);

n_frames = 200;
samples_per_frame = 1;
cell_positions = generate_cell_positions(flowmap, mask, ...
                                         pp, widths, ...
                                         n_cells, n_frames+1, ...
                                         [], samples_per_frame);
    
return

