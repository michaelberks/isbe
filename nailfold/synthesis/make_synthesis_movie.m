function make_synthesis_movie()

% base parameters
p = struct( ...
    'imgroot', 'U:/projects/nailfold/synthesis/showcase/', ...
    'n_points', 200, ...
    'sigma_x', 0.02, ...
    'sigma_y', 0.01, ...
    'w0', 2.0, ...
    'sigma_w', 0.005, ...
    'im_width', 160, ...
    'max_flow', 8.0, ... % in pixel/frame
    'cell_sz', 11, ...
    'n_cells', 800, ...
    'n_frames', 5000, ...
    'frame_skip', 1, ...
    'burn_in', 500, ...
    'sigma_f', 0.025, ...
    'jitter_sigma', 0.0, ...
    'brightness0', 240, ...
    'brightness_sigma', 0, ...
    'contrast0', 64, ...
    'contrast_sigma', 0, ...
    'bg_weight', 0, ...
    'mask_weight', 0, ...
    'sigma_n', 0, ...
    'f_export', false, ...
    'f_halfsize', false ...
);

imgpath = fullfile(p.imgroot, 'rendered/ideal');
cell_positions = func(p, imgpath);

p.contrast0 = 24;
p.brightness0 = 145;
imgpath = fullfile(p.imgroot, 'rendered/realbc');
func(p, imgpath, cell_positions);

p.bg_weight = 24;
imgpath = [imgpath, '+sb'];
func(p, imgpath, cell_positions);

p.mask_weight = 1.0;
imgpath = [imgpath, '+sc'];
func(p, imgpath, cell_positions);

p.brightness_sigma = 4;
imgpath = [imgpath, '+tb'];
func(p, imgpath, cell_positions);

p.contrast_sigma = 0.05;
imgpath = [imgpath, '+tc'];
func(p, imgpath, cell_positions);

p.jitter_sigma = 2;
imgpath = [imgpath, '+j'];
func(p, imgpath, cell_positions);

p.sigma_n = 4;
imgpath = [imgpath, '+n'];
func(p, imgpath, cell_positions);


%% Generate the images
function cell_positions = func(p, imgpath, cell_positions)
% Make a movie of the synthesis process.

s = [ 3145040330; 1496042031 ];
randn('state', s);
rand('twister', s(1)); 

pts = generate_path(p.n_points, p.sigma_x, p.sigma_y);

[widths, dirs, inner, outer] = generate_edges(pts, p.sigma_w, p.w0);

[flowmap, mask, xx, yy] = create_flowmap(pts, widths, dirs, ...
                                         [inner; outer], p.im_width);

save(fullfile(p.imgroot,'flowmap.mat'), ...
     'pts', 'widths', 'dirs', 'inner', 'outer', ...
     'flowmap', 'mask', 'xx', 'yy');
                                     
if (p.brightness0 > 200)
    make_path_movie(p, xx, yy, pts);
    make_edge_movie(p, xx, yy, pts, inner, outer);
end

% make_flow_movie(); % ???

units_per_pixel = xx(1,2) - xx(1,1);
flowmap = flowmap * units_per_pixel; % max(flowmap) = 1 pixel/frame
flowmap = flowmap * p.max_flow; % max(flowmap) = max_flow pixel/frame

p.sigma_f = p.sigma_f * p.max_flow;

if ~exist('cell_positions','var')
    [cell_positions] = generate_cell_positions(xx, yy, ...
                                               flowmap / p.frame_skip, ...
                                               mask, ...
                                               pts, widths, ...
                                               p.n_cells, ...
                                               p.n_frames * p.frame_skip, ...
                                               p.burn_in * p.frame_skip, ...
                                               p.sigma_f / p.frame_skip);
    cell_positions = cell_positions(:,:,1:p.frame_skip:end);
end

p.n_frames = size(cell_positions,3);
flowmap = flowmap / units_per_pixel; % flowmap now in pixels/frame

if (p.brightness0 > 200)
    make_cell_movie(p, xx, yy, inner, outer, cell_positions);
end

jitter = p.jitter_sigma * randn(p.n_frames+2, 2);
jitter = conv2(jitter, [1; 2; 1]/4, 'valid');

b_add = p.brightness_sigma * randn(p.n_frames+2, 1);
b_add = conv2(b_add, [1; 2; 1]/4, 'valid');

c_mult = 1 + (p.contrast_sigma * randn(p.n_frames+2, 1));
c_mult = conv2(c_mult, [1; 2; 1]/4, 'valid');

bg_sz = size(xx) + mod(size(xx),2);
cloud_bg = noiseonf(max(bg_sz), 1.5);
cloud_bg = p.bg_weight * ...
           normim(cloud_bg(1:size(xx,1), 1:size(xx,2)), 'stretch_fixed');

cloud_mult = noiseonf(max(bg_sz), 1.5);
cloud_mult = p.mask_weight * ...
             normim(cloud_mult(1:size(xx,1), 1:size(xx,2)), 'stretch_fixed');
cloud_mult = 1 - cloud_mult;

if ~exist(imgpath, 'dir')
    mkdir(imgpath);
else
    delete(fullfile(imgpath,'frame*.png'));
end

for i = 1:p.n_frames
    img0 = make_frame(xx, yy, cell_positions(:,:,i), ...
                      p.cell_sz, mask, ...
                      (p.brightness0 + b_add(i)) - cloud_bg, ...
                      (p.contrast0 * c_mult(i)) * cloud_mult, ...
                      p.sigma_n);
                  
    % Consider adding jitter in make_frame to ensure that the noise SD is
    % correct
	if (p.jitter_sigma > 0)
        [x2,y2] = meshgrid(1:size(xx,2), 1:size(xx,1));
        x2 = x2 + jitter(i,1);
        y2 = y2 + jitter(i,2);

        img = interp2(img0, x2, y2);
    else
        img = img0;
    end
    
    % Rotate the image
    img = img(end:-1:1, end:-1:1);
    
    % Write image to disk
    filename = sprintf('frame_%04d.png', i);
    imwrite(uint8(img), fullfile(imgpath, filename));
    
    if (p.contrast_sigma > 0)
        filename = sprintf('sc_%04d.png', i);
        imwrite(uint8(12*(p.contrast0*c_mult(i)) * cloud_mult), fullfile(imgpath, filename));
    elseif (p.brightness_sigma > 0)
        filename = sprintf('sb_%04d.png', i);
        imwrite(uint8((p.brightness0+b_add(i))-cloud_bg), fullfile(imgpath, filename));
    elseif (p.mask_weight > 0)
        filename = sprintf('sc.png');
        imwrite(uint8(12*(p.contrast0*c_mult(i)) * cloud_mult), fullfile(imgpath, filename));
    elseif (p.bg_weight > 0)
        filename = sprintf('sb.png');
        imwrite(uint8((p.brightness0+b_add(i))-cloud_bg), fullfile(imgpath, filename));
    end
end


%% Scale a colour toward white
function sc = scaled_colour(base_colour, lightness)
sc = base_colour + lightness*(1-base_colour);


%% Save the PNG file at a given width
function savepng(filename, width)
set(gca,'position', [0,0,1,1]);

ar = get(gca,'PlotBoxAspectRatio');
ar = ar(2)/ar(1);

w = width;
set(gcf, 'paperposition', [0,0, w,ar*w], ...
         'position', [1000,200, w,ar*w]);

res = sprintf('-r%i', round(72));
print('-dpng', res, filename);


%% Generate frames of the path being created
function make_path_movie(p, xx, yy, pts)

clc;

% Rotate through 180 deg
pts = -pts;
xx = -xx(end:-1:1, end:-1:1);
yy = -yy(end:-1:1, end:-1:1);

for m = 1:8
    imgpath = fullfile(p.imgroot, sprintf('generate_path/movie%02d', m));
    if ~exist(imgpath, 'dir')
        mkdir(imgpath);
    else
        delete(fullfile(imgpath, 'frame*.*'));
    end

    f = 1;
    scl = 0;
    for i = p.n_points:-8:1
        figure(1); cla;
            inds = i:size(pts,1)-i+1;

            plot(pts(inds,1), pts(inds,2), '.-', ...
                 'color', scaled_colour([0,0,1], scl), ...
                 'markersize', 8);
            axis('ij','equal');
            xlim(xx(1,[1,end]));
            ylim(yy([1,end],1));
            set(gca,'xtick',[], 'ytick',[], 'box','on');
        drawnow;

        filename = sprintf('frame_%04d.png', f);
        savepng(fullfile(imgpath, filename), size(xx,2));

        f = f + 1;
    end
    
    scl = linspace(0, 1, 25);
    for i = 1:length(scl)
        figure(1); cla;
            plot(pts(inds,1), pts(inds,2), '.-', ...
                 'color', scaled_colour([0,0,1], scl(i)), ...
                 'markersize', 8);
            axis('ij','equal');
            xlim(xx(1,[1,end]));
            ylim(yy([1,end],1));
            set(gca,'xtick',[], 'ytick',[], 'box','on');
        drawnow;

        filename = sprintf('frame_%04d.png', f);
        savepng(fullfile(imgpath, filename), size(xx,2));

        f = f + 1;
    end

    % Now generate new (random) points and repeate
    pts = generate_path(p.n_points, p.sigma_x, p.sigma_y);
    
    % Rotate through 180 deg
    pts = -pts;
end


%% Generate frames of the edges being created
function make_edge_movie(p, xx, yy, pts, inner, outer)

clc;

% Rotate through 180 deg
pts = -pts;
inner = -inner;
outer = -outer;

xx = -xx(end:-1:1, end:-1:1);
yy = -yy(end:-1:1, end:-1:1);

for m = 1:8
    imgpath = fullfile(p.imgroot, sprintf('generate_edges/movie%02d', m));
    if ~exist(imgpath, 'dir')
        mkdir(imgpath);
    else
        delete(fullfile(imgpath, 'frame*.*'));
    end

    f = 1;
    scl = 0;
    for i = p.n_points:-8:1
        figure(1); cla; hold on;
            inds = i:size(pts,1)-i+1;

            plot(pts(:,1), pts(:,2), '-', ...
                 'color', scaled_colour(0.8*[1,1,1], scl));
            plot(inner(inds,1), inner(inds,2), '.-', ...
                 'color', scaled_colour([0,0,0], scl), ...
                 'markersize', 8);
            plot(outer(inds,1), outer(inds,2), '.-', ...
                 'color', scaled_colour([0,0,0], scl), ...
                 'markersize', 8);
                        axis('ij','equal');
            xlim(xx(1,[1,end]));
            ylim(yy([1,end],1));
            set(gca,'xtick',[], 'ytick',[], 'box','on');
        drawnow;

        filename = sprintf('frame_%04d.png', f);
        savepng(fullfile(imgpath, filename), size(xx,2));

        f = f + 1;
    end
    
    % Fade out the image
    scl = linspace(0, 1, 25);
    for i = 1:length(scl)
        figure(1); cla; hold on;
            plot(pts(:,1), pts(:,2), '-', ...
                 'color', scaled_colour(0.8*[1,1,1], scl(i)));
            plot(inner(inds,1), inner(inds,2), '.-', ...
                 'color', scaled_colour([0,0,0], scl(i)), ...
                 'markersize', 8);
            plot(outer(inds,1), outer(inds,2), '.-', ...
                 'color', scaled_colour([0,0,0], scl(i)), ...
                 'markersize', 8);
            axis('ij','equal');
            xlim(xx(1,[1,end]));
            ylim(yy([1,end],1));
            set(gca,'xtick',[], 'ytick',[], 'box','on');
        drawnow;

        filename = sprintf('frame_%04d.png', f);
        savepng(fullfile(imgpath, filename), size(xx,2));

        f = f + 1;
    end
    
    % Now generate new (random) points and repeate
    [widths, dirs, inner, outer] = generate_edges(-pts, p.sigma_w, p.w0);

    % Rotate through 180 deg
    inner = -inner;
    outer = -outer;
end


%% Generate schematics of cells in motion
function make_cell_movie(p, xx, yy, inner, outer, cell_positions)
clc;

% Rotate through 180 deg
inner = -inner;
outer = -outer;
cell_positions = -cell_positions;

xx = -xx(end:-1:1, end:-1:1);
yy = -yy(end:-1:1, end:-1:1);

m = 1;
imgpath = fullfile(p.imgroot, sprintf('generate_cells/movie%02d', m));
if ~exist(imgpath, 'dir')
    mkdir(imgpath);
else
    delete(fullfile(imgpath, 'frame*.*'));
end

f = 1;
nf = size(cell_positions,3);
figure(1); clf;
    set(1, 'position', [1000,200, size(xx,2),size(xx,1)]);
for i = 1:nf
    figure(1); cla; hold on;
        plot(cell_positions(:,1,f), cell_positions(:,2,f), 'r.', ...
             'markersize',24);
        plot(inner(:,1), inner(:,2), 'k-');
        plot(outer(:,1), outer(:,2), 'k-');
        
        axis('ij','equal');
        xlim(xx(1,[1,end]));
        ylim(yy([1,end],1));
        set(gca,'xtick',[], 'ytick',[], 'box','on');
    drawnow;

    filename = sprintf('frame_%04d.png', f);
    savepng(fullfile(imgpath, filename), size(xx,2));

    f = f + 1;
end
