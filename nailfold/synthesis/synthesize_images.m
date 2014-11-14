function imgpath = synthesize_images(varargin)

if (nargin==0 && nargout==0), test_script(); return; end

p = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    ...
    'imgroot', 'U:/projects/nailfold/synthesis', ...
    'path_points', [], ...
    'flowmap', [], ...
    ...
    'n_points', 300, ...
    'sigma_x', 0.0, ...
    'roundness', 1.0, ...
    'w0', 2.0, ...
    'width_ratio', 1.0, ...
    'sigma_w', 0.0, ...
    'im_width', 160, ...
    'trim_borders', [], ...
    'trim_top', [], ...
    'max_flow', 1.0, ... % in pixel/frame
    'n_frames', 180, ...
    'n_cells', [], ...
    'cell_sz', [], ...
    'frames_per_exposure', 1, ...
    'n_exposed', inf, ...
    'burn_in', 50, ...
    'sigma_f', 0, ...
    'jitter_sigma', 0.0, ...
    'brightness0', 200, ...
    'brightness_sigma', 0, ...
    'contrast0', 64, ...
    'contrast_sigma', 0, ...
    'bg_weight', 0, ...
    'mask_weight', 0, ...
    'sigma_n', 0, ...
    'f_export', false ...
);

clc;

% s = [ 211418064; 2094606752 ];
% s = [ 3145040330; 1496042031 ];
s = [ 3556622558; 3168652225 ]; % Good hairpin

randn('state', s);
rand('twister', s(1));

% Limit n_exposed to range [1,frames_per_exposure]
p.n_exposed = min(max(p.n_exposed, 1), p.frames_per_exposure);



if isempty(p.flowmap) && isempty(p.path_points)
    fprintf('Generating vessel centreline...');
    pts = generate_path(p.n_points, p.sigma_x, p.roundness);
    fprintf('done\n');

    fprintf('Generating vessel edges...');
    [widths, dirs, inner, outer] = ...
        generate_edges(pts, p.sigma_w, p.w0, p.width_ratio);
    fprintf('done\n');
else
    pts = p.path_points;
end

% Convert from arbitrary units to pixels
[pts, widths, inner, outer, imsz] = ...
    scale_vessel(pts, widths, inner, outer, p.im_width);

% If no cell size specified then set it to 
% 0.8 times the width of the vessel at its narrowest point.
if isempty(p.cell_sz)
    p.cell_sz = min(widths) * 0.8;
end
% Trim 20% of the height off the top of the image if none defined
if isempty(p.trim_top)
    p.trim_top = 8 * ceil(0.2 * imsz(1) / 8);
end
% Trim 5% of the width off every side.
if isempty(p.trim_borders)
    p.trim_borders = 8 * ceil(0.05 * imsz(2) / 8);
end

if isempty(p.flowmap)
    fprintf('Creating flow map...');
    % Flowmap is in pixels here. Max flow = 1 pixel/frame.
    [flowmap, mask] = create_flowmap(pts, widths, dirs, imsz);
    % Scale flowmap as required
    flowmap = flowmap * p.max_flow; % max(flowmap) = max_flow pixel/frame
    fprintf('done\n');
else
    pts = [];
    widths = [];
    flowmap = p.flowmap;
    imsz = size(p.flowmap);
    mask = ones(imsz);
end

% Choose a number of cells based on the size of the proposed vessel space
if isempty(p.n_cells)
    n_mask_pixels = sum(mask(:)>0);
    n_pixels_per_cell = p.cell_sz^2;
    
    % Use a density of n cells per space
    n = 8;
    p.n_cells = ceil(n * n_mask_pixels / n_pixels_per_cell);
end

% Create an extra frame so that we have displacements for all frames
% (sigma_f is scaled to be in proportion to the max_flow)
fprintf('Generating cell positions...');
sigma_f = p.sigma_f*p.max_flow;
cell_positions = generate_cell_positions(flowmap, mask, ...
                                         pts, widths, ...
                                         p.n_cells, p.n_frames+1, ...
                                         sigma_f, ...
                                         p.frames_per_exposure);
fprintf('done\n');

imgroot = p.imgroot;
imgpath = fullfile(imgroot, datestr(now, 'yyyymmddTHHMMSS'));
if ~exist(imgpath, 'dir')
    mkdir(imgpath);
end

%% Export figures if requested
if p.f_export
    export_figures(pts, inner, outer, flowmap, cell_positions, imgpath);
end

%% Generate frames
fprintf('Creating images...');
% Create camera 'jitter'
jitter = p.jitter_sigma * randn(p.n_frames+2, 2);
% Smooth it so it's not quite so random and normalize so that displacement
% at frame 1 is zero (as a reference).
jitter = conv2(jitter, [1; 2; 1]/4, 'valid');
jitter = jitter - jitter(ones(p.n_frames,1),:);

% jitter = zeros(size(jitter));
% % jitter(2:end,1) = 2.0 * (2:p.n_frames) / p.n_frames;
% % jitter(:,2) = jitter(:,1);
% jitter(2:end, :) = 0.5;

% Modulate brightness randomly
b_add = p.brightness_sigma * randn(p.n_frames+2, 1);
b_add = conv2(b_add, [1; 2; 1]/4, 'valid');

% Modulate contrast randomly
c_mult = 1 + (p.contrast_sigma * randn(p.n_frames+2, 1));
c_mult = conv2(c_mult, [1; 2; 1]/4, 'valid');

% Create spatially correlated 'cloud noise' image (from Peter Kovesi)
cloud_bg = noiseonf(max(imsz), 1.5);
cloud_bg = p.bg_weight * ...
           normim(cloud_bg(1:imsz(1),1:imsz(2)), 'stretch_fixed');

% Multiply the image by a 'mottled' texture
% p.mask_weight = 0 means no weighting
cloud_mult = noiseonf(max(imsz), 1.5);
cloud_mult = p.mask_weight * ...
             normim(cloud_mult(1:imsz(1),1:imsz(2)), 'stretch_fixed');
cloud_mult = 1 - cloud_mult;
       
% This is written in such a way that brightness is approximately the value
% of the brightest pixel in the image. Everything is subtracted from this
% maximum value.
flowStack = complex( zeros([imsz, p.n_frames]), ...
                     zeros([imsz, p.n_frames]) );

flowStackSum = flowStack(:,:,1);
flowStackCount = real(flowStackSum);

tb = p.trim_borders;
tt = p.trim_top;

f_rng = 1:1+p.frames_per_exposure;
for i = 1:p.n_frames
    brightness = (p.brightness0 + b_add(i)) - cloud_bg;
    contrast = (p.contrast0 * c_mult(i)) * cloud_mult;
    
    [img0, flow0] = make_frame( cell_positions(:,:,f_rng), ...
                                imsz, p.cell_sz, mask, ...
                                brightness, contrast, ...
                                p.n_exposed);
    f_rng = f_rng + p.frames_per_exposure;
    
	if (i > 1)
        [x2,y2] = meshgrid(1:imsz(2), 1:imsz(1));
        % New image is a copy of the first, displaced by jitter(i,:)
        % Therefore, we need to sample points at (x,y) - jitter(i,:) rather
        % than (x,y) + jitter(i,:).
        x2 = x2 - jitter(i,1);
        y2 = y2 - jitter(i,2);

        img = interp2(img0, x2, y2);
        flowimg = cat(3, interp2(flow0(:,:,1), x2, y2, '*linear'), ...
                         interp2(flow0(:,:,2), x2, y2, '*linear'));
    else
        img = img0;
        flowimg = flow0;
    end

    % Add noise to the image
    img = img + p.sigma_n*randn(size(img));
    
    % Trim image
    if (tt > 0)
        img = img(tt+1:end,:);
    end
    if (tb > 0)
        img = img(1+tb:end-tb,1+tb:end-tb);
    end
    
    % Add to big stack.
    flowStackSum = flowStackSum + flowimg(:,:,1);
    flowStackCount = flowStackCount + flowimg(:,:,2);
    flowStack(:,:,i) = flowimg(:,:,1) ./ flowimg(:,:,2);
    
    % Write image to disk
    filename = sprintf('frame_%04d.png', i);
    imwrite(uint8(img), fullfile(imgpath, filename));
end

flowStackMean = flowStackSum ./ flowStackCount;

% Average cell positions over exposed frames
% FIXME: Doesn't yet account for death/rebirth of cells
cp_new = cell_positions(:,:,1:p.frames_per_exposure:end);
for i = 2:p.n_exposed
    cp_new = cp_new + cell_positions(:,:,i:p.frames_per_exposure:end);
end
cell_positions = cp_new(:,:,1:p.n_frames) / p.n_exposed;

% Trim ground truth values
if (tt > 0)
    flowmap = flowmap(1+tt:end,:);
    flowStack = flowStack(1+tt:end,:,:);
    flowStackMean = flowStackMean(1+tt:end,:);
    
    mask = mask(1+tt:end,:);
    
    pts(:,2) = pts(:,2) - tt;
    cell_positions(:,2,:) = cell_positions(:,2,:) - tt;
end
if (tb > 0)
    flowmap = flowmap(1+tb:end-tb,1+tb:end-tb);
    flowStack = flowStack(1+tb:end-tb,1+tb:end-tb,:);
    flowStackMean = flowStackMean(1+tb:end-tb,1+tb:end-tb);
    
    mask = mask(1+tb:end-tb,1+tb:end-tb);
    
    pts = pts - tb;
    cell_positions = cell_positions - tb;
end

parameters = p;

save(fullfile(imgpath,'_ground_truth.mat'), ...
     'parameters', ...
     'pts', 'widths', 'cell_positions', ...
     'flowmap', 'mask', ...
     'flowStack', 'flowStackMean', ...
     'jitter', 'b_add', 'c_mult');
 
fid = fopen(fullfile(imgpath,'_parameters.txt'), 'w');
if (fid ~= 0)
    fprintf(fid,'%s\n',evalc('parameters'));
    fclose(fid);
end

fprintf('done\n');


%% Save figures to disk
function export_figures(pts, inner, outer, flowmap, cell_positions, ...
                        imgpath)

imsz = size(flowmap);
figure(1); clf; hold on;
    plot(pts(:,1), pts(:,2), 'b-');
    axis('equal', [1,imsz(2),1,imsz(1)]);
    set(gca,'box','on','xticklabel',[],'yticklabel',[]);
    exportfig(fullfile(imgpath,'path_points'));

    plot(inner(:,1), inner(:,2), 'r-');
    plot(outer(:,1), outer(:,2), 'r-');
    axis('equal', [1,imsz(2),1,imsz(1)]);
    exportfig(fullfile(imgpath,'edge_points'));

% 	figure(2); clf;
%         cmap = [1 1 1; redgreen(255)];
%         uu = normim(-real(flowmap(end:-1:1,end:-1:1)),'stretch_fixed');
%         subplot(1,2,1); image(uint8(1+uu*254)); axis('xy'); colormap(cmap);
%         imwrite(uint8(1+uu*254), cmap, fullfile(imgpath,'u_gt.png'));
%         vv = normim(-imag(flowmap(end:-1:1,end:-1:1)),'stretch_fixed');
%         subplot(1,2,2); image(uint8(1+vv*254)); axis('xy'); colormap(cmap);
%         imwrite(uint8(1+vv*254), cmap, fullfile(imgpath,'v_gt.png'));

flow_rgb = show_flow_as('rgb', flowmap);
imwrite(uint8(flow_rgb), fullfile(imgpath,'flow_gt.png'));
    
figure(2); clf; hold on;
    plot(cell_positions(:,1,end), cell_positions(:,2,end), 'r.', ...
         'markersize', 16);
    plot(inner(:,1), inner(:,2), 'b-');
    plot(outer(:,1), outer(:,2), 'b-');
    axis('equal', [1,imsz(2),1,imsz(1)]);
    set(gca,'box','on','xticklabel',[],'yticklabel',[]);
    exportfig(fullfile(imgpath,'cells_positions'));
    

%% Test script
function test_script()
clc;

image_type = 'test';

f_profile = false;

if f_profile, profile clear; profile on; end

switch image_type
    case 'ideal',
        % Ideal image parameters
        imgpath = synthesize_images( ...
            'n_frames', 180, 'frames_per_exposure', 4, ...
            'sigma_x', 0.00, 'roundness', 1.0, ...
            'w0', 2.0, 'sigma_w', 0.00, ...
            'n_cells', 640, 'cell_sz', 11, ...
            'max_flow', 1.0, 'sigma_f', 0.0, ...
            'brightness0', 200, 'brightness_sigma', 0, ...
            'contrast0', 64, 'contrast_sigma', 0.0, ...
            'jitter_sigma', 0.0, ...
            'bg_weight', 0, 'mask_weight', 0, 'sigma_n', 0, ...
            'f_export', false, ...
            'n_points', 200 ...
        );
    
    case 'test',
        % Ideal image parameters
        imgpath = synthesize_images( ...
            'im_width', 128, 'trim_borders', [], ...
            'n_cells', [], 'cell_sz', [], ...
            'n_frames', 60, 'frames_per_exposure', 8, 'n_exposed', inf, ...
            'sigma_x', 0.015, 'roundness', 1.0, ...
            'w0', 12.0, 'sigma_w', 0.005, 'width_ratio', 2.0, ...
            'max_flow', 4.0, 'sigma_f', 0.0, ...
            'brightness0', 140, 'brightness_sigma', 0, 'bg_weight', 0, ...
            'contrast0', 24, 'contrast_sigma', 0.0, 'mask_weight', 0.0, ...
            'jitter_sigma', 0.0, 'sigma_n', 4.0, ...
            'f_export', false, ...
            'n_points', 120 ...
        );

    case 'real',
        % Real image parameters
        imgpath = synthesize_images( ...
            'n_frames', 180, 'frames_per_exposure', 4, ...
            'sigma_x', 0.02, 'roundness', 1.0, ...
            'w0', 2.0, 'sigma_w', 0.0005, ...
            'n_cells', 1600, 'cell_sz', 11, ...
            'max_flow', 8.0, 'sigma_f', 0.05, ...
            'brightness0', 145, 'brightness_sigma', 0, 'bg_weight', 8, ...
            'contrast0', 24, 'contrast_sigma', 0.0, 'mask_weight', 1.0, ...
            'jitter_sigma', 1, 'sigma_n', 4, ...
            'f_export', false, ...
            'n_points', 200 ...
        );
end

if f_profile, profile off; profile report; end

d = dir(fullfile(imgpath,'*.png'));
figure(1); clf; hold off; colormap(gray(256));
for i = 1:length(d)
    img = imread(fullfile(imgpath, d(i).name));
    image(uint8(img));
    axis('image');
    pause(1/30);
end
