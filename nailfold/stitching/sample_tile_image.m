function [tgt_tiles, tgt_mask] = sample_tile_image(src_tiles, weights, ...
                                                   t2s_transform, ...
                                                   tgt_sz, src_mask)
if (nargin==0 && nargout==0), test(); return; end

% Deal with missing parameters
if ~exist('src_mask','var'), src_mask = []; end

if isempty(weights)
    weights = ones(size(src_tiles{1}));
end

% Return a cell array of sampled images, given a cell array of input images
% and their location in the larger mosaic (specified by the compound
% transforms).

tgt_height = tgt_sz(1);
tgt_width = tgt_sz(2);

% Destination sample points
[x_tgt, y_tgt] = meshgrid(1:tgt_width, 1:tgt_height);

if isempty(src_mask)
    src_mask = true(tgt_sz);
end
xy_tgt = [ x_tgt(src_mask) y_tgt(src_mask) ];


% Source sample points
xy_src = (t2s_transform * [xy_tgt ones(size(xy_tgt,1),1)]')';

if 1
    tgt_mask = zeros(tgt_sz);
    tgt_mask(src_mask) = interp2(weights, xy_src(:,1), xy_src(:,2), '*linear');
    
    tgt_tiles = cell(size(src_tiles));
    for i = 1:numel(src_tiles)
        tgt_tiles{i} = zeros(tgt_sz);
        tgt_tiles{i}(src_mask) = interp2(src_tiles{i}, xy_src(:,1), xy_src(:,2), ...
                                  '*linear');
        tgt_tiles{i} = tgt_mask .* tgt_tiles{i};
    end
    
    return
end

% Limits of source image samples
src_sz = size(src_tiles{1});
src_width = src_sz(2);
src_height = src_sz(1);

if isempty(weights)
    weights = ones(src_sz);
end

% Find those that actually map to a point in the source image
valid_pts = (1 <= xy_src(:,1)) & (xy_src(:,1) <= src_width) & ...
            (1 <= xy_src(:,2)) & (xy_src(:,2) <= src_height);

x0 = floor(xy_src(valid_pts,1));
x1 = ceil(xy_src(valid_pts,1));
y0 = floor(xy_src(valid_pts,2));
y1 = ceil(xy_src(valid_pts,2));

% Compute the interpolation weights for each point
x1_weight = xy_src(valid_pts,1) - x0;
y1_weight = xy_src(valid_pts,2) - y0;
x0_weight = 1 - x1_weight;
y0_weight = 1 - y1_weight;

%Compute the tile weights multplied by the interpolation weights
% if isempty(src_mask)
%     curr_weight = weights(valid_pts);
% else
%     curr_weight = weights(src_mask(valid_pts));
% end

% Emulate ind2sub()
inds_11 = ((x1-1) * src_sz(1)) + y1;
inds_01 = ((x0-1) * src_sz(1)) + y1;
inds_10 = ((x1-1) * src_sz(1)) + y0;
inds_00 = ((x0-1) * src_sz(1)) + y0;

% Assign weights
pixel_weights = zeros(src_sz(1)*src_sz(2), 4);
pixel_weights(inds_11, 1) = x1_weight .* y1_weight;
pixel_weights(inds_01, 2) = x0_weight .* y1_weight;
pixel_weights(inds_10, 3) = x1_weight .* y0_weight;
pixel_weights(inds_00, 4) = x0_weight .* y0_weight;

x_tgt = xy_tgt(valid_pts, 1);
y_tgt = xy_tgt(valid_pts, 2);

% Add pixel (later removed) for interpolation.
tgt_width = tgt_width + 1;
tgt_height = tgt_height + 1;

for i = 1:numel(src_tiles)
    % Add the tile intensities * interpolation weights to the mosaic

    tile = double(src_tiles{i});
    if isempty(src_mask)
        src_pixels = tile(:);
    else
        src_pixels = tile(src_mask(src_inds));
    end
        
    tgt_pixels = src_pixels .* sum(pixel_weights, 2);

    tgt_tiles{i} = ...
        sparse(y_tgt, x_tgt, tgt_pixels(valid_pts), tgt_height, tgt_width);
        
%     tgt_tiles{i} = ...
%         sparse(y_tgt, x_tgt, tile .* weights_00, tgt_height, tgt_width) + ...
%         sparse(y_tgt, x_tgt, tile .* weights_01, tgt_height, tgt_width) + ...
%         sparse(y_tgt, x_tgt, tile .* weights_10, tgt_height, tgt_width) + ...
%         sparse(y_tgt, x_tgt, tile .* weights_11, tgt_height, tgt_width);
    
    % Remove pixel that was added earlier
    tgt_tiles{i} = full(tgt_tiles{i}(1:end-1, 1:end-1));
end

if (nargout > 1)
    tgt_mask = ...
        sparse(y_pix, x_pix,     weights_00, height, width) + ...
        sparse(y_pix, x_pix+1,   weights_01, height, width) + ...
        sparse(y_pix+1, x_pix,   weights_10, height, width) + ...
        sparse(y_pix+1, x_pix+1, weights_11, height, width);

    % Remove pixel that was added earlier
    tgt_mask = full(tgt_mask(1:end-1, 1:end-1));
end


%% Test script
function test()
clc;

% Get built-in test image
figure(1); colormap(gray(64));
    subplot(2,2,1);
    image(); axis('image', 'ij', 'off');
    img = get(get(gca, 'children'), 'cdata');
    img = {img{1}};
    
theta = 15 * pi/180;
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
t = 32*[1; 1];
t2s_transform = [R -t; 0 0 1];

% Warp forwards
[img_out] = func(img, [], t2s_transform, 4*size(img{1}));
subplot(2,2,3);
    image(img_out{1}); axis('image', 'ij', 'off');

% Warp back (should give something close to original input image).
[img_in] = func(img_out, [], inv(t2s_transform), 2*[64, 64]);
subplot(2,2,4);
    image(img_in{1}); axis('image', 'ij', 'off');
    
return
