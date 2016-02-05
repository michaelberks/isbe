function [mosaic, mosaic_weights, t2m_transforms, cleaned_tiles, edge_mask, g_lims] = ...
    create_mosaic(tiles, t2m_transforms, tile_weight_type, tile_masks, mosaic, mosaic_weights)

if ~exist('tile_weight_type','var') || isempty(tile_weight_type), tile_weight_type = 'rect'; end
if ~exist('tile_masks','var'), tile_masks = []; end
if ~exist('mosaic','var'), mosaic = []; end
if ~exist('mosaic_weights','var'), mosaic_weights = []; end

tiles_loaded = isnumeric(tiles);


%Get tile size and number of tiles
if tiles_loaded
    tile_sz = size(tiles(:,:,1));
    num_tiles = size(tiles,3);
else
    num_tiles = length(tiles);
    im_header = imfinfo(tiles{1});
    tile_sz(1) = im_header.Height;
    tile_sz(2) = im_header.Width;
end

if isempty(tile_masks)
    tile_masks = zeros(0, 0, num_tiles);
end

clean_tiles = false;
if nargout >= 4
    clean_tiles = true;
end

if isempty(mosaic)
    % Get size of mosaic
    [mosaic_sz, t2m_transforms] = ...
        mosaic_limits(tile_sz, t2m_transforms);
    
    %Always create even dimensions for mosaics
    if rem(mosaic_sz(1), 2)
        mosaic_sz(1) = mosaic_sz(1) + 1;
    end
    if rem(mosaic_sz(2), 2)
        mosaic_sz(2) = mosaic_sz(2) + 1;
    end
    
    % Pre-allocate space for the nailfold mosaic and the sum of tile weights
    mosaic = zeros(mosaic_sz);
    mosaic_weights = zeros(mosaic_sz);
else
    mosaic_sz = size(mosaic);
    if rem(mosaic_sz(1), 2)
        mosaic_sz(1) = mosaic_sz(1) + 1;
        mosaic(end+1,:) = 0;
        mosaic_weights(end+1,:) = 0;
    end
    if rem(mosaic_sz(2), 2)
        mosaic_sz(2) = mosaic_sz(2) + 1;
        mosaic(:,end+1) = 0;
        mosaic_weights(:,end+1) = 0;
    end
    mosaic = mosaic.*mosaic_weights;
    mosaic(isnan(mosaic))=0;
end

if clean_tiles
    cleaned_tiles = zeros(mosaic_sz(1), mosaic_sz(2), num_tiles);
    edge_mask = false(mosaic_sz);
end

t_weights = tile_weights(tile_sz, tile_weight_type);

%Loop the the tiles again adding each one to mosaic and recording how much
%weight it adds at each pixel
if exist('tb', 'var')
    timebar(tb, 'title', 'Creating mosaic');
    timebar(tb, 'reset');
end

for i_tile = 1:num_tiles
    % Transform the points using the compounded transforms

    % (For i_tile==1, the transform is the identity so this could be more efficient,
    % but it's easier to read this way.)

    %Get current tile
    if tiles_loaded
        tile_curr = double(tiles(:,:,i_tile));
    else
        tile_curr = double(imread(tiles{i_tile}));
    end
    
    % The tiles to sample i.e. the image and the weight map.       
    tiles_in = {tile_curr, ones(tile_sz)};

    tiles_out = sample_tile_image(tiles_in, t_weights, ...
                                  inv(t2m_transforms(:,:,i_tile)), ...
                                  mosaic_sz, tile_masks(:,:,i_tile));
                              
    if clean_tiles
        clean_tile_i = tiles_out{1} ./ tiles_out{2};
        edge_mask = edge_mask | isnan(clean_tile_i);
        %clean_tile_i(edge_mask) = 0;
        cleaned_tiles(:,:,i_tile) = clean_tile_i;
    end
        
    valid = ~isnan(tiles_out{1});
    mosaic(valid) = mosaic(valid) + tiles_out{1}(valid);
    mosaic_weights(valid) = mosaic_weights(valid) + tiles_out{2}(valid);

    if exist('tb', 'var')
        timebar(tb, 'advance');
    end
end

if exist('tb', 'var')
    timebar(tb, 'close');
end

% Point divide the nailfold intensities by the weights to produce the final
% nailfold mosaic
mosaic = mosaic ./ mosaic_weights;
mosaic(mosaic_weights == 0) = NaN;

if clean_tiles
    g_lims = [inf, -inf];
    if all(edge_mask(:))
        cleaned_tiles(:) = 0;
        return;
    end
    for i_tile = 1:num_tiles
        clean_tile_i = cleaned_tiles(:,:,i_tile);
        g_lims(1) = min(g_lims(1), min(clean_tile_i(~edge_mask)));
        g_lims(2) = max(g_lims(2), max(clean_tile_i(~edge_mask)));
        cleaned_tiles(:,:,i_tile) = clean_tile_i;
    end
    for i_tile = 1:num_tiles
        clean_tile_i = cleaned_tiles(:,:,i_tile);
        clean_tile_i = (clean_tile_i - g_lims(1)) / (g_lims(2)-g_lims(1));
        clean_tile_i(edge_mask) = 0;
        cleaned_tiles(:,:,i_tile) = clean_tile_i;
    end
end
