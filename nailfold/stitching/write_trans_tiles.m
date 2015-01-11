function [diff_img seq_mask g_lims] = write_trans_tiles(tiles, compound_transforms, ...
                           trans_dir, trans_name, g_lims, mosaic, ...
                           varargin)
%WRITE_TRANS_TILES *Insert a one line summary here*
%   [nailfold_mosaic mosaic_weights offsets thetas] = write_trans_tiles(tiles, offsets, thetas, trans_dir, trans_name, g_lims)
%
% Inputs:
%      tiles - *Insert description of input variable here*
%
% Outputs:
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 23-Aug-2011
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'match_grey_levels', 1,...
    'diff_image', [],...
    'tile_masks', [],...
    'format', 'png',...
    'debug', 0);
clear varargin;

no_write = false;
if isnan(trans_dir)
    no_write = true;
end

if trans_dir(end) ~= filesep
    trans_dir = [trans_dir filesep];
end
if ~isdir(trans_dir)
    mkdir(trans_dir);
end

%Get dimensions of the input tiles
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

%--------------------------------------------------
%Loop each transform computing the size the mosaic needs to be
[mosaic_sz, compound_transforms] = ...
    mosaic_limits(tile_sz, compound_transforms);

if strcmp(get_username(), 'ptresadern')
    tb = timebar('title', 'Writing images', ...
                 'limit', num_tiles);
end

%Loop the the tiles again adding each one to mosaic and recording how much
%weight it adds at each pixel

diff_img = zeros(tile_sz);
for i_tile = 1:num_tiles
    
    display(['Correcting tile ' num2str(i_tile) ' of ' num2str(num_tiles)]);
    
    % Sample the mosaic into the camera coordinate frame (where the vessel
    % moves and dirt on the lens stays still).
    [mosaic_cam] = sample_tile_image({mosaic}, ones(mosaic_sz), ...
                                     compound_transforms(:,:,i_tile), ...
                                     tile_sz, []);
	mosaic_cam = full(mosaic_cam{1});
    
    %Get current tile
    if tiles_loaded
        tile_curr = double(tiles(:,:,i_tile));
    else
        tile_curr = double(imread(tiles{i_tile}));
    end
    
    if ~isempty(args.diff_image)
        tile_curr = tile_curr - args.diff_image;
    end
    
    diff_img_i = tile_curr - mosaic_cam;
    diff_img_i(isnan(diff_img_i)) = 0;
    diff_img = diff_img + diff_img_i;
    
    %Flag to say we don't want to write out any tiles - just skip to next
    %iteration and keep computing the diff image
    if no_write
        continue;
    end
        

    % Sample the tile into the vessel/mosaic coordinate frame (where the 
    % vessel is static).
    [trans_tile, mask] = sample_tile_image({tile_curr}, ones(tile_sz), ...
                                           inv(compound_transforms(:,:,i_tile)), ...
                                           mosaic_sz, []);
    mask = full(mask > 0.5);
    trans_tile = full(trans_tile{1});
    trans_tile(~mask) = NaN;
    
    % Pad image if odd rows/cols - movies don't like odd numbers
    if (rem(size(trans_tile,1), 2) == 1)
        trans_tile(end+1,:) = NaN;
        mask(end+1,:) = false;
    end
    if (rem(size(trans_tile,2), 2) == 1)
        trans_tile(:,end+1) = NaN;
        mask(:,end+1) = false;
    end

    if (i_tile == 1)
        % Store copy of first image and its mask
        img0 = trans_tile;
        mask0 = mask;
        
        g_min = min(trans_tile(mask));
        g_max = max(trans_tile(mask));
        
        seq_mask = mask;
    else
        % Find overlapping pixels between this image and the first - also
        % make sure we only use non NaN pixels
        pixels = mask0 & mask & ~isnan(img0) & ~isnan(trans_tile);
        
        if args.match_grey_levels
            % Match greyvalues by computing a 1D affine transformation.
            src_pixels = trans_tile(pixels);
            tgt_pixels = img0(pixels);
            gain_bias = [src_pixels ones(size(src_pixels))] \ tgt_pixels;
            trans_tile = (gain_bias(1) * trans_tile) + gain_bias(2);
        end
        
        min_i = min(trans_tile(pixels));
        max_i = max(trans_tile(pixels));
        
        g_min = min(min_i, g_min);
        g_max = max(max_i, g_max);

        seq_mask = seq_mask & mask;
    end
    
    % Write out the image
    filename = fullfile(trans_dir, ...
                        sprintf('%s%04d.%s',trans_name, i_tile, args.format));
                    
    if exist(filename,'file')
        % Overwriting images is slooooow.
        delete(filename);
    end
    
    trans_tile = uint8(255 * (trans_tile - g_lims(1)) / diff(g_lims));
    imwrite(trans_tile, filename);

    if exist('tb', 'var')
        timebar(tb, 'advance');
    end
end

g_lims = [g_min g_max];

% Get average difference between a frame and the mosaic (where camera-based
% artefacts have, with luck, been removed).
diff_img = diff_img / num_tiles;

if ~no_write
    filename = fullfile(trans_dir, 'seq_mask.png');
    imwrite(uint8(255*seq_mask), filename);
end

if exist('tb', 'var')
    timebar(tb, 'close');
end

% The sequence mask can later be applied to the outputted images using
% Imagemagick's <convert> or <mogrify> commands, e.g.
%   mogrify -clip-mask seq_mask.png -threshold 101% frame*.png

% Frames can be built into a video using ffmpeg.
% For a 1200kbs (lossy), 30 fps, mpeg movie (no sound):
%   ffmpeg -b 1200k -r 30 -i frame_%04d.png -an movie.mpg
% For a lossless, 30 fps, avi movie (no sound):
%   ffmpeg -i frame_%04d.png -vcodec huffyuv -r 30 -an movie.avi
