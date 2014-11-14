function [diff_img seq_mask] = dual_camera_register(frame1, frame2, compound_transforms, ...
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
    'tile_masks', [],...
    'format', 'png',...
    'debug', 0);
clear varargin;

if trans_dir(end) ~= filesep
    trans_dir = [trans_dir filesep];
end
if ~isdir(trans_dir)
    mkdir(trans_dir);
end

%Get dimensions of the input tiles
[rows cols num_tiles] = size(tiles);
tile_sz = [rows, cols];

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
for tt = 1:num_tiles
    % Sample the mosaic into the camera coordinate frame (where the vessel
    % moves and dirt on the lens stays still).
    [mosaic_cam] = sample_tile_image({mosaic}, ones(mosaic_sz), ...
                                     compound_transforms(:,:,tt), ...
                                     tile_sz, []);
	mosaic_cam = full(mosaic_cam{1});
    
    img = double(tiles(:,:,tt));
    diff_img = diff_img + (img - mosaic_cam);

    % Sample the tile into the vessel/mosaic coordinate frame (where the 
    % vessel is static).
    [trans_tile, mask] = sample_tile_image({img}, ones(tile_sz), ...
                                           inv(compound_transforms(:,:,tt)), ...
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

    if (tt == 1)
        % Store copy of first image and its mask
        img0 = trans_tile;
        mask0 = mask;
        
        seq_mask = mask;
    else
        % Find overlapping pixels between this image and the first
        pixels = (mask0 & mask);
        
        % Match greyvalues by computing a 1D affine transformation.
        src_pixels = trans_tile(pixels);
        tgt_pixels = img0(pixels);
        gain_bias = [src_pixels ones(size(src_pixels))] \ tgt_pixels;
        trans_tile = (gain_bias(1) * trans_tile) + gain_bias(2);

        seq_mask = seq_mask & mask;
    end
    
    % Write out the image
    filename = fullfile(trans_dir, ...
                        sprintf('%s%04d.%s',trans_name, tt, args.format));
                    
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

% Get average difference between a frame and the mosaic (where camera-based
% artefacts have, with luck, been removed).
diff_img = diff_img / num_tiles;

filename = fullfile(trans_dir, 'seq_mask.png');
imwrite(uint8(255*seq_mask), filename);

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
