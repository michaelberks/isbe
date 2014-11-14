function [diff_imgs] = dual_camera_translate(frames, compound_transforms, mosaic)
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

%Get dimensions of the input tiles
[rows cols num_tiles] = size(frames);
tile_sz = [rows, cols];
diff_imgs = zeros(rows, cols, num_tiles);

%--------------------------------------------------
%Loop each transform computing the size the mosaic needs to be
[mosaic_sz, compound_transforms] = ...
    mosaic_limits(tile_sz, compound_transforms);

for i_tile = 1:2

    % Sample the mosaic into the camera coordinate frame (where the vessel
    % moves and dirt on the lens stays still).
    [mosaic_cam] = sample_tile_image({mosaic}, ones(mosaic_sz), ...
                                     compound_transforms(:,:,i_tile), ...
                                     tile_sz, []);
    mosaic_cam = full(mosaic_cam{1});

    frame = double(frames(:,:,i_tile));
    diff_imgs(:,:,i_tile) = frame - mosaic_cam;
end

