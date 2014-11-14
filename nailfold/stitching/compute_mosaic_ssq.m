function [nailfold_ssq ssq_count] = compute_mosaic_ssq(tiles, offsets, thetas, varargin)
%REGISTER_FRAME *Insert a one line summary here*
%   [nailfold_mosaic mosaic_weights offsets thetas] = register_tile(tiles, tile_masks, offset_lim, theta_lim, debug)
%
% Inputs:
%      tiles - *Insert description of input variable here*
%
%      offset_lim - *Insert description of input variable here*
%
%      theta_lim - *Insert description of input variable here*
%
%
% Outputs:
%      nailfold_mosaic - *Insert description of input variable here*
%
%      offsets - *Insert description of input variable here*
%
%      thetas - *Insert description of input variable here*
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
    'offset_lim', [],...
    'theta_lim', [],...
    'weights', 'rect',...
    'debug', 0);
clear varargin;

%Get dimensions of the input tiles
[rows cols num_tiles] = size(tiles);

%Initialise the pts that keep track of the corners of transformed points
x_min = 0;
x_max = cols;
y_min = 0;
y_max = rows;

%--------------------------------------------------
%Loop each transform computing the size the mosaic needs to be
for tt = 1:num_tiles-1
    
    %-------------------------------------------------------------------
    %Transform the points for this tile to workout where they lie in the
    %full mosaic
    
    %Compund the offsets and rotations for all tiles up to this tile
    compound_offset = sum(offsets(1:tt,:), 1);
    compund_theta = sum(pi*thetas(1:tt)/180);
    R = [cos(compund_theta) -sin(compund_theta); sin(compund_theta) cos(compund_theta)];
    
    %Transform the corner points using the compunded transforms to update
    %the max size of the final mosaic
    xyt = [-cols/2 -rows/2; cols/2 rows/2]*R;
    x_corner = xyt(:,1) + compound_offset(1) + cols/2;
    y_corner = xyt(:,2) + compound_offset(2) + rows/2;
    x_min = floor(min(x_min, min(x_corner)));
    x_max = ceil(max(x_max, max(x_corner)));
    y_min = floor(min(y_min, min(y_corner)));
    y_max = ceil(max(y_max, max(y_corner)));
  
end

%Shift max x,y by min values so indexing starts ats 1, and add 1 as a
%buffer needed in the interpolation steps
x_max = x_max - x_min + 2;
y_max = y_max - y_min + 2;

%pre-allocate space for the nailfold mosaic and the sum of tile weights
nailfold_ssq = zeros(y_max, x_max);
ssq_count = zeros(y_max, x_max);

%Make some containers of x,y coordinates the same size as each tile
x_tile = repmat(1:cols, rows, 1);
y_tile = repmat((1:rows)', 1, cols);

%Loop the the tiles again adding each one to mosaic and recording how much
%weight it adds at each pixel
for tt = 1:num_tiles
     
    %Get transformed x,y points for this tile
    if tt == 1
        %For the first tile this is just x_tile,y_tile shifted by the
        %minimum x,y vals
        if isempty(args.tile_masks)
            x_pts = x_tile(:);
            y_pts = y_tile(:);
        else
            x_pts = x_tile(args.tile_masks(:,:,tt));
            y_pts = y_tile(args.tile_masks(:,:,tt));
        end
    else
        %For frames 2+, workout the transform required
        %Compund the offsets and rotations for all tiles up to this tile
        compound_offset = sum(offsets(1:tt-1,:),1);
        compund_theta = sum(pi*thetas(1:tt-1)/180);
        R = [cos(compund_theta) -sin(compund_theta); sin(compund_theta) cos(compund_theta)];

        %Transform the points using the compunded transforms
        if isempty(args.tile_masks)
            xyt = [x_tile(:)-cols/2 y_tile(:)-rows/2]*R;
        else
            xyt = [...
                x_tile(args.tile_masks(:,:,tt))-cols/2 ...
                y_tile(args.tile_masks(:,:,tt))-rows/2]*R;
        end
        x_pts = xyt(:,1) + compound_offset(1) + cols/2;
        y_pts = xyt(:,2) + compound_offset(2) + rows/2;
    end
    
    %Shift the x,y pts by the minimum so indexing starts at 1
    x_pts = x_pts  - x_min + 1;
    y_pts = y_pts  - y_min + 1;
    
    %Workout which pixels in the mosaic the transformed points contribute
    x_pix = floor(x_pts);
    y_pix = floor(y_pts);
    
    %Compute the interpolation weights for each point
    x_weight = x_pts - x_pix;
    y_weight = y_pts - y_pix;
    
    %Compute the the interpolation weights
    weights_00 = (1-x_weight) .* (1-y_weight);
    weights_01 = x_weight .* (1-y_weight);
    weights_10 = (1-x_weight) .* y_weight;
    weights_11 = x_weight .* y_weight;
    
    %Add the tile intensities * interpolation weights to the mosiac
    curr_tile = double(tiles(:,:,tt));
    if isempty(args.tile_masks)
        curr_tile = curr_tile(:);
    else
        curr_tile = curr_tile(args.tile_masks(:,:,tt));
    end
    
    tile_ssq = zeros(y_max, x_max);
    ssq_weights = zeros(y_max, x_max);
    
    tile_ssq = tile_ssq +...
        sparse(y_pix, x_pix, curr_tile .* weights_00, y_max, x_max);
        
    tile_ssq = tile_ssq +...
        sparse(y_pix, x_pix + 1, curr_tile .* weights_01, y_max, x_max);
        
    tile_ssq = tile_ssq +...
        sparse(y_pix + 1, x_pix, curr_tile .* weights_10, y_max, x_max);
        
    tile_ssq = tile_ssq +...
        sparse(y_pix + 1, x_pix + 1, curr_tile .* weights_11, y_max, x_max);
    
    %Add the tile weights * interpolation weights to the mosaic weights
    ssq_weights = ssq_weights +...
        sparse(y_pix, x_pix, weights_00, y_max, x_max);
        
    ssq_weights = ssq_weights +...
        sparse(y_pix, x_pix + 1, weights_01, y_max, x_max);
        
    ssq_weights = ssq_weights +...
        sparse(y_pix + 1, x_pix, weights_10, y_max, x_max);
        
    ssq_weights = ssq_weights +...
        sparse(y_pix + 1, x_pix + 1, weights_11, y_max, x_max);
    
    %Get a mask of vaild pixels
    ssq_mask = ssq_weights > (1 - 1e-6);
    
    %Add the square of tile to main mosiac sum of squares
    nailfold_ssq(ssq_mask) = nailfold_ssq(ssq_mask) + tile_ssq(ssq_mask).^2;
    ssq_count(ssq_mask) = ssq_count(ssq_mask) + 1;

end

