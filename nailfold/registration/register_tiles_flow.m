function [nailfold_mosaic mosaic_weights compound_transforms] = ...
    register_tiles_flow(tiles, varargin)
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

if (nargin==0 && nargout==0), test_script(); return; end
[nailfold_mosaic mosaic_weights compound_transforms] = func(tiles, varargin{:});

function [nailfold_mosaic mosaic_weights compound_transforms] = ...
    func(tiles, varargin)

args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'tile_masks', [],...
    'sigma', [2],...
    'offset_lim', [],...
    'theta_range', [],...
    'max_pts', inf,...
    'weights', 'rect',...
    'debug', 0);
clear varargin;

%% Parse arguments
%Get dimensions of the input tiles
[rows cols num_tiles] = size(tiles);
tile_sz = [rows, cols];
cx = (1+cols)/2; 
cy = (1+rows)/2;

%--------------------------------------------------
%Generate offset and theta ranges to be used when registering each tile
if isempty(args.tile_masks)
    args.tile_masks = zeros(0, 0, num_tiles);
end

%% Pre-allocate space for output offsets and thetas
offsets = zeros(num_tiles - 1, 2);
thetas = zeros(num_tiles - 1, 1);
transforms = zeros(3, 3, num_tiles);
compound_transforms = zeros(3, 3, num_tiles);

transforms(:,:,1) = eye(3);
compound_transforms(:,:,1) = transforms(:,:,1);

%% Loop through each of the remaining tiles, registering it to the previous 
%  tile and saving the offset and theta
if strcmp(username(), 'ptresadern')
    tb = timebar('title', 'Registering frames', ...
                 'limit', num_tiles-1);
end

for tt = 1:num_tiles-1
    % Get offset between consecutive tiles
    transform = eye(3);
    tile0 = double(tiles(:,:,tt));
    dp = [inf, inf];
    n_iterations = 0;
    while (norm(dp) > 1e-3) && ...
          (n_iterations < 10)
        tile1 = double(tiles(:,:,tt+1));
        t_weights = ones(size(tile1));
        
        tile1 = sample_tile_image({tile1}, t_weights, transform, ...
                                   tile_sz, tile_sz);
        tile1 = tile1{1};
                              
        output = pt_flow_derivs(tile0, tile1);
        dx = median(output.dx(~isnan(output.dx)));
        dy = median(output.dy(~isnan(output.dy)));
        dp = [dx, dy];

        transform = [ 1 0 -dx;
                      0 1 -dy;
                      0 0 1 ] * transform;

        n_iterations = n_iterations + 1;
    end
    
    transforms(:,:,tt+1) = transform;
    
    %-------------------------------------------------------------------
    %Transform the points for this tile to workout where they lie in the
    %full mosaic
    
    %Compound the offsets and rotations for all tiles up to this tile
	compound_transforms(:,:,tt+1) = ...
        transforms(:,:,tt+1) * compound_transforms(:,:,tt);
    
    if exist('tb', 'var')
        timebar(tb, 'advance');
    end
end

% Get size of mosaic
[mosaic_sz, compound_transforms] = ...
    mosaic_limits(tile_sz, compound_transforms);

% Pre-allocate space for the nailfold mosaic and the sum of tile weights
nailfold_mosaic = zeros(mosaic_sz);
mosaic_weights = zeros(mosaic_sz);

t_weights = tile_weights(tile_sz, args.weights);

%Loop the the tiles again adding each one to mosaic and recording how much
%weight it adds at each pixel
if exist('tb', 'var')
    timebar(tb, 'title', 'Creating mosaic');
    timebar(tb, 'reset');
end
for tt = 1:num_tiles
    % Transform the points using the compounded transforms
    
    % (For tt==1, the transform is the identity so this could be more efficient,
    % but it's easier to read this way.)
    
    % The tiles to sample i.e. the image and the weight map.
    tiles_in = {double(tiles(:,:,tt)), ...
                ones(size(tiles(:,:,tt)))};
    
    tiles_out = sample_tile_image(tiles_in, t_weights, ...
                                  compound_transforms(:,:,tt), ...
                                  tile_sz, mosaic_sz, ...
                                  args.tile_masks(:,:,tt));

    nailfold_mosaic = nailfold_mosaic + tiles_out{1};
    mosaic_weights = mosaic_weights + tiles_out{2};
    
    if exist('tb', 'var')
        timebar(tb, 'advance');
    end
end
if exist('tb', 'var')
    timebar(tb, 'close');
end

%Point divide the nailfold intenities by the weights to produce the final
%nailfold mosaic
nailfold_mosaic = nailfold_mosaic ./ mosaic_weights;
nailfold_mosaic(mosaic_weights == 0) = NaN;

if args.debug
    figure; imagesc(mosaic_weights); axis image;
    figure; imagesc(nailfold_mosaic); colormap(gray(256)); axis image;
end


%% Test script
function test_script()
clc; clear; close all;
timebar('closeall');

imgroot = 'U:\projects\nailfold\capture\2012_10_22\Left\Digit4\x300\';
imgpath = fullfile(imgroot, 'seq2\preprocessed\registered_g1d\masked\tmp');

d = dir(fullfile(imgpath,'frame_*.png'));

for i = 2:length(d)
    img1 = mean(imread(fullfile(imgpath, d(i-1).name)), 3);
    dx = 2;
    img2 = img1(:, 1:(end-dx));
    img1 = img1(:, (dx+1):end);
end

tiles = cat(3, img1, img2);

[mosaic mosaic_weights compound_transforms] = func(tiles);



