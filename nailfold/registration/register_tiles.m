function [nailfold_mosaic mosaic_weights compound_transforms] = ...
    register_tiles(tiles, varargin)
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

% offset_lim = [ ox_min ox_max 
%                oy_min oy_max ]
if isempty(args.offset_lim)
    %Use default +/-20 in x and y if range not user specified
    offset_lim = [0 cols; -ceil(rows/2) ceil(rows/2)];
else
    offset_lim = args.offset_lim;
    if size(offset_lim, 1) == 1
        % i.e. offset_lim = [o_min o_max]
        offset_lim = [offset_lim; offset_lim];
    end
    if size(offset_lim, 2) == 1
        % i.e. offset_lim = [ox; oy]
        offset_lim = [-offset_lim offset_lim];
    end
end
offset_sz = offset_lim(:,2) - offset_lim(:,1) + 1;

if isempty(args.theta_range)
    theta_range = -10:10; %Use default 10 if range not user specified
else
    theta_range = args.theta_range;
end
theta_range = theta_range * pi/180;
theta_sz = length(theta_range);

%% Pre-allocate space for output offsets and thetas
offsets = zeros(num_tiles - 1, 2);
thetas = zeros(num_tiles - 1, 1);
transforms = zeros(3, 3, num_tiles);
compound_transforms = zeros(3, 3, num_tiles);

transforms(:,:,1) = eye(3);
compound_transforms(:,:,1) = transforms(:,:,1);

%% Compute line points in first tile as the initial target
[x_tgt, y_tgt] = process_tile(tiles(:,:,1), args.sigma, ...
                              args.tile_masks(:,:,1));
                          
if args.debug
    figure(1); clf;
    all_vessels = axes; 
    set(all_vessels, ...
        'Ydir', 'reverse',...
        'DataAspectRatio', [1 1 1],...
        'NextPlot', 'add');
    setappdata(all_vessels,'PlotHoldStyle',true);
    plot(all_vessels, x_tgt, y_tgt, '.', 'markersize', 4);
end

%% Loop through each of the remaining tiles, registering it to the previous
%  tile and saving the offset and theta
if strcmp(username(), 'ptresadern')
    tb = timebar('title', 'Registering frames', ...
                 'limit', num_tiles);
    timebar(tb, 'advance');
end

for tt = 1:num_tiles-1
    % Get feature points from the image
    [x_src, y_src] = process_tile(double(tiles(:,:,tt+1)), args.sigma, ...
                                  args.tile_masks(:,:,tt+1));
    
    %--------------------------------------------------
    %Loop through thetas finding rotation that highest offset max

    max_count = 0;
    for th = 1:theta_sz
        theta = theta_range(th);

        [best_offset_theta count_theta] = ...
            evaluate_theta(theta, [x_src y_src], [x_tgt y_tgt], ...
                           offset_lim, args.max_pts);
                       
        % If this vote count is larger than the existing maximum record the
        % offset index and theta
        if count_theta > max_count
            max_count = count_theta;
            
            offsets(tt,:) = best_offset_theta;
            thetas(tt) = theta;
            
            transforms(:,:,tt+1) = ...
                [ cos(theta) -sin(theta) best_offset_theta(1);
                  sin(theta)  cos(theta) best_offset_theta(2);
                  0           0           1 ];
        end
    end
    
    %-------------------------------------------------------------------
    %Transform the points for this tile to workout where they lie in the
    %full mosaic
    
    %Compound the offsets and rotations for all tiles up to this tile
% 	compound_transforms(:,:,tt+1) = ...
%         transforms(:,:,tt+1) * compound_transforms(:,:,tt);
	compound_transforms(:,:,tt+1) = transforms(:,:,tt+1);
    
    if args.debug
        % Display the transformed vessel points with respect to the first
        % frame
        xy_src_t = (compound_transforms(:,:,tt+1) * ...
                      [x_src y_src ones(size(x_src,1),1)]')';
        plot(all_vessels, xy_src_t(:,1), xy_src_t(:,2), '.', 'markersize', 4);

        % Display the two tile and the vessel points found on each
        figure(2); clf; colormap(gray(256));
        subplot(1,2,2);
            imagesc(tiles(:,:,tt+1)); axis image; hold on;
            plot(cx+x_src, cy+y_src, 'g.', 'markersize', 2);
        subplot(1,2,1);
            imagesc(tiles(:,:,tt)); axis image; hold on;
            xy_src_t = (transforms(:,:,tt+1) * ...
                          [x_src y_src ones(size(x_src,1),1)]')';
            % Display the vessel points from the src tile on the target tile
            plot(cx+x_tgt, cy+y_tgt, 'r.', 'markersize', 2);
            plot(cx+xy_src_t(:,1), cy+xy_src_t(:,2), 'g.', 'markersize', 2);
    end
    
    %Update the target (x,y) coordinates for the next iteration
%     x_tgt = x_src;
%     y_tgt = y_src;
    
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
                                  mosaic_sz, ...tile_sz, ...
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


function [x, y] = ...
    process_tile(tile, sigma, tile_mask)

feature_type = 'g1d';

switch feature_type
    case {'g2d'},
        [line_strength, line_orientation] = gaussian_clover_line(tile, sigma);
    case {'g1d'},
        [line_strength, line_orientation] = ...
            gaussian_1st_derivative_gradient(tile, sigma);
    otherwise
        error(['Feature type ',feature_type,' not recognized']);
end

%Apply non-maximal suppression to skeletonise the line strength map
line_nms = mb_non_maximal_supp(line_strength, line_orientation);

%Discard points from edges of map
if ~isempty(tile_mask)
    line_nms(~tile_mask) = 0;
end

%Apply hysterisis to select suitable lines from the NMS map
[line_mask] = hysterisis(line_nms);

%Extract (x,y) coordinates of the remaining lines
[y x] = find(line_mask);

%Make points relative to centre
[rows cols] = size(tile);
cx = (cols + 1) / 2;
cy = (rows + 1) / 2;
y = y - cy;
x = x - cx;


function [best_offset_theta count_theta] = ...
    evaluate_theta(theta, xy_src, xy_tgt, offset_lim, max_pts)

%Generate rotation matrix for this theta
R = [cos(theta) -sin(theta); 
     sin(theta) cos(theta)];

% Rotate line points extracted from tile 1. It is assumed that the points
% have been normalized with respect to the centre of the image.
xy_src = xy_src * R;

%Set up hough matrix to count the votes for each offset
offset_sz = offset_lim(:,2) - offset_lim(:,1) + 1;
offset_counts = zeros(offset_sz(2), offset_sz(1));

%Loop through each line point from tile 1 and compute the offset to all
%points in tile 2
for ii = 1:min(max_pts, size(xy_src,1));
    %Computes offset and decentre
    x_offset = round(xy_tgt(:,1) - xy_src(ii,1));
    y_offset = round(xy_tgt(:,2) - xy_src(ii,2));

    %Only count votes less than the offset limit
    keep = x_offset >= offset_lim(1,1) & ...
           x_offset <= offset_lim(1,2) & ...
           y_offset >= offset_lim(2,1) & ...
           y_offset <= offset_lim(2,2);

    %Use sparse matrix trick to counts votes
    if any(keep)
        offset_counts = offset_counts + ...
            sparse(y_offset(keep)-offset_lim(2,1)+1, ...
                   x_offset(keep)-offset_lim(1,1)+1, ...
                   1, offset_sz(2), offset_sz(1));
    end
end

%Workout the offset with maximum vote for this theta
[count_theta max_idx_theta] = max(offset_counts(:));

%Convert the maximum offset ID
[max_row max_col] = ind2sub([offset_sz(2) offset_sz(1)], max_idx_theta);
best_offset_theta = [max_col max_row] + offset_lim(:,1)' - 1;


