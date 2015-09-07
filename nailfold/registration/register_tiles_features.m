function [compound_transforms, match_counts] = register_tiles_features(tiles, varargin)
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
[compound_transforms, match_counts] = func(tiles, varargin{:});


%% The function
function [compound_transforms, match_counts] = func(tiles, varargin)
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'ref_type', 'previous',...
    'tile_masks', [],...
    'sigma', [2],...
    'feature_type', 'g1d',...
    'offset_lim', [],...
    'theta_range', [],...
    'max_pts', inf,...
    'weights', 'rect',...
    'debug', false, ...
    'compound_transforms', [], ...
    'offset_centres', [],...
    'mosaic', [],....
    'max_iterations', 1);

clear varargin;


%% Parse arguments
% Get dimensions of the input tiles

%Check if tiles are already in memory, in which case they should be an r x
%c x n_tiles array. Otherwise they should be cell of strings, giving the
%full filepath to each tile
tiles_loaded = isnumeric(tiles);

%Get tile size and number of tiles
if tiles_loaded
    [rows cols num_tiles] = size(tiles);
else
    num_tiles = length(tiles);
    im_header = imfinfo(tiles{1});
    rows = im_header.Height;
    cols = im_header.Width;
end

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

if (size(args.offset_centres,1)==num_tiles && size(args.offset_centres,2)==2)
    offset_centres = args.offset_centres;
else
    offset_centres = zeros(num_tiles,2);
end

if isempty(args.theta_range)
    theta_range = -10:10; %Use default 10 if range not user specified
else
    theta_range = args.theta_range;
end
theta_range = theta_range * pi/180;
theta_sz = length(theta_range);

match_counts = zeros(num_tiles, 1);

%% Pre-allocate space for compound transforms
if strcmpi(args.ref_type, 'previous')
    
    max_iterations = 1;
    
    compound_transforms = zeros(3, 3, num_tiles);
    for i = 1:num_tiles
        compound_transforms(:,:,i) = eye(3);
    end
    if tiles_loaded
        tile0 = double(tiles(:,:,1));
    else
        tile0 = double(imread(tiles{1}));
    end
    tile_range = 2:num_tiles;
else %args.ref_type = 'mosaic';
    
    max_iterations = args.max_iterations;
    compound_transforms = args.compound_transforms;
    
    if isempty(args.mosaic)      
        [tile0, ~, compound_transforms] = ...
            create_mosaic(tiles, compound_transforms, args.weights, args.tile_masks);
    else
        tile0 = args.mosaic;
        if isempty(compound_transforms)
            compound_transforms = zeros(3, 3, num_tiles);
            for i = 1:num_tiles
                compound_transforms(:,:,i) = eye(3);
            end
        end
    end

    tile_range = 1:num_tiles;
end

% Register images using feature point locations first
registration_method = 'features';

%% Loop through each of the remaining tiles, registering it to the previous
%  tile and saving the offset and theta
if strcmp(username(), 'ptresadern')
    ttl = ['Registering to ', args.ref_type, ' via features'];
    tb = timebar('title', ttl, ...
                 'limit', length(tile_range) * max_iterations);
end

it = 0;
while (it < max_iterations)
    if strcmp(registration_method, 'features')
        % Define the origin to be the centre of the reference tile.
        cx = (1 + size(tile0,2)) / 2;
        cy = (1 + size(tile0,1)) / 2;
       
        % Compute line points in reference tile as the initial target
        [x_tgt, y_tgt] = get_line_points(tile0, args.sigma, ...
                                         args.tile_masks(:,:,1), args.feature_type);
        x_tgt = x_tgt - cx - offset_centres(1,1);
        y_tgt = y_tgt - cy - offset_centres(1,2);
        
        n_updated = 0;
    end

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

    for i_tile = tile_range
        
        if tiles_loaded
            tile_curr = double(tiles(:,:,i_tile));
        else
            tile_curr = double(imread(tiles{i_tile}));
        end
        
        display(['Registering tile ' num2str(i_tile-1) ' of ' num2str(length(tile_range))]);
            
        if strcmp(registration_method, 'features')
            % Get feature-based estimate of motion between current frame and
            % reference.
            
            
            [x_src, y_src] = get_line_points(tile_curr, args.sigma, ...
                                             args.tile_masks(:,:,i_tile), args.feature_type);
            x_src = x_src - cx;
            y_src = y_src - cy;

            if strcmp(args.ref_type,'mosaic')
                % Transform to reference tile frame
                pts = compound_transforms(:,:,i_tile) * ...
                      [x_src(:)'; y_src(:)'; ones(1,numel(x_src))];
                x_src = pts(1,:)';
                y_src = pts(2,:)';
            end

            max_count = -inf;
            for th = 1:theta_sz
                theta = theta_range(th);
                [offset, count] = get_best_offset(theta, ...
                                                  [x_src y_src], [x_tgt y_tgt], ...
                                                  offset_lim, args.max_pts);

                % If this vote count is larger than the existing maximum record the
                % vote count and corresponding transformation
                if count > max_count
                    max_count = count;
                    
                    transform = [ cos(theta) -sin(theta) 0;
                                  sin(theta)  cos(theta) 0;
                                  0           0          1 ];
                              
                    offset_t = [cx cy] - (transform(1:2, 1:2)*[cx; cy])';
                          
                    if strcmp(args.ref_type,'mosaic')
                        transform(1,3) = offset(1) + offset_t(1);
                        transform(2,3) = offset(2) + offset_t(2);
                    else
                        transform(1,3) = offset(1) + offset_t(1) + offset_centres(i_tile-1,1);
                        transform(2,3) = offset(2) + offset_t(2) + offset_centres(i_tile-1,2);
                    end
                end
                match_counts(i_tile) = max_count;
            end

            transform_diff = transform - eye(3);
            if any(transform_diff(:) ~= 0)
                n_updated = n_updated + 1;
            end
        else
%             % Register using flow
%             transform = register_flow(double(tiles(:,:,i_tile+1)), ...
%                                       tile0, ...
%                                       compound_transforms(:,:,i_tile+1), 0);

            % Register using Levenberg-Marquadt
            ct0 = compound_transforms(:,:,i_tile);
            opts = optimset('lsqnonlin');
            opts = optimset(opts, ...
                            'display', 'iter', ...
                            'diffminchange', 1.0);
            dt = lsqnonlin(...
                    @(dt) lmerrfun(dt, tile_curr, tile0, ct0), ...
                    [0; 0], [], [], opts);
                
            transform = [1 0 dt(1); 0 1 dt(2); 0 0 1];
            disp(transform);
        end
        
        %-------------------------------------------------------------------
        % Transform the points for this tile to work out where they lie in the
        % full mosaic
       
        % Compound the offsets and rotations for all tiles up to this tile
        if strcmp(args.ref_type, 'previous')
            % Update compound transform by chaining existing sequence
            compound_transforms(:,:,i_tile) = transform * ...
                                          compound_transforms(:,:,i_tile-1);
        else
            % Update current tile->mosaic transformation
            compound_transforms(:,:,i_tile) = transform * ...
                                          compound_transforms(:,:,i_tile);
        end            

        if args.debug
            
            if strcmp(args.ref_type, 'previous') 
                if i_tile == 2
                    tile_prev = tile0;
                end
            
                % Display the transformed vessel points with respect to the first
                % frame
                xy_src_t = (compound_transforms(:,:,i_tile) * ...
                              [x_src y_src ones(size(x_src,1),1)]')';
                plot(all_vessels, xy_src_t(:,1), xy_src_t(:,2), '.', 'markersize', 4);

                % Display the two tile and the vessel points found on each
                figure(2); clf; colormap(gray(256));
                subplot(1,2,2);
                    imagesc(tile_curr); axis image; hold on;
                    plot(cx+x_src, cy+y_src, 'g.', 'markersize', 2);
                subplot(1,2,1);
                    imagesc(tile_prev); axis image; hold on;
                    plot(cx+x_tgt, cy+y_tgt, 'r.', 'markersize', 2);
                    % Display the vessel points from the src tile on the target tile
                    xy_src_t = (transform * [x_src y_src ones(size(x_src,1),1)]')';
                    plot(cx+xy_src_t(:,1), cy+xy_src_t(:,2), 'g.', 'markersize', 2);

                display(transform);
            
                %Save the previous tile to display next time...
                tile_prev = tile_curr; 
                
            else
                
                
                xy_src_t = (transform * [x_src y_src ones(size(x_src,1),1)]')';
                xy_src_o = (compound_transforms(:,:,i_tile) \ [xy_src_t(:,1:2) ones(size(x_src,1),1)]')';
                
                figure(2); clf;
                subplot(1,2,1);
                    imgray(tile_curr);
                    plot(cx+xy_src_o(:,1), cy+xy_src_o(:,2), 'g.', 'markersize', 2);
                    
                subplot(1,2,2); hold on;
                    plot(x_tgt, y_tgt, 'r.', 'markersize', 2);
                    plot(xy_src_t(:,1), xy_src_t(:,2), 'g.', 'markersize', 2); 
                    axis equal ij;
                    axis([min(xy_src_t(:,1)) max(xy_src_t(:,1)) min(xy_src_t(:,2)) max(xy_src_t(:,2))]);
                  
            end
                
        end

        % Update target points if we are to match consecutive frames
        if strcmp(args.ref_type, 'previous')
            x_tgt = x_src - offset_centres(i_tile,1);
            y_tgt = y_src - offset_centres(i_tile,2);
        end

        if exist('tb', 'var')
            timebar(tb, 'advance');
        end
    end

    if 0%strcmp(args.ref_type,'mosaic')
        % Redefine the reference tile
        [tile0, mosaic_weights, compound_transforms] = ...
            create_mosaic(tiles, compound_transforms, args.weights, args.tile_masks);
        
        if args.debug
            figure; imagesc(mosaic_weights); axis image;
            figure; imagesc(tile0); colormap(gray(256)); axis image;
        end
    end
    
    it = it + 1;
    
%     if (strcmp(args.ref_type, 'mosaic') && ...
%         strcmp(registration_method, 'features') && ...
%         ((it == max_iterations) || (n_updated == 0)))
%         % Switch to flow-based registration and restart iterations if no 
%         % further feature-based corrections were applied.
%         registration_method = 'flow';
%         
%         ttl = ['Registering to ', args.ref_type, ' via flow'];
%         timebar(tb, 'reset', 'title', ttl, ...
%                     'limit', length(tile_range) * max_iterations);
%         it = 0;
%         continue;
%     end
end

if exist('tb', 'var')
    timebar(tb, 'close');
end


%% Find line structures in the image
function [x, y] = ...
    get_line_points(tile, sigma, tile_mask, feature_type)

switch feature_type
    case {'g2d'},
        [line_strength, line_orientation] = gaussian_2nd_derivative_line(tile, sigma);
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


%% Determine best offset between sets of points for a given rotation
function [best_offset_theta count_theta] = ...
    get_best_offset(theta, xy_src, xy_tgt, offset_lim, max_pts)

%Generate rotation matrix for this theta
R = [cos(theta) sin(theta); 
     -sin(theta) cos(theta)];

% Rotate line points extracted from tile 1. It is assumed that the points
% have been normalized with respect to the centre of the image.
xy_src = xy_src * R;

%Set up hough matrix to count the votes for each offset
offset_sz = offset_lim(:,2) - offset_lim(:,1) + 1;
offset_counts = zeros(offset_sz(2), offset_sz(1));

%Loop through each line point from tile 1 and compute the offset to all
%points in tile 2
for ii = 1:min(max_pts, size(xy_src,1));
    % Compute offset and decentre
    x_offset = round(xy_tgt(:,1) - xy_src(ii,1));
    y_offset = round(xy_tgt(:,2) - xy_src(ii,2));

    % Count only those votes less than the offset limit
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

offset_counts = imfilter(offset_counts, fspecial('gaussian', 7, 2), 'replicate');

%Workout the offset with maximum vote for this theta
[count_theta max_idx_theta] = max(offset_counts(:));

%Convert the maximum offset ID
[max_row max_col] = ind2sub([offset_sz(2) offset_sz(1)], max_idx_theta);
best_offset_theta = [max_col max_row] + offset_lim(:,1)' - 1;


%% Use optical flow to register two images
function [transform_out] = register_flow(tile_src, tile_tgt, ...
                                         transform_in, ...
                                         max_iterations, tolerance)

if ~exist('transform_in', 'var'), transform_in = eye(3); end
if ~exist('max_iterations', 'var'), max_iterations = 1; end
if ~exist('tolerance', 'var'), tolerance = 1e-2; end

% Returns the identity transform by default, indicating no relative
% movement between the frames.
transform_out = eye(3);

n_iterations = 0;
dp_prev = [0, 0];
delta_rel = tolerance + 1.0;
Ix = []; Iy = [];
while (n_iterations < max_iterations) && ...
      (delta_rel > tolerance)
      
    % Sample pixels from the reference (tgt) tile that roughly correspond to
    % those in the source (src) tile.
    t_weights = ones(size(tile_tgt));
    tile_tgt2 = sample_tile_image({tile_tgt}, t_weights, ...
                                 inv(transform_out * transform_in), ...
                                 size(tile_tgt), size(tile_src));
    tile_tgt2 = tile_tgt2{1};

    output = pt_flow_derivs(tile_src, tile_tgt2, Ix, Iy);
    Ix = output.Ix;
    Iy = output.Iy;

    dx = median(output.dx(~isnan(output.dx)));
    dy = median(output.dy(~isnan(output.dy)));
    transform_out = [ 1 0 dx;
                      0 1 dy;
                      0 0 1 ] * transform_out;

    % Compute relative change in dp
    dp = [dx, dy];
    delta = (dp - dp_prev);
    delta_rel = norm(delta) / norm(dp_prev);
    dp_prev = dp;

    n_iterations = n_iterations + 1;
end
disp(n_iterations);


%% Error function to be used in LM optimization
function err = lmerrfun(dt, tile_src, tile_tgt, transform0)
transform = [1 0 dt(1); 0 1 dt(2); 0 0 1] * transform0;
tile_tgt2 = sample_tile_image({tile_tgt}, [], inv(transform), size(tile_src));
tile_tgt2 = tile_tgt2{1};

err = tile_tgt2(:) - tile_src(:);

% Replace NaNs with values that maintain the Mean Squared Error (but which 
% will affect the Sum of Squared Error).
err2 = err(~isnan(err));
mse = err2'*err2 / length(err2);
err(isnan(err)) = sqrt(mse);

sse = err'*err;
fprintf('%g %g %f\n', dt(1), dt(2), sse);

ignore = 0;


%% Evaluate whether displacements are a local minimum to the nearest pixel
function test_disps(tile_tgt, tiles, transforms)

hw = 1;
for i_tile = 1:size(tiles,3)
    tile_src = double(tiles(:,:,i_tile));

    emat = nan(3,3);
    for dx = -hw:hw
        for dy = -hw:hw
            dt = transforms(1:2, 3, i_tile);
            xrng = (1:size(tile_src,2)) + dt(1) + dx;
            yrng = (1:size(tile_src,1)) + dt(2) + dy;
            
            validx = (1 <= xrng) & (xrng <= size(tile_tgt,2));
            validy = (1 <= yrng) & (yrng <= size(tile_tgt,1));
            
            samples_tgt = tile_tgt(yrng(validy),xrng(validx));
            samples_src = tile_src(validy,validx);
            
            valid = ~isnan(samples_src) & ~isnan(samples_tgt);
            s = samples_src(valid);
            t = samples_tgt(valid);
            
            if 1
                % Mean squared error
                errvec = s-t;
                err = (errvec'*errvec) / length(errvec);
            else
                % Normalized correlation
                s = s - mean(s);
                t = t - mean(t);
                err = (s'*t) / (sqrt(s'*s)*sqrt(t'*t));
            end
                
            emat(dy+hw+1, dx+hw+1) = err;
        end
    end
    
    disp(emat);
end


%% Test script
function test_script()
clc; close all;
timebar('closeall');

imgroot = 'U:\projects\nailfold\capture\';
imgpath = fullfile(imgroot,'2012_10_22\Left\Digit4\x300\seq2\corrected');
d = dir(fullfile(imgpath,'frame_*.png'));
d = d(1:18:end);
% d = d([1,end]);
n_frames = length(d);

outroot = 'U:\tmp\ncm\registration';
if ~exist(outroot, 'dir')
    mkdir(outroot);
end

frame1 = imread(fullfile(imgpath,d(1).name));
frame1 = mean(frame1,3);

frames = zeros([size(frame1), n_frames], 'uint8');
frames(:,:,1) = uint8(frame1);
for i = 2:n_frames
    frame = imread(fullfile(imgpath,d(i).name));
    frames(:,:,i) = uint8(mean(frame,3));
end

% Register every frame to its predecessor
[t2m_transforms0] = func(frames, 'theta_range', [0], ...
                                 'offset_lim', 40);

% xt = compound_transforms0(1,3,:);
% yt = compound_transforms0(2,3,:);
% figure(); clf;
%     plot(xt(:), yt(:), 'b-');
                                  
[mosaic, m_weights, t2m_transforms1] = ...
    create_mosaic(frames, t2m_transforms0);
figure(); clf; colormap(gray(256));
    imagesc(mosaic); axis('image');

outpath = fullfile(outroot, 'round1');
if ~exist(outpath, 'dir')
    mkdir(outpath);
else
    delete(fullfile(outpath,'*.png'));
end
for i = 1:size(frames,3)
    sampled = sample_tile_image({double(frames(:,:,i))}, [], ...
                                inv(t2m_transforms1(:,:,i)), ...
                                size(mosaic));
    filename = sprintf('frame_%03d.png', i);
    imwrite(uint8(full(sampled{1})), fullfile(outpath, filename));
end

% Now register all frames to the mosaic until convergence
% Note that the offset limit can be reduced now because they have already been
% approximately registered.
[t2m_transforms2] = func(frames, 'theta_range', [0], ...
                                 'offset_lim', 2, ...
                                 'compound_transforms', t2m_transforms1, ...
                                 'max_iterations', 1);

[mosaic, m_weights, t2m_transforms3] = ...
    create_mosaic(frames, t2m_transforms2);
figure(); clf; colormap(gray(256));
    imagesc(mosaic); axis('image');

outpath = fullfile(outroot, 'round3');
if ~exist(outpath, 'dir')
    mkdir(outpath);
else
    delete(fullfile(outpath,'*.png'));
end
for i = 1:size(frames,3)
    sampled = sample_tile_image({double(frames(:,:,i))}, [], ...
                                inv(t2m_transforms3(:,:,i)), ...
                                size(mosaic));
    filename = sprintf('frame_%03d.png', i);
    imwrite(uint8(full(sampled{1})), fullfile(outpath, filename));
end
