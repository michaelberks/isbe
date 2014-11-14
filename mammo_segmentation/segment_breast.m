function [segmentation normal] = segment_breast(varargin)
%SEGMENT_BREAST Applies breast segementation to a mammogram
%   [failed_list] = segment_breast(varargin)
%
% SEGMENT_BREAST uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Required Arguments:
%
% 'mammo'
%   - mammogram to segment
%
% 'right'
%   - {0,1} indicates if breast is right (1) or left (0)
%
% 'mlo'
%   - {0,1} indicates if breast is ML/MLO (1) or CC (0)
%
% Optional Arguments:
%
% 'offset [50 150]'
%   - Pixel offset to sample normal profiles inside and outside the
%   approximate breast edge
%
% 'alpha(0.5),beta(0.5)'
%   - Weighting arguments for snake algorithm (MB_SNAKE_NORMAL), should work for 0-255
%   (8-bit) grayscale images with height 1024 pixels but may need to
%   experiment for a given dataset (see Mike Berks' thesis for more
%   details)
%
% 'gamma (3)'
%   - Weighting function to find breast edge - may need varying depending
%   for different datasets (see Mike Berks' thesis for more
%   details)
%
% 'resolution (1)' 
%   - pixel resolution of search region in MB_SNAKE_NORMAL
%
% 'search_width_fraction (0.1)' 
%   - Set watch percentage of the normal
%  profile we'll search in 1 iteration (e.g. setting this to 1 will search
%  the whole profile in a single iteration and return the global optimum but
%  will prob very slow)
%
%  'iterate (1)' 
%   - if true, run iteration of snake fitting
%
%  'nipple (0)'
%   - if true, try and find the nipple (this doesn't really work - it was never properly tested or
%   used, but may be experimented with further)
%
% Outputs:
%  segmentation - structure containing fields:
%       - breast_border, [Nx2] array of x,y coordinates giving the final
%       border of the whole breast region, completely encloses the breast
%
%      - breast_air, [nx1] array of indices that specify the points of
%      breast_border that lie on the skin/air boundary
%
%      - pectoral_edge (MLO only) [2x2] giving x,y points on top and side of image
%        defining straight line fitted to pectoral muscle
%
%      - nipple, xy coords of nipple, or empty if not found (the default)
%
%     - size of image, stored so segmentation can resized for any change in
%     image size (see SEGMENT_BREAST_RESIZE)
%
%  normal - structure containing normal profiles, useful for further debugging 
%
%
% Example: [segmentation] = segment_breast('mammo', mam, 'mlo', mlo, 'right', right, 'plot', 0);
%
% Notes: Assumes mammograms are named containing either their viewcode
% either: LML, LCC, RML, RCC. Will attempt to correctly orient the
% mammogram - right breasts will then be flipped in SEGMENT_BREAST
%
% See also: SEGMENT_BREAST_BATCH, SEGMENT_BREAST_RESIZE, MB_SNAKE_NORMAL
%
% Created: 09-May-2006
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 5125 
% Copyright: (C) University of Manchester

args = u_packargs(varargin,... % the user's input
			 '0', ... % strict mode
			 {'mammo',...
             'mlo',...
             'right'}, ... % the mandatory arguments
			 'alpha', 0.5,...
             'beta', 0.5,...
             'gamma', 3,...
             'border_erosion_size', 75,...
             'offset', [50 50], ...
             'resolution', 1,...
             'search_width_fraction', 0.1,...
             'iterate', true,...
             'nipple', false,...
             'plot', 1,...
             'quiet', 0);

%Prepare the mammogram, it should have a height of 1024 pixels, an
%intensity range of 0-256, and have the chestwall on the left edge of the
%image for left breasts and right edge for right breasts. (i.e. the breast
%should be upright). The anonymised label should be white (or at least
%intensity)
if args.right
    
    %anonymised label in top left corner
    anon_r = find(args.mammo(:,1), 1) - 1;
    anon_c = find(args.mammo(1,:), 1) - 1;
    
    args.mammo(1:anon_r, 1:anon_c) = 255;
    
    %Now flip the image
    args.mammo = fliplr(args.mammo);
    
    
else
    %anonymised label in bottom right corner
    anon_r = find(args.mammo(end:-1:1,end), 1) - 2;
    anon_c = find(args.mammo(end,end:-1:1), 1) - 2;
    
    args.mammo(end-anon_r:end, end-anon_c:end) = 255;
end

%Resize of mammogram to 1024 rows if not already
if size(args.mammo,1) ~= 1024;
    args.mammo = imresize(args.mammo, [1024 NaN], 'bilinear');
end
args.mammo = double(args.mammo);

%Get size of resized mammogram
[rows, cols] = size(args.mammo);

%**************************************************************************
% The first stage of the segmentation algorithm is to find a reasonable
% inner border for the breast. We select a suitable threshold to compute a
% binary mask, and then perform a series of morphological operations on
% this mask to removes artefacts not belonging to the breast region
%**************************************************************************

%Compute histogram of the mammogram
[counts, x_range] = hist(args.mammo(:), 128);

%Throw away the counts of pixels with zero intensity from the first bin is
%these are due to the anonymisation
counts(1) = counts(1) - sum(~args.mammo(:));

%Find the index of the maximal bin in the first half of the intensity range
[dummy max_idx] = max(counts(1:64));

%As the intensity increases, the bin counts decrease - find the point at
%which they reach a minimum before increasing again, and use this as an
%upper threshold
counts(1:max_idx-1) = [];
x_range(1:max_idx-1) = [];
min_idx = find(diff(counts) >= 0, 1);
upper_threshold = ceil(x_range(min_idx));

%Also for later, find out a threshold to throw-out the really high stuff
% In fact, fixing one at 240 seems fine, although we could try something
% adaptive in the future
high_threshold = 240;

%Get binary masks based on the threshold
upper_mask = args.mammo > upper_threshold;

%Throw away all but the main connected region in the mask, and fill the
%region in
[ll, no_objects] = bwlabel(upper_mask, 4); %#ok
[nn] = hist(ll(:), no_objects + 1);
[mm main_label] = max(nn(2:end));
upper_mask = ll == main_label;

%**************************debugging**********************************
% figure; 
% subplot(1,2,1); imagesc(upper_mask); axis image; hold on;
%*************************************************************************

%************************************************
%Close the mask (to fill in any holes inside the main region)
upper_mask = imdilate(upper_mask, strel('disk', 11));
upper_mask = imerode(upper_mask, strel('disk', 11));

%Open the mask to remove the straight line border artifacts
upper_mask = imerode(upper_mask, strel('disk', args.border_erosion_size));
upper_mask = imdilate(upper_mask, strel('disk', args.border_erosion_size));

%We may still have artefacts across the top and bottom - find which rows
%are all 1's
full_rows = all(upper_mask, 2);

%Start at the top and find the first non full row and then the first 0
%pixel, set everything above and to the right of this point to zero
rr = find(~full_rows, 1, 'first');
cc = find(~upper_mask(rr,:), 1, 'first');
upper_mask(1:rr, cc:cols) = 0;

%Do the same at the bottom
rr = find(~full_rows, 1, 'last');
cc = find(~upper_mask(rr,:), 1, 'first');
upper_mask(rr:rows, cc:cols) = 0;

%************************************************

%Again, throw away any regions we have now disconnected - note the breast
%region may no longer be the largest region, so we assume the point
%half-way up the chestwall will always belong to the main breast region
[ll] = bwlabel(upper_mask, 4);
main_label = ll(round(rows/2), 1);
upper_mask = ll == main_label;

%Finally get  a mask of all pixels above the threshold, not in the main
%region and connected to one of the four corners - set the pixels to the
%mean of all the back ground not in these regions
if ~args.mlo %we don't need to do this for MLO
    discard_mask = bwselect((args.mammo > upper_threshold) & ~upper_mask,...
        [1 1 cols cols], [1 rows rows 1]);
    %args.mammo(discard_mask) = mean(args.mammo(~upper_mask & ~discard_mask));
    fill_value = mean(args.mammo(~upper_mask & ~discard_mask));
    args.mammo(discard_mask) = fill_value; %nan;
end

%**************************debugging**********************************
% subplot(1,2,2); imagesc(upper_mask); axis image; hold on;
%*************************************************************************

%Now find the two points (1 in the top half,
% 1 in the bottom) of minimum width in the mask
widths = inf(rows, 1);
for ii = 1:rows
    w = find(~upper_mask(ii,:), 1);
    if ~isempty(w)
        widths(ii) = w;
    end
end

min_w = min(widths(ceil(rows/2):end)) + 1;
y_end = find(widths(ceil(rows/2):end) <= min_w, 1, 'first') + ceil(rows/2) - 2;
x_end = widths(y_end+1);

%Use these points as cutoffs for the main breast region
if args.mlo
    %For args.mlo go vertically upwards at top and horizontally outwards    
    y_start = 1;
    x_start = find(upper_mask(1,:), 1, 'last');
    
else    
    min_w = min(widths(1:floor(rows/2))) + 1;
    y_start = find(widths(1:floor(rows/2)) <= min_w, 1, 'last') + 1; 
    x_start = widths(y_start-1);
end

%************************** debugging**********************************
% plot([x_start x_end], [y_start y_end], 'gx', 'MarkerSize', 10);
%*************************************************************************

%Trace the border
inner_border = bwtraceboundary(upper_mask, [y_start, x_start], 'E');

%Make this border x-y, not rc
inner_border = fliplr(inner_border);

%************************** debugging**********************************
% plot(inner_border(:,1), inner_border(:,2), 'b.');
%*************************************************************************

%Throw away the boundary between the start and end points we defined
end_idx = find(inner_border(:,1) == x_end & inner_border(:,2) == y_end, 1);
inner_border(end_idx+1:end,:) = [];

%************************** debugging**********************************
% plot(inner_border(:,1), inner_border(:,2), 'g.');
%*************************************************************************

if args.mlo
    %Throw away points that appear to lie on a vertical edge at start of the 
    % border for MLOs
    go_on = true;
    while go_on
        if inner_border(1,1) ~= inner_border(2,1)
            go_on = false;
        else
            %plot(inner_border(1,1),inner_border(1,2), 'ro');
            inner_border(1,:) = [];                 
        end
    end
    go_on = true;
    while go_on
        if any(inner_border(1,:) ~= (inner_border(2,:) + [1 -1]))
            go_on = false;
        else
            %plot(inner_border(1,1),inner_border(1,2), 'ro');
            inner_border(1,:) = [];        
        end
    end
    
    %Throw away points that appear to lie on a horizontal edge at end of the 
    % border for MLOs
    go_on = true;
    while go_on
        if inner_border(end,2) ~= inner_border(end-1,2)
            go_on = false;
        else
            %plot(inner_border(end,1),inner_border(end,2), 'ro');
            inner_border(end,:) = [];                 
        end
    end
    go_on = true;
    while go_on
        if any(inner_border(end,:) ~= (inner_border(end-1,:) + [-1 1]))
            go_on = false;
        else
            %plot(inner_border(end,1),inner_border(end,2), 'ro');
            inner_border(end,:) = [];       
        end
    end
end

%************************** debugging**********************************
% plot(inner_border(:,1), inner_border(:,2), 'y.');
%*************************************************************************

%Now take a smoother average of this border
if args.mlo
    % For MLO, take a moving average of the border
    inner_border = imfilter(inner_border(1:5:end,:), ones(5,1)/5, 'replicate', 'same');
    
    %Compute the normal vectors at each point on the inner border
    [fx, fy] = gradient(inner_border); clear fx;
    border_angle = atan2(fy(:,1),fy(:,2));

    %Discard points at the start on an overhanging horizontal edge
    start_pt = find(border_angle > -pi/4, 1, 'first');
    inner_border = inner_border(start_pt:end,:);
else
    %For CC, taking the convex hull of the intial border works well
    %Add points to the start and end so border joins up with edge of image
    inner_border = [1 inner_border(1,2); inner_border; 1 inner_border(end,2)];
    
    %take the convex hull of this border
    discard_idx = setdiff(1:size(inner_border,1), convhull(inner_border(:,1), inner_border(:,2)));
    inner_border(discard_idx,:) = [];

    %Discard the first point and last points
    inner_border([1 end],:) = [];
end

%************************** debugging**********************************
% plot(inner_border(:,1), inner_border(:,2), 'w');
%*************************************************************************

%************************** debugging**********************************
% plot(inner_border(:,1), inner_border(:,2), 'g');
%*************************************************************************

if 1
%Linearly interpolate the points to generate a border with evenly spaced
%points
len = 50; %why 50?
dists = cumsum([0; sqrt(sum(diff(inner_border).^2, 2))]);
smooth_border = interp1(dists, inner_border, linspace(0, dists(end), len), 'linear');

%**************************************************************************
% Having obtained a reasonable inner border for the breast, we now attempt
% to find the actual skin-air boundary. To achieve this we take normal
% profiles of the mammogram at a set points sampled on the inner border to
% form an unwrapped image of the border region. The skin-air border is then
% an approximately vertical edge within this image, with a sample point on
% each row. It remains to find the appropriate x co-ordinate along each
% row, which we do by fitting a snake
%**************************************************************************

%median filter the mammogram - why, justify this?
filtered_image = medfilt2(args.mammo, [11 11]);

%The width of the normal profiles is determined  the
%start and end pts relative to the inner border in which to search for the
%skin-air border, as defined in the optional argument offsets
norm_width = sum(args.offset);

% Pre-allocate containers for the normal profiles and the associated
% sampling points
normal_p = zeros(len, norm_width);
normal_x = zeros(len, norm_width);
normal_y = zeros(len, norm_width);

%Compute the normal vectors at each point on the inner border
[fx, fy] = gradient(smooth_border); clear fx;

%normalise fy
fy = fy ./ [sqrt(sum(fy.^2, 2)), sqrt(sum(fy.^2, 2))];

%Compute normal profiles of the image at every point
for ii = 1:length(fy(:,1)) %= number of rows in skin_air
    
    n1_x = smooth_border(ii, 1) - args.offset(1)*fy(ii, 2);
    n1_y = smooth_border(ii, 2) + args.offset(1)*fy(ii, 1);
    n2_x = smooth_border(ii, 1) + args.offset(2)*fy(ii, 2);
    n2_y = smooth_border(ii, 2) - args.offset(2)*fy(ii, 1);

    [cx, cy, cp] = improfile(filtered_image, [n1_x, n2_x], [n1_y, n2_y], norm_width);
    normal_p(ii, :) = cp';
    normal_x(ii, :) = cx';
    normal_y(ii, :) = cy';

end

%Having constructed the unwrapped border image comprised of normal
%profiles, compute the edge map. We are looking for an edge right to left,
%such that the intensity on the right is minimal. The optional argument
%gamma compromises between edge strength, and the minimality of intensity
%to the right of the edge
norm_edge = imfilter(normal_p, [1 2 4 -args.gamma -4 -2 -1], 'replicate');

%Define what percentage of the normal profile we will search in each
%iteration of the snake
search_width = round(args.search_width_fraction*norm_width);

%Compute the maximal edge strength on each row and use as the intial points
%for the snake
[dummy max_idx] = max(norm_edge,[],2);
initial_line = [max_idx (1:len)'];

%Run first iteration of the snake
[snake_pnts,e] = mb_snake_normal(initial_line, args.alpha, args.beta, search_width, args.resolution, norm_edge, normal_x, normal_y);

%Repeat iterations of the snake until the overall energy increases again
%(note this may find a local minima)
if args.iterate
    go_on = true;
    ii = 2;
    while go_on
        [snake_pnts,e(ii)] = mb_snake_normal(snake_pnts, args.alpha, args.beta, search_width, args.resolution, norm_edge, normal_x, normal_y);
        go_on = e(ii) < e(ii-1);
        if args.quiet
            display(['Iteration: ', num2str(ii)]);
        end
        ii = ii+1;
    end
end

%Having found the skin-air boundary co-ordinates in the unwrapped border
%image, convert these to co-ordinates in the original mammogram
breast_border = zeros(len,2);
for ii = 1:len
    breast_border(ii,1) = normal_x(ii, snake_pnts(ii,1));
    breast_border(ii,2) = normal_y(ii, snake_pnts(ii,1));
end

%Choose where we extrapolate the first and last points depending on image
%type
if args.mlo
    p_start = polyfit(breast_border(1:5,2),breast_border(1:5,1),1);
    p_end = polyfit(breast_border(end-4:end,1),breast_border(end-4:end,2),1);
    breast_border = [1 1; sum(p_start) 1; breast_border; 1 sum(p_end)];
    
else
    p_start = polyfit(breast_border(1:5,1),breast_border(1:5,2),1);
    p_end = polyfit(breast_border(end-4:end,1),breast_border(end-4:end,2),1);
    breast_border = [1 sum(p_start); breast_border; 1 sum(p_end)];
end

%Make a mask of the image using the breast border
breast_mask = roipoly(args.mammo, breast_border(:,1), breast_border(:,2));

%Now throw out the bright white edges near the chest-wall edges
%filtered_image(discard_mask) = fill_value;
breast_mask = breast_mask & (filtered_image < high_threshold);

%Throw away all but the main connected region in the mask, and fill the
%region in
[ll, no_objects] = bwlabel(breast_mask, 4); %#ok
[nn] = hist(ll(:), no_objects + 1);
[mm ind] = max(nn(2:end));
breast_mask = ll == ind;

%Take the perimeter of the breast mask
[y_start, x_start] = find(breast_mask, 1);
breast_perim = fliplr(bwtraceboundary(breast_mask, [y_start, x_start], 'E'));
    
%Finally take every 10th point along the border - why 10 should this be an
%argument?
final_border = breast_perim(1:10:end,:);

%Use the previous breast_border as the breast/air border, but project each
%point onto the nearest point in the final border, then take unique rows
if args.mlo
    breast_border(1,:) = [];
end

dists = (final_border(:,1)-breast_border(1,1)).^2+ ...
        (final_border(:,2)-breast_border(1,2)).^2;
[dummy start_idx] = min(dists);

dists = (final_border(:,1)-breast_border(end,1)).^2+ ...
        (final_border(:,2)-breast_border(end,2)).^2;
[dummy end_idx] = min(dists);

if start_idx > end_idx
    end_idx = length(final_border(:,1));
end
breast_air = (start_idx:end_idx)';

%--------------------------------------------------------------------------
%Find the nipple (optional - doesn't always work very well)
%--------------------------------------------------------------------------
if args.nipple
    nba = length(breast_air);
    si = round(nba/4);
    ei = 3*si;
    border_sums = zeros(nba,1);
    
    [xx yy] = meshgrid(-20:20, -20:20);
    [circ_idx] = find(xx.^2 + yy.^2 <= 20^2);
    for ii = si:ei
        cx = round(final_border(breast_air(ii),1));
        cy = round(final_border(breast_air(ii),2));
        idx = sub2ind(size(args.mammo), cy+yy(circ_idx), cx+xx(circ_idx));
        border_sums(ii) = mean(args.mammo(idx));
    end
    [dummy nipple_idx] = max(border_sums);
    nipple = final_border(breast_air(nipple_idx),:);
else
    nipple = [];
end

%--------------------------------------------------------------------------
% If mlo do pectoral segmentation
%--------------------------------------------------------------------------
if args.mlo
    [yy xx] = find(breast_mask);
    breast_centre = [mean(xx) mean(yy)];
    
    %Extract pectoral region as the box bounded by the start of the breast
    %air border and the centorid of the breast
    pectoral = double(args.mammo(1:round(breast_centre(2)), 1:round(final_border(breast_air(1),1))));
    
    %Filter the pectoral region and discard high/low pixels
    pectoral_med = medfilt2(pectoral, [9 9]);
    pectoral_mask = pectoral_med > 240 | pectoral_med < 2;
    pectoral_med(pectoral_mask) = 0;
    
    %Perform canny edge detection on pectoral region
    [pectoral_canny ori] = canny_edge(pectoral_med, [], 1, 0.98, .6);
    
    %Discard edges that aren't approximatelt aligned on an upperleft
    %diagonal
    ignore_ori = ori > pi/4 | ori < pi/24;
    pectoral_ori = pectoral_canny;
    pectoral_ori(ignore_ori) = 0;
    
    %Take hough transformation of edges to locate most likely staright line
    [line_scores rhos thetas] = hough_line(pectoral_ori);
    [max_rhos rho_idx] = max(line_scores);
    [dummy theta_idx] = max(max_rhos);
    
    %Convert from rho,theta parametisation to y = mx + c;
    theta = thetas(theta_idx);
    rho = rhos(rho_idx(theta_idx));
    m = tan(theta); 
    c = rho / cos(theta);
    
    %Work out points on pectoral edge at top and bottom of image
    y1 = 1; x1 = (y1 - c)/m;
    y2 = rows; x2 = (y2 - c)/m;
    pectoral_edge = [x1 y1; x2 y2];
end

%--------------------------------------------------------------------------
% Display visual results of the algorithm - show the breast and the initial
% and final segmented borders, also show the profile images used to find
% the skin air border
%--------------------------------------------------------------------------
if args.plot
    figure;
    subplot(1,6,1:2); imagesc(filtered_image); axis image; colormap(jet(256));
    hold on;
    plot(inner_border(:,1), inner_border(:,2), 'y:');
    plot(final_border(:,1), final_border(:,2), 'y', 'linewidth', 2);
    plot(final_border(breast_air,1), final_border(breast_air,2), 'g:', 'linewidth', 2);
    if args.mlo
        plot([1 breast_centre(1)], [breast_centre(2) breast_centre(2)], 'c');
        plot([breast_centre(1) breast_centre(1)], [1 breast_centre(2)], 'c');
        plot(pectoral_edge(:,1), pectoral_edge(:,2), 'm');
    end
    
    if args.nipple
        plot(nipple(1), nipple(2), 'r*', 'MarkerSize', 10);
    end
    
    subplot(1,6,3:4); imagesc(upper_mask); axis image; colormap(gray(256));
    hold on;
    plot(inner_border(:,1), inner_border(:,2), 'y:');
    plot(final_border(:,1), final_border(:,2), 'y', 'linewidth', 2);
    plot(final_border(breast_air,1), final_border(breast_air,2), 'g:', 'linewidth', 2);
    if args.mlo
        plot([1 breast_centre(1)], [breast_centre(2) breast_centre(2)], 'c');
        plot([breast_centre(1) breast_centre(1)], [1 breast_centre(2)], 'c');
        plot(pectoral_edge(:,1), pectoral_edge(:,2), 'm');
    end
    
    if args.nipple
        plot(nipple(1), nipple(2), 'r*', 'MarkerSize', 10);
    end

    subplot(1,6,5); imagesc(normal_p); axis ij; colormap(jet(256));
    hold on;
    plot(snake_pnts(:,1), snake_pnts(:,2), 'y');
    subplot(1,6,6); imagesc(norm_edge); axis ij; colormap(jet(256));
    hold on;
    plot(snake_pnts(:,1), snake_pnts(:,2), 'y');
    plot(initial_line(:,1), initial_line(:,2), 'w');
    plot(initial_line(:,1), initial_line(:,2), 'wx');
end

%Finally, if right breast, reflect x-coords of final border
if args.right
    final_border(:,1) = cols - final_border(:,1) + 1;
    if args.mlo
        pectoral_edge(:,1) = cols - pectoral_edge(:,1) + 1;
    end
    if args.nipple
        nipple(1) = cols - nipple(1) + 1;
    end
end

segmentation.breast_border = final_border;
segmentation.breast_air = breast_air;
segmentation.size = [rows, cols];
segmentation.nipple = nipple;
if args.mlo
    segmentation.pectoral_edge = pectoral_edge;
end
if nargout > 1
    normal.profiles = normal_p;
    normal.x = normal_x;
    normal.y = normal_y;
    normal.edge = norm_edge;
    normal.gamma = args.gamma;
end

%*****************End of function******************************************
end



