function [grad_map edge_map] = segment_breast_snake(image_in, mlo)
% segment_breast
% created: 09/05/2006 16:34
% author: Michael Berks
% function: segment the breast border using Ferrari et al (2004)
%           method

%Get size of mammogram
[rows, cols] = size(image_in); %#ok

%Compute histogram of the mammogram
[counts, x_range] = hist(image_in(:), 128);

%Find the index of the maximal bin in the first half of the intensity range
[dummy max_idx] = max(counts(1:64));

%As the intensity increases, the bin counts decrease - find the point at
%which they reach a minimum before increasing again, and use this as an
%upper threshold
counts(1:max_idx-1) = [];
x_range(1:max_idx-1) = [];
min_idx = find(diff(counts) >= 0, 1);
upper_threshold = ceil(x_range(min_idx));

%%Also for later, find out a threshold to throw-out the really high stuff
%counts = counts(end:-1:1);
%x_range = x_range(end:-1:1);
%min_idx = find(diff(counts) >= 0, 1);

%Get binary masks based on the threshold
upper_mask = image_in > upper_threshold;

%Throw away all but the main connected region in the mask, and fill the
%region in
[ll, no_objects] = bwlabel(upper_mask, 4);
[nn] = hist(ll(:), no_objects + 1);
[mm ind] = max(nn(2:end));
upper_mask = upper_mask == ind;
upper_mask = bwmorph(upper_mask, 'fill');

%Now in the using the upper_mask find the two points (1 in the top half,
% 1 in the bottom) of minimum width
widths = inf(rows, 1);
for ii = 1:rows
    w = find(~upper_mask(ii,:), 1);
    if ~isempty(w)
        widths(ii) = w;
    end
end

[x_top, y_top] = min(widths(1:floor(rows/2)));
[x_bot, y_bot] = min(widths(ceil(rows/2):end));
y_bot = y_bot + ceil(rows/2) - 1;

%Use these points as cutoffs for the main breast region
if mlo
    %For MLO go vertically upwards
    upper_mask(1:y_top, x_top+1) = 0;
    upper_mask(y_bot, :) = 0;
    
    y_start = 1;
else
    %For CC go horizontally inwards
    upper_mask([y_top y_bot], :) = 0;
    y_start = y_top+1;
end

%label the new mask and throw away outside the cuttoffs
ll = bwlabel(upper_mask, 4);
upper_mask = ll == ll(y_start, 1);

%Now trace the boundary of the main breast
x_start = find(upper_mask(y_start,:), 1, 'last');

y_end = y_bot-1;
x_end = find(upper_mask(y_end,:), 1, 'last');

inner_border = bwtraceboundary(upper_mask, [y_start, x_start], 'E');

%Throw away the boundary between the start and end points we defined
end_idx = find(inner_border(:,1) == y_end & inner_border(:,2) == x_end, 1);
inner_border(end_idx+1:end,:) = [];

%figure;
%imagesc(upper_mask); axis image; colormap gray; hold on;
%plot([x_start x_end], [y_start y_end], 'go');
    
%Throw away points that appear to lie on a vertical edge at start and end
%of border (not at start for MLO)
go_on = true;
ii = 1;
while go_on
    if std(inner_border(ii:ii+9,2)) > 2
        go_on = false;
    else
        if mlo
            ii = ii+1;
        else
            inner_border(ii,:) = [];
        end
        %plot(inner_border(ii,2),inner_border(ii,1), 'rx');        
    end
end

go_on = true;
while go_on
    if std(inner_border(end-9:end,2)) > 2
        go_on = false;
    else
        %plot(inner_border(end,2),inner_border(end,1), 'rx');
        inner_border(end,:) = [];
    end
end

%Add points to the start and end so border joins up with edge of image
inner_border = [inner_border(1,1) 1; inner_border; inner_border(end,1) 1];

%Now take the convex hull of this border
if mlo
    discard_idx1 = setdiff(1:ii, convhull(inner_border(1:ii,2), inner_border(1:ii,1)));
    discard_idx2 = setdiff(ii+1:size(inner_border,1), convhull(inner_border(ii+1:end,2), inner_border(ii+1:end,1))+ii);
    discard_idx = [discard_idx1 discard_idx2];
    %plot(inner_border(discard_idx1,2), inner_border(discard_idx1,1), 'bx');
    %plot(inner_border(discard_idx2,2), inner_border(discard_idx2,1), 'rx');
else
    discard_idx = setdiff(1:size(inner_border,1), convhull(inner_border(:,2), inner_border(:,1)));
end

inner_border(discard_idx,:) = [];

%Discard the first point for MLOs
if mlo
    inner_border(1,:) = [];
end

%Linearly interpolate the points to generate a border with evenly spaced
%points
len = 200;
dists = cumsum([0; sqrt(sum(diff(inner_border).^2, 2))]);
smooth_border = interp1(dists, inner_border, linspace(0, dists(end), len), 'linear');

%Make border x-y not [r c]
smooth_border = fliplr(smooth_border);

%Compute normal profiles of the image at every point
[dummy, fxy] = gradient(smooth_border);

%normalise fxy
fxy = fxy ./ [sqrt(sum(fxy.^2, 2)), sqrt(sum(fxy.^2, 2))];
fxy = fliplr(fxy);
fxy(:,2) = -fxy(:,2);
start_pts = smooth_border + 20*fxy;

%Now create a gradient map so that each pixel is assigned the normal
%gradient of the border pixel it is nearest too
xx = repmat(1:cols, rows, 1);
yy = repmat((1:rows)', 1, cols);

angles = atan(-fxy(:,2) ./ fxy(:,1));
grad_map = griddata(smooth_border(:,1), smooth_border(:,2), angles, xx, yy, 'nearest');

%Find edges using gaussian deriative filters
sigma = 1;

% Magic numbers
GaussianDieOff = .0001;

% Design the filters - a gaussian and its derivative

pw = 1:30; % possible widths
ssq = sigma^2;
width = find(exp(-(pw.*pw)/(2*ssq))>GaussianDieOff,1,'last');
if isempty(width)
    width = 1;  % the user entered a really small sigma
end

t = (-width:width);
gau = exp(-(t.*t)/(2*ssq))/(2*pi*ssq);     % the gaussian 1D filter

% Find the directional derivative of 2D Gaussian (along X-axis)
% Since the result is symmetric along X, we can get the derivative along
% Y-axis simply by transposing the result for X direction.
[x,y]=meshgrid(-width:width,-width:width);
dgau2D=-x.*exp(-(x.*x+y.*y)/(2*ssq))/(pi*ssq);

% Convolve the filters with the image in each direction

%smooth the image out
aSmooth=imfilter(image_in,gau,'conv','replicate');   % run the filter accross rows
aSmooth=imfilter(aSmooth,gau','conv','replicate'); % and then accross columns

%apply directional derivatives
ax = imfilter(aSmooth, dgau2D, 'conv','replicate');
ay = imfilter(aSmooth, dgau2D', 'conv','replicate');

%gradient and normal alignment
alignment_map = angle(complex(cos(grad_map), -sin(grad_map)).*conj(complex(ax,ay)));

edge_map = sqrt((ax.*ax) + (ay.*ay));
edge_map = edge_map - 0.05*image_in;

edge_map(abs(alignment_map) < 2.5) = min(edge_map(:));

nms_map = mb_non_maximal_supp(edge_map, grad_map, 0, 1, min(edge_map(:)));
go_on = true;
ii = 2;
[snake_pts,e] = mb_snake_normal(start_pts, 0, .1, 20, 1, fxy(:,1), fxy(:,2), nms_map);
while go_on
    [snake_pts,e(ii)] = mb_snake_normal(snake_pts, 0, .1, 20, 1, fxy(:,1), fxy(:,2), nms_map);
    go_on = e(ii) < e(ii-1);
    display(['Iteration: ', num2str(ii)]);
    ii = ii+1;
end


breast_border = snake_pts;


%Finally smooth the breast border
short_len = 50;
dists = cumsum([0; sqrt(sum(diff(breast_border).^2, 2))]);
breast_border = interp1(dists, breast_border, linspace(0, dists(end), short_len), 'linear');
dists = cumsum([0; sqrt(sum(diff(breast_border).^2, 2))]);
breast_border = interp1(dists, breast_border, linspace(0, dists(end), len), 'spline');

figure;
subplot(1,4,1:2); imagesc(image_in); axis image; colormap(jet(256));
hold on;
plot(breast_border(:,1), breast_border(:,2), 'y:');

% subplot(1,4,1:2); imagesc(alignment_map); axis image; colormap(jet(256)); 
% hold on;
% plot(snake_pts(:,1), snake_pts(:,2), 'rx');

subplot(1,4,3:4); imagesc(nms_map); axis image; colormap(jet(256)); caxis([-1 1]);
hold on;
plot(snake_pts(:,1), snake_pts(:,2), 'rx');



