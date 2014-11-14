function c_rgb = two_map_colour(map1, map2, lims1, lims2)
% Map a complex number to the rgb colorspace using a transform into the hsv
% colorspace, such that complex phase is represented as hue, saturation is
% 1, and magnitude is represented as brightness
if nargin < 4
    lims1 = [];
    lims2 = [];
end
if isempty(lims1)
    lims1 = [min(map1(:)) max(map1(:))];
end
if isempty(lims2)
    lims2 = [min(map2(:)) max(map2(:))];
end

%Get size of image;
[row col] = size(map1);

%pre-allocate c_rgb
c_rgb = zeros(row, col, 3);

%transform both maps to scale to [0 1]
map1 = (map1 - lims1(1)) / (lims1(2) - lims1(1));
map1(map1 < 0) = 0;
map1(map1 > 1) = 1;
map2 = (map2 - lims2(1)) / (lims2(2) - lims2(1));
map2(map2 < 0) = 0;
map2(map2 > 1) = 1;

%Scale map2 to lie between 0 and pi/2
map2 = map2*pi/2;

c_rgb(:,:,1) = map1 .* sin(map2);
c_rgb(:,:,2) = map1 .* cos(map2);
c_rgb(:,:,3) = 0;