function write_im_from_colormap(image_in, file_name, cmap, clims, nan_color)
%make double
image_in = double(image_in);
nan_idx = isnan(image_in);

if nargin < 4 || isempty(clims);
    clims = [min(image_in(:)) max(image_in(:))];
end

if nargin < 5
    nan_color = [0 0 0];
end

%scale image from 0 to 1
image_in = (image_in - clims(1)) / (clims(2) - clims(1));

%scale to 0-255 and round
num_colors = size(cmap,1);
image_in = uint8((num_colors-1)*image_in);

if nargin > 2
    %convert image to colormap
    
    image_in = ind2rgb(image_in, cmap);

    r = image_in(:,:,1);
    g = image_in(:,:,2);
    b = image_in(:,:,3);
    
    r(nan_idx) = nan_color(1);
    g(nan_idx) = nan_color(2);
    b(nan_idx) = nan_color(3);
    
    image_in = cat(3, r, g, b);
    
end

%save image
imwrite(image_in, file_name);