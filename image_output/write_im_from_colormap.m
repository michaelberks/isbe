function image_out = write_im_from_colormap(image_in, file_name, cmap, clims, nan_color)
%make double
image_out = double(image_in);
nan_idx = isnan(image_out);

if nargin < 4 || isempty(clims);
    clims = [min(image_out(:)) max(image_out(:))];
end

if nargin < 5
    nan_color = [0 0 0];
end

%scale image from 0 to 1
if clims(2) > clims(1)
    image_out = (image_out - clims(1)) / (clims(2) - clims(1));
else
    warning(['Contrast limits [' num2str(clims) '] not valid, no contrast scaling applied.']);
end

%scale to 0-255 and round
num_colors = size(cmap,1);
image_out = uint8((num_colors-1)*image_out);

if nargin > 2
    %convert image to colormap
    
    image_out = ind2rgb(image_out, cmap);

    r = image_out(:,:,1);
    g = image_out(:,:,2);
    b = image_out(:,:,3);
    
    r(nan_idx) = nan_color(1);
    g(nan_idx) = nan_color(2);
    b(nan_idx) = nan_color(3);
    
    image_out = cat(3, r, g, b);
    
end

%save image
if ~isempty(file_name)
    imwrite(image_out, file_name);
end