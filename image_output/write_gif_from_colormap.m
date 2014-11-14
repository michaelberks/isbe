function write_gif_from_colormap(image_in, file_name, cmap, append, clims, delay)
%make double
image_in = double(image_in);

if nargin < 4
    append = 'append';
end
if nargin < 5
    clims = [min(image_in(:)) max(image_in(:))];
end
if nargin < 6
    delay = 0.5;
end

%scale image from 0 to 1
image_in = image_in - clims(1);
if (clims(2) ~= clims(1))
    image_in =  image_in / (clims(2) - clims(1));
end

%scale to 0-255 and round (assume cmap is 256 element)
image_in = uint8(255*image_in);

%save image
if strcmpi(append, 'overwrite')
    imwrite(image_in, cmap, file_name, 'gif', 'WriteMode', append, 'LoopCount', Inf, 'DelayTime', delay);
else
    imwrite(image_in, cmap, file_name, 'gif', 'WriteMode', append, 'DelayTime', delay);
end