function img_out = imscale(img_in)
% Return a copy of img_in, scaled to range [0..255]

img_out = img_in - min(img_in(:));
img_out = img_out / (max(img_out(:))/255);
