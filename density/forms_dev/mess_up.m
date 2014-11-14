function [image_out] = mess_up(image_in)

%Function to apply some random transformations of an image - imitating the
%effects of manually scanning

image_out = double(image_in);
image_out = (image_out - min(image_out(:))) / (max(image_out(:)) - min(image_out(:)));

% Choose whether or not to flip image;
if rand > 0.5
    image_out = rot90(image_out, 2);
end

%scale image by a factor between 0.9 and 1.1 (i.e. +/- 10%)
scale_factor = 1.1 - (rand/5);
image_out = imresize(image_out, scale_factor, 'bilinear');

%rotate by an angle between +/- 10 degrees (do something clever so we don't
%have the effect of zeros round the edge
theta = 10 - (20*rand);

add_mask = ~imrotate(ones(size(image_out)), theta, 'nearest');
image_out = imrotate(image_out, theta, 'nearest') + add_mask;


%Now add some noise - gaussian blur then speckle;
speckle = rand / 100;
g_var = rand / 100;
image_out = imnoise(image_out, 'gaussian', 0, g_var);
image_out = imnoise(image_out, 'speckle', speckle);


