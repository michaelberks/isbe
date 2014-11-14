clc;

halfwidth = 0.5;
contrast = 128;
orientation = 30 * pi/180;
row = 32; col = 32;
centre_x = 15; centre_y = 15;

[image_out, label] = create_rect_bar(halfwidth, contrast, orientation * 180/pi, ...
                                     row, col, centre_x, centre_y);
                                 
figure(1); clf; hold on; colormap(gray(256));
imagesc(image_out); axis('image','ij');

n1 = 6; n2 = 6;
t1 = linspace(0,pi,n1+1);
image_out = [];
for i = 1:length(t1)-1
    t2 = linspace(t1(i),t1(i+1),n2+1);
    bars = [];
    for j = 1:length(t2)-1
        bar_image = create_bar_image(halfwidth, centre_x, t2(j));
        bars = [bars bar_image zeros(size(bar_image,1),1)];
    end
    image_out = [image_out; bars; zeros(1,size(bar_image,2))];
end

figure(2); clf; hold on; colormap(gray(256));
imagesc(image_out); axis('image','ij');

bg = 255 * ones(size(image_out));
fg = 255 * image_out;
alpha = 0.5;
final_image = alpha * fg + ...
              (1 - alpha) * bg;
figure(2); image(final_image);

