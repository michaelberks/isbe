num = xlsread('C:\isbe\density\view_distribution.xls', 'm11:m130');
num = kron(num', ones(40, 5));

bar_image_r = zeros(size(num));
bar_image_g = zeros(size(num));
bar_image_b = zeros(size(num));

%do the blue
bar_image_r(num == 1) = 0;
bar_image_g(num == 1) = 0;
bar_image_b(num == 1) = 1;

%do the yellow
bar_image_r(num == 2) = 1;
bar_image_g(num == 2) = 1;
bar_image_b(num == 2) = 0;

%do the green
bar_image_r(num == 3) = 0;
bar_image_g(num == 3) = 1;
bar_image_b(num == 3) = 0;

%do the red
bar_image_r(num == 4) = 1;
bar_image_g(num == 4) = 0;
bar_image_b(num == 4) = 0;

%Write the image
imwrite(cat(3, bar_image_r, bar_image_g, bar_image_b), 'C:\isbe\density\breast_volume_bar.bmp');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%