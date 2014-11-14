clc;

imsz = 1024;
im = round(rand(imsz,imsz)*255);

tic;
dt = compute_dual_tree(im, 5);
toc;

tic;
ft = dt_to_pixel_subset(dt, 1:100, 1:100, [1], 'linear');
toc;

return

im = 127*ones(8,8); im(:,4:5) = 255;
% im = zeros(16,16); im(:,7:10) = 255;

for i = 1:size(im,1)
    for j = 1:size(im,2)
        im(i,j) = i+j-2;
    end
end

dt = compute_dual_tree(im, 2);

ft = dt_to_pixel_subset(dt, [4], [4], [1], 'linear');
% ft = dt_to_pixel_subset(dt, [8], [8], [2], 'linear');

display(ft(:))
