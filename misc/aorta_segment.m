%Load image
im = imread('P:\isbe\misc\aifDetection.png');

%Take raw image to be green channel - although would be better if had
%original image without marks (so the top and bottom of the snake contours
%aren't distorted by the marks)
img = im(:,:,2);

%Workout the location of the marked lines
marked_lines = bwlabel((im(:,:,1) - img) > 250);
[line1_y, line1_x] = find(marked_lines==1);
[line2_y, line2_x] = find(marked_lines==2);

%Take a rectangular ROI from the bounding box provided by the lines
x1 = min([line1_x; line2_x]); 
x2 = max([line1_x; line2_x]); 
y1 = min([line1_y; line2_y]); 
y2 = max([line1_y; line2_y]);

aorta_roi = img(y1:y2, x1:x2);
nc = x2-x1+1;
nr = y2-y1+1;

%
%Compute feature image as simple horizontal edges. For the left edge, negate
%the image
feature_im_r = imfilter(double(aorta_roi), fspecial('sobel')', 'replicate');
feature_im_r = feature_im_r / max(abs(feature_im_r(:))); %Scale between [-1 1]

%Initialise the snakes to have one control point per row. Start the left on
%the left edge and right on the right edge, although as we'll search
%exhaustively it doesn't matter what x coordinates we use
snake_l = [ones(nr,1) (1:nr)']; 
snake_r = [nc*ones(nr,1) (1:nr)'];
%    
alpha = 0; %Don't really need alpha, but can try positive values
beta = 0.1; %Experiment with beta in powers of 10 to get rough idea, then fine tune
%If you set beta to zero, the edge will lock onto the strongest edge point in
%each row, regardless of how zigzaggy this makes the snake. If you set beta
%very large, it will be straight line, regardless of the edge strength

max_delta_x = nc-1; %Allow each x-coordinate to move anywhere in ROI
resol_x = 1;
max_delta_y = 0; %Fix each y-coordinate to stay in its row
resol_y = 1;

%Fit the snakes
[snake_l, e_l] = mb_snake(snake_l, alpha, beta, max_delta_x, resol_x, ...
	max_delta_y, resol_y, -feature_im_r);
[snake_r, e_r] = mb_snake(snake_r, alpha, beta, max_delta_x, resol_x, ...
	max_delta_y, resol_y, feature_im_r);

%Check they've locked on to the appropriate edges
figure; 
subplot(1,2,1); imgray(abs(feature_im_r));
plot(snake_l(:,1), snake_l(:,2), 'g');
plot(snake_r(:,1), snake_r(:,2), 'r');
subplot(1,2,2); imgray(aorta_roi);
plot(snake_l(:,1), snake_l(:,2), 'g');
plot(snake_r(:,1), snake_r(:,2), 'r');

if any(snake_l(:,1) > snake_r(:,1))
    display('Snakes have crossed... bugger!');
end

%Flip the right edge upside down to form a continous contour with the left
%edge and make an image mask
aorta_boundary = [snake_l; flipud(snake_r)];
aorta_boundary(:,1) = aorta_boundary(:,1) + x1 - 1;
aorta_boundary(:,2) = aorta_boundary(:,2) + y1 - 1;
aorta_mask = poly2mask(aorta_boundary(:,1), aorta_boundary(:,2), size(img,1), size(img,2));

figure; 
subplot(1,2,1); imgray(im);
plot([x1 x2 x2 x1 x1], [y1 y1 y2 y2 y1], 'b--');
plot(aorta_boundary(:,1), aorta_boundary(:,2), 'g');
subplot(1,2,2); imgray(aorta_mask);