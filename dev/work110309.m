mammo_list = dir('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\*.mat');

[mammo_data] = sample_mammo_dt_data(...
    'num_samples', 2e3,...
    'image_dir', 'C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\',...
    'mask_dir', 'C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\',...
    'win_size', 1,...
    'num_levels', 3,...
    'feature_type', 'conj',...
    'feature_shape', 'rect',...
    'rotate', 0,...
    'do_max', 0,...
    'image_list', mammo_list(1:10), ...
    'resize_image', 0.25,...
    'use_nag', 0);
%
[background_data] = sample_mammo_dt_data(...
    'num_samples', 2e3,...
    'image_dir', 'C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\',...
    'mask_dir', 'C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\',...
    'win_size', 1,...
    'num_levels', 3,...
    'feature_type', 'conj',...
    'feature_shape', 'rect',...
    'rotate', 0,...
    'do_max', 0,...
    'image_list', mammo_list(1:10), ...
    'resize_image', 0.25,...
    'invert_mask', 1,...
    'use_nag', 0);
%
y = [true(2e3,1); false(2e3,1)];
[tree] = mb_tree_class_train([mammo_data; background_data], y);
%%
[training_data training_labels] = sample_mammo_segment_data(...
    'num_samples', 2e3,... % the mandatory arguments
    'mammo_dir', 'C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\',...
    'mask_dir', 'C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\',...
    'num_mammos', 10,...
    'view', [],...
    'win_size', 1,...
    'num_levels', 5,...
    'feature_type', 'all',...
    'feature_shape', 'rect',...
    'rotate', 0,...
    'do_max', 0,...
    'image_type', '.mat',...
    'resize', 0.25,...
    'use_nag', 0,...
    'save_path', []);
%%
ori_list = dir('Z:\data\mammograms\2004_screening\abnormals\results\306079\*.mat');
for ii = 1:20
    try
        ori_map = load_uint8(['Z:\data\mammograms\2004_screening\abnormals\results\306079\' ori_list(ii).name]);
        display(num2str(ii));
    catch
        display('She''s too big, cap''n');
        continue;
    end
    if strcmpi(ori_list(ii).name(4), 'R')
        ori_map = fliplr(ori_map);
    end
    edge_map = [true(size(ori_map,1),100) false(size(ori_map) - [0 100])] &...
        (mod(angle(ori_map),pi) > 11*pi/24) &...
        (mod(angle(ori_map),pi) < 13*pi/24);
    top_map = [true(100,size(ori_map,2)); false(size(ori_map) - [100 0])] &...
        ((mod(angle(ori_map),pi) > 23*pi/24) | (mod(angle(ori_map),pi) < pi/24));
    
    bottom_map = [false(size(ori_map) - [100 0]); true(100,size(ori_map,2))] &...
        ((mod(angle(ori_map),pi) > 23*pi/24) | (mod(angle(ori_map),pi) < pi/24));
    
    diagonal_map = [true(size(ori_map,1),100) false(size(ori_map) - [0 100])] &...
        [true(100,size(ori_map,2)); false(size(ori_map) - [100 0])] &...
        (mod(angle(ori_map),pi) > 3*pi/24) &...
        (mod(angle(ori_map),pi) < 9*pi/24);
 
    line_map = abs(ori_map);
    line_map(top_map | edge_map | bottom_map | diagonal_map) = 0;
    figure; 
    subplot(1,2,1); imagesc(abs(ori_map)); axis image;
    subplot(1,2,2); imagesc(line_map); axis image;
    
    clear *_map
end
%%
ori_map = load_uint8('Z:\data\mammograms\2004_screening\abnormals\results\306079\024RCC_class.mat');
ori_map = imresize(ori_map, 0.5, 'bilinear');
seg = u_load('C:\isbe\asymmetry_project\data\segmentations\2004_screening\abnormals\024RCC_segmentation.mat');

[rows cols] = size(ori_map);
xy = segment_breast_resize([rows cols], seg);

num_pts = size(xy,1);

% figure; 
% image(zeros(segmentation.size)); axis image; hold on; colormap(gray(256));
% plot(segmentation.breast_border(:,1),segmentation.breast_border(:,2)); axis image;
%

uv = zeros(size(xy));
for ii = 1:num_pts
    i1 = max(1, ii-3);
    i2 = min(num_pts, ii+3);
    uv(ii,:) = mean(diff(xy(i1:i2,:)));
    uv(ii,:) = uv(ii,:) / sqrt(sum(uv(ii,:).^2));
%     quiver(xy(ii,1), xy(ii,2), uv(ii,1), uv(ii,2), 20, 'r');
end
%
outer_mask = poly2mask(xy(:,1),xy(:,2),rows, cols);
outer_mask(:,1) = 0;
outer_mask(1,:) = 0;
outer_mask(:,end) = 0;
outer_mask(end,:) = 0;
inner_mask = imerode(outer_mask, strel('disk', 50));
edge_mask = outer_mask & ~inner_mask;
%
% figure; imagesc(edge_mask); axis image
%
[yi xi] = find(edge_mask);

uvi = zeros(length(yi),2);
uvi(:,1) = griddata(xy(:,1),xy(:,2),uv(:,1),xi,yi, 'nearest');
uvi(:,2) = griddata(xy(:,1),xy(:,2),uv(:,2),xi,yi, 'nearest');
%quiver(xi, yi, uvi(:,1), uvi(:,2), 'c');

edge_idx = find(edge_mask);
edge_ori = zeros(rows, cols);
edge_ori(edge_idx) = complex(uvi(:,1), -uvi(:,2)) ./ sqrt(sum(uvi.^2,2));
%%
figure; 
subplot(1,2,1); image(complex2rgb(ori_map.^2)); axis image;
subplot(1,2,2); image(complex2rgb(edge_ori.^2)); axis image;
%%
edge_mask = edge_mask & abs(angle((ori_map .* conj(edge_ori)).^2)) < pi/12;
%%
pad = 50;
left_mask = [true(size(ori_map,1),50) false(size(ori_map) - [0 50])] &...
    (mod(angle(ori_map),pi) > 11*pi/24) &...
    (mod(angle(ori_map),pi) < 13*pi/24);
right_mask = [false(size(ori_map) - [0 50]) true(size(ori_map,1),50)] &...
    (mod(angle(ori_map),pi) > 11*pi/24) &...
    (mod(angle(ori_map),pi) < 13*pi/24);
top_mask = [true(50,size(ori_map,2)); false(size(ori_map) - [50 0])] &...
    ((mod(angle(ori_map),pi) > 23*pi/24) | (mod(angle(ori_map),pi) < pi/24));
bottom_mask = [false(size(ori_map) - [50 0]); true(50,size(ori_map,2))] &...
    ((mod(angle(ori_map),pi) > 23*pi/24) | (mod(angle(ori_map),pi) < pi/24));
%%
discard_mask = edge_mask | top_mask | bottom_mask | left_mask | right_mask;  
figure; imagesc(discard_mask); axis image;
%%
line_map = abs(ori_map);
line_map(discard_mask) = 0;
figure; 
subplot(1,2,1); imagesc(abs(ori_map)); axis image;
subplot(1,2,2); imagesc(line_map); axis image;
%%
mammo = imresize(u_load('C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\abnormals\024RCC.mat'), 0.5, 'bilinear');
%
[grad_strength, grad_orientation] = gaussian_1st_derivative_gradient(double(mammo), 10);
%
grad_complex = grad_strength .* exp(i*(2*grad_orientation+pi));
figure; image(complex2rgb(grad_complex)); axis image;

grad_mask = abs(angle(ori_map.^2 .* conj(grad_complex))) < pi/12;
figure; 
subplot(1,2,1); image(complex2rgb(ori_map.^2)); axis image;
subplot(1,2,2); image(complex2rgb(grad_complex)); axis image;
%%
line_map = abs(ori_map);
line_map(grad_mask) = 0;
figure; 
subplot(1,2,1); imagesc(abs(ori_map)); axis image;
subplot(1,2,2); imagesc(line_map); axis image;
%%
discard_mask = (edge_mask & grad_mask) | top_mask | bottom_mask | left_mask | right_mask;  
figure; imagesc(discard_mask); axis image;