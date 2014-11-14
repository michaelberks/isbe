    
for ii = 1:16
    f1 =  imread(['I:\nailfold\Visit1\Lhand\Digit4\X300\Scene1\2V1LD4X3S1F' num2str(ii) '.bmp']);
    mask = f1 < 100;
    f1(mask) = mean(f1(~mask));

    [nailfold_strength, nailfold_orientation] = ...
                gaussian_2nd_derivative_line2(f1, [4 8]);
    nailfold_nms = mb_non_maximal_supp(nailfold_strength, nailfold_orientation);

    [nailfold_h] = hysterisis(nailfold_nms);
    nailfold_bao = bwareaopen(nailfold_nms > 0, 10);

    [y_b x_b] = find(nailfold_bao);
    [y_h x_h] = find(nailfold_h);

    figure; imagesc(f1); axis image; colormap(gray(256)); hold on;
    plot(x_b, y_b, 'r.', 'markersize', 2);
    plot(x_h, y_h, 'g.', 'markersize', 2);
end
%%
for ff = 1:3
    f1 =  imread(['C:\isbe\nailfold\ncm\Visit1\Lhand\Digit4\X300\Scene1\2V1LD4X3S1F' num2str(ff) '.bmp']);
    f2 =  imread(['C:\isbe\nailfold\ncm\Visit1\Lhand\Digit4\X300\Scene1\2V1LD4X3S1F' num2str(ff+1) '.bmp']);
    mask1 = f1 < 100;
    mask2 = f2 < 100;
    f1(mask1) = mean(f1(~mask1));
    f2(mask2) = mean(f2(~mask2));

    [nailfold_strength1, nailfold_orientation1] = ...
        gaussian_2nd_derivative_line2(f1, [4 8]);
    nailfold_nms1 = mb_non_maximal_supp(nailfold_strength1, nailfold_orientation1);
    [nailfold_h1] = hysterisis(nailfold_nms1);
    [y_h1 x_h1] = find(nailfold_h1);


    [nailfold_strength2, nailfold_orientation2] = ...
        gaussian_2nd_derivative_line2(f2, [4 8]);
    nailfold_nms2 = mb_non_maximal_supp(nailfold_strength2, nailfold_orientation2);
    [nailfold_h2] = hysterisis(nailfold_nms2);
    [y_h2 x_h2] = find(nailfold_h2);

    tic;
    range = -11:11;
    t_sz = length(range) - 2;
    max_val = 0;
    thetas = pi*(-15:15)/180;
    r_sz = length(thetas);
    
    xc = size(f1,2)/2;
    yc = size(f1,1)/2;
    max_data = zeros(r_sz,3);
    for tt = 1:r_sz
        theta = thetas(tt);
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        xy_r = [x_h1-xc y_h1-yc]*R;

        translations = zeros(t_sz);
        for ii = 1:length(x_h1);
            x_offset = x_h2 - xy_r(ii,1) - xc;
            y_offset = y_h2 - xy_r(ii,2) - yc;
            hough_counts = hist3([y_offset x_offset], {range, range});

            translations = translations + hough_counts(2:end-1, 2:end-1);
        end
        [max_val_theta max_id_theta] = max(translations(:));
%         [max_row_theta max_col_theta] = ind2sub([t_sz t_sz], max_id_theta);
%         max_x_theta = range(max_col_theta);
%         max_y_theta = range(max_row_theta);

        if max_val_theta > max_val
            max_val = max_val_theta;
            max_id = max_id_theta;
            max_theta = theta;
        end
%         max_data(tt,:) = [max_val_theta max_x_theta max_y_theta];
    end
    toc;

    [max_row max_col] = ind2sub([t_sz t_sz], max_id);
    max_x = range(max_col);
    max_y = range(max_row);
    display([ff max_x max_y round(180*max_theta/pi)]);
end
%% 

theta = 0;
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
xy_r = [x_h1-xc y_h1-yc]*R;

tic;
translations = zeros(t_sz);
for ii = 1:length(x_h1);
    x_offset = x_h2 - xy_r(ii,1) - xc;
    y_offset = y_h2 - xy_r(ii,2) - yc;
    hough_counts = hist3([y_offset x_offset], {range, range});

    translations = translations + hough_counts(2:end-1, 2:end-1);
end
toc;
tic;
translations2 = zeros(t_sz);
for ii = 1:length(x_h1);
    x_offset = round(x_h2 - xy_r(ii,1) - xc);
    y_offset = round(y_h2 - xy_r(ii,2) - yc);
    
    keep = abs(x_offset) <= 10 & abs(y_offset) <= 10;

    translations2 = translations2 + sparse(y_offset(keep)+11, x_offset(keep)+11, 1, t_sz, t_sz);
end
toc;
    
%%
for ff = 1:3
    f1 =  imread(['C:\isbe\nailfold\ncm\Visit1\Lhand\Digit4\X300\Scene1\2V1LD4X3S1F' num2str(ff) '.bmp']);
    f2 =  imread(['C:\isbe\nailfold\ncm\Visit1\Lhand\Digit4\X300\Scene1\2V1LD4X3S1F' num2str(ff+1) '.bmp']);
    mask1 = f1 < 100;
    mask2 = f2 < 100;
    f1(mask1) = mean(f1(~mask1));
    f2(mask2) = mean(f2(~mask2));
    
    tic;
    [offset, theta] = register_frame(f1, f2);
    toc;
    
    display([ff offset theta]);
end

f1a = zeros(size(f1));
f1a(3:end, 4:end) =  f1(1:end-2, 1:end-3);
mask1a = f1a < 100;
f1a(mask1) = mean(f1a(~mask1));
[offset, theta] = register_frame(f1, f1a);
display([offset, theta]);
%%
f1b = imrotate(f1, 2, 'bicubic', 'crop');
mask1b = f1b < 100;
f1b(mask1) = mean(f1b(~mask1));
[offset, theta] = register_frame(f1, f1b);
display([offset, theta]);
%%
f1c = imrotate(f1a, 2, 'bicubic', 'crop');
mask1c = f1c < 100;
f1c(mask1) = mean(f1c(~mask1));
[offset, theta] = register_frame(f1, f1c);
display([offset, theta]);
%%
offsets = zeros(15,2);
thetas = zeros(15,1);
for ff = 1:15
    f1 =  imread(['C:\isbe\nailfold\ncm\Visit1\Lhand\Digit4\X300\Scene1\2V1LD4X3S1F' num2str(ff) '.bmp']);
    f2 =  imread(['C:\isbe\nailfold\ncm\Visit1\Lhand\Digit4\X300\Scene1\2V1LD4X3S1F' num2str(ff+1) '.bmp']);
    mask1 = f1 < 100;
    mask2 = f2 < 100;
    f1(mask1) = mean(f1(~mask1));
    f2(mask2) = mean(f2(~mask2));
    
    [offsets(ff,:) thetas(ff)] = register_frame(f1, f2);
end
f1 =  imread(['C:\isbe\nailfold\ncm\Visit1\Lhand\Digit4\X300\Scene1\2V1LD4X3S1F' num2str(1) '.bmp']);
mask1 = f1 < 100;
f1(mask1) = mean(f1(~mask1));
[offsets1 thetas1] = register_frame(f1, f2);
%%
f1 =  imread(['C:\isbe\nailfold\ncm\Visit1\Lhand\Digit4\X300\Scene1\2V1LD4X3S1F' num2str(1) '.bmp']);
f2 =  imread(['C:\isbe\nailfold\ncm\Visit1\Lhand\Digit4\X300\Scene1\2V1LD4X3S1F' num2str(2) '.bmp']);
mask1 = f1 < 100;
mask2 = f2 < 100;
f1(mask1) = mean(f1(~mask1));
f2(mask2) = mean(f2(~mask2));

[offset theta] = register_frame(f1, f2);
frames = add_frame(f1, f2, offset, theta);    
%%

frame_tgt =  imread(['C:\isbe\nailfold\ncm\Visit1\Lhand\Digit4\X300\Scene1\2V1LD4X3S1F' num2str(1) '.bmp']);    

offsets2 = zeros(15,2);
thetas2 = zeros(15,1);
for ff = 1:15
    
    frame_src =  imread(['C:\isbe\nailfold\ncm\Visit1\Lhand\Digit4\X300\Scene1\2V1LD4X3S1F' num2str(ff+1) '.bmp']);    
    [offsets2(ff,:) thetas2(ff)] = register_frame(frame_src, frame_tgt);
end
%%
for ss = 3:9
    str = num2str(ss);
    ztr = zerostr(ss,2);
    %mkdir(['C:\isbe\nailfold\ncm\Visit1\Lhand\Digit4\X300\Scene' ztr]);
    for ff = 10:16
        movefile(...
            ['C:\isbe\nailfold\ncm\Visit1\Lhand\Digit4\X300\Scene' ztr '\2V1LD4X3S' str 'F' zerostr(ff,1) '.bmp'],...
            ['C:\isbe\nailfold\ncm\Visit1\Lhand\Digit4\X300\Scene' ztr '\2V1LD4X3S' ztr 'F' zerostr(ff,2) '.bmp']);    
    end
%     movefile(...
%         ['C:\isbe\nailfold\ncm\Visit1\Lhand\Digit4\X300\Scene' str '\2V1LD4X3S' str 'Fcomp.bmp'],...
%         ['C:\isbe\nailfold\ncm\Visit1\Lhand\Digit4\X300\Scene' ztr '\2V1LD4X3S' ztr 'Fcomp.bmp']);  
end
%%
for ss = 3:11;
    str = zerostr(ss,2);
    for ff = 1:16
        frame =  imread(['C:\isbe\nailfold\ncm\Visit1\Lhand\Digit4\X300\Scene' str '\2V1LD4X3S' str 'F' zerostr(ff,2) '.bmp']);
        frames(:,:,ff) =  frame(2:end, 4:end);
    end

    [frames_reg offsets thetas] = register_frames2(frames, 50, 5);
    discard = any(frames_reg < 100, 3);
    frames_reg_std = std(double(frames_reg),0,3);
    frames_reg_std(discard) = 0;
    
    frame_reg = mean(frames_reg,3) - frames_reg_std;

    discard = any(frames < 100, 3);
    frames_std = std(double(frames),0,3);
    frames_std(discard) = 0;
    frame = mean(frames,3) - frames_std;

    clims = [0 max(max(frames_std(:)), max(frames_reg_std(:)))];

    figure('Name', ['Scene ' num2str(ss)]); 
    subplot(2,2,1); imagesc(frames_std); axis image; caxis(clims);
    subplot(2,2,2); imagesc(frame); axis image;
    subplot(2,2,3); imagesc(frames_reg_std); axis image; caxis(clims);
    subplot(2,2,4); imagesc(frame_reg); axis image;
    
    display(['Scene ' num2str(ss)]);
    display([offsets thetas]);
    
    save(['C:\isbe\nailfold\ncm\Visit1\Lhand\Digit4\X300\mb_reg\frames' str '.mat'], 'frames_reg');
    save(['C:\isbe\nailfold\ncm\Visit1\Lhand\Digit4\X300\mb_reg\trans' str '.mat'], 'offsets', 'thetas');
end
%%
for ii = 1:16
    figure; imagesc(frames(2:end,4:end,ii)); axis image; colorbar;
end
%%
for ii = 2:16
    figure; 
    subplot(2,1,1); imagesc((double(frames(2:end,4:end,1)) + double(frames(2:end,4:end,ii)))/2); axis image; colorbar;
    subplot(2,1,2); imagesc((double(frames(2:end,4:end,ii-1)) + double(frames(2:end,4:end,ii)))/2); axis image; colorbar;
end
%%
figure; imagesc(mean(frames_reg,3)); axis image;
figure; imagesc(frame_reg); axis image;
%%
for ii = 1:8
    [line_strength, line_orientation] = ...
        gaussian_2nd_derivative_line2(frames(:,:,1), ii);
    
    figure; imagesc(line_strength); axis image;
    caxis([min(min(line_strength(10:end-9,10:end-9))) max(max(line_strength(10:end-9,10:end-9)))]); colorbar;
end
%%

%Apply non-maximal suppression to skeletonise the line strength map
line_nms = mb_non_maximal_supp(line_strength, line_orientation);

%Discard points from edges of map
line_nms([1:32 end-32:end], :) = 0;
line_nms(:, [1:32 end-32:end]) = 0;

%Apply hysterisis to select suitable lines from the NMS map
[line_mask] = hysterisis(line_nms);

%Extract (x,y) coordinates of the remaining lines
[y_h x_h] = find(line_mask);
[y_n x_n] = find(line_nms);

figure; imagesc(frames(:,:,1)); axis image; colormap(gray(256)); hold on; caxis([100 255]);
plot(x_n, y_n, 'r.', 'markersize', 2);
plot(x_h, y_h, 'g.', 'markersize', 2);
%%
for ii = 0:15:165
    
    frame = create_gauss_bar(4, -1, ii, 256, 256, 128, 128);
    [line_strength, line_orientation] = ...
        gaussian_2nd_derivative_line2(frame, 8);

    %Apply non-maximal suppression to skeletonise the line strength map
    line_nms = mb_non_maximal_supp(line_strength, line_orientation);

    %Discard points from edges of map
    line_nms([1:32 end-32:end], :) = 0;
    line_nms(:, [1:32 end-32:end]) = 0;

    %Apply hysterisis to select suitable lines from the NMS map
    [line_mask] = hysterisis(line_nms);

    %Extract (x,y) coordinates of the remaining lines
    [y_h x_h] = find(line_mask);
    [y_n x_n] = find(line_nms);

    figure; imagesc(frame); axis image; colormap(gray(256)); hold on;
    plot(x_n, y_n, 'r.', 'markersize', 2);
    plot(x_h, y_h, 'g.', 'markersize', 2);
end
%%
frame = ones(256, 32);
old_c = 0.5;
for ii = 1:8
    frame = [frame zeros(256,ii) ones(256,32)]; %#ok
    line_centres(ii) = old_c + 32 + ii/2 + (ii-1)/2; %#ok
    old_c = line_centres(ii);
end
[rows cols] = size(frame);
sample_points = sub2ind([rows cols], 128*ones(1,8), line_centres);
frame = ones(rows, cols);
for ii = 1:8
    frame = frame +...
        create_ellipse_bar(ii, -1, 90, rows, cols, line_centres(ii), 128);
        %create_gauss_bar(ii, -1, 90, rows, cols, line_centres(ii), 128);
end
%
[line_strength, line_orientation line_scale] = ...
    gaussian_2nd_derivative_line(frame, 1:8);
    
% [line_strength line_orientation line_scale] = ...
%     line_operator_conv(1-frame, 8, 1, 9,'radians');
%%
%Apply non-maximal suppression to skeletonise the line strength map
line_nms = mb_non_maximal_supp(line_strength, line_orientation);

%Apply hysterisis to select suitable lines from the NMS map
[line_mask] = hysterisis(line_nms);

%Extract (x,y) coordinates of the remaining lines
[y_h x_h] = find(line_mask);
[y_n x_n] = find(line_nms);

figure; imagesc(frame); axis image; colormap(gray(256)); hold on;
plot([line_centres; line_centres], [zeros(1,8); 256*ones(1,8)], 'b');
%plot(x_n, y_n, 'r.', 'markersize', 2);
%plot(x_h, y_h, 'g.', 'markersize', 2);

figure; imagesc(line_strength); axis image; colormap(gray(256)); hold on;
plot([line_centres; line_centres], [zeros(1,8); 256*ones(1,8)], 'b');
figure; plot(1:cols, line_strength(128,:), 'b'); hold on;
plot([line_centres; line_centres], [min(line_strength(128,:))*ones(1,8); max(line_strength(128,:))*ones(1,8)], 'r--');

figure; imagesc(line_orientation); axis image; colormap(hsv(256)); hold on; caxis([-pi/2 pi/2]);
plot([line_centres; line_centres], [zeros(1,8); 256*ones(1,8)], 'b');

% figure; imagesc(round(log(line_scale)/log(2))); axis image; caxis([0 3]); 
% colormap([1 0 0; 0 1 0; 0 0 1; 0 0 0]); hold on; colorbar;
% plot([line_centres; line_centres], [zeros(1,8); 256*ones(1,8)], 'y');

figure; imagesc(line_scale); axis image;
colormap(jet(length(unique(line_scale(128,:))))); hold on; colorbar;
plot([line_centres; line_centres], [zeros(1,8); 256*ones(1,8)], 'k--');
%%
figure; hold all;
max_line = zeros(1,cols);
scales = zeros(1,cols);
for ii = 1:8
    line_strength = ...
        gaussian_2nd_derivative_line2(frame, ii);
    line_strength = line_strength(128,:);
    
    subplot(2,4,ii); hold on;
    plot([line_centres; line_centres],...
        [min(line_strength)*ones(1,8);...
        max(line_strength)*ones(1,8)], 'r--');
    plot(1:cols, line_strength);
    axis([0 300 0 0.4]);
    
    max_vals(ii) = max(line_strength);
    swap = line_strength > max_line;
    max_line(swap) = line_strength(swap);
    scales(swap) = ii;
end
figure; plot(1:cols, max_line); hold on;
plot([line_centres; line_centres],...
        [min(max_line)*ones(1,8);...
        max(max_line)*ones(1,8)], 'r--');
    
figure; plot(1:cols, scales); hold on;
plot([line_centres; line_centres],[zeros(1,8); 8*ones(1,8)], 'r--');
%%
r = linspace(0, 5*pi, 1000);
figure; plot(r.*cos(r), r.*sin(r));

spiral_x = 10*r.*cos(r);
spiral_y = 10*r.*sin(r);
spiral_cols = round(spiral_x) - min(round(spiral_x)) + 16;
spiral_rows = round(spiral_y) - min(round(spiral_y)) + 16;
spiral_bw = false(max(spiral_rows)+16, max(spiral_cols)+16);
spiral_bw(sub2ind(size(spiral_bw), spiral_rows, spiral_cols)) = 1;
figure; imagesc(spiral_bw); axis image;

halfwidth = 4;
dx = bwdist(spiral_bw);
sigma2 = (halfwidth^2) / log(2);
spiral_bar = 1 - exp(-(dx.^2 / sigma2));

ex = halfwidth*halfwidth - dx.*dx;
ex(ex<=0) = 0;

image_out = (contrast/halfwidth)*sqrt(ex);
%%
for ii = 3:8
    
    %test_im = spiral_bar;
    test_im = frames(:,:,1);
%     rng = 101:105;
%     [str1, ori1] = ...
%         gaussian_2nd_derivative_line(test_im, ii);
%     [str2, ori2] = ...
%         gaussian_clover_line(test_im, ii);
%     
%     disp('-----');
%     str1(rng,rng)
%     str2(rng,rng)
%     ori1(rng,rng)-ori2(rng,rng)
%     return
    
    [r c] = size(test_im);
    [line_strength, line_orientation] = ...
        gaussian_2nd_derivative_line2(test_im, ii);
%         gaussian_2nd_derivative_line(test_im, ii);
    
    figure; imagesc(line_strength); axis image;
    caxis([min(min(line_strength(40:end-39,40:end-39))) max(max(line_strength(40:end-39,40:end-39)))]); colorbar;
    hold on;
    quiver(32:4:c-31, 32:4:r-31,...
        line_strength(32:4:r-31,32:4:c-31).*cos(line_orientation(32:4:r-31,32:4:c-31)),...
       -line_strength(32:4:r-31,32:4:c-31).*sin(line_orientation(32:4:r-31,32:4:c-31)), 'w');
    
end
%%
load(['C:\isbe\nailfold\ncm\Visit1\Lhand\Digit4\X300\mb_reg\frames' zerostr(1,2) '.mat'], 'frames_reg');

discard = any(frames_reg < 10, 3);
frames_reg(discard) = 0;
frames_reg_std = std(double(frames_reg),0,3);
frames_reg_std(discard) = 0;

mosaic_frames = mean(frames_reg,3) - frames_reg_std; clear frames;
mosaic_masks = ~discard;
for ss = 2:11
    load(['C:\isbe\nailfold\ncm\Visit1\Lhand\Digit4\X300\mb_reg\frames' zerostr(ss,2) '.mat'], 'frames_reg');

    %figure; imagesc(all(frames_reg,3)); axis image;
    
    discard = any(frames_reg < 10, 3);
    frames_reg(discard) = 0;
    frames_reg_std = std(double(frames_reg),0,3);
    frames_reg_std(discard) = 0;

    frame_reg = mean(frames_reg,3) - frames_reg_std; clear frames;
    
    mosaic_frames = add_frame(frame_reg, mosaic_frames, [0 0], 0);
    mosaic_masks = add_frame(~discard, mosaic_masks, [0 0], 0);
    
    
end
mosaic_frames = flipdim(mosaic_frames, 3);
mosaic_masks = flipdim(mosaic_masks, 3);
clear discard frame*
%
figure;
for ii = 1:6
    frame = mosaic_frames(:,:,ii);
    subplot(3,2,ii); imagesc(frame); axis image; colormap(gray(256)); hold on;
    caxis([min(frame(mosaic_masks(:,:,ii))) max(frame(mosaic_masks(:,:,ii)))]);
end
figure;
for ii = 7:11
    frame = mosaic_frames(:,:,ii);
    subplot(3,2,ii-6); imagesc(frame); axis image; colormap(gray(256)); hold on;
    caxis([min(frame(mosaic_masks(:,:,ii))) max(frame(mosaic_masks(:,:,ii)))]);
end
    
%%
profile on;
[mosaic offsets thetas] = register_frames2(mosaic_frames, max(size(mosaic_frames)), 5, 1, 0);
profile viewer;
display([offsets thetas]);
%%

offset_lim = [0 size(mosaic_frames,2); -floor(size(mosaic_frames,1)/2) floor(size(mosaic_frames,1)/2)];
%offset_lim = max(size(mosaic_frames));

for ii = 1:10
    [mosaic offsets thetas] = register_tiles(mosaic_frames(:,:,ii:ii+1), mosaic_masks(:,:,ii:ii+1), offset_lim, 0, 0, 1);
    display([offsets thetas]);
end
%%
for ii = 1:10
    frame = sum(mosaic(:,:,ii:ii+1), 3) ./ sum(mosaic(:,:,ii:ii+1) > 50, 3);
    [r c] = size(frame);
    sr = find(isnan(frame(1:r/2,c/2)), 1, 'last');
    er = find(isnan(frame(r/2+1:end,c/2)), 1, 'first') + r/2 - 1;
    sc = find(isnan(frame(r/2,1:c/2)), 1, 'last');
    ec = find(isnan(frame(r/2,c/2+1:end)), 1, 'first') + c/2 - 1;
    frame = frame(sr:er, sc:ec);
    
    figure; 
    subplot(3,2,1); imagesc(mosaic_frames(:,:,ii)); axis image;
    subplot(3,2,2); imagesc(mosaic_frames(:,:,ii+1)); axis image;
    subplot(3,2,3:6); imagesc(frame); axis image;
end
%%
figure;
for ii = 1:5
    subplot(5,1,ii); imagesc(mosaic(:,:,ii)); axis image;
end
figure;
for ii = 6:11
    subplot(6,1,ii-5); imagesc(mosaic(:,:,ii)); axis image;
end
%%

for ss = 1:11
    [frames] = load_frames(['C:\isbe\nailfold\ncm\Visit1\Lhand\Digit4\X300\Scene' zerostr(ss, 2) '\']);
    [frames_reg offsets thetas] = register_frames2(frames, 50, 5); %#ok
    
    save(['C:\isbe\nailfold\ncm\Visit1\Lhand\Digit4\X300\mb_reg\frames' zerostr(ss, 2) '.mat'], 'frames_reg');
    save(['C:\isbe\nailfold\ncm\Visit1\Lhand\Digit4\X300\mb_reg\trans' zerostr(ss, 2) '.mat'], 'offsets', 'thetas');
    
    discard = any(frames_reg < 10, 3);
    frames_reg(discard) = 0;
    frames_reg_std = std(double(frames_reg),0,3);
    frames_reg_std(discard) = 0;

    frame_reg = mean(frames_reg,3) - frames_reg_std; clear frames;
    
    if ss == 1
        mosaic_frames = frame_reg;
        mosaic_masks = ~discard;
    else
        mosaic_frames = add_frame(frame_reg, mosaic_frames, [0 0], 0);
        mosaic_masks = add_frame(~discard, mosaic_masks, [0 0], 0);
    end 
end
tiles = flipdim(mosaic_frames, 3);
tile_masks = flipdim(mosaic_masks, 3);
clear discard frame*

[offsets thetas] = register_tiles(tiles, tile_masks, [], 0, 0, 1);
%%
for ss = 1:11
    load(['C:\isbe\nailfold\ncm\Visit1\Lhand\Digit4\X300\mb_reg\frames' zerostr(ss, 2) '.mat'], 'frames_reg');
    load(['C:\isbe\nailfold\ncm\Visit1\Lhand\Digit4\X300\mb_reg\trans' zerostr(ss, 2) '.mat'], 'offsets', 'thetas');
    
    discard = any(frames_reg < 10, 3);
    frames_reg(discard) = 0;
    frames_reg_std = std(double(frames_reg),0,3);
    frames_reg_std(discard) = 0;

    frame_reg = mean(frames_reg,3) - frames_reg_std; clear frames;
    
    if ss == 1
        mosaic_frames = frame_reg;
        mosaic_masks = ~discard;
    else
        mosaic_frames = add_frame(frame_reg, mosaic_frames, [0 0], 0);
        mosaic_masks = add_frame(~discard, mosaic_masks, [0 0], 0);
    end 
end
tiles = flipdim(mosaic_frames, 3);
tile_masks = flipdim(mosaic_masks, 3);
clear discard frame*

[offsets thetas] = register_tiles(tiles, tile_masks, [], 0, 0, 1);
 
    
    