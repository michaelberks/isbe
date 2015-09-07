num_segs = length(segments_s);
discard = false(num_segs,1);
for i_seg = 1:num_segs    
    discard(i_seg) = length(segments_s{i_seg}) < 120;
end
segments_s(discard) = [];
num_segs = length(segments_s);
%%
sequence_dir = 'N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\002wellcome\2015_02_27\R4_10_32_15\';
frame_r_sum = zeros(480,640);
frame_r_total = 0;

dirt_images = zeros(480, 640, i_seg);
for i_seg = 1:num_segs-1
    
    for i_f = segments_s{i_seg}(end)+1:segments_s{i_seg+1}(1)-1
        frame_r_sum = frame_r_sum + double(rot90(imread([sequence_dir 'frame' zerostr(i_f, 5) '.bmp']),2));
        frame_r_total = frame_r_total+1;
    end
    
    if i_seg == 1
        dirt_images(:,:,i_seg) = frame_r_sum / frame_r_total;
    else
        dirt_images(:,:,i_seg) = (frame_r_sum + frame_l_sum) / (frame_r_total + frame_l_total);
    end
    
    frame_l_sum = frame_r_sum;
    frame_l_total = frame_r_total;
    
    frame_r_sum = zeros(480,640);
    frame_r_total = 0;
    
    figure; imgray(dirt_images(:,:,i_seg));
end
dirt_images(:,:,num_segs) = frame_l_sum / frame_l_total;
figure; imgray(dirt_images(:,:,num_segs));
%%
i_seg = 4;
frame_idx = s.focused_idx{i_seg};
for i_f = 1:10
    frame = double(rot90(imread([sequence_dir 'frame' zerostr(frame_idx(i_f), 5) '.bmp']),2));
    frame_r = double(imread([sequence_dir 'registered' zerostr(i_seg, 2) '\frame' zerostr(i_f, 4) '.png']));
    figure;
    subplot(2,2,1); imgray(frame);
    subplot(2,2,2); imgray(frame - dirt_images(:,:,i_seg));
    subplot(2,2,3); imgray(frame - s.dirt_image);
    subplot(2,2,4); imgray(frame_r);
end
%%
white_lines = find(sharpness > 20 & sharpness < 46.2366);
%%
for i_f = 21:40
    frame = double(rot90(imread([sequence_dir 'frame' zerostr(white_lines(i_f), 5) '.bmp']),2));
    figure; imgray(frame);
end
%%
i_seg = 4;
    
frame_idx = segments_s{i_seg};
num_frames = length(frame_idx);

seg_x = motor_x(frame_idx);
seg_y = motor_y(frame_idx);
seg_s = sharpness(frame_idx);
seg_s_smoothed = conv(seg_s, ones(1,30)/30, 'valid');   
seg_s_cleaned = seg_s;
seg_s_cleaned(seg_s > 20) = 0;

figure; 
hold all;
plot(1:length(seg_s), seg_s', 'r');
plot(1:length(seg_s_smoothed), seg_s_smoothed', 'b');
plot(1:length(seg_s_cleaned), seg_s_cleaned', 'g');
[~, best_frames] = sort(seg_s_cleaned, 'descend');

frame = imread([sequence_dir 'frame' zerostr(frame_idx(best_frames(1)), 5) '.bmp']);
%%
frame_idx = segments_s{i_seg};
num_frames = length(frame_idx);
s1 = zeros(num_frames,1);
s2 = zeros(num_frames,1);
m1 = zeros(num_frames,1);
for i_f = 1:length(frame_idx)
    frame = double(imread([sequence_dir 'frame' zerostr(frame_idx(i_f), 5) '.bmp']));
    
    roi = frame(121:360, 161:480);
    [s1(i_f)] = image_stats(roi, 31);   
    g = gaussian_1st_derivative_gradient(roi, 2);   
    m1(i_f) = max(roi(:));
    s2(i_f) = median(g(:));
end
%%
frame_idx = segments_s{i_seg};
num_frames = length(frame_idx);
reject = false(num_frames,1);
for i_f = 1:length(frame_idx)
    frame = double(imread([sequence_dir 'frame' zerostr(frame_idx(i_f), 5) '.bmp']));
    reject(i_f) = sum(frame(:) < 100 | frame(:) > 150) > thresh;
end
%%
seg_s_cleaned = seg_s;
seg_s_cleaned(seg_s > 20 | reject) = 0;
[~, best_frames] = sort(seg_s_cleaned, 'descend');
for i_f = 1:30
    frame = double(imread([sequence_dir 'frame' zerostr(frame_idx(best_frames(i_f)), 5) '.bmp']));
    figure; imgray(frame);
end
%%
figure;
subplot(1,2,1); plot(s1, s2, 'r.');
subplot(1,2,2); plot(s1, s2 ./ m1, 'r.');
%% 
[~, best_frames] = sort(s2 ./ m1, 'descend');
thresh = 64*480;
for i_f = 1:30
    frame = double(imread([sequence_dir 'frame' zerostr(frame_idx(best_frames(i_f)), 5) '.bmp']));
    
    figure; 
    subplot(2,2,1); imgray(frame); 
    if sum(frame(:) < 100 | frame(:) > 150) > thresh
        title('Reject!!!');
    else
        title('Ok');
    end
    subplot(2,2,2); hist(frame(:), 0:16:256); 
    
    frame = double(imread([sequence_dir 'frame' zerostr(frame_idx(end-i_f), 5) '.bmp']));   
    subplot(2,2,3); imgray(frame);    
    if sum(frame(:) < 100 | frame(:) > 150) > thresh
        title('Reject!!!');
    else
        title('Ok');
    end
    subplot(2,2,4); hist(frame(:), 0:16:256); 
end
%%
for i_f = 1:10
    frame = double(imread([sequence_dir 'frame' zerostr(frame_idx(end-i_f), 5) '.bmp']));
    
    figure; 
    subplot(1,2,1); imgray(frame);    
    subplot(1,2,2); hist(frame(:), 0:16:256); 
end
%%
frame1 = double(imread([sequence_dir 'frame' zerostr(frame_idx(best_frames(9)), 5) '.bmp']));
frame2 = double(imread([sequence_dir 'frame' zerostr(frame_idx(best_frames(10)), 5) '.bmp']));

roi = frame1(121:360, 161:480);
[s1(i_f) m1(i_f)] = image_stats(roi, 31);

sobel_x = conv2([1 2 1]/4, [1 0 -1]/2, roi, 'valid');
sobel_y = conv2([1 0 -1]/2, [1 2 1]/4, roi, 'valid');
sobel_mag1 = sqrt(sobel_x.^2 + sobel_y.^2);

roi = frame2(121:360, 161:480);
[s1(i_f) m1(i_f)] = image_stats(roi, 31);

sobel_x = conv2([1 2 1]/4, [1 0 -1]/2, roi, 'valid');
sobel_y = conv2([1 0 -1]/2, [1 2 1]/4, roi, 'valid');
sobel_mag2 = sqrt(sobel_x.^2 + sobel_y.^2);

line_strength1 = gaussian_1st_derivative_gradient(frame1, 1);
line_strength2 = gaussian_1st_derivative_gradient(frame2, 1);

figure;
subplot(2,2,1); imgray(frame1);
subplot(2,2,2); imgray(frame2);
subplot(2,2,3); imgray(line_strength1);
subplot(2,2,4); imgray(line_strength2);