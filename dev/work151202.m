frames_hist = zeros(10000, 256);
bins = 0:255;

curr_frame = 1;
for i_seq = 1:1000
    
    for i_frame = 1:10
        frame = imread([sequence_dirs{i_seq} 'frame' zerostr(i_frame,5) '.bmp']);
        frame([1:120 361:480],:) = [];
        frames_hist(curr_frame,:) = hist(frame(:), bins);
        curr_frame = curr_frame + 1;
    end
end
%%
for i_im = 1:5
    figure; 
    subplot(2,1,1); imgray(frames_hist((i_im-1)*1000 + (1:1000),:)');
    subplot(2,1,2); imgray(frames_hist((i_im+4)*1000 + (1:1000),:)');
end
%%
frames_hist_sum_l = cumsum(frames_hist,2);
frames_hist_sum_r = fliplr(cumsum(fliplr(frames_hist),2));
num_p = 240*640;
%%
pcts = [0.005 0.01 0.02 0.08 0.16];
for i_p = 1:5
    figure; 
    subplot(1,2,1);
    plot(sum(frames_hist_sum_l > pcts(i_p)*num_p) / 100, 'r.');
    subplot(1,2,2);
    plot(sum(frames_hist_sum_r > pcts(i_p)*num_p) / 100, 'r.');
end
%%
thresh = 200;
pct = 0.01;

for i_seq = 1:1000
    if frames_hist_sum_r(i_seq*10, thresh) > pct*num_p
        frame = imread([sequence_dirs{i_seq} 'frame00010.bmp']);
        figure; imgray(frame);
        title(['Bad pixels = ' num2str(frames_hist_sum_r(i_seq*10, thresh)/num_p, 2)]);
    end
end
%%
thresh = 64;
pct = 0.16;

for i_seq = 1:1000
    if frames_hist_sum_l(i_seq*10, thresh) > pct*num_p
        frame = imread([sequence_dirs{i_seq} 'frame00010.bmp']);
        figure; imgray(frame);
        title(['Bad pixels = ' num2str(frames_hist_sum_r(i_seq*10, thresh)/num_p, 2)]);
    end
end