for i_im = 1:20
    frame = ...
        double(rot90(imread([args.sequence_dir 'frame' zerostr(focused_idx_i(i_im), 5) '.bmp']),2));
        figure; imgray(frame);

end
%%
[~, best_frames] = sort(seg_s_cleaned, 'descend');
thresh = 480*640*0.2;
target_frames = min(30,round(num_frames/4));

frames = zeros(args.frame_h,args.frame_w,target_frames);
frame_count = 0;
i_im = 1;
while frame_count < target_frames && i_im <= num_frames
    frame = double(rot90(imread([args.sequence_dir 'frame' zerostr(frame_idx(best_frames(i_im)), 5) '.bmp']),2));       
    if sum(frame(:) < 100 | frame(:) > 150) < thresh;
        frame_count = frame_count + 1;
        frames(:,:,frame_count) = frame  - args.dirt_image;       
        figure; imgray(frame);
        
    elseif args.debug
        %figure; imgray(frame);
    end
    %Increment frame index to load in next frame
    i_im = i_im + 1;
end
%%
root_dir = 'N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\';

session_dir = [root_dir '002wellcome\2015_02_27\'];
sequences_s = dir(session_dir);

num_seqs = length(sequences_s) - 2;
sequence_dir = cell(num_seqs,1);
for i_seq = 1:num_seqs
    sequence_dir{i_seq} = [session_dir sequences_s(i_seq+2).name '\'];
end

thresh = 1000;%480*640*0.01;
lower_g_lim = 0.25;
upper_g_lim = 0.75;

good_dir = ('C:\isbe\nailfold\good_frames\');
bad_dir = ('C:\isbe\nailfold\bad_frames\');
    
create_folder(good_dir);
create_folder(bad_dir);
g_lims = [128 128];
g_range = g_lims(2) - g_lims(1);
%%
for i_seq = 1:10
    sequence = u_load([sequence_dir{i_seq} 'sequences.mat']);
    [segments_s, segments_ns, motor_x, motor_y, motor_z, sharpness] = ...
        get_stationary_segments(sequence, 120);
    num_segs = length(segments_s);

    for i_seg = 1:num_segs

        frame_idx = segments_s{i_seg};
        seg_s = sharpness(frame_idx);
        seg_s(seg_s > 50) = 0;
        [~, max_i] = max(seg_s);
        frame_name = ['frame' zerostr(frame_idx(max_i), 5) '.bmp'];
        frame_u = imread([sequence_dir{i_seq} frame_name]);  
        frame = double(frame_u);
        %lims_f(1) = min(frame(:));
        %lims_f(2) = max(frame(:));
        
        %g_lims(1) = min(g_lims(1), lims_f(1));
        %g_lims(2) = max(g_lims(2), lims_f(2));

%         frame_n = (frame - g_lims(1)) / g_range; 
% 
%         %is_valid = sum(frame_n(:) < lower_g_lim | frame_n(:) > upper_g_lim) < thresh;
%         is_valid = sum(frame_n(:) > upper_g_lim) < thresh;

        %frame_r = imresize(frame, 1/16);
        %[dip, p] = HartigansDipSignifTest(frame_r(:), 500);
        %is_valid = p > 0.05;
        is_valid = sum(frame(:) > 200) < thresh;
        if is_valid
            frame_r = imresize(frame, 1/16);
            [dip, p] = HartigansDipSignifTest(frame_r(:), 500);
            is_valid = p > 0.0001;
        end
        
        if is_valid
            imwrite(frame_u, [good_dir zerostr(i_seq,2) frame_name]);        
        else
            imwrite(frame_u, [bad_dir zerostr(i_seq,2) frame_name]); 
        end
    end
end
%%
frame_list = dir([bad_dir '*.bmp']);
lower_g_lim = 100;
upper_g_lim = 200;
n = 480*640;
for i_f = 1:length(frame_list)
    frame_u = imread([bad_dir frame_list(i_f).name]);  
    frame = double(frame_u);
    %lims_f(1) = min(frame(:));
    %lims_f(2) = max(frame(:));

    %g_lims(1) = min(g_lims(1), lims_f(1));
    %g_lims(2) = max(g_lims(2), lims_f(2));

    %frame_n = (frame - g_lims(1)) / g_range; 

    %is_valid = sum(frame_n(:) < lower_g_lim | frame_n(:) > upper_g_lim) < thresh;

    lower_p = sum(frame(:) < lower_g_lim) / n;
    upper_p = sum(frame(:) > upper_g_lim) / n;
    figure; 
    subplot(1,2,1); imgray(frame);
    title(['L = ' num2str(lower_p,3) ', U = ' num2str(upper_p,3) ', T = ' num2str(lower_p+upper_p,3)]);
    %subplot(1,2,2); hist(frame_n(:), linspace(0, 1, 20));
    subplot(1,2,2); hist(frame(:), 0:8:256);
end
%%
frame_u = imread('C:\isbe\nailfold\good_frames\09frame07927.bmp');
frame = double(frame_u);
counts = hist(frame(:), 0:8:256);
[counts bin_i] = sort(counts, 'descend');

