im_folder = 'C:\isbe\nailfold\camera_capture\testing_new_cameras\bar\2014_11_26\L4_17_05_42\mosaic\';
new_folder = 'C:\isbe\nailfold\camera_capture\testing_new_cameras\bar\2014_11_26\L4_17_05_42\mosaic2\';

frame_list = dir([im_folder '*.bmp']);

num_frames = length(frame_list);

frames = zeros(128,128,num_frames);
for i_im = 1:num_frames
    frames(:,:,i_im) = imread([im_folder frame_list(i_im).name]);
end
%%
create_folder(new_folder);
min_g = min(frames(:));
max_g = max(frames(:));
for i_im = 1:num_frames
    write_im_from_colormap(frames(:,:,i_im), [new_folder frame_list(i_im).name], gray(256), [min_g max_g]);
end
%%
%session_dir = 'C:\isbe\nailfold\camera_capture\testing_new_cameras\bar\';
%sequence_dir = {'2014_12_03\14_47_54','2014_12_03\14_50_15','2014_12_03\15_58_27', '2014_12_03\16_16_10', '2014_12_10\11_11_51\' };

session_dir = 'N:\Nailfold Capillaroscopy\camera_capture\testing_new_camera\an\2014_12_10\';
sequences = dir(session_dir);

num_seqs = length(sequences) - 2;
sequence_dir = cell(num_seqs,1);
for i_s = 1:num_seqs
    sequence_dir{i_s} = sequences(i_s+2).name;
end
%
num_seqs = length(sequence_dir);
sequence = cell(num_seqs,1);
for i_s = 1:num_seqs
    [sequence{i_s}] = read_sequence_from([session_dir sequence_dir{i_s} '\sequence_properties.txt']);
end
%%
load([session_dir 'sequences.mat'], 'sequence');
num_seqs = length(sequence);

motor_x = cell(num_seqs,1);
motor_y = cell(num_seqs,1);
motor_z = cell(num_seqs,1);
sharpness = cell(num_seqs,1);

for i_s = 1:num_seqs
    num_frames = sequence{i_s}.num_frames;
    motor_x{i_s} = zeros(num_frames,1);
    motor_y{i_s} = zeros(num_frames,1);
    motor_z{i_s} = zeros(num_frames,1);
    sharpness{i_s} = zeros(num_frames,1);

    for i_f = 1:num_frames
        motor_x{i_s}(i_f) = sequence{i_s}.frames(i_f).motor_x;
        motor_y{i_s}(i_f) = sequence{i_s}.frames(i_f).motor_y;
        motor_z{i_s}(i_f) = sequence{i_s}.frames(i_f).motor_z;
        sharpness{i_s}(i_f) = sequence{i_s}.frames(i_f).sharpness;
        
        %frame = imread([session_dir sequence_dir{i_s} '\' sequence{i_s}.frames(i_f).frame_name]);
        %brightness{i_s}(i_f) = mean(frame(:));
    end
end
save([session_dir 'sequences.mat'], 'sequence');
%%
for i_s = 1:num_seqs

    num_frames = sequence{i_s}.num_frames;
    
    xy_speed = [0; sqrt(diff(motor_x{i_s}).^2 + diff(motor_y{i_s}).^2)];
    z_speed = [0; diff(motor_z{i_s})];
    
    figure; 
    
    subplot(2,3,1); hold all; plot(motor_x{i_s}, sharpness{i_s}, 'r.');
    xlabel('X motor position - mm');
    ylabel('Frame sharpness');
    set(gca, 'ylim', [15 40]);
    
    subplot(2,3,2); hold all; plot(motor_y{i_s}, sharpness{i_s}, 'r.');
    xlabel('Y motor position - mm');
    ylabel('Frame sharpness');
    set(gca, 'ylim', [15 40]);
    
    subplot(2,3,3); hold all; plot(motor_z{i_s}, sharpness{i_s}, 'r.');
    xlabel('Z motor position - mm');
    ylabel('Frame sharpness');
    set(gca, 'ylim', [15 40]);
    
    subplot(2,3,4); hold all; plot(1:num_frames, sharpness{i_s}, 'r.');
    xlabel('Time');
    ylabel('Frame sharpness');
    set(gca, 'ylim', [15 40]);
    
    subplot(2,3,5); hold all; plot(1:num_frames, motor_z{i_s}, 'r.');
    xlabel('Time');
    ylabel('Z motor position - mm');
    
    subplot(2,3,6); hold all; 
    plot(1:num_frames, xy_speed, 'b.');
    plot(1:num_frames, z_speed, 'r.');
    xlabel('Time');
    ylabel('Motor movements');
    legend({'(x,y) movement', 'z movement'});
    
%     brightness_sm = medfilt1(brightness{i_s}, 12);
%     figure; 
%     subplot(2,1,1); 
%     plot(1:24, brightness(1:24));
%     plot(1:24, brightness_sm(1:24), 'r'); set(gca, 'ylim', [118 132]);
%     subplot(2,1,2); hold on;
%     plot(1:num_frames, brightness);
%     plot(1:num_frames, brightness_sm, 'r'); set(gca, 'ylim', [118 132]);
    
end
%%
diff_rng = 3;
sm_win = 15;
half_win = (sm_win-1)/2;
for i_s = 1:num_seqs
    xy_speed = [0; sqrt(diff(motor_x{i_s}).^2 + diff(motor_y{i_s}).^2)];
    xy0 = bwareaopen(~xy_speed, 12);

    z_speed = [0; diff(motor_z{i_s})];
    z0 = bwareaopen(~z_speed, 12);
    
    %xy_speed = conv(xy_speed, ones(sm_win,1)/sm_win, 'valid');
    %az_speed = conv(abs(z_speed), ones(sm_win,1)/sm_win, 'valid');
    %sharpness_sm = conv(sharpness{i_s}, ones(sm_win,1)/sm_win, 'valid');
    
    frame_nums = 1:length(xy_speed);
    
    idx1 = xy0 & z0;
    idx2 = ~xy0;
    idx3 = ~z0;
    figure;
    subplot(3,1,1);hold on; 
    plot(xy_speed, 'b');
    plot(abs(z_speed), 'r');
    xlabel('Time');
    ylabel('Motor speed');
    legend({'|XY| speed', '|Z| speed'});
    
    %plot(sharpness_change / (max(sharpness_change)-min(sharpness_change)), 'g--');
    subplot(3,1,2); hold on;
    plot(motor_z{i_s}(idx1), 'g.');
    plot(motor_z{i_s}(idx2), 'b.');
    plot(motor_z{i_s}(idx3), 'r.');
    xlabel('Time');
    ylabel('Z position');
    legend({'Stationary', '|XY| movement', 'Z movement'});
    
    subplot(3,1,3); hold on;
    plot(frame_nums(idx1), sharpness{i_s}(idx1), 'g.');
    plot(frame_nums(idx2), sharpness{i_s}(idx2), 'b.');
    plot(frame_nums(idx3), sharpness{i_s}(idx3), 'r.');
    xlabel('Time');
    ylabel('Sharpness');
    legend({'Stationary', 'XY movement', 'Z movement'});
    
    
    
end
%%
%%

i_s = 3;
xy_speed = [0; sqrt(diff(motor_x{i_s}).^2 + diff(motor_y{i_s}).^2)];
xy0 = bwareaopen(~xy_speed, 12);
xy_lab = bwlabel(xy0);

num_segments = max(xy_lab);

for i_seg = 1:num_segments
    
    seg_idx = xy_lab == i_seg;
    
    if any(diff(motor_z{i_s}(seg_idx)))
%         figure;
%         plot(motor_z{i_s}(seg_idx), sharpness{i_s}(seg_idx), 'r.');  
%         ylim = get(gca, 'ylim');
%         ylim(2) = ylim(1)+4;
%         set(gca, 'ylim', ylim);
        focus_expt2(d, seg_idx, motor_z{i_s});
    end

    
end
%%
for i_seg = 1:num_segments
    
    seg_idx = xy_lab == i_seg;
    
    if any(diff(motor_z{i_s}(seg_idx)))
        
        brightness = zeros(sum(seg_idx),1);
        
        di = d(seg_idx);
        for i_f = 1:length(di)
            frame = double(imread(fullfile([session_dir sequence_dir{i_s}],d(i_f).name)));
            brightness(i_f) = mean(frame(:));
        end
        figure;
        plot(motor_z{i_s}(seg_idx), brightness / max(brightness), 'r.');  
    end

    
end

%%
lims = 1339:1735;
motor_xi = motor_x{15}(lims);
motor_yi = motor_y{15}(lims);
motor_zi = motor_z{15}(lims);
sharpness_i = sharpness{15}(lims);

num_frames = lims(end)-lims(1)+1;

figure; 
    
subplot(2,3,1); hold all; plot(motor_xi, sharpness_i, 'r.');
xlabel('X motor position - mm');
ylabel('Frame sharpness');
%set(gca, 'ylim', [15 40]);

subplot(2,3,2); hold all; plot(motor_yi, sharpness_i, 'r.');
xlabel('Y motor position - mm');
ylabel('Frame sharpness');
% set(gca, 'ylim', [15 40]);

subplot(2,3,3); hold all; plot(motor_zi, sharpness_i, 'r.');
xlabel('Z motor position - mm');
ylabel('Frame sharpness');
% set(gca, 'ylim', [15 40]);

subplot(2,3,4); hold all; plot(1:num_frames, sharpness_i, 'r.');
xlabel('Time');
ylabel('Frame sharpness');
% set(gca, 'ylim', [15 40]);

subplot(2,3,5); hold all; plot(1:num_frames, motor_zi, 'r.');
xlabel('Time');
ylabel('Z position');
% set(gca, 'ylim', [15 40]);


%%
im_folder = 'C:\isbe\nailfold\camera_capture\testing_new_cameras\bar\2014_11_26\L4_17_05_42\';
reg_folder = [im_folder 'reg_g1d\'];

frame_list = dir([im_folder 'frame*.bmp']);

num_frames = length(frame_list);

frames = zeros(480,640,num_frames);
for i_im = 1:num_frames
    frames(:,:,i_im) = imread([im_folder frame_list(i_im).name]);
end
[compound_transforms] = ...
        register_tiles_features(frames, ...
                                'theta_range', [0], ...
                                'offset_lim', 40, ...
                                'debug', false);

[nailfold_mosaic, mosaic_weights, compound_transforms] = ...
    create_mosaic(frames, compound_transforms);
%%
[diff_img, ~, g_lims] = write_trans_tiles(frames, compound_transforms, ...
                             reg_folder, 'frame_', [0 255], nailfold_mosaic);
                         
figure; imgray(diff_img);
%%
corrected_folder = [im_folder 'reg_g1d\corrected\'];
create_folder(corrected_folder);
for f = 1:num_frames
    frame = double(frames(:,:,f)) - diff_img;
    filename = sprintf('frame_%04d.png', f);
    imwrite(uint8(frame), fullfile(corrected_folder, filename));
end
%%
corrected_frames = bsxfun(@minus, double(frames), diff_img);
[compound_transforms2] = ...
        register_tiles_features( corrected_frames, ...
                                'theta_range', [0], ...
                                'offset_lim', 40, ...
                                'debug', false);
[nailfold_mosaic2, ~, compound_transforms2] = ...
    create_mosaic(corrected_frames, compound_transforms2);
%%
reg_corrected_folder = [im_folder 'reg_g1d\reg_corrected\'];
diff_img2 = write_trans_tiles(corrected_frames, ct, ...
                             reg_corrected_folder, 'frame_', g_lims, nailfold_mosaic2, 'match_grey_levels', 1, 'diff_image', []);
                         
figure; imgray(diff_img2);
%%
[compound_transforms2] = register_tiles_features( corrected_frames(:,:,1:10:end), ...
    'theta_range', [0], ...
    'offset_lim', 40, ...
    'debug', false);

[nailfold_mosaic2, ~, ct] = ...
    create_mosaic(corrected_frames, ct);
                            
%%
register_tiles_features( corrected_frames(:,:,[1 end]), ...
    'theta_range', [0], ...
    'offset_lim', 150, ...
    'debug', true);
%%
fid = fopen(fullfile(reg_folder,'make_movie.bat'),'w');
    fprintf(fid,'@echo off\r\n');
    
    fprintf(fid,'mkdir masked\r\n');
    fprintf(fid,'del masked\\*.png\r\n');
    fprintf(fid,'copy frame*.png masked\r\n');
    
    fprintf(fid,'mogrify -clip-mask seq_mask.png -threshold 101%% masked\\frame*.png\r\n');
    fprintf(fid,'ffmpeg -b 1200k -i masked\\frame_%%04d.png -an masked\\movie.mpg\r\n');
    
    fprintf(fid,'mkdir corrected_masked\r\n');
    fprintf(fid,'del corrected_masked\\*.png\r\n');
    fprintf(fid,'copy reg_corrected\\frame*.png corrected_masked\r\n');
        
    fprintf(fid,'mogrify -clip-mask seq_mask.png -threshold 101%% corrected_masked\\frame*.png\r\n');
    fprintf(fid,'ffmpeg -b 1200k -i corrected_masked\\frame_%%04d.png -an corrected_masked\\movie.mpg\r\n');
fclose(fid);
%%