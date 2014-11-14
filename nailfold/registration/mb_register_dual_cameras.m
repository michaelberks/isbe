%Script experimenting with registering frames captured on the dual camera
%system


%%-------------------------------------------------------------------------
% 1) Get list of frame names
frames_dir = 'C:\isbe\nailfold\data\dual_wavelength\2013_07_12\';

orig_dir = [frames_dir 'original\'];
rotated_dir = [frames_dir 'rotated\']; create_folder(rotated_dir);
reg_dir = [frames_dir 'registered_frames\']; create_folder(reg_dir);
first_corrected_dir = [frames_dir 'first_corrections\']; create_folder(first_corrected_dir);
final_corrected_dir = [frames_dir 'final_corrections\']; create_folder(final_corrected_dir);

camera1_frame_names = dir_to_file_list([orig_dir '*_C_1_*.png'], []);
camera2_frame_names = dir_to_file_list([orig_dir '*_C_2_*.png'], []);

num_frames1 = length(camera1_frame_names);
num_frames2 = length(camera2_frame_names);

if num_frames1 ~= num_frames2
    error('Number of frames do not match');
end
num_frames = num_frames1; clear num_frames1 num_frames2;
%%
%--------------------------------------------------------------------------
% Flip all the frames from camera 2, and rotate all the frames from both
% cameras
for i_frame = 1:num_frames
    
    frame1 = imread([orig_dir camera1_frame_names{i_frame}]);    
    frame2 = imread([orig_dir camera2_frame_names{i_frame}]);
    
    frame1 = rot90(frame1, 1);    
    frame2 = fliplr(rot90(frame2, 1));
    
    imwrite(frame1, [rotated_dir camera1_frame_names{i_frame}]);    
    imwrite(frame2, [rotated_dir camera2_frame_names{i_frame}]);   
end
%%
%--------------------------------------------------------------------------
%Lets look at some frames
for i_frame = 500+(1:20)
    frame1 = imread([rotated_dir camera1_frame_names{i_frame}]);
    frame2 = imread([rotated_dir camera2_frame_names{i_frame}]);
    
    figure;
    subplot(1,2,1); imgray(frame1);
    subplot(1,2,2); imgray(frame2);
    
end
%%
%--------------------------------------------------------------------------
% 3) Get list of their times relative to the first frame of camera 1
mm = str2num(camera1_frame_names{1}(28:29));
ss = str2num(camera1_frame_names{1}(31:32));
ms = str2num(camera1_frame_names{1}(34:36));
t1_base = 60*mm + ss + ms/1000;

times1 = zeros(num_frames, 1);
times2 = zeros(num_frames, 1);

for i_frame = 1:num_frames
    
    mm = str2num(camera1_frame_names{i_frame}(28:29));
    ss = str2num(camera1_frame_names{i_frame}(31:32));
    ms = str2num(camera1_frame_names{i_frame}(34:36));
    
    times1(i_frame) = (60*mm + ss + ms/1000) - t1_base;
    
    mm = str2num(camera2_frame_names{i_frame}(28:29));
    ss = str2num(camera2_frame_names{i_frame}(31:32));
    ms = str2num(camera2_frame_names{i_frame}(34:36));
    
    times2(i_frame) = (60*mm + ss + ms/1000) - t1_base;
end

for ii = 0:100:500
    figure;
    plot([times1(ii+(1:100))'; times2(ii+(1:100))'], [0 1]);
    title( 'Time recorded for each camera 2 frame')
    xlabel( 'Time recorded for each camera 1 frame')
end
%--------------------------------------------------------------------------
%%
% Instead, we try this approach:
% 1) Register a selection of camera 2 frames, spread throughout the
% sequence to the first frame from camera 1
% 2) Use the resulting transforms to make a mosaic of the camera 2 frames,
% and use this to compute the difference image which estimates the camera 2
% dirt
% 3) Correct all the camera 2 frames, a re-register the entire sequence in
% the standard consecutive frames manner
% 4) Recompute the diff image, and update the corrected frames
%
% 5) Repeat 1-4 but swapping the two cameras
% 6) Pair of individual frames or time interval mosaics from camera 1 and
% camera 2 can then be registered, and a difference image, (hopefully
% showing signs of oxidisation!) computed. Note that at this stage we need
% to flip one of the cameras

%Load in first frame from camera 1
frame1 = imread([rotated_dir camera1_frame_names{1}]);

%Pre-allocate storage for frames from camera 2 and transforms 
[nrows ncols] = size(frame1);
max_frames = 20;
frames2 = zeros(nrows, ncols, max_frames);
compound_transforms2_to_1 = zeros(3,3,max_frames);

%Loop through camera 2 frames, registering each one to frame1
spacing = floor(num_frames / max_frames);
for i_frame = 1:max_frames
    
    frame2 = imread([rotated_dir camera2_frame_names{spacing*i_frame}]);

    frame12 = cat(3,frame1,frame2);
    frames2(:,:,i_frame) = frame2;

    [compound_transforms_i] = ...
        register_tiles_features(frame12, ...
                                'theta_range', -15:3:15, ...
                                'offset_lim', 40, ...
                                'debug', false);
    compound_transforms2_to_1(:,:,i_frame) = compound_transforms_i(:,:,2);
end

%Create mosaic from the registered camera 2 frames
[nailfold_mosaic2, dummy, mosaic_transforms2] = ...
    create_mosaic(frames2, compound_transforms2_to_1);
clear dummy;

%Compute the mean diff image of each frame from the mosaic
[diff_img2] = write_trans_tiles(frames2, mosaic_transforms2, NaN, ...
    [], [], nailfold_mosaic2);
figure; imgray(diff_img2);
%
%Now correct all the camera 2 frames with this image and save the output
for i_frame = 1:num_frames
    frame2 = imread([rotated_dir camera2_frame_names{i_frame}]);    
    corrected_frame = uint8(double(frame2) - diff_img2);
    imwrite(corrected_frame, [first_corrected_dir camera2_frame_names{i_frame}]);
    if ~rem(i_frame, 30)
        figure; 
        subplot(1,2,1); imgray(frame2);
        subplot(1,2,2); imgray(corrected_frame);
    end
end

%Now we can register all the corrected frames sequentially and create a
%mosaic
[compound_transforms2] = ...
        register_tiles_features(strcat(first_corrected_dir, camera2_frame_names(1:60)), ...
                                'theta_range', [], ...
                                'offset_lim', 40, ...
                                'debug', false);
[nailfold_mosaic2, dummy, mosaic_transforms2] = ...
    create_mosaic(strcat(first_corrected_dir, camera2_frame_names(1:60)), compound_transforms2);
clear dummy;

%save these transforms
save ...
    
%Create a new diff image
[diff_img2] = write_trans_tiles(strcat(first_corrected_dir, camera2_frame_names(1:60)),...
    mosaic_transforms2, NaN, [], [], nailfold_mosaic2);

%Write out the new corrected frames
for i_frame = 1:60%num_frames
    frame2 = imread([first_corrected_dir camera2_frame_names{i_frame}]);    
    corrected_frame = uint8(double(frame2) - diff_img2);
    imwrite(corrected_frame, [final_corrected_dir camera2_frame_names{i_frame}]);
    if ~rem(i_frame, 30)
        figure; 
        subplot(1,2,1); imgray(frame2);
        subplot(1,2,2); imgray(corrected_frame);
    end
end
%%
figure; imgray(diff_img1);

for i_frame = 1:600
    frame1 = imread([frames_dir camera1_frame_names{i_frame}]);
    frame1 = rot90(frame1, 1);
    
    corrected_frame = double(frame1) - diff_img1;
    save([corrected_dir 'frame1' zerostr(i_frame,4) '.mat'], 'corrected_frame');
    if ~rem(i_frame, 20)
        figure; 
        subplot(1,2,1); imgray(frame1);
        subplot(1,2,2); imgray(corrected_frame);
    end
end
                         
                         

figure; imgray(diff_img2);

%%-------------------------------------------------------------------------
%Match them up by their time stamps

%Co register them

%Optionally register pairs of frames in time in time
%%

%%
for i_frame = 1:600
    frame2 = imread([frames_dir camera2_frame_names{i_frame}]);
    frame2 = fliplr(rot90(frame2, 1));
    
    corrected_frame = double(frame2) - diff_img2;
    save([corrected_dir 'frame2' zerostr(i_frame,4) '.mat'], 'corrected_frame');
    if ~rem(i_frame, 20)
        figure; 
        subplot(1,2,1); imgray(frame2);
        subplot(1,2,2); imgray(corrected_frame);
    end
end
%%

%Now repeat, registering frame 1 to frame 2
frame2 = frames2(:,:,1) - diff_img2;

frames1 = zeros(nrows, ncols, max_frames);
compound_transforms1_to_2 = zeros(3,3,max_frames);

for i_frame = 1:max_frames
    
    frame1 = imread([frames_dir camera1_frame_names(30*i_frame).name]);
    frame1 = rot90(frame1, 1);
    frame21 = cat(3,frame2,frame1);
    frames1(:,:,i_frame) = frame1;

    [compound_transforms_i] = ...
        register_tiles_features(frame21, ...
                                'theta_range', -15:3:15, ...
                                'offset_lim', 40, ...
                                'debug', false);
    compound_transforms1_to_2(:,:,i_frame) = compound_transforms_i(:,:,2);
end

[nailfold_mosaic1, mosaic_weights1, mosaic_transforms1] = ...
    create_mosaic(frames1, compound_transforms1_to_2);

[diff_img1] = write_trans_tiles(frames1, mosaic_transforms1, ...
                             reg_dir, 'frame_', [0 255], nailfold_mosaic1);

%%
[compound_transforms1] = ...
        register_tiles_features(bsxfun(@minus, frames1, diff_img1), ...
                                'theta_range', -15:3:15, ...
                                'offset_lim', 40, ...
                                'debug', true);
                            
[nailfold_mosaic1, mosaic_weights1, mosaic_transforms1] = ...
    create_mosaic(bsxfun(@minus, frames1, diff_img1), compound_transforms1);

[diff_img1a] = write_trans_tiles(bsxfun(@minus, frames1, diff_img1), mosaic_transforms1, ...
                             reg_dir, 'frame_', [0 255], nailfold_mosaic1);
                         
[nailfold_mosaic1a] = ...
    create_mosaic(bsxfun(@minus, frames1, diff_img1+diff_img1a), compound_transforms1);                         
                         
%%
for i_frame = 1:600
    frame1 = imread([frames_dir camera1_frame_names{i_frame}]);
    frame1 = rot90(frame1, 1);
    
    corrected_frame = double(frame1) - (diff_img1 + diff_img1a);
    save([corrected_dir 'frame1a' zerostr(i_frame,4) '.mat'], 'corrected_frame');
    if ~rem(i_frame, 20)
        figure; 
        subplot(1,2,1); imgray(frame1);
        subplot(1,2,2); imgray(corrected_frame);
    end
end
%%

c1_transforms = u_load('C:\isbe\nailfold\data\dual_wavelength\2013_07_12\original\transforms\camera1_final_reg_transforms.mat');
c2_transforms = u_load('C:\isbe\nailfold\data\dual_wavelength\2013_07_12\original\transforms\camera2_final_reg_transforms.mat');
frames_dir = 'C:\isbe\nailfold\data\dual_wavelength\2013_07_12\original\final_corrections\';

register_dual_camera_frames(frames_dir, c1_transforms, c2_transforms, [1 60], 'frames');





























%%
%%
%--------------------------------------------------------------------------
% THIS DOESN'T WORK!!!!!
% Register the tiles together (i.e. compute the transformation from one to
% the next over the sequence) ...
% ... because there is too much dirt in the images, that is registered
% instead of the vessels!

[compound_transforms_orig] = ...
    register_tiles_features(frames1, ...
                            'theta_range', [], ...
                            'offset_lim', 40, ...
                            'debug', false);

%                        
[nailfold_mosaic1, mosaic_weights, compound_transforms] = ...
    create_mosaic(frames1, compound_transforms_orig);

[diff_img] = write_trans_tiles(frames1, compound_transforms, ...
                             reg_dir, 'frame_', [0 255], nailfold_mosaic1);
%
frames1_corrected = bsxfun(@minus, frames1, diff_img);
[nailfold_mosaic1_corrected] = ...
    create_mosaic(frames1_corrected, compound_transforms_orig);                         

[diff_img_corrected seq_mask] = write_trans_tiles(frames1_corrected, compound_transforms, ...
                             corrected_dir, 'frame_', [0 255], nailfold_mosaic1_corrected);                                                  
figure; imgray(diff_img);
figure; imgray(diff_img_corrected);
%--------------------------------------------------------------------------