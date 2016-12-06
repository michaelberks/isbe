function [] = prepare_dual_camera_frames(frames_dir, varargin)
%REGISTER_DUAL_CAMERA_FRAMES *Insert a one line summary here*
%   [] = register_dual_camera_frames(varargin)
%
% REGISTER_DUAL_CAMERA_FRAMES uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Optional Arguments:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 14-Aug-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, '0', ...
    'camera1_ext', '_C_1_',...
    'camera2_ext', '_C_2_',...
    'image_format', 'bmp',...
    'rotated_dir', [], ...
    'rotation1', 1,...
    'rotation2', 1,...
    'flip1', 0,...
    'flip2', 1,...
    'transforms_dir', [], ...
    'first_corrected_dir', [], ...
    'final_corrected_dir', [], ...
    'reg_dir', [], ...
    'display_output', 1); %Should probably add options for the registration theta/offset limits
clear varargin;

%%-------------------------------------------------------------------------
% Create folder to store the transformed images
if isempty(args.rotated_dir)
    rotated_dir = [frames_dir 'rotated\'];
else
    rotated_dir = args.rotated_dir;
end

if isempty(args.transforms_dir)
    transforms_dir = [frames_dir 'transforms\'];
else
    transforms_dir = args.transforms_dir;
end

if isempty(args.first_corrected_dir)
    first_corrected_dir = [frames_dir 'first_corrections\'];
else
    first_corrected_dir = args.first_corrected_dir;
end

if isempty(args.rotated_dir)
    final_corrected_dir = [frames_dir 'final_corrections\'];
else
    final_corrected_dir = args.final_corrected_dir;
end

if isempty(args.reg_dir)
    reg_dir = [frames_dir 'registered_frames\'];
else
    reg_dir = args.reg_dir;
end

create_folder(rotated_dir);
create_folder(transforms_dir);
create_folder(first_corrected_dir);
create_folder(final_corrected_dir);
create_folder(reg_dir);

%%-------------------------------------------------------------------------
% Get list of frame names
camera1_frame_names = dir_to_file_list([frames_dir '*' args.camera1_ext '*.' args.image_format]);
camera2_frame_names = dir_to_file_list([frames_dir '*' args.camera2_ext '*.' args.image_format]);

num_frames1 = length(camera1_frame_names);
num_frames2 = length(camera2_frame_names);

if num_frames1 > num_frames2
    excess = num_frames1 - num_frames2;
    
    display('Number of frames do not match');
    display(['Camera 1 has ' num2str(num_frames1) ' frames, camera 2 has ' num2str(num_frames2) ' frames']);
    display(['The last ' num2str(excess) ' frames from camera 1 will not be used']);
    
    num_frames1 = num_frames2;
    camera1_frame_names(num_frames1+1:end) = [];
    
elseif num_frames2 > num_frames1
    excess = num_frames2 - num_frames1;
    
    display('Number of frames do not match');
    display(['Camera 1 has ' num2str(num_frames1) ' frames, camera 2 has ' num2str(num_frames2) ' frames']);
    display(['The last ' num2str(excess) ' frames from camera 2 will not be used']);
    
    num_frames2 = num_frames1;
    camera2_frame_names(num_frames2+1:end) = [];
else
    display(['Processing ' num2str(num_frames1) ' from cameras 1 and 2']);
end
    
num_frames = num_frames1; clear num_frames1 num_frames2;

%%
%--------------------------------------------------------------------------
% Flip all the frames from camera 2, and rotate all the frames from both
% cameras
for i_frame = 1:num_frames
    
    frame1 = imread([frames_dir camera1_frame_names{i_frame}]);    
    frame2 = imread([frames_dir camera2_frame_names{i_frame}]);
    
    if args.flip1
        frame1 = fliplr(frame1);
    end
    if args.flip2
        frame2 = fliplr(frame2);
    end
    if args.rotation1
        frame1 = rot90(frame1, args.rotation1);  
    end
    if args.rotation2
        frame2 = rot90(frame2, args.rotation2);
    end
    
    imwrite(frame1, [rotated_dir camera1_frame_names{i_frame}]);    
    imwrite(frame2, [rotated_dir camera2_frame_names{i_frame}]);   
end
%%
%--------------------------------------------------------------------------
%%
% Now prepare the frames from each camera
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

%
prepare_frames(camera2_frame_names, camera1_frame_names, ...
    rotated_dir, first_corrected_dir, final_corrected_dir, transforms_dir, ...
    'camera1', args.display_output)

prepare_frames(camera1_frame_names, camera2_frame_names, ...
    rotated_dir, first_corrected_dir, final_corrected_dir, transforms_dir, ...
    'camera2', args.display_output)

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function prepare_frames(camera1_frame_names, camera2_frame_names, ...
    rotated_dir, first_corrected_dir, final_corrected_dir, transforms_dir, ...
    camera_name, display_output)
%The preparation steps are
% 1) Register a selection of camera 2 frames, spread throughout the
% sequence to the first frame from camera 1
% 2) Use the resulting transforms to make a mosaic of the camera 2 frames,
% and use this to compute the difference image which estimates the camera 2
% dirt
% 3) Correct all the camera 2 frames, a re-register the entire sequence in
% the standard consecutive frames manner
% 4) Recompute the diff image, and update the corrected frames

%Load in first frame from camera 1
frame1 = imread([rotated_dir camera1_frame_names{1}]);

%Pre-allocate storage for frames from camera 2 and transforms 
[nrows ncols] = size(frame1);
max_frames = 20;
num_frames = length(camera1_frame_names);
frames2 = zeros(nrows, ncols, max_frames);
compound_transforms2_to_1 = zeros(3,3,max_frames);

%Loop through camera 2 frames, registering each one to frame1
spacing = floor(num_frames / max_frames);
for i_frame = 1:max_frames
    
    frame2 = imread([rotated_dir camera2_frame_names{spacing*i_frame}]);

    frame12 = cat(3,frame1,frame2);
    frames2(:,:,i_frame) = frame2;

    %At this stage the transforms may be quite large, so need large theta
    %and offset ranges
    [compound_transforms_i] = ...
        register_tiles_features(frame12, ...
                                'theta_range', -15:3:15, ...
                                'offset_lim', 100, ...
                                'debug', false);
    compound_transforms2_to_1(:,:,i_frame) = compound_transforms_i(:,:,2);
end

%Create mosaic from the registered camera 2 frames
[nailfold_mosaic2to1, dummy, mosaic_transforms2to1] = ...
    create_mosaic(frames2, compound_transforms2_to_1);
clear dummy;

%Compute the mean diff image of each frame from the mosaic
[first_dirt_img] = write_trans_tiles(frames2, mosaic_transforms2to1, NaN, ...
    [], [], nailfold_mosaic2to1);
save([transforms_dir camera_name '_first_correction_img.mat'], 'first_dirt_img');

if display_output
    figure; imgray(first_dirt_img);
    title('First estimate of correction image (e.g. dirt on sensor)');
end

%Now correct all the camera 2 frames with this image and save the output
for i_frame = 1:num_frames
    frame2 = imread([rotated_dir camera2_frame_names{i_frame}]);    
    corrected_frame = uint8(double(frame2) - first_dirt_img);
    imwrite(corrected_frame, [first_corrected_dir camera2_frame_names{i_frame}]);
    
    if ~rem(i_frame, 60) && display_output
        figure; 
        subplot(1,2,1); imgray(frame2); title('Original frame'); 
        subplot(1,2,2); imgray(corrected_frame); title('Frame after first correction'); 
    end
end

%Now we can register all the corrected frames sequentially and create a
%mosaic
[compound_transforms] = ...
        register_tiles_features(strcat(first_corrected_dir, camera2_frame_names), ...
                                'theta_range', 0, ...
                                'offset_lim', 40, ...
                                'debug', false);
[nailfold_mosaic, dummy, mosaic_transforms] = ...
    create_mosaic(strcat(first_corrected_dir, camera2_frame_names), compound_transforms);
clear dummy;

if display_output
    figure; imgray(nailfold_mosaic);
    title('Compound mosaic of all consecutively registered frames');
end

%save these transforms
save([transforms_dir camera_name '_final_reg_transforms.mat'], 'mosaic_transforms');
    
%Create a new diff image
[final_dirt_img] = write_trans_tiles(strcat(first_corrected_dir, camera2_frame_names),...
    mosaic_transforms, NaN, [], [], nailfold_mosaic);
save([transforms_dir camera_name '_final_correction_img.mat'], 'final_dirt_img');

if display_output
    figure; imgray(first_dirt_img);
    title('Final estimate of correction image (e.g. dirt on sensor)');
end

%Write out the new corrected frames
for i_frame = 1:num_frames
    frame2 = imread([first_corrected_dir camera2_frame_names{i_frame}]);    
    corrected_frame = uint8(double(frame2) - final_dirt_img);
    imwrite(corrected_frame, [final_corrected_dir camera2_frame_names{i_frame}]);
    
    if ~rem(i_frame, 60) && display_output
        figure; 
        subplot(1,2,1); imgray(frame2); title('Frame after first correction'); 
        subplot(1,2,2); imgray(corrected_frame); title('Frame after final correction'); 
    end
end