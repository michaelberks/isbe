function [] = register_dual_camera_frames(frames_dir, camera1_transforms, camera2_transforms, ranges, range_type, varargin)
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
    'save_images', 1,...
    'save_dir', [], ...
    'save_name', [], ...
    'theta_range', -15:3:15,...
    'offset_range', 100,...
    'display_output', 1);
clear varargin;

if args.save_images && isempty(args.save_dir)
    args.save_dir = [frames_dir 'difference_images\'];
end
if args.save_images
    create_folder(args.save_dir);
end

%%-------------------------------------------------------------------------
% Get list of frame names
camera1_frame_names = dir_to_file_list([frames_dir '*_C_1_*.png']);
camera2_frame_names = dir_to_file_list([frames_dir '*_C_2_*.png']);

num_frames1 = length(camera1_frame_names);
num_frames2 = length(camera2_frame_names);

if num_frames1 ~= num_frames2
    error('Number of frames do not match');
end
num_frames = num_frames1; clear num_frames1 num_frames2;

if strcmp(range_type, 'time')
    %Get list of their times relative to the first frame of camera 1
    [dummy, frame_name] = fileparts(camera1_frame_names{1});
    clear dummy;
    
    mm = str2num(frame_name(28:29));
    ss = str2num(frame_name(31:32));
    ms = str2num(frame_name(34:36));
    t1_base = 60*mm + ss + ms/1000;

    times1 = zeros(num_frames, 1);
    times2 = zeros(num_frames, 1);

    for i_frame = 1:num_frames

        [dummy, frame_name] = fileparts(camera1_frame_names{i_frame});
        clear dummy;
        
        mm = str2num(frame_name(28:29));
        ss = str2num(frame_name(31:32));
        ms = str2num(frame_name(34:36));

        times1(i_frame) = (60*mm + ss + ms/1000) - t1_base;

        [dummy, frame_name] = fileparts(camera2_frame_names{i_frame});
        clear dummy;
        
        mm = str2num(frame_name(28:29));
        ss = str2num(frame_name(31:32));
        ms = str2num(frame_name(34:36));

        times2(i_frame) = (60*mm + ss + ms/1000) - t1_base;
    end
end

num_ranges = size(ranges, 1);
for i_rng = 1:num_ranges
    if strcmp(range_type, 'time')
        %Get frames from each camera in this time ranges
        include_frames1 = find(...
            (times1 >= ranges(i_rng,1)) & (times1 < ranges(i_rng,2)) );
        include_frames2 = find(...
            (times2 >= ranges(i_rng,1)) & (times2 < ranges(i_rng,2)) );

    else
        %Otherwise the ranges should just be the ith row of ranges
        include_frames1 = ranges(i_rng,1):ranges(i_rng,2);
        include_frames2 = ranges(i_rng,1):ranges(i_rng,2);
    end
    
    %Create mosaic for each camera
    [mosaic1, weights1] = create_mosaic(...
        camera1_frame_names(include_frames1), camera1_transforms(:,:,include_frames1), 'uniform');    
    tile_mask1 = weights1 == length(include_frames1);
    
    [mosaic2, weights2] = create_mosaic(...
        camera2_frame_names(include_frames2), camera2_transforms(:,:,include_frames2), 'uniform');
    tile_mask2 = weights2 == length(include_frames2);
    
    %Get the relative sizes of the two mosaics, increase the size of the
    %smaller
    [nrows1 ncols1] = size(mosaic1);
    [nrows2 ncols2] = size(mosaic2);
    
    if nrows1 < nrows2
        mosaic1(nrows2, :) = 0;
        tile_mask1(nrows2, :) = 0;
    elseif nrows1 > nrows2
        mosaic2(nrows1, :) = 0;
        tile_mask2(nrows1, :) = 0;
    end
    
    if ncols1 < ncols2
        mosaic1(:, ncols2) = 0;
        tile_mask1(:, ncols2) = 0;
    elseif ncols1 > ncols2
        mosaic2(:, ncols1) = 0;
        tile_mask2(:, ncols1) = 0;
    end
    tile_sz = size(mosaic1);
    
    %Register the two mosaics
    mosaic12 = cat(3, mosaic1, mosaic2);
    tile_mask12 = cat(3, tile_mask1, tile_mask2);
    [frame_transforms] = register_tiles_features(mosaic12, ...
        'tile_masks', tile_mask12,...
        'theta_range', args.theta_range, ...
        'offset_lim', args.offset_range, ...
        'debug', false);
    
    % Get size of the combined reference frame
    [mosaic_sz, frame_transforms] = mosaic_limits(tile_sz, frame_transforms);
    
    registered_mosaics = zeros(mosaic_sz(1), mosaic_sz(2), 2);
    for i_cam = 1:2
        
        [tile_out, mask] = sample_tile_image(...
            {mosaic12(:,:,i_cam)}, ones(tile_sz), ...
            inv(frame_transforms(:,:,i_cam)), ...
            mosaic_sz, []);
        
        mask = full(mask > 0.5);
        tile_out = full(tile_out{1});
        tile_out(~mask) = NaN;
        
        registered_mosaics(:,:,i_cam) = tile_out;

    end
    
    registered_difference = registered_mosaics(:,:,1) - registered_mosaics(:,:,2);
    
    if args.save_images
        save([args.save_dir args.save_name 'difference_image' zerostr(i_rng, 3) '.mat'],...
            'registered_difference', 'registered_mosaics', 'include_frames1', 'include_frames2');
    end
    
    if args.display_output
    
        figure;
        subplot(1,3,1); imgray(registered_mosaics(:,:,1));
        title('Registered compound frame from camera 1');
        subplot(1,3,2); imgray(registered_mosaics(:,:,2));
        title('Registered compound frame from camera 2');
        subplot(1,3,3); imgray(registered_difference);
        title('Difference between compound frames (1 - 2)');
    end
end
                            
                            

