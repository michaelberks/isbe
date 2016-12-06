function [] = prepare_sequential_camera_frames(frames_dir, varargin)
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
    'time1_ext', '_t0_',...
    'time2_ext', '_t1_',...
    'image_format', 'bmp',...
    'dirt_image', [],...
    'theta_range', 0,...
    'offset_lim', 120,...
    'sigma', 8,...
    'transforms_dir', [], ...
    'corrected_dir', [], ...
    'make_videos', false,...
    'display_output', 1);
clear varargin;

%%-------------------------------------------------------------------------
% Create folder to store the transformed images
if isempty(args.transforms_dir)
    transforms_dir = [frames_dir 'transforms\'];
else
    transforms_dir = args.transforms_dir;
end

if isempty(args.corrected_dir)
    corrected_dir = [frames_dir 'corrections\'];
else
    corrected_dir = args.corrected_dir;
end

create_folder(transforms_dir);
create_folder(corrected_dir);

%%-------------------------------------------------------------------------
% Get list of frame names
time1_frame_names = dir([frames_dir '*' args.time1_ext '*.' args.image_format]);
time2_frame_names = dir([frames_dir '*' args.time2_ext '*.' args.image_format]);

time1_frame_names = {time1_frame_names(:).name}';
time2_frame_names = {time2_frame_names(:).name}';

num_frames1 = length(time1_frame_names);
num_frames2 = length(time2_frame_names);

if ~num_frames1
    display('No frames found for time period 1.');
    display(['Is the pattern ' frames_dir '*' args.time1_ext '*.' args.image_format ' correct?']);
    return;
    
elseif ~num_frames2
    display('No frames found for time period 2');
    display(['Is the pattern ' frames_dir '*' args.time2_ext '*.' args.image_format ' correct?']);
    return;  

elseif num_frames1 > num_frames2
    excess = num_frames1 - num_frames2;
    
    display('Number of frames do not match');
    display(['Time 1 has ' num2str(num_frames1) ' frames, time 2 has ' num2str(num_frames2) ' frames']);
    display(['The last ' num2str(excess) ' frames from time 1 will not be used']);
    
    num_frames1 = num_frames2;
    time1_frame_names(num_frames1+1:end) = [];
    
elseif num_frames2 > num_frames1
    excess = num_frames2 - num_frames1;
    
    display('Number of frames do not match');
    display(['Time 1 has ' num2str(num_frames1) ' frames, time 2 has ' num2str(num_frames2) ' frames']);
    display(['The last ' num2str(excess) ' frames from time 2 will not be used']);
    
    num_frames2 = num_frames1;
    time2_frame_names(num_frames2+1:end) = [];
else
    display(['Processing ' num2str(num_frames1) ' from time periods 1 and 2']);
end

%%
%--------------------------------------------------------------------------
% Don't need to worry about flipping or rotating frames now - they're from
% the same camera
%--------------------------------------------------------------------------
%%
% Now prepare the frames from each time period
% 1) Load in each frame, subtracting the dirt image if supplied
% 2) Register sequential frames
% 3) Make mosaic and write out the cleaned frames, saving the tranforms
reg_params.theta_range = args.theta_range;
reg_params.offset_lim = args.offset_lim;
reg_params.sigma = args.sigma;

prepare_frames(time1_frame_names, frames_dir, ...
    args.dirt_image, corrected_dir, transforms_dir, ...
    reg_params,...
    'time1', args.display_output, args.make_videos);

prepare_frames(time2_frame_names, frames_dir, ...
    args.dirt_image, corrected_dir, transforms_dir, ...
    reg_params,...
    'time2', args.display_output, args.make_videos);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function prepare_frames(frame_names, frames_dir, ...
    dirt_image, corrected_dir, transforms_dir, ...
    reg_params,...
    display_name, display_output, make_video)
%The preparation steps are
% 1) Load in each frame, subtracting the dirt image if supplied
% 2) Register sequential frames
% 3) Make mosaic and write out the cleaned frames, saving the tranforms

%Load in frame subtracting dirt image
num_frames = length(frame_names);

for i_frame = 1:num_frames
    
    frame = double(imread([frames_dir frame_names{i_frame}]));

    if ~isempty(dirt_image)
        frame = frame - dirt_image;
        imwrite(uint8(frame), [corrected_dir frame_names{i_frame}]);
    end
    
    if i_frame == 1
        [nrows, ncols] = size(frame);
        frames = zeros(nrows, ncols, num_frames);
    end
    frames(:,:,i_frame) = frame;
    
end

%At this stage the transforms may be quite large, so need large theta
%and offset ranges
[compound_transforms] = ...
    register_tiles_features(frames, ...
                            'dirt_image', dirt_image,...
                            'theta_range', reg_params.theta_range, ...
                            'offset_lim', reg_params.offset_lim, ...
                            'sigma', reg_params.sigma,...
                            'debug', display_output > 1);


[nailfold_mosaic, ~, mosaic_transforms] = ...
    create_mosaic(frames, compound_transforms); %#ok

if display_output
    figure; imgray(nailfold_mosaic);
    title(['Compound mosaic of all consecutively registered frames for' display_name]);
end

%save these transforms
save([transforms_dir display_name '_reg_transforms.mat'], 'mosaic_transforms');

if make_video
    reg_folder = [frames_dir 'registered\'];
    g_lims = [min(nailfold_mosaic(:)) max(nailfold_mosaic(:))];

    delete([reg_folder 'frame*']);

    write_trans_tiles(frames, mosaic_transforms, ...
                             reg_folder, 'frame', g_lims, nailfold_mosaic);

    cmd = ['ffmpeg -y -r 15 -i "' reg_folder 'frame%04d.png" -c:v libx264 -preset slow -crf 18 -an "' reg_folder '_' display_name '_movie.mp4"'];
    system(cmd);

end