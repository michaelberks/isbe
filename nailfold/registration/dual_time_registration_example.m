%Example script for preparing and the registering frames captured from the
%dual camera system

%Set folder containing all the frames
frames_dir = 'C:\isbe\nailfold\data\dual_wavelength\550 nm t_compare_t0_t2\';

%Make dirt image - I don't have frames for this, so just reusing the above.
%You want to use the (out-of-focus) frames you captured especially for this
dirty_dir = 'C:\isbe\nailfold\data\dual_wavelength\550 nm t_compare_t0_t2\'; %
[dirt_image] = make_dirt_image(dirty_dir);
figure; imgray(dirt_image); colorbar;

%Prepare the frames - I've included all the optional arguments below with
%their default values, you shouldn't need to change these, but you might
%want to:
% - Increase offset_lim if it is too small to capture the motion between
% frames (the value is the number of pixels in each direction successive
% frames can move and still be matched)
% - Add a range of theta's if you think you are seeing rotational movement
% between frames
% - Play with different values of sigma - this will change the feature
% points used to match vessels between frames. Higher values of sigma pick
% large scale smoother vessels edges. A value about half the width (in
% pixels) of the vessel is probably about right
% - Change what output you display 1 = summary figures, 2 = full debug
% matching during registration. I suggest running with 2 to start with, and
% watching the registration in real-time to check frames are being matched
% appropriately
prepare_sequential_camera_frames(frames_dir, ... 
    'time1_ext', '_t0_',...
    'time2_ext', '_t2_',...
    'image_format', 'bmp',...
    'dirt_image', dirt_image,...
    'theta_range', 0,...
    'offset_lim', 120,...
    'sigma', 8,...
    'corrections_dir', [], ...
    'transforms_dir', [], ...
    'make_videos', false,... %Switch to true if you have ffmpeg installed and included on your system path and want to make videos of vessel flow
    'display_output', 0);
%%

%Load in the transforms from the above registrations
t1_transforms = u_load([frames_dir 'transforms\time1_reg_transforms.mat']);
t2_transforms = u_load([frames_dir 'transforms\time2_reg_transforms.mat']);

%Register pairs of compund frames from each camera - define the range used
%in each compound frame either by frame numbers (which will be the same for
%both cameras) or by time interval in seconds over the sequence (in which
%case different frames may be used for each camera because the times may
%not match up)
register_dual_camera_frames([frames_dir 'corrections\'],...
    t1_transforms, t2_transforms, [1 60], 'frames',...
    'camera1_ext', '_t0_',...
    'camera2_ext', '_t2_',...
    'image_format', 'bmp',...
    'theta_range', 0,...
    'offset_range', 240,...
    'save_images', 1, ...,
    'save_dir', [frames_dir 'difference_images_f\']);
%%
