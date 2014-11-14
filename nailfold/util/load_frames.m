function [frames] = load_frames(frames_dir, crop_size, max_frames, image_type)
%LOAD_FRAMES *Insert a one line summary here*
%    frames = load_frames()
%
% Inputs:
%
%
% Outputs:
%      frames - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 23-Aug-2011
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
if nargin < 2
    crop_size = 4;
end
if nargin < 3
    max_frames = 16;
end
if nargin < 4
    image_type = 'bmp';
end

%Check directory is / terminated
if ~frames_dir(end)==filesep
    frames_dir = [frames_dir filesep];
end

%Get list of frames from dir
frame_list = dir([frames_dir '*.' image_type]); 

%Workout number of frames
num_frames = min(length(frame_list), max_frames);

%Loop through list loading each frame
for ff = 1:num_frames
    frame =  double(imread([frames_dir frame_list(ff).name]));
    
    %crop the frame
    frame = frame(crop_size+1:end-crop_size,crop_size+1:end-crop_size);
    
    if ff == 1
        %Get size of frame, we're assuming all subsequent frames are the
        %same size
        [r c] = size(frame);
        
        %Workout number of rows following interleaving
        r2 = 2*r-1;
        
        %Pre-allocate space for all frames
        frames = zeros(r2, c, num_frames);
    end
    
    %Take the odd rows from the loaded frame
    frames(1:2:r2,:,ff) = frame;
    
    %Compute the even rows as the average of the neighbouring odd rows
    frames(2:2:r2-1,:,ff) = conv2(frame, [.5 .5]', 'valid');
end