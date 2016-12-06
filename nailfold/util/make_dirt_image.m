function [dirt_image] = make_dirt_image(frames_dir, image_format)
%MAKE_DIRT_IMAGE *Insert a one line summary here*
%   [dirt_image] = make_dirt_image(frame_dir)
%
% Inputs:
%      frames_dir - Folder of frames used to make dirt image
%
%      image_format (optional): Default 'bmp'
%
%
% Outputs:
%      dirt_image
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 05-Dec-2016
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

if ~exist('image_format', 'var') || isempty(image_format)
    image_format = 'bmp';
end

%Get list of frames
frames_list = dir([frames_dir '*.' image_format]);

num_frames = length(frames_list);
for i_frame = 1:num_frames
    frame = double(imread([frames_dir frames_list(i_frame).name]));
    
    if i_frame == 1
        dirt_image = frame;
    else
        dirt_image = dirt_image + frame;
    end
end

dirt_image = dirt_image / num_frames;
dirt_image = dirt_image - mean(dirt_image(:));