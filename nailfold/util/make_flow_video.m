function [] = make_flow_video(frames, video_path, frame_rate)
%MAKE_FLOW_VIDEO *Insert a one line summary here*
%   [] = make_flow_video(frames, video_path)
%
% Inputs:
%      frames - *Insert description of input variable here*
%
%      video_path - *Insert description of input variable here*
%
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 12-Feb-2016
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
num_frames = size(frames,3);
g_lims = prctile(double(frames(:)), [1 99]);
g_range = g_lims(2) - g_lims(1);
frames = 255*(double(frames)-g_lims(1))/g_range;

rr = rem(size(frames, 1),4);
rc = rem(size(frames, 2),4);
frames(1:rr,:,:) = [];
frames(:,1:rc,:) = [];

temp_frames_dir = tempname;
mkdir(temp_frames_dir);

for i_fr = 1:num_frames;
    imwrite(uint8(frames(:,:,i_fr)),...
        [temp_frames_dir '\frame' zerostr(i_fr,4) '.bmp']);
end
cmd = ['ffmpeg -y -r ' num2str(frame_rate) ' -i "' temp_frames_dir '\frame%04d.bmp" -c:v libx264 -preset slow -crf 18 -an "' video_path '"'];
system(cmd);

delete([temp_frames_dir '\*']);
rmdir(temp_frames_dir);
