function [frames_reg offsets thetas] = ...
    register_frames((frames, offset_lim, theta_lim, first_frame)
%REGISTER_FRAME *Insert a one line summary here*
%   [offset, theta] = register_frame(frame1, frame2, offset_lim, theta_lim)
%
% Inputs:
%      frames - *Insert description of input variable here*
%
%      offset_lim - *Insert description of input variable here*
%
%      theta_lim - *Insert description of input variable here*
%
%
% Outputs:
%      frames_reg - *Insert description of input variable here*
%
%      offsets - *Insert description of input variable here*
%
%      thetas - *Insert description of input variable here*
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

if (nargin==0 && nargout==0), test_script(); return; end

if ~exist('offset_lim','var'), offset_lim = 10; end
if ~exist('theta_lim','var'), theta_lim = 10; end
if ~exist('first_frame','var'), first_frame = 1; end

%Get dimensions of the input frames
[rows cols num_frames] = size(frames);

%--------------------------------------------------
%Compute line points in first frame as the target

%Apply Gaussian 2nd derivatives
[line_strength, line_orientation] = ...
    gaussian_clover_line(frames(:,:,first_frame), [4 8]);%gaussian_2nd_derivative_line2

%Apply non-maximal suppression to skeletonise the line strength map
line_nms = mb_non_maximal_supp(line_strength, line_orientation);

%Discard points from edges of map
line_nms([1:16 end-15:end], :) = 0;
line_nms(:, [1:16 end-15:end]) = 0;

%Apply hysterisis to select suitable lines from the NMS map
[line_mask] = hysterisis(line_nms);

%Extract (x,y) coordinates of the remaining lines
[y_tgt x_tgt] = find(line_mask);

%Make points relative to centre
y_tgt = y_tgt - rows/2;
x_tgt = x_tgt - cols/2;

%--------------------------------------------------
%Generate offset and theta ranges to be used when registering each frame
offset_sz = 2*offset_lim + 1;

theta_range = pi*(-theta_lim:theta_lim)/180;
theta_sz = length(theta_range);

%Use first frame ito initial the target frames and pre-allocate space for
%output offsets and thetas
frames_reg = frames(:,:,1);
offsets = zeros(num_frames - 1, 2);
thetas = zeros(num_frames - 1, 1);

%--------------------------------------------------
%Loop through each of the remaining frames, registering it to the first
%frame and adding to output frames
frame_idx = [1:first_frame-1 first_frame+1:num_frames];
for ff = 1:num_frames-1
    %Apply Gaussian 2nd derivatives
    [line_strength, line_orientation] = ...
        gaussian_clover_line(frames(:,:,frame_idx(ff)), [4 8]);%gaussian_2nd_derivative_line2

    %Apply non-maximal suppression to skeletonise the line strength map
    line_nms = mb_non_maximal_supp(line_strength, line_orientation);

    %Discard points from edges of map
    line_nms([1:16 end-15:end], :) = 0;
    line_nms(:, [1:16 end-15:end]) = 0;

    %Apply hysterisis to select suitable lines from the NMS map
    [line_mask] = hysterisis(line_nms);

    %Extract (x,y) coordinates of the remaining lines
    [y_src x_src] = find(line_mask);
    
    %Make points relative to centre
    y_src = y_src - rows/2;
    x_src = x_src - cols/2;
    
    %--------------------------------------------------
    %Loop through thetas finding rotation that highest offset max

    max_count = 0;
    for tt = 1:theta_sz
        theta = theta_range(tt);

        %Generate rotation matrix for this theta
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

        %Rotate line points extracted from frame 1 - need to centre by (xc, yc)
        %first
        xyt = [x_src y_src]*R;

        %Set up hough matrix to count the votes for each offset
        offset_counts = zeros(offset_sz);

        %Loop through each line point from frame 1 and compute the offset to all
        %points in frame 2
        for ii = 1:length(x_src);

            %Computes offset and decentre
            x_offset = round(x_tgt - xyt(ii,1));
            y_offset = round(y_tgt - xyt(ii,2));

            %Only count votes less than the offset limit
            keep = (abs(x_offset) <= offset_lim & ...
                    abs(y_offset) <= offset_lim);

            %Use sparse matrix trick to counts votes
            if any(keep)
                offset_counts = offset_counts + ...
                    sparse(y_offset(keep)+offset_lim+1, ...
                           x_offset(keep)+offset_lim+1, ...
                           1, offset_sz, offset_sz);
            end
        end

        %Workout the offset with maximum vote for this theta
        [max_count_theta max_idx_theta] = max(offset_counts(:));

        %If this vote count is larger than the existing maximum record the
        %offset index and theta
        if max_count_theta > max_count
            max_count = max_count_theta;
            max_idx = max_idx_theta;
            max_theta = theta;
        end

    end

    %Convert the maximum offset ID
    [max_row max_col] = ind2sub([offset_sz offset_sz], max_idx);
    offsets(ff,:) = [max_col max_row] - offset_lim - 1;

    %Contvert theta to degrees
    thetas(ff) = 180*max_theta/pi;
    
    %Add the source frame to the target frames
    frames_reg = add_frame(frames(:,:,frame_idx(ff)), frames_reg, ...
                           offsets(ff,:), thetas(ff,:));
end


%% Test script
function test_script()
clc;

% Read input
imgroot = 'C:\isbe\nailfold\playground\dmk\capture_software\';
imgpath = fullfile(imgroot,'2012_10_16\Left\Digit4\x300\seq5\');
d = dir(fullfile(imgpath,'*.png'));

n_frames = length(d);
if (n_frames == 0), return; end

frame1 = imread(fullfile(imgpath,d(1).name));
frame1 = mean(frame1,3);

frames = zeros([size(frame1),n_frames]);
frames(:,:,1) = frame1;

for i = 2:n_frames
    frame = imread(fullfile(imgpath,d(i).name));
    frame = mean(frame,3);
    frames(:,:,i) = frame;
end

% Register frames
[frames_reg offsets thetas] = register_frames(frames, 20, 0);

% Write output
outpath = fullfile(imgpath,'registered');
if ~exist(outpath,'dir')
    mkdir(outpath);
end

for i = 1:n_frames
   filename = fullfile(outpath, sprintf('frame_%04d.png', i));
   imwrite(uint8(frames_reg(:,:,i)), filename);
end
