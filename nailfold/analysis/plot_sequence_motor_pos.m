function [] = plot_sequence_motor_pos(varargin)
%REGISTER_SEQUENCE *Insert a one line summary here*
%   [] = register_sequence()
%
% Inputs:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 15-Jun-2015
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
args = u_packargs(varargin, 0, ... % the user's input
    'sequence',                 [],...
    'sequence_data_path',       [],...
    'sequence_dir',             [],...
    'pixels_per_mm',            1000,...
    'min_stationary_frames',    120);
clear varargin;

%Check if we've already got the sequence structure, if not load it in or
%compute from sequence properties (latter will be slooowwwww)
if isempty(args.sequence)
    if isempty(args.sequence_data_path)
         args.sequence = read_sequence_from([args.sequence_dir '\sequence_properties.txt']);
    else
        args.sequence = u_load(args.sequence_data_path);
    end
end

%--------------------------------------------------------------------------
%Compute stationary segments, discard any that were stationary for less
%than 120 frames (1 sec)
[segments_s, segments_ns, motor_x, motor_y, motor_z, sharpness] = ...
    get_stationary_segments(args.sequence, args.min_stationary_frames);

plot_segment_traces(args.pixels_per_mm, segments_s, segments_ns, motor_x, motor_y, motor_z, sharpness);

%--------------------------------------------------------------------------
%**************************************************************************
%**************************************************************************
function plot_segment_traces(pixels_per_mm, segments_s, segments_ns, motor_x, motor_y, motor_z, sharpness)

xy_speed = [0; sqrt(diff(motor_x).^2 + diff(motor_y).^2)];
z_speed = [0; diff(motor_z)];   
frame_nums = 1:length(xy_speed);

figure;
a1 = subplot(3,1,1);hold on; 
plot(xy_speed, 'b');
plot(abs(z_speed), 'r');
xlabel('Time');
ylabel('Motor speed');
legend({'|XY| speed', '|Z| speed'});
set(a1, 'ylim', [ 0 0.04]);

%plot(sharpness_change / (max(sharpness_change)-min(sharpness_change)), 'g--');
a2 = subplot(3,1,2); hold on;
xlabel('Time');
ylabel('Z position');
plot(0,0,'r-');
plot(0,0,'b-');
legend(a2, {'Stationary/moving in Z', '|XY| movement'});

a3 = subplot(3,1,3); hold on;
xlabel('Time');
ylabel('Sharpness');
plot(0,0,'r-');
plot(0,0,'b-');
legend(a3, {'Stationary/moving in Z', 'XY movement'});

for i_sg = 1:length(segments_s)
    plot(a2, frame_nums(segments_s{i_sg}), motor_z(segments_s{i_sg}), 'r.');
    plot(a3, frame_nums(segments_s{i_sg}), sharpness(segments_s{i_sg}), 'r.');
end

for i_sg = 1:length(segments_ns)
    plot(a2, frame_nums(segments_ns{i_sg}), motor_z(segments_ns{i_sg}), 'b.');
    plot(a3, frame_nums(segments_ns{i_sg}), sharpness(segments_ns{i_sg}), 'b.');
end
    
hw = 320 / pixels_per_mm;
hh = 240 / pixels_per_mm;
frame_rect_x = [-hw hw hw -hw -hw];
frame_rect_y = [-hh -hh hh hh -hh];

figure; hold all; axis equal; a1 = gca;
plot(a1, motor_x, motor_y, 'b.', 'markersize', 4);
%
for i_seg = 1:length(segments_s)

    frame_idx_i = segments_s{i_seg};

    seg_x = motor_x(frame_idx_i);
    seg_y = motor_y(frame_idx_i);

    cx = median(seg_x);
    cy = median(seg_y);

    plot(a1, cx, cy, 'rx');
    plot(a1, cx+frame_rect_x, cy+frame_rect_y, 'r', 'linewidth', 2);  
end