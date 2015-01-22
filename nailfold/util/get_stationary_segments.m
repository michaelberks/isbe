function [segments_s, segments_ns, motor_x, motor_y, motor_z, sharpness] = ...
    get_stationary_segments(sequence)
%GET_STATIONARY_SEGMENTS
%   [segments_s, segments_ns] = get_stationary_segments(sequence)
% Inputs:
%      sequence - *Insert description of input variable here*
%
%
% Outputs:
%      segments_s - *Insert description of input variable here*
%
%      segments_ns - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 20-Jan-2015
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

%Get the x,y motor positions from the sequence structure
num_frames = sequence.num_frames;
if nargout > 2
    motor_x = zeros(num_frames,1);
    motor_y = zeros(num_frames,1);
    motor_z = zeros(num_frames,1);
    sharpness = zeros(num_frames,1);
end



%Note that the below is a very non-matlab way of coding this (usually I'd
%just use bwareaopen and bwlabel on the complete list of frame positions)
%however I'm writing this as a test bed for a C++ implementation so doing it
% the 'old-fashioned' loopy way!
is_moving = true;
min_stat_frames = 20;
stat_count = 0;

segments_s = cell(0,1);
segments_ns = cell(0,1);

mx0 = sequence.frames(1).motor_x;
my0 = sequence.frames(1).motor_y;

for i_f = 1:num_frames
    
    mx = sequence.frames(i_f).motor_x;
    my = sequence.frames(i_f).motor_y;
    
    moving_frame = abs(mx-mx0) || abs(my-my0);
    mx0 = mx;
    my0 = my;
    
    if (~moving_frame)
        
        if (is_moving)
            %We were in moving segment, so this is potentially the first
            %frame in a new stationary segment
            stat_count = 1;
            is_moving = false;
            segments_s(end+1,1) = {i_f}; %#ok
  
        else
            %We weren't moving so add this to stationary segment
            segments_s{end}(end+1,1) = i_f;
            
            %If we're still in the probationary period, add to the
            %stationary count
            if (stat_count < min_stat_frames)
                stat_count = stat_count + 1;          
            end           
        end
        
    else
        if (is_moving)
            %We were already moving so add this to the current moving segment
            segments_ns{end}(end+1,1) = i_f;
            
        else
            is_moving = true;
            if (stat_count < min_stat_frames)
                %We're weren't moving, but were in the probationary period -
                %move all the frames from the temporary stationary segment into
                %the previous moving segment, and append this frame
                
                if isempty(segments_ns)
                    segments_ns(end+1,1) = segments_s(end); %#ok
                else
                    segments_ns{end} = [segments_ns{end}; segments_s{end}];
                end
                segments_s(end) = [];
                segments_ns{end}(end+1,1) = i_f;
                
            else
                %We were in a real stationary segment, so end that and
                %start a new non_stationary segment
                segments_ns(end+1,1) = {i_f}; %#ok
            end
        end               
    end
    
    if nargout > 2
        motor_x(i_f) = mx;
        motor_y(i_f) = my;
        motor_z(i_f) = sequence.frames(i_f).motor_z;
        sharpness(i_f) = sequence.frames(i_f).sharpness;
    end
end
    
