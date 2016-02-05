function [is_connected] = is_mosaic_connected(frame_centres, frame_w, frame_h)
%IS_MOSAIC_CONNECTED *Insert a one line summary here*
%   [is_connected] = is_mosaic_connected(frame_centres, frame_w, frame_h)
%
% Inputs:
%      frame_centres - *Insert description of input variable here*
%
%      frame_w - *Insert description of input variable here*
%
%      frame_h - *Insert description of input variable here*
%
%
% Outputs:
%      is_connected - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 12-Jan-2016
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
num_segs = size(frame_centres,1);
x_diffs = zeros(num_segs);
y_diffs = zeros(num_segs);
for i_seg = 1:num_segs
    x_diffs(:,i_seg) = abs(frame_centres(:,1) - frame_centres(i_seg,1));
    y_diffs(:,i_seg) = abs(frame_centres(:,2) - frame_centres(i_seg,2));
end

num_components = 0;
adjacency_matrix = (x_diffs < frame_w) & (y_diffs < frame_h);

%Result array.
marked = zeros(num_segs, 1);

%Enumerate all vertices, 
for i_vi = 1:num_segs

    %if for vertex number i, marks[i] == 0 then
    if ~marked(i_vi)
        
        %Increment components
        num_components = num_components+1;

        %Put this vertex into queue, and 
        queue = i_vi;
        
        while (~isempty(queue)) 
            %Pop current vertex from queue
        	current = queue(end);
            queue(end) = [];

            %Add all adjacent and not currently marked vertices to the
            %queue
            for i_vj = 1:num_segs
                if (adjacency_matrix(current, i_vj) && ~marked(i_vj))
                    marked(i_vj) = num_components;
                    queue(end+1) = i_vj; %#ok
                end
            end
        end
    end
end

is_connected = num_components == 1;
