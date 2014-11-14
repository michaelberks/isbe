function [streams_count_matrix streams_sum_matrix ] = ...
    stream_connect_segments(prob_map, orientation_map, segment_map, step_length, num_streams, max_steps, segment_di)
%PROB_TRACK *Insert a one line summary here*
%   [path_map] = prob_track(I_mag, I_ori_theta, I_ori_D, xs, xy)
%
% Inputs:
%      I_prob - Map of structure (e.g. vessel) probabilites, most likely
%       returned from a classifier
%
%      I_ori_theta - Map of structures orientations - direction (i.e.
%      angle) component (this should NOT be doubled)
%
%      I_ori_D - Map of structures orientations - dispersion (i.e magnitude) component
%
%      x_orig - x coord of start point
%
%      y_orig - y coord of start point
%
%      step_length - does what it says on the tin
%
%      junction_map - (optional) binary map indicating whether there is
%       structure junction (i.e. bifurcation or crossing) present
%
%
% Outputs:
%      path_map - Map traced by path
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 23-Feb-2012
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

ori_theta = angle(orientation_map) / 2;
ori_D = abs(orientation_map);

neighbours_r = [-1 -1 -1; 0 0 0; 1 1 1];
neighbours_c = [-1 0 1; -1 0 1; -1 0 1];

max_dist = ceil(max_steps * step_length);
patch_sz = 2*max_dist + 1;
y_orig = max_dist + 1;
x_orig = max_dist + 1;

%Find segment ends
segment_labels = bwlabel(segment_map, 8);
segment_ends = segment_labels;
segment_ends(~bwmorph(segment_map, 'endpoints')) = 0;
segment_labels_dilated = imdilate(segment_labels, strel('disk', segment_di));
num_segs = max(segment_labels(:));

%Get number of ends
[ends_y ends_x] = find(segment_ends);
num_ends = length(ends_y);

streams_sum_matrix = zeros(num_segs, num_segs);
streams_count_matrix = zeros(num_segs, num_segs);
total_stream_sums = zeros(num_segs,1);

for i_end = 1:num_ends
    
    ex_i = ends_x(i_end);
    ey_i = ends_y(i_end);
    
    %Get patch of end nodes about this point
    labels_patch = sample_window(segment_labels_dilated, patch_sz, ey_i, ex_i, 0);
    
    %Mask out this end, from the ends patch
    seg_num = labels_patch(y_orig, x_orig);
    labels_patch(labels_patch == seg_num) = 0;
    
    %If we've got no other ends in the vicinity, we can give up now!
    if ~any(labels_patch(:))
        continue;
    end
    
    %Otherwise, get the other patches we need
    prob_patch = sample_window(prob_map, patch_sz, ey_i, ex_i, 0);
    ori_theta_patch = sample_window(ori_theta, patch_sz, ey_i, ex_i, 0);
    ori_D_patch = sample_window(ori_D, patch_sz, ey_i, ex_i, 1e-6);  

    %Workout initial direction from connected neighbours
    neighbours = sample_window(segment_map, 3, ey_i, ex_i, 0);
    xs_orig = mean(neighbours_c(neighbours));
    ys_orig = mean(neighbours_r(neighbours));

    l0 = sqrt(xs_orig^2 + ys_orig^2);
    xs_orig = -xs_orig/l0;
    ys_orig = -ys_orig/l0;
    
    %preallocate path map
    path_map = zeros(patch_sz, patch_sz);

    %get orientation at initisl point
    for i_stream = 1:num_streams
        %Set up intitial points and driection
        xi = x_orig;
        yi = y_orig;
        xii = round(xi);
        yii = round(yi);
        
        xs0 = xs_orig;
        ys0 = ys_orig;

        go_on = true;
        step = 0;
        stream_sum = 0;
        while go_on
            %Increment the step count
            step = step + 1;

            %sample from orientation distribution (defined by direction and dispersion) at current location
            theta = ori_theta_patch(yii, xii);
            rho = ori_D_patch(yii, xii);
            new_theta = theta + wrapped_normal_sample(0, rho, 1)/2;

            %Convert direction into x,y steps
            xs1 = cos(new_theta);
            ys1 = -sin(new_theta);

            %Check we're not going back on ourselves *by more than 45
            %degrees*
            if xs1*xs0 + ys1*ys0 < 0%-0.5
                xs1 = -xs1;
                ys1 = -ys1;
            end

            %Step in the chosen direction
            xi = xi + step_length*xs1;
            yi = yi + step_length*ys1;
            xs0 = xs1;
            ys0 = ys1;

            %Check path hasn't gone off the edge of the image
            xii = round(xi);
            yii = round(yi);

            %Add the accumulated probability at this location
            stream_sum = stream_sum + prob_patch(yii, xii);
            path_map(yii, xii) = path_map(yii, xii) + 1;

            %Stop if we've reached max steps or hit another segments end
            if (step == max_steps) || labels_patch(yii, xii)
                go_on = false;
            end        
        end
        
        %Divide the stream sum by the stream length
        stream_sum = stream_sum / step;
        
        %Update the total streams
        total_stream_sums(seg_num) = total_stream_sums(seg_num) + stream_sum;
        
        %If we hit an end, update the sum for that end 
        if labels_patch(yii, xii)
            streams_sum_matrix(seg_num, labels_patch(yii, xii)) = ...
                streams_sum_matrix(seg_num, labels_patch(yii, xii)) + stream_sum;
            
            streams_count_matrix(seg_num, labels_patch(yii, xii)) = ...
                streams_count_matrix(seg_num, labels_patch(yii, xii)) + 1;
        end

    end
end

0;
            
            
    
    





