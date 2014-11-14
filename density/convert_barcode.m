function [barcode barcode_warning] = convert_barcode(in_pts, barcode_length)

curr_state = 0;
curr_pos = 1;
curr_length = 1;
segment_lengths = [];

num_pts = length(in_pts);
step_length = num_pts / barcode_length;
barcode = zeros(1, barcode_length);
for ii = 2:num_pts+1
    if ii <= num_pts && in_pts(ii) == curr_state
        curr_length = curr_length + 1;
    else
        num_bits = round(curr_length / step_length);
        segment_lengths = [segment_lengths; curr_length]; %#ok
        barcode(curr_pos:curr_pos+num_bits-1) = curr_state;
        
        curr_pos = curr_pos + num_bits;
        curr_state = abs(1 - curr_state);
        curr_length = 1;
    end
end

barcode_warning = ...
    (curr_pos ~= barcode_length + 1) ||... %Wrong number of bits in total
    (sum(abs(.5 - rem(segment_lengths / step_length, 1) ) < .2) > 1); %At least two segments of ambiguous number of bits
    