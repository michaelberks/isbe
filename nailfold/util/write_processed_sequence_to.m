function [] = write_processed_sequence_to(filename, sequence_data)
% Parse the *_markup.txt file output by the ncm_qmarkup.exe application
% into a structure whose fields have corresponding names.

%Get string stream from file
fid = fopen(filename, 'wt');

%Find property tags and extract values
fprintf(fid, '%s %s\n', 'Sequence_name:', sequence_data.sequence_name);
fprintf(fid, '%s %d\n', 'Num_frames:', sequence_data.num_frames);
fprintf(fid, '%s %d\n', 'Num_segments:', sequence_data.num_segments);
if sequence_data.num_segments < 2
    fclose(fid);
    return;
end

fprintf(fid, '%s %d\n', 'Sequence_is_valid:', sequence_data.sequence_is_valid);
fprintf(fid, '%s %d\n', 'Num_required_segments:', sequence_data.num_required_segments);
fprintf(fid, '%s\n', 'Segment_alignment: {');
for i_seg = 1:sequence_data.num_segments
    fprintf(fid, '%d %d %d %.4f %.4f %.4f %.4f %d\n', ...
        i_seg,...
        sequence_data.segment_rejected(i_seg),...
        sequence_data.segment_required(i_seg),...
        925*sequence_data.segment_motor_centres(i_seg,1)/1000,...
        -925*sequence_data.segment_motor_centres(i_seg,2)/1000,...
        sequence_data.segment_displacements(i_seg,1),...
        sequence_data.segment_displacements(i_seg,2),...
        sequence_data.segment_matched(i_seg,1));
end
fprintf(fid, '%s\n', '} //Segment_alignment');
fprintf(fid, '%s %d %d\n', 'Final_mosaic_size:', sequence_data.final_mosaic_size(1), sequence_data.final_mosaic_size(2));
fprintf(fid, '%s %d %d\n', 'Frame_dims:', sequence_data.frame_dims(1), sequence_data.frame_dims(2));
fprintf(fid, '%s %d %d\n', 'Inner_frame_dims:', sequence_data.inner_frame_dims(1), sequence_data.inner_frame_dims(2));
fprintf(fid, '%s %.4f\n', 'Pixels_per_mm:', sequence_data.pixels_per_mm);
fprintf(fid, '%s %.4f\n', 'Intra_segment_offset:', sequence_data.intra_segment_offset);
fprintf(fid, '%s %.4f\n', 'Inter_segment_offset:', sequence_data.inter_segment_offset);
fprintf(fid, '%s %d\n', 'Min_segment_size:', sequence_data.min_segment_size);
fprintf(fid, '%s %d\n', 'Min_frames_for_flow:', sequence_data.min_frames_for_flow);
fprintf(fid, '%s %d\n', 'Min_frames_for_mosaic:', sequence_data.min_frames_for_mosaic);
fprintf(fid, '%s %d\n', 'Frame_border:', sequence_data.frame_border);

fclose(fid);