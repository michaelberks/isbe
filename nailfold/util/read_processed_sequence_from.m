function [sequence_data] = read_processed_sequence_from(filename)
% Parse the *_markup.txt file output by the ncm_qmarkup.exe application
% into a structure whose fields have corresponding names.

%Get string stream from file
fid = fopen(filename);
frewind(fid);
s = textscan(fid, '%s', 'commentstyle', '//'); 
fclose(fid);
s = s{1};

%Find property tags and extract values
sequence_data = [];
sequence_data.sequence_name = s{find(strcmpi(s, 'Sequence_name:'), 1) + 1};
sequence_data.num_frames = str2double(s{find(strcmpi(s, 'Num_frames:'), 1) + 1});
sequence_data.num_segments = str2double(s{find(strcmpi(s, 'Num_segments:'), 1) + 1});
if sequence_data.num_segments < 2
    return;
end
sequence_data.sequence_is_valid = str2double(s{find(strcmpi(s, 'Sequence_is_valid:'), 1) + 1}) > 0;
sequence_data.num_required_segments = str2double(s{find(strcmpi(s, 'Num_required_segments:'), 1) + 1});
fi = find(strcmpi(s, 'Final_mosaic_size:'), 1);
sequence_data.final_mosaic_size = [str2double(s{fi + 1}) str2double(s{fi + 2})];
fi = find(strcmpi(s, 'Frame_dims:'), 1);
sequence_data.frame_dims = [str2double(s{fi + 1}) str2double(s{fi + 2})];
fi = find(strcmpi(s, 'Inner_frame_dims:'), 1);
sequence_data.inner_frame_dims = [str2double(s{fi + 1}) str2double(s{fi + 2})];
sequence_data.pixels_per_mm = str2double(s{find(strcmpi(s, 'Pixels_per_mm:'), 1) + 1});
sequence_data.intra_segment_offset = str2double(s{find(strcmpi(s, 'Intra_segment_offset:'), 1) + 1});
sequence_data.inter_segment_offset = str2double(s{find(strcmpi(s, 'Inter_segment_offset:'), 1) + 1});
sequence_data.min_segment_size = str2double(s{find(strcmpi(s, 'Min_segment_size:'), 1) + 1});
sequence_data.min_frames_for_flow = str2double(s{find(strcmpi(s, 'Min_frames_for_flow:'), 1) + 1});
sequence_data.min_frames_for_mosaic = str2double(s{find(strcmpi(s, 'Min_frames_for_mosaic:'), 1) + 1});
sequence_data.frame_border = str2double(s{find(strcmpi(s, 'Frame_border:'), 1) + 1});

%Find braces and extract alignment data
s1 = find(strcmp(s, '{'));
s2 = find(strcmp(s, '}'));
ss = s(s1+1:s2-1);

%Current version is missing a entry for the first row, fill this in as a
%zero
n = length(ss);
ssn = zeros(n,1);
offset = 0;
for i_s = 1:n; 
    str_i = ss{i_s};
    if str_i(end) == ';'
        str_i(end) = [];
        d = str2double(str_i);
        if d == -10
            ssn(i_s+offset) = -1;
            ssn(i_s+offset+1) = 0;
            offset = offset+1;
        else
            ssn(i_s+offset) = d;
        end
    else
        ssn(i_s+offset) = str2double(str_i);
    end
     
end
ssn = reshape(ssn, 8, [])';

sequence_data.segment_rejected = ssn(:,2) > 0;
sequence_data.segment_required = ssn(:,3) > 0;
sequence_data.segment_matched = ssn(:,8) > 0;
sequence_data.segment_motor_centres = ssn(:,4:5);
sequence_data.segment_displacements = ssn(:,6:7);
% sequence_data.segment_required(1) = 1;
% sequence_data.segment_matched(1) = 1;
