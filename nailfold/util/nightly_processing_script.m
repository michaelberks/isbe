study_dir = 'N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\';

subject_dirs = [
    dir([study_dir '0*'])
    dir([study_dir '1*'])];

num_subs = length(subject_dirs);
recruitment_dates = zeros(num_subs,1);

sequence_dirs = cell(0,1);
for i_sub = 1:num_subs
    session_dirs = dir([study_dir subject_dirs(i_sub).name '\2015*']);
    
    for i_ses = 1:length(session_dirs)
        sequence_dirs_i = [
            dir([study_dir subject_dirs(i_sub).name '\' session_dirs(i_ses).name '\L*'])
            dir([study_dir subject_dirs(i_sub).name '\' session_dirs(i_ses).name '\R*'])];
        
        for i_seq = 1:length(sequence_dirs_i)
            sequence_dirs{end+1} = ...
                [study_dir subject_dirs(i_sub).name '\' session_dirs(i_ses).name '\' sequence_dirs_i(i_seq).name '\']; %#ok
        end
    end
end

num_seqs = length(sequence_dirs) - 2;
%%
last_sequence = u_load('C:\isbe\nailfold\data\wellcome_study\last_sequence.mat');
for i_seq = last_sequence:num_seqs
    clc;
    display(['processing sequence ' num2str(i_seq) ' of ' num2str(num_seqs)]);
    display(['Sequence ' sequence_dirs{i_seq}]);
    register_sequence(... % the user's input
        'sequence',                 [],...
        'sequence_data_path',       [sequence_dirs{i_seq} '\sequence_frames_data.dat'],...
        'sequence_dir',             sequence_dirs{i_seq},...
        'dirt_image',               [],...
        'pixels_per_mm',            1000,...
        'min_stationary_frames',    120,...
        'mosaic_lims_x',            -310:309,...
        'mosaic_lims_y',            -230:229,...
        'frame_h',                  480,...
        'frame_w',                  640,...
        'sigma',                    6,...
        'intra_segment_offset',     120,...
        'inter_segment_offset',     240,...
        'plot',                     0,...
        'debug',                    0);
    last_sequence = last_sequence+1;
    save('C:\isbe\nailfold\data\wellcome_study\last_sequence.mat','last_sequence');
end