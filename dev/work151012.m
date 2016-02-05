base_dir = 'N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\';
seq_dir = [base_dir '035wellcome\2015_03_23\R2_12_24_58\'];
frames_dir = [seq_dir 'registered07\vessel03\'];
frame_list = dir([frames_dir 'frame*.png']);
n_frames = length(frame_list);
%%
for i_f = 4:5%[2 3 4 5]
    frames_dir = [seq_dir 'registered07\vessel' zerostr(i_f,2) '\'];
    [flowPyramidEst, flowConfidence] = ...
        estimate_flow_multilevel(frames_dir, n_frames, [frames_dir 'flow'], 1:3);
end