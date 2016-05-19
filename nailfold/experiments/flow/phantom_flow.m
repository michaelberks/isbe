root_dir = 'O:\flow_project\camera_frames\500um_videos\';
%%
seq_dirs = {...
    [root_dir '\4ml\h_500_1\2016_05_03\'];
    [root_dir '\5ml\h_500_1\2016_05_03\'];
    [root_dir '\5.5ml\h_500_1\2016_05_03\'];
    [root_dir '\6ml\h_500_1\2016_05_03\'];
    [root_dir '\6.5ml\h_500_1\2016_05_03\'];
    [root_dir '\7ml\h_500_1\2016_05_03\'];
    [root_dir '\7.5ml\h_500_1\2016_05_03\'];
    [root_dir '\8ml\h_500_1\2016_05_03\'];
    [root_dir '\8.5ml\h_500_1\2016_05_03\'];
    [root_dir '\9ml\h_500_1\2016_05_03\'];
    [root_dir '\9.5ml\h_500_1\2016_05_03\'];
    [root_dir '\10ml\h_500_1\2016_05_03\'];
    [root_dir '\70ml\h_500_1\2016_05_03\']};

%%
for i_s = 9:length(seq_dirs)
    seq_dir = seq_dirs{i_s};
    frames = dir([seq_dir '*.bmp']);

    num_repeats = length(frames) / 180;
    for i_rpt = 1:num_repeats
        mkdir([seq_dir 'v' zerostr(i_rpt,2)]);
        for i_f = (i_rpt-1)*180 + (1:180)
            copyfile([seq_dir frames(i_f).name],...
                [seq_dir 'v' zerostr(i_rpt,2) '\frame' zerostr(i_f - 180*(i_rpt-1),4) '.bmp']);
        end
        cmd = ['ffmpeg -y -r 60 -i "' seq_dir 'v' zerostr(i_rpt,2) '\frame%04d.bmp" -c:v libx264 -preset slow -crf 18 -an "' seq_dir 'v' zerostr(i_rpt,2) '\flow.mp4"'];
        system(cmd);
        delete([seq_dir 'v' zerostr(i_rpt,2) '\*.bmp']);
    end
end
%%
for i_s = 8%1:length(seq_dirs)
    seq_dir = seq_dirs{i_s};
    frames_list = dir([seq_dir '*.bmp']);

    num_repeats = length(frames_list) / 180;
    for i_rpt = 2%1:num_repeats
        frames = zeros(480, 160, 180);
        
        for i_f = (i_rpt-1)*180 + (1:180)
            frame = imread([seq_dir frames_list(i_f).name]);
            frames(:,:,i_f - (i_rpt-1)*180) = double(frame(:,240 + (1:160)));
        end
        
        %Contrast normalise the frames
        g_lims = prctile(frames(:), [1 99]);
        g_range = g_lims(2) - g_lims(1);

        %Compute flow results
        flow_results = [];
        [flow_results.flowPyramidEst, flow_results.flowConfidence] = ...
            estimate_flow_multilevel(255*(frames-g_lims(1))/g_range, [],[],1:4);   
        flow_results.g_lims = g_lims;   
        save([seq_dir 'v' zerostr(i_rpt,2) '\flow_results.mat'], 'flow_results');
        
    end
end

