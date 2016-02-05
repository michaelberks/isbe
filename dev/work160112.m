%Get list of sequences to copy
args.study_dir = 'N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\';
fid = fopen([args.study_dir 'subject_summary.txt']);
    frewind(fid);
    s = textscan(fid, '%s'); 
    fclose(fid);
    s = s{1};
    subject_summary = reshape(s, 3, [])';
    [~,s_idx] = sort(subject_summary(:,1));
    subject_summary = subject_summary(s_idx,:);
%%
    num_subs = size(subject_summary,1);
    sequence_names = cell(num_subs, 10);
    
    for i_sub = 1:num_subs
        session_dirs = dir([args.study_dir subject_summary{i_sub,1} '\2015*']);

        num_sessions = length(session_dirs);
        if ~num_sessions
            display(['Subject ' num2str(i_sub) ' has no sessions available.']);
            continue;
        else
            if num_sessions > 1
                display(['Subject ' num2str(i_sub) ' has ' num2str(num_sessions) ' sessions']);
            end
        end
        for i_ses = 1:num_sessions
            curr_seq = 1;
            session_dir = [args.study_dir subject_summary{i_sub,1} '\' session_dirs(i_ses).name '\'];
            for hand = 'LR'
                for digit = '12345'
                    sequence_dirs = ...
                        dir([session_dir hand digit '*']);

                    num_sequences = length(sequence_dirs);
                    if ~num_sequences
                        display(['Subject ' num2str(i_sub) ' has no sequences for ' hand digit]);
                    else
                        sequence_dir = sequence_dirs(end).name;
                        if num_sequences > 1
                            display(['Subject ' num2str(i_sub) ' has ' num2str(num_sequences) ' sequences for ' hand digit '. Using the latest:' sequence_dir]);
                        end
                        if ~exist([session_dir sequence_dir '\sequence_frames_data.dat'], 'file')
                            display(['Sequence ' sequence_dir ' for subject ' num2str(i_sub) ' has no sequence frames data file']);
                        else
                            sequence_names{i_sub,curr_seq} = [session_dir sequence_dir];
                        end
                    end
                    curr_seq = curr_seq+1;
                end
            end
        end
    end
%%
for i_sub = 1:num_subs
    for i_dig = 1:10
        if isempty(sequence_names{i_sub,i_dig})
            display(['Subject ' num2str(i_sub) ' has no sequences for digit ' num2str(i_dig)]);
        end
    end
end
%%
missing_sequence_frames = cell(0,1);
missing_sequence_mosaic = cell(0,1);
missing_vessel_preds = cell(0,1);
missing_apex_preds = cell(0,1);
present_and_correct = cell(0,1);

missing_sequence_frames_id = cell(0,1);
missing_sequence_mosaic_id = cell(0,1);
missing_vessel_preds_id = cell(0,1);
missing_apex_preds_id = cell(0,1);
present_and_correct_id = cell(0,1);
%
for i_seq = 1:num_seqs
    
    seq_dir = sequence_dirs{i_seq};
    if ~exist([seq_dir 'sequence_frames_data.dat'], 'file')
        missing_sequence_frames{end+1,:} = seq_dir; %#ok
        missing_sequence_frames_id{end+1,:} = i_seq; %#ok
    elseif ~exist([seq_dir 'sequence_data\full_mosaic.png'], 'file')
        missing_sequence_mosaic{end+1,:} = seq_dir; %#ok
        missing_sequence_mosaic_id{end+1,:} = i_seq; %#ok
    elseif ~exist([seq_dir 'capillary_data\vessels_v_pred.txt'], 'file')
        missing_vessel_preds{end+1,:} = seq_dir; %#ok
        missing_vessel_preds_id{end+1,:} = i_seq; %#ok
    elseif ~exist([seq_dir 'capillary_data\apex_candidates.txt'], 'file')
        missing_apex_preds{end+1,:} = seq_dir; %#ok
        missing_apex_preds_id{end+1,:} = i_seq; %#ok
    else
        present_and_correct{end+1,:} = seq_dir; %#ok
        present_and_correct_id{end+1,:} = i_seq; %#ok
    end

end
%%
filename = 'N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\subject_summary.txt';
fid = fopen(filename);
frewind(fid);
s = textscan(fid, '%s'); 
fclose(fid);
s = s{1};
subject_summary = reshape(s, 3, [])';
[~,s_idx] = sort(subject_summary(:,1));
subject_summary = subject_summary(s_idx,:);
%%
group_codes = str2num(cell2mat(subject_summary(:,3))); %#ok
recruitment_im = repmat(group_codes', 50, 1);
figure; imgray(recruitment_im); colormap(lines(4));

write_im_from_colormap(recruitment_im, 'C:\isbe\nailfold\data\wellcome_study\recruitment_im.bmp', lines(4));
%%
sequence_name = [sequence_names{112, 10} '\'];
plot_sequence_motor_pos(...
    'sequence',                 [],...
    'sequence_data_path',       [sequence_name 'sequence_frames_data.dat'],...
    'sequence_dir',             sequence_name,...
    'pixels_per_mm',            925,...
    'min_stationary_frames',    30);
plot_sequence_motor_pos(...
    'sequence',                 [],...
    'sequence_data_path',       [sequence_name 'sequence_frames_data.dat'],...
    'sequence_dir',             sequence_name,...
    'pixels_per_mm',            925,...
    'min_stationary_frames',    120);
%%
pixels_per_mm = 925;
frame_w = 640 / pixels_per_mm;
frame_h = 480 / pixels_per_mm;

sequences_that_change = [];
for i_sub = 1:num_subs
    display(['Checking subject ' num2str(i_sub)]);
    for i_dig = 1:10
        if isempty(sequence_names{i_sub,i_dig})
            continue;
        end
            
        sequence = load([sequence_names{i_sub,i_dig} '\sequence_frames_data.dat']);
        [segments_s, segments_ns, motor_x, motor_y] = get_stationary_segments(sequence, 30);
        num_segs = length(segments_s);
        frame_centres = zeros(num_segs,2);

        kept_segs = true(num_segs,1);
        for i_seg = 1:num_segs
            frame_centres(i_seg,1) = motor_x(segments_s{i_seg}(1));
            frame_centres(i_seg,2) = motor_y(segments_s{i_seg}(1));
            if size(segments_s{i_seg},1) < 120
                kept_segs(i_seg) = 0;
            end
        end
        kept_segs([1 end]) = 1;
        [is_connected1] = is_mosaic_connected(frame_centres, frame_w, frame_h);
        [is_connected2] = is_mosaic_connected(frame_centres(kept_segs,:), frame_w, frame_h);
    
        if is_connected1 && ~is_connected2
            sequences_that_change{end+1,1} = sequence_names{i_sub,i_dig}; %#ok
        end
    end
end
%%
copy_dir = 'C:\isbe\nailfold\data\wellcome_study\copy_scripts\';

copy_to_file_name = [copy_dir 'copy_to_rds_dregs.sh'];
fid_to = fopen(copy_to_file_name, 'wt');

copy_from_file_name = [copy_dir 'copy_from_rds_dregs.sh'];
fid_from = fopen(copy_from_file_name, 'wt');

for i_seq = 1:length(sequences_that_change)

    seq_name = sequences_that_change{i_seq}(67:end);
    dividers = seq_name == '\';
    seq_name(dividers) = '/';
    subject_name = seq_name(1:12);
    rds_name = [subject_name seq_name(24:end)];
    make_cmd = ['mkdir rds-scratch/nailfold/' subject_name];
    fprintf(fid_to, '%s\n', make_cmd);
    
    copy_cmd = ['rsync -ra --update --exclude-from="rds-scratch/ignore_file.txt" ' ...
        'mhsrfs01/epi-musculo/Nailfold\ Capillaroscopy/camera_capture/wellcome_nailfold_study/' ...
        seq_name ' rds-scratch/nailfold/' subject_name];
    display(copy_cmd);
    
    fprintf(fid_to, '%s\n', copy_cmd);

    copy_cmd = ['rsync -ra --update --exclude-from="rds-scratch/ignore_file_from.txt" ' ...
        'rds-scratch/nailfold/' rds_name ...
        ' mhsrfs01/epi-musculo/Nailfold\ Capillaroscopy/camera_capture/wellcome_nailfold_study/' seq_name];
    display(copy_cmd);
    fprintf(fid_from, '%s\n', copy_cmd);

end
fclose(fid_to);
fclose(fid_from);
%%
fid = fopen('C:\isbe\nailfold\data\wellcome_study\sequences_to_process.txt', 'wt');
for i_seq = 1:length(sequences_that_change)    
    seq_name = sequences_that_change{i_seq}([67:77 89:end]);
    dividers = seq_name == '\';
    seq_name(dividers) = '/';
    fprintf(fid, '-I rds-scratch/nailfold/%s/sequence_frames_data.dat -p rds-scratch/models/csf_capillary_detection_w3_parameters.txt\n', seq_name);
end
fclose(fid);
%%
for i_seq = 1:length(sequences_that_change)
    display(['Moving sequence ' num2str(i_seq)]);
    seq_dir = sequences_that_change{i_seq};
    movefile([seq_dir '/sequence_data'], [seq_dir '/old_sequence_data']);
    movefile([seq_dir '/capillary_data'], [seq_dir '/old_capillary_data']);
    mkdir([seq_dir '/old_reg']);
    movefile([seq_dir '/registered*'], [seq_dir '/old_reg/']);
end
%%
for i_seq = [4 8 11 16 24 29 30 32 34:37] %1:length(sequences_that_change)
    display(['Moving sequence ' num2str(i_seq)]);
    seq_dir = sequences_that_change{i_seq};
%     movefile([seq_dir seq_dir(89:end) '\*'], [seq_dir '\']);
%     rmdir([seq_dir seq_dir(89:end)]);
    
    if exist([seq_dir '\capillary_data\vessels_v_pred.txt'], 'file')
        im = imread([seq_dir '\capillary_data\vessels_v_pred.txt']);
        imwrite(im, [seq_dir '\capillary_data\vessels_v_pred.png']);

        im = imread([seq_dir '\capillary_data\vessels_o_pred.txt']);
        imwrite(im, [seq_dir '\capillary_data\vessels_o_pred.png']);

        im = imread([seq_dir '\capillary_data\vessels_w_pred.txt']);
        imwrite(im, [seq_dir '\capillary_data\vessels_w_pred.png']);
    else
        display(['Missing capillary data: ' seq_dir(67:end)]);
    end
end
    
%%
register_sequence(... % the user's input
    'sequence',                 [],...
    'sequence_data_path',       [sequence_name 'sequence_frames_data.dat'],...
    'sequence_dir',             sequence_name,...
    'selected_segments',        [],...[3 6 7 11 14]
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
    'plot',                     1,...
    'debug',                    0);
%%
%tidy frames
sequence_names2 = sequence_names';
for i_seq = [347 504]%numel(sequence_names)
    
    seq_dir = [sequence_names2{i_seq} '\'];
    display(['Tidying sequence ' num2str(i_seq) ': ' seq_dir]);
    
    reg_dirs = dir([seq_dir 'registered*']);
    for i_reg = 1:length(reg_dirs)
        f1 = dir([seq_dir reg_dirs(i_reg).name '\frame*.png']);
        f1 = {f1.name}';
        matches = regexp(f1, 'frame\d\d\d\d\d.png');
        num_frames = length(f1);
        remaining = num_frames;
        for i_frame = 1:num_frames
            if isempty(matches{i_frame})
                delete([seq_dir reg_dirs(i_reg).name '\' f1{i_frame}]);
                remaining = remaining - 1;
            else
                break;
            end
        end
        if exist([seq_dir reg_dirs(i_reg).name '\seq_mask.png'], 'file')
            delete([seq_dir reg_dirs(i_reg).name '\seq_mask.png']);
        end
        if exist([seq_dir reg_dirs(i_reg).name '\movie.mp4'], 'file')
            delete([seq_dir reg_dirs(i_reg).name '\movie.mp4']);
        end
        if ~remaining
            rmdir([seq_dir reg_dirs(i_reg).name '\'], 's');
        end
    end
end
%%
capillary_data_dir = 'C:\isbe\nailfold\data\wellcome_study\capillary_data\';

capillary_type = {'distal', 'nondistal'};
fields = {...
    'width_at_apex' ...
	'mean_width' ...
	'median_width' ...
	'max_width' ...
	'min_width' ...
	'std_width' ...
	'mean_weighted_width' ...
	'mean_weighted_prob' ...
	'total_prob' ...
	'orientation_hist' ...
	'base_orientation' ...
	'initial_connected_orientation' ...
	'connected_orientation' ...
	'weighted_orientation' ...
	'apex_xy'};
num_fields = length(fields);

for i_seq = 1:numel(sequence_names)
    
    if isempty(sequence_names2{i_seq})
        continue;
    end
    
    seq_dir = [sequence_names2{i_seq} '\'];%sequence_dirs{i_seq};
    seq_name = seq_dir;
    dividers = seq_name == '\';
    pos = find(dividers, 4, 'last');
    seq_name(dividers) = '_';
    seq_name = seq_name(pos(1)+1:end-1);

    load([capillary_data_dir seq_name '_capillary_data.mat'], 'apex_measures');
    if isempty(apex_measures.distal) && isempty(apex_measures.nondistal)
        continue;
    end
    modified = false;
    for i_type = 1:2
        
        discard_apices = ~(apex_measures.(capillary_type{i_type}).mean_weighted_width > 0);
        
        if any(discard_apices)
            display(['Discarding ' num2str(sum(discard_apices)) ' ' capillary_type{i_type} ' apexes from ' seq_name]);
            
%             for i_field = 1:num_fields
%                 apex_measures.(capillary_type{i_type}).(fields{i_field})(discard_apices,:) = nan;
%             end
            apex_measures.(capillary_type{i_type}).bounding_box(:,:,discard_apices) = nan;
            modified = true;
        end
    end
    if modified
        save([capillary_data_dir seq_name '_capillary_data.mat'], 'apex_measures', '-append');
    end
    
end
                
%%
load('C:\isbe\nailfold\data\wellcome_study\capillary_data\001wellcome_2015_02_27_L1_09_38_22_capillary_data.mat');
vessel_pred = imread([sequence_names{1} '\capillary_data\vessels_v_pred.png']);

v_thresh = 0.5;

vessel_mask = false(size(vessel_pred));
for i_cap = 1:length(apex_measures.distal.width_at_apex)
    if isnan(apex_measures.distal.width_at_apex(i_cap))
        continue;
    end
    rows = apex_measures.distal.bounding_box(1,2,i_cap):apex_measures.distal.bounding_box(2,2,i_cap);
    cols = apex_measures.distal.bounding_box(1,1,i_cap):apex_measures.distal.bounding_box(3,1,i_cap);
    v_patch = vessel_pred(rows, cols);
    cx = apex_measures.distal.apex_xy(i_cap,1) - cols(1);
    cy = apex_measures.distal.apex_xy(i_cap,2) - rows(1);
    v_mask = bwselect(v_patch > v_thresh, cx, cy, 4);
    
    vessel_mask(rows,cols) = vessel_mask(rows,cols) | v_mask;
end
figure; 
subplot(2,1,1); imgray(vessel_pred);
subplot(2,1,2); imgray(vessel_mask);
%%
load('C:\isbe\nailfold\data\wellcome_study\results\auto_stats.mat');
figure; plot(auto_stats.total_scores(:), auto_stats.mean_median_flow(:), 'rx')
[r c] = find(auto_stats.total_scores > 40);

for i_im = 1:length(r)
    nailfold_mosaic = imread([sequence_names{r(i_im), c(i_im)} '\sequence_data\full_mosaic.png']);
    figure; imgray(nailfold_mosaic);
end
%%
flow_seqs = dir('C:\isbe\nailfold\data\wellcome_study\flow_results\035wellcome_2015_03_23_R2_12_24_58_*.mat');
num_seqs = length(flow_seqs);

for i_seq = 1:num_seqs
    frames = u_load(['N:\Nailfold Capillaroscopy\wellcome\flow_data\' flow_seqs(i_seq).name]);
    load(['C:\isbe\nailfold\data\wellcome_study\flow_results\' flow_seqs(i_seq).name]);
    
    mask = full(flow_results.flowConfidence{1}) < flow_prctiles(i_seq,3);
    flow = flow_results.flowPyramidEst{1};
    flow(mask) = complex(0,0);
    
    figure; 
    subplot(1,3,1); show_flow_as('rgb', flow);
    subplot(1,3,2); imgray(flow_results.flowConfidence{1});
    caxis(flow_prctiles(i_seq,[2 4]));
    subplot(1,3,3); imgray(mean(frames,3));   
       
end
%%
%%
flow_seqs = dir('C:\isbe\nailfold\data\wellcome_study\flow_results\*v01*.mat');
%%
num_seqs = length(flow_seqs);
flow_prctiles = zeros(num_seqs, 5);
confidence_prctiles = zeros(num_seqs, 5);
num_frames = zeros(num_seqs,1);

for i_seq = 1:num_seqs
    display(['processing ' num2str(i_seq)]);
    frames = u_load(['N:\Nailfold Capillaroscopy\wellcome\flow_data\' flow_seqs(i_seq).name]);
    load(['C:\isbe\nailfold\data\wellcome_study\flow_results\' flow_seqs(i_seq).name]);
    
    num_frames(i_seq) = size(frames,3);
    confidence_prctiles(i_seq,:) = prctile(flow_results.flowConfidence{1}(:), [0 25 50 75 100]);
    
    mask = full(flow_results.flowConfidence{1}) > confidence_prctiles(i_seq,3);
    flow_prctiles(i_seq,:) = prctile(abs(flow_results.flowPyramidEst{1}(mask)), [0 25 50 75 100]);
         
end
results_path = 'C:\isbe\nailfold\data\wellcome_study\results\flow_confidence_exp.mat';    
save(results_path, 'flow_prctiles', 'num_frames');
%%
hi_flow = find(flow_prctiles(:,3) > 1);
%%
flow_data_dir = 'N:\Nailfold Capillaroscopy\wellcome\flow_data\';
temp_frames_dir = 'C:\isbe\nailfold\data\wellcome_study\temp_frames\';
fast_videos_dir = 'C:\isbe\nailfold\data\wellcome_study\flow_videos\fast\';
create_folder(temp_frames_dir);
create_folder(fast_videos_dir);
%%
for i_seq = 11:length(hi_flow)
    seq_name = flow_seqs(hi_flow(i_seq)).name;
    frames = u_load([flow_data_dir seq_name]);
    num_frames = size(frames,3);
    g_lims = prctile(double(frames(:)), [1 99]);
    g_range = g_lims(2) - g_lims(1);

    delete([temp_frames_dir '*.bmp']);

    for i_fr = 1:num_frames;
        imwrite(uint8(255*(double(frames(:,:,i_fr))-g_lims(1))/g_range ),...
            [temp_frames_dir 'frame' zerostr(i_fr,4) '.bmp']);
    end
    cmd = ['ffmpeg -y -r 60 -i "' temp_frames_dir 'frame%04d.bmp" -c:v libx264 -preset slow -crf 18 -an "' fast_videos_dir seq_name '.mp4"'];
    system(cmd);
end
%%
lo_flow = find(adjusted_confidence > 6e5);
flow_videos_dir = 'C:\isbe\nailfold\data\wellcome_study\flow_videos\slow\';
create_folder(flow_videos_dir);
%
for i_seq = 1:length(lo_flow)
    seq_name = flow_seqs(lo_flow(i_seq)).name;
    frames = u_load([flow_data_dir seq_name]);
    num_frames = size(frames,3);
    g_lims = prctile(double(frames(:)), [1 99]);
    g_range = g_lims(2) - g_lims(1);

    delete([temp_frames_dir '*.bmp']);

    for i_fr = 1:num_frames;
        imwrite(uint8(255*(double(frames(:,:,i_fr))-g_lims(1))/g_range ),...
            [temp_frames_dir 'frame' zerostr(i_fr,4) '.bmp']);
    end
    cmd = ['ffmpeg -y -r 60 -i "' temp_frames_dir 'frame%04d.bmp" -c:v libx264 -preset slow -crf 18 -an "' flow_videos_dir seq_name '.mp4"'];
    system(cmd);
end
%%
for i_seq = 1:length(lo_flow)
    seq_name = flow_seqs(lo_flow(i_seq)).name;
    frames = u_load(['N:\Nailfold Capillaroscopy\wellcome\flow_data\' seq_name]);
    load(['C:\isbe\nailfold\data\wellcome_study\flow_results\' seq_name]);
    
    %mask = full(flow_results.flowConfidence{1}) < confidence_prctiles(lo_flow(i_seq),3);
    flow = flow_results.flowPyramidEst{1};
    %flow(mask) = complex(0,0);
    
    figure; 
    subplot(1,3,1); show_flow_as('rgb', flow);
    subplot(1,3,2); imgray(flow_results.flowConfidence{1});
    caxis(confidence_prctiles(lo_flow(i_seq),[2 4]));
    subplot(1,3,3); imgray(mean(frames,3));   
       
end
%%
base_dir = 'N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\';
match_flow_to_vessel_orientation(sequence_names2(1));
    
%%
flow_metrics_dir = 'C:\isbe\nailfold\data\wellcome_study\flow_metrics\';
flow_list = dir([flow_metrics_dir '*.mat']);
num_vessels = length(flow_list);
mean_errors = zeros(num_vessels,1);
mean_weighted_errors = zeros(num_vessels,1);
weighted_flow_rates = zeros(num_vessels,1);
total_vessel_probs = zeros(num_vessels,1);
mean_vessel_widths = zeros(num_vessels,1);
%
for i_ve = 1:num_vessels
    f = load([flow_metrics_dir flow_list(i_ve).name],...
        'mean_error', 'mean_weighted_error', 'weighted_flow_rate', 'total_vessel_prob', 'weighted_width');
    mean_errors(i_ve,1) = f.mean_error;
    mean_weighted_errors(i_ve,1) = f.mean_weighted_error;
    weighted_flow_rates(i_ve,1) = f.weighted_flow_rate;
    total_vessel_probs(i_ve,1) = f.total_vessel_prob;
    mean_vessel_widths(i_ve,1) = f.weighted_width;
end
%%
figure;
hist(mean_weighted_errors, 100);
title('Histogram of weighted flow errors');

figure;
hist(weighted_flow_rates, 100);
title('Histogram of weighted flow rates');

figure;
hist(total_vessel_probs, 100);
title('Histogram of total vessel probabilities');

figure;
subplot(1,2,1);
plot(mean_weighted_errors, weighted_flow_rates,'rx');
xlabel('Weighted flow error');
ylabel('Weighted flow rate');

subplot(1,2,2);
plot(mean_weighted_errors, total_vessel_probs,'rx');
xlabel('Weighted flow error');
ylabel('Total vessel probability');

figure;
plot(mean_vessel_widths, weighted_flow_rates,'rx');
xlabel('Vessel width');
ylabel('Weighted flow rate');

%%
flow_data_dir = 'N:\Nailfold Capillaroscopy\wellcome\flow_data\';
slow_idx = find(mean_weighted_errors < 0.5 & weighted_flow_rates < 0.08 & total_vessel_probs > 1e3);
flow_videos_dir = 'C:\isbe\nailfold\data\wellcome_study\flow_videos\slow\';
temp_frames_dir = 'C:\isbe\nailfold\data\wellcome_study\temp_frames\';
create_folder(flow_videos_dir);

for i_ve = 1:length(slow_idx)
    seq_name = flow_list(slow_idx(i_ve)).name;
    frames = u_load([flow_data_dir seq_name]);
    num_frames = size(frames,3);
    g_lims = prctile(double(frames(:)), [1 99]);
    g_range = g_lims(2) - g_lims(1);

    delete([temp_frames_dir '*.bmp']);

    for i_fr = 1:num_frames;
        imwrite(uint8(255*(double(frames(:,:,i_fr))-g_lims(1))/g_range ),...
            [temp_frames_dir 'frame' zerostr(i_fr,4) '.bmp']);
    end
    cmd = ['ffmpeg -y -r 60 -i "' temp_frames_dir 'frame%04d.bmp" -c:v libx264 -preset slow -crf 18 -an "' flow_videos_dir seq_name '.mp4"'];
    system(cmd);   
       
end
%%
fast_idx = find(mean_weighted_errors < 0.5 & weighted_flow_rates > 2 & total_vessel_probs > 2e3);
flow_videos_dir = 'C:\isbe\nailfold\data\wellcome_study\flow_videos\fast\';
temp_frames_dir = 'C:\isbe\nailfold\data\wellcome_study\temp_frames\';
create_folder(flow_videos_dir);

for i_ve = 1:10%length(fast_idx)
    seq_name = flow_list(fast_idx(i_ve)).name;
    frames = u_load([flow_data_dir seq_name]);
    num_frames = size(frames,3);
    g_lims = prctile(double(frames(:)), [1 99]);
    g_range = g_lims(2) - g_lims(1);

    delete([temp_frames_dir '*.bmp']);

    for i_fr = 1:num_frames;
        imwrite(uint8(255*(double(frames(:,:,i_fr))-g_lims(1))/g_range ),...
            [temp_frames_dir 'frame' zerostr(i_fr,4) '.bmp']);
    end
    cmd = ['ffmpeg -y -r 60 -i "' temp_frames_dir 'frame%04d.bmp" -c:v libx264 -preset slow -crf 18 -an "' flow_videos_dir seq_name '.mp4"'];
    system(cmd);   
       
end
%%
wide_idx = find(mean_weighted_errors < 0.5 & weighted_flow_rates > 0.1 & total_vessel_probs > 2e3 & mean_vessel_widths > 15 & weighted_flow_rates < 0.5);
flow_videos_dir = 'C:\isbe\nailfold\data\wellcome_study\flow_videos\wide\';
temp_frames_dir = 'C:\isbe\nailfold\data\wellcome_study\temp_frames\';
create_folder(flow_videos_dir);

for i_ve = 1:length(fast_idx)
    seq_name = flow_list(wide_idx(i_ve)).name;
    frames = u_load([flow_data_dir seq_name]);
    num_frames = size(frames,3);
    g_lims = prctile(double(frames(:)), [1 99]);
    g_range = g_lims(2) - g_lims(1);

    delete([temp_frames_dir '*.bmp']);

    for i_fr = 1:num_frames;
        imwrite(uint8(255*(double(frames(:,:,i_fr))-g_lims(1))/g_range ),...
            [temp_frames_dir 'frame' zerostr(i_fr,4) '.bmp']);
    end
    cmd = ['ffmpeg -y -r 60 -i "' temp_frames_dir 'frame%04d.bmp" -c:v libx264 -preset slow -crf 18 -an "' flow_videos_dir seq_name '.mp4"'];
    system(cmd);   
       
end
%%
valid_idx = find(mean_weighted_errors < 0.5 & total_vessel_probs > 2e3);
valid_flow = weighted_flow_rates(valid_idx);
min_flow = min(valid_flow);
max_flow = max(valid_flow);
num_steps = 20;
flow_steps = linspace(min_flow, max_flow, num_steps);

flow_videos_dir = 'C:\isbe\nailfold\data\wellcome_study\flow_videos\flow_steps\';
temp_frames_dir = 'C:\isbe\nailfold\data\wellcome_study\temp_frames\';
create_folder(flow_videos_dir);

for i_ve = 1:num_steps
    [~, idx] = min(abs(valid_flow - flow_steps(i_ve)));
    
    seq_name = flow_list(valid_idx(idx)).name;
    frames = u_load([flow_data_dir seq_name]);
    num_frames = size(frames,3);
    g_lims = prctile(double(frames(:)), [1 99]);
    g_range = g_lims(2) - g_lims(1);

    delete([temp_frames_dir '*.bmp']);

    for i_fr = 1:num_frames;
        imwrite(uint8(255*(double(frames(:,:,i_fr))-g_lims(1))/g_range ),...
            [temp_frames_dir 'frame' zerostr(i_fr,4) '.bmp']);
    end
    cmd = ['ffmpeg -y -r 30 -i "' temp_frames_dir 'frame%04d.bmp" -c:v libx264 -preset slow -crf 18 -an "' flow_videos_dir zerostr(i_ve,2) '_' seq_name '.mp4"'];
    system(cmd);   
       
end
%%
valid_idx = find(mean_weighted_errors < 0.5 & total_vessel_probs > 2e3 & weighted_flow_rates > 0.1);
valid_widths = mean_vessel_widths(valid_idx);
min_width = min(valid_widths);
max_width = max(valid_widths);
num_steps = 20;
width_steps = linspace(min_width, max_width, num_steps);

flow_videos_dir = 'C:\isbe\nailfold\data\wellcome_study\flow_videos\width_steps\';
temp_frames_dir = 'C:\isbe\nailfold\data\wellcome_study\temp_frames\';
create_folder(flow_videos_dir);

for i_ve = 1:num_steps
    [~, idx] = min(abs(valid_widths - width_steps(i_ve)));
    
    seq_name = flow_list(valid_idx(idx)).name;
    frames = u_load([flow_data_dir seq_name]);
    num_frames = size(frames,3);
    g_lims = prctile(double(frames(:)), [1 99]);
    g_range = g_lims(2) - g_lims(1);

    delete([temp_frames_dir '*.bmp']);

    for i_fr = 1:num_frames;
        imwrite(uint8(255*(double(frames(:,:,i_fr))-g_lims(1))/g_range ),...
            [temp_frames_dir 'frame' zerostr(i_fr,4) '.bmp']);
    end
    cmd = ['ffmpeg -y -r 30 -i "' temp_frames_dir 'frame%04d.bmp" -c:v libx264 -preset slow -crf 18 -an "' flow_videos_dir zerostr(i_ve,2) '_' seq_name '.mp4"'];
    system(cmd);   
       
end
%%
flow_videos_dir = 'C:\isbe\nailfold\data\wellcome_study\flow_videos\misc\'; 
seq_name = '113wellcome_2015_08_26_R3_15_11_52_s07_v02';
frames = load([flow_data_dir seq_name], 'frames_i');

figure;
subplot(1,2,1); imgray(mean(frames,3));
subplot(1,2,2); imgray(frames(:,:,1));
%%
rows = any(frames.frames_i(:,:,1),2);
cols = any(frames.frames_i(:,:,1),1);
cropped_frames = frames.frames_i(rows, cols, :);

[vessel_transforms] = ...
        register_tiles_features(cropped_frames, ...
                            'ref_type', 'mosaic',...
                            'theta_range', 0, ...
                            'offset_lim', 20, ...
                            'mosaic', mean(cropped_frames,3),...
                            'sigma', 6,...
                            'tile_masks', [],...
                            'debug', 1);
[frames_mosaic, ~, ~, cleaned_frames, edge_mask] = ...
        create_mosaic(cropped_frames, vessel_transforms);   
%    
rows = ~all(edge_mask,2);
cols = ~all(edge_mask);
cropped_frames = cleaned_frames(rows, cols, :);

g_lims = prctile(double(cropped_frames(:)), [1 99]);
g_range = g_lims(2) - g_lims(1);


[flow_results.flowPyramidEst, flow_results.flowConfidence] = ...
    estimate_flow_multilevel(255*(cropped_frames-g_lims(1))/g_range, [], [], 1:3);

%%
flow_results_dir = 'C:\isbe\nailfold\data\wellcome_study\flow_results\';

flow_results_old = u_load([flow_results_dir seq_name]);

figure;
subplot(1,2,1); show_flow_as('rgb', flow_results_old.flowPyramidEst{1});
subplot(1,2,2); show_flow_as('rgb', flow_results.flowPyramidEst{1});

figure;
subplot(1,2,1); imgray(abs(flow_results_old.flowPyramidEst{1})); colorbar;
subplot(1,2,2); imgray(abs(flow_results.flowPyramidEst{1})); colorbar;

figure;
subplot(1,2,1); imgray(flow_results_old.flowConfidence{1}(11:end-10, 11:end-10)); colorbar;
subplot(1,2,2); imgray(flow_results.flowConfidence{1}(11:end-10, 11:end-10)); colorbar;
%%
figure;
subplot(1,2,1); imgray(mean(frames,3));
subplot(1,2,2); imgray(mean_clean_machine);
caxis([min(mean_clean_machine(~edge_mask)) max(mean_clean_machine(~edge_mask))]);
%%
g_lims = prctile(double(f.cropped_frames(:)), [1 99]);
g_range = g_lims(2) - g_lims(1);
if rem(size(f.cropped_frames,1),2)
    f.cropped_frames(end,:,:) = [];
end
if rem(size(f.cropped_frames,2),2)
    f.cropped_frames(:,end,:) = [];
end

delete([temp_frames_dir '*.bmp']);

for i_fr = 1:size(f.cropped_frames,3);
    imwrite(uint8(255*(f.cropped_frames(:,:,i_fr)-g_lims(1))/g_range ),...
        [temp_frames_dir 'frame' zerostr(i_fr,4) '.bmp']);
end
cmd = ['ffmpeg -y -r 30 -i "' temp_frames_dir 'frame%04d.bmp" -c:v libx264 -preset slow -crf 18 -an "' flow_videos_dir 'cleaned_' seq_name '.mp4"'];
system(cmd);
    