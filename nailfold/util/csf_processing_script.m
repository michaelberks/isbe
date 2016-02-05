%Get list of sequences to copy
study_dir = 'N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\';
rds_dir = 'rds-scratch/nailfold/';
subject_dirs = [
    dir([study_dir '0*'])
    dir([study_dir '1*'])];

num_subs = length(subject_dirs);
recruitment_dates = zeros(num_subs,1);
seq_nums = zeros(num_subs, 2);

sequence_dirs = cell(0,1);
sequence_rds_dirs = cell(0,1);
sequence_names = cell(0,1);
seq_count = 1;
for i_sub = 1:num_subs
    session_dirs = dir([study_dir subject_dirs(i_sub).name '\2015*']);
    
    seq_nums(i_sub,1) = seq_count;
    if i_sub > 1
        seq_nums(i_sub-1,2) = seq_count-1;
    end
    for i_ses = 1:length(session_dirs)
        sequence_dirs_i = [
            dir([study_dir subject_dirs(i_sub).name '\' session_dirs(i_ses).name '\L*'])
            dir([study_dir subject_dirs(i_sub).name '\' session_dirs(i_ses).name '\R*'])];
        
        for i_seq = 1:length(sequence_dirs_i)
            sequence_dirs{end+1} = ...
                [study_dir subject_dirs(i_sub).name '\' session_dirs(i_ses).name '\' sequence_dirs_i(i_seq).name '\']; %#ok
            sequence_rds_dirs{end+1} = ...
                [rds_dir subject_dirs(i_sub).name '/' session_dirs(i_ses).name '/' sequence_dirs_i(i_seq).name '/']; %#ok
            sequence_names{end+1} = sequence_dirs_i(i_seq).name; %#ok
            
            seq_count = seq_count + 1;
        end
    end
end

num_seqs = length(sequence_dirs);
%%
%Make batches of copy scripts
copy_dir = 'C:\isbe\nailfold\data\wellcome_study\copy_scripts\';
create_folder(copy_dir);

for i_batch = 1:12
    seqs = (1:10) + (i_batch-1)*10;
    copy_to_file_name = [copy_dir 'copy_to_rds' zerostr(i_batch,2) '.sh'];
    fid_to = fopen(copy_to_file_name, 'wt');
    
    copy_from_file_name = [copy_dir 'copy_from_rds' zerostr(i_batch,2) '.sh'];
    fid_from = fopen(copy_from_file_name, 'wt');
    
    for i_seq = 1:10
        
        if seqs(i_seq) > 112; break; end
        copy_cmd = ['rsync -ra --update --exclude-from="rds-scratch/ignore_file.txt" ' ...
            'mhsrfs01/epi-musculo/Nailfold\ Capillaroscopy/camera_capture/wellcome_nailfold_study/' ...
            zerostr(seqs(i_seq),3) 'wellcome rds-scratch/nailfold/'];
        %display(copy_cmd);
        fprintf(fid_to, '%s\n', copy_cmd);
        fprintf(fid_to, 'echo %s copy complete $(date) \n', [zerostr(seqs(i_seq),3) 'wellcome']);
        
        copy_cmd = ['rsync -ra --update --exclude-from="rds-scratch/ignore_file_from.txt" ' ...
            'rds-scratch/nailfold/' zerostr(seqs(i_seq),3) 'wellcome/ ' ...
            'mhsrfs01/epi-musculo/Nailfold\ Capillaroscopy/camera_capture/wellcome_nailfold_study/' zerostr(seqs(i_seq),3) 'wellcome/'];
        %display(copy_cmd);
        fprintf(fid_from, '%s\n', copy_cmd);
        fprintf(fid_from, 'echo %s copy complete $(date) \n', [zerostr(seqs(i_seq),3) 'wellcome']);
        
    end
    fclose(fid_to);
    fclose(fid_from);
end
%%
%make_list sequence parameter file inputs
fid = fopen('C:\isbe\nailfold\data\wellcome_study\sequences_to_process.txt', 'wt');
for i_seq = 1:num_seqs
    fprintf(fid, '-i %ssequence_frames_data.dat -p rds-scratch/models/csf_capillary_detection_w3_parameters.txt\n', sequence_rds_dirs{i_seq});
end
fclose(fid);
%%
%%
%make_list sequence parameter file inputs
fid = fopen('C:\isbe\nailfold\data\wellcome_study\sequences_to_copy.txt', 'wt');
for i_seq = 1:num_seqs
    seq_dir = sequence_dirs{i_seq};
    if exist([seq_dir 'capillary_data\apex_candidates.txt'], 'file')
        dividers = seq_dir == '\';
        seq_dir(dividers) = '/';
        fprintf(fid, '%s\n', seq_dir);
    end
end
fclose(fid);
%%
%Check which sequences failed to process, and if so, at what stage:
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
for i_sub = 35
    for i_digit = 1:10
        seq_dir = [sequence_names{i_sub, i_digit} '\sequence_data\'];
        cap_dir = [sequence_names{i_sub, i_digit} '\capillary_data\'];
        if ~exist([cap_dir 'apex_candidates.txt'], 'file')
            display([cap_dir ' full mosaic missing']);
            continue;
        end
        vessel_pred = imread([cap_dir 'vessels_v_pred.png']);
        ac = load([cap_dir 'apex_candidates.txt']);

        if isempty(ac)
            figure;
            subplot(2,1,1); 
            imgray(nailfold_im); caxis([gmin gmax]);
            title(seq_dir);
            subplot(2,1,2);
            imgray(vessel_pred);
            continue;
        end
        selected_distal = ac(:,11)>0;
        selected_non_distal = ac(:,12)>0;
        dy = ac(:,2) - ac(:,8);
        dxy = sortrows([ac(:,1) dy]);

        nailfold_im = imread([seq_dir 'full_mosaic.png']);
        count_map = imread([seq_dir 'count_map.png']);
        mask = count_map > 0;
        mask([1 end],:) = 0;
        mask(:, [1 end]) = 0;
        mask = imerode(mask, strel('disk', 2));

        gmin = min(nailfold_im(mask));
        gmax = max(nailfold_im(mask));
        figure;
        subplot(2,1,1); 
        imgray(nailfold_im); caxis([gmin gmax]);
        title(seq_dir);
        subplot(2,1,2);
        imgray(vessel_pred);

        plot(ac(selected_distal,1), ac(selected_distal,2), 'rx');
        plot(ac(selected_non_distal,1), ac(selected_non_distal,2), 'b+');
        plot(dxy(:,1), dxy(:,2), 'g--');

    end
end
%%
capillary_data_dir = 'C:\isbe\nailfold\data\wellcome_study\capillary_data\';
create_folder(capillary_data_dir);

apex_args.prob_sigma = 2;
apex_args.ori_sigma = 0;
apex_args.width_sigma = 2;
apex_args.num_c_pts = 20;
apex_args.max_dist = 120;
apex_args.border_sz = 16;
apex_args.num_ori_bins = 36;
apex_args.connect_thresh = 0.5;
apex_args.base_width = 20;
apex_args.width_predictor = [];
apex_args.plot = 0;

do_capillary_type = [true true];
capillary_type = {'distal', 'nondistal'};
resize_factor = 0.6220;

%
for i_seq = 347%1:numel(sequence_names)
    
    apex_measures.distal = [];
    apex_measures.nondistal = [];
    
    seq_dir = [sequence_names2{i_seq} '\'];%sequence_dirs{i_seq};
    seq_name = seq_dir;
    dividers = seq_name == '\';
    pos = find(dividers, 4, 'last');
    seq_name(dividers) = '_';
    seq_name = seq_name(pos(1)+1:end-1);
    
    if ~exist([seq_dir 'capillary_data\apex_candidates.txt'], 'file')
        display(['Missing sequence ' num2str(i_seq) ': ' seq_name]);
        continue;
    end
    
    apex_data = load([seq_dir 'capillary_data\apex_candidates.txt']);
    
    if isempty(apex_data)
        display(['saving empty sequence ' num2str(i_seq) ': ' seq_name]);
        save([capillary_data_dir seq_name '_capillary_data.mat'],...
            'resize_factor', 'apex_data', 'apex_measures');
        continue;
    end
    display(['processing sequence ' num2str(i_seq) ': ' seq_name]);
        
    %vessel_im = imread([seq_dir 'sequence_data\full_mosaic.png']);
    vessel_prob = imread([seq_dir 'capillary_data\vessels_v_pred.png']);
    vessel_prob = double(vessel_prob)/100;
    vessel_ori = imread([seq_dir 'capillary_data\vessels_o_pred.png']);
    vessel_ori = rgb2complex(vessel_ori, [], 1, [], 0);
    vessel_width = imread([seq_dir 'capillary_data\vessels_w_pred.png']);
    vessel_width = double(vessel_width);   
    
    candidate_xy = apex_data(:,1:2);
    selected_distal = apex_data(:,11) > 0;
    selected_non_distal = apex_data(:,12) > 0;
    
    for i_type = 1:2
    
        if ~do_capillary_type(i_type); continue; end
        
        if i_type == 1
            selected_idx = selected_distal;
        else
            selected_idx = selected_non_distal;
        end
    
        [measures_struc] = extract_apex_measures_newest(...
            vessel_prob, vessel_ori, vessel_width, vessel_prob, candidate_xy(selected_idx,:),... 
            apex_args);
        
        %Copying the fields in this way means we don't overwrite existing
        %data
        fnames = fieldnames(measures_struc);
        for i_f = 1:length(fnames)
            apex_measures.(capillary_type{i_type}).(fnames{i_f}) = measures_struc.(fnames{i_f});
        end
        
        apex_measures.(capillary_type{i_type}).candidate_scores = apex_data(selected_idx,7);
        apex_measures.(capillary_type{i_type}).candidate_displacements = apex_data(selected_idx,8);
        apex_measures.(capillary_type{i_type}).candidate_class_probs = apex_data(selected_idx,10);
        
    end

    %vessel_predictions = cat(3, vessel_prob, vessel_ori, vessel_width);
    save([capillary_data_dir seq_name '_capillary_data.mat'],...
        'resize_factor', 'apex_data', 'apex_measures');%vessel_predictions
end
%%
capillary_data_dir = 'C:\isbe\nailfold\data\wellcome_study\capillary_data\';
flow_data_dir = 'N:\Nailfold Capillaroscopy\Wellcome\flow_data.old\';
create_folder(flow_data_dir);

load('C:\isbe\nailfold\data\wellcome_study\sequence_names.mat');
sequence_names2 = sequence_names';
%
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
    display(['processing sequence ' num2str(i_seq) ': ' seq_name]);
    
    vessel_data = load([capillary_data_dir seq_name '_capillary_data.mat'],...
        'resize_factor', 'apex_measures');
    
    if ~isempty(vessel_data.apex_measures.distal)
        save_path = [flow_data_dir seq_name];
        save_mosaic_flow_patches(seq_dir, vessel_data, save_path, 'plot', 0);
    end
end
%%
flow_data_dir = 'N:\Nailfold Capillaroscopy\Wellcome\flow_data.old\';
flow_results_dir = 'N:\Nailfold Capillaroscopy\Wellcome\flow_results\';
flow_mosaic_dir = 'N:\Nailfold Capillaroscopy\Wellcome\flow_mosaics\';
capillary_data_dir = 'C:\isbe\nailfold\data\wellcome_study\capillary_data\';
flow_metrics_dir = 'C:\isbe\nailfold\data\wellcome_study\flow_metrics\';
create_folder(flow_mosaic_dir);

load('C:\isbe\nailfold\data\wellcome_study\sequence_names.mat');
sequence_names2 = sequence_names';
%%
for i_seq = 1%72:numel(sequence_names)
    
    if isempty(sequence_names2{i_seq})
        continue;
    end
    
    seq_dir = [sequence_names2{i_seq} '\'];%sequence_dirs{i_seq};
    seq_name = seq_dir;
    dividers = seq_name == '\';
    pos = find(dividers, 4, 'last');
    seq_name(dividers) = '_';
    seq_name = seq_name(pos(1)+1:end-1);
    display(['processing sequence ' num2str(i_seq) ': ' seq_name]);   
    
    save_path = [flow_mosaic_dir seq_name '_flow_mosaic.mat'];
    [mosaic_flow] = ...
         compute_mosaic_flow_c(seq_dir, capillary_data_dir, flow_data_dir, flow_results_dir,...
         flow_metrics_dir, seq_name, 'plot', 2);
     save(save_path, 'mosaic_flow');
    
end
%%
capillary_data_dir = 'C:\isbe\nailfold\data\wellcome_study\capillary_data\';
flow_mosaic_dir = 'C:\isbe\nailfold\data\wellcome_study\flow_mosaics\';
flow_mask_dir = 'C:\isbe\nailfold\data\wellcome_study\flow_masks\';
create_folder(flow_mask_dir);
%%
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
    display(['processing sequence ' num2str(i_seq) ': ' seq_name]);
    
    vessel_data = load([capillary_data_dir seq_name '_capillary_data.mat'],...
        'resize_factor', 'apex_measures');
    
    if ~isempty(vessel_data.apex_measures.distal)
        mosaic_flow = u_load([flow_mosaic_dir seq_name '_flow_mosaic.mat']);
         [mosaic_flow_mask, vessel_data] = ...
             compute_mosaic_flow_mask(seq_dir, vessel_data, mosaic_flow, 'plot', 0);
         save([flow_mask_dir seq_name '_mask.mat'], 'mosaic_mask');
         apex_measures = vessel_data.apex_measures;
         save([capillary_data_dir seq_name '_capillary_data.mat'], 'apex_measures', '-append');
    end
end
%%
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
    display(['processing sequence ' num2str(i_seq) ': ' seq_name]);
    
    vessel_data = load([capillary_data_dir seq_name '_capillary_data.mat'],...
        'resize_factor', 'apex_measures');
    
    if ~isempty(vessel_data.apex_measures.distal)
        mosaic_flow = u_load([flow_mosaic_dir seq_name '_flow_mosaic.mat']);
        mosaic_flow_mask = u_load([flow_mask_dir seq_name '_mask.mat']);
         [~, vessel_data] = ...
             compute_mosaic_flow_mask(mosaic_flow_mask, vessel_data, mosaic_flow, 'plot', 0);
         
         apex_measures = vessel_data.apex_measures;
         save([capillary_data_dir seq_name '_capillary_data.mat'], 'apex_measures', '-append');
         save([flow_mask_dir seq_name '_mask.mat'], 'mosaic_flow_mask');
    end
end
%%
load('C:\isbe\nailfold\data\wellcome_study\sequence_names.mat');
load('C:\isbe\nailfold\data\wellcome_study\subject_summary.mat');
%%
[auto_stats] = analyse_wellcome_apex_measures(...
    'auto_stats',       [],...
    'subject_summary',  subject_summary,...
    'sequence_names',   sequence_names,...
    'study_dir',        'N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\',...
    'capillary_dir',    'C:\isbe\nailfold\data\wellcome_study\capillary_data\',...
    'selected_features', [],...
    'do_xls',           0,...
    'do_auto_stats',    1,...
    'do_image_plots',   0, ...
    'do_people_plots',  0,...
    'fig_dir',          [],...
    'save_dir',         'C:\isbe\nailfold\data\wellcome_study\results\',...
    'um_per_pix',       1.6077);
%%
selected_features = {...
    'num_distal_vessels',...
    'num_giant_vessels',...
    'num_enlarged_vessels',...
    'vessel_density',...
    'mean_mean_width',...
    'max_mean_width',...
    'mean_connected_orientation_dispersion',...
    'dispersion_connected_orientation',...
    'mean_mean_flow',...
    'mean_median_flow'};
    
feature_display_names = {...
    'Num distal caps',...
    'Num giant caps',...
    'Num enlarged caps',...
    'Capillary density',...
    'Mean width',...
    'Max width',...
    'Shape score',...
    'Derangement score',...
    'Mean flow',...
    'Median flow'};
%%
analyse_wellcome_apex_measures(...
    'auto_stats',       auto_stats,...
    'subject_summary',  subject_summary,...
    'sequence_names',   sequence_names,...
    'study_dir',        'N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\',...
    'capillary_dir',    'C:\isbe\nailfold\data\wellcome_study\capillary_data\',...
    'selected_features', selected_features,...
    'feature_display_names', feature_display_names,...
    'do_xls',           0,...
    'do_auto_stats',    0,...
    'do_image_plots',   1, ...
    'do_people_plots',  0,...
    'fig_dir',          'C:\isbe\nailfold\data\wellcome_study\results\figs\',...
    'save_dir',         'C:\isbe\nailfold\data\wellcome_study\results\',...
    'um_per_pix',       1.6077);
%%
analyse_wellcome_apex_measures(...
    'auto_stats',       auto_stats,...
    'subject_summary',  subject_summary,...
    'sequence_names',   sequence_names,...
    'study_dir',        'N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\',...
    'capillary_dir',    'C:\isbe\nailfold\data\wellcome_study\capillary_data\',...
    'selected_features', selected_features,...
    'feature_display_names', feature_display_names,...
    'do_xls',           0,...
    'do_auto_stats',    0,...
    'do_image_plots',   0, ...
    'do_people_plots',  1,...
    'fig_dir',          'C:\isbe\nailfold\data\wellcome_study\results\figs\',...
    'save_dir',         'C:\isbe\nailfold\data\wellcome_study\results\',...
    'um_per_pix',       1.6077);
%%
selected_features = {...
    'num_distal_vessels',...
    'num_giant_vessels',...
    'num_enlarged_vessels',...
    'vessel_density',...
    'mean_mean_width',...
    'max_mean_width',...
    'mean_connected_orientation_dispersion',...
    'dispersion_connected_orientation'};
    
feature_display_names = {...
    'Num distal caps',...
    'Num giant caps',...
    'Num enlarged caps',...
    'Capillary density',...
    'Mean width',...
    'Max width',...
    'Shape score',...
    'Derangement score'};

analyse_wellcome_apex_measures(...
    'auto_stats',       auto_stats,...
    'subject_summary',  subject_summary,...
    'sequence_names',   sequence_names,...
    'study_dir',        'N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\',...
    'capillary_dir',    'C:\isbe\nailfold\data\wellcome_study\capillary_data\',...
    'selected_features', selected_features,...
    'feature_display_names', feature_display_names,...
    'do_xls',           1,...
    'do_auto_stats',    0,...
    'do_image_plots',   0, ...
    'do_people_plots',  0,...
    'fig_dir',          [],...
    'save_dir',         'C:\isbe\nailfold\data\wellcome_study\results\',...
    'um_per_pix',       1.6077,...
    'xls_wide_format',  1,...
    'xls_filename_subjects',     'subject_auto_stats.xls');
%%
selected_features = {...
    'vessel_density',...
    'mean_mean_width',...
    'mean_connected_orientation_dispersion',...
    'dispersion_connected_orientation',...
    'mean_median_flow',...
    'mean_scores'};

feature_display_names = {...
    'Capillary density',...
    'Mean width',...
    'Shape score',...
    'Derangement score',...
    'Median flow',...
    'Mean score'};

num_features = length(selected_features);
all_stats = zeros(numel(sequence_names), num_features);

for i_f = 1:num_features
    all_stats(:,i_f) = auto_stats.(selected_features{i_f})(:);
end
discard_ims = any(isnan(all_stats),2);
all_stats(discard_ims,:) = [];
%%
for i_f1 = 1:num_features
    for i_f2 = 1:i_f1-1
        figure;
        plot(all_stats(:,i_f1), all_stats(:,i_f2), 'rx');
        xlabel(feature_display_names{i_f1});
        ylabel(feature_display_names{i_f2});
        [rho pval] = corr(all_stats(:,i_f1), all_stats(:,i_f2));
        title(['R = ' num2str(rho,2) ', p = ' num2str(pval, 2)]);
    end
end
        
%%
%%
selected_sequences = sequence_names(auto_stats.num_distal_vessels > 50);

for i_seq = 1:length(selected_sequences)

    seq_dir = [selected_sequences{i_seq} '\sequence_data\'];
    cap_dir = [selected_sequences{i_seq} '\capillary_data\'];
    if ~exist([cap_dir 'apex_candidates.txt'], 'file')
        display([cap_dir ' full mosaic missing']);
        continue;
    end
    vessel_pred = imread([cap_dir 'vessels_v_pred.png']);
    ac = load([cap_dir 'apex_candidates.txt']);

    if isempty(ac)
        figure;
        subplot(2,1,1); 
        imgray(nailfold_im); caxis([gmin gmax]);
        title(sequence_names{i_seq});
        subplot(2,1,2);
        imgray(vessel_pred);
        continue;
    end
    selected_distal = ac(:,11)>0;
    selected_non_distal = ac(:,12)>0;
    dy = ac(:,2) - ac(:,8);
    dxy = sortrows([ac(:,1) dy]);

    nailfold_im = imread([seq_dir 'full_mosaic.png']);
    count_map = imread([seq_dir 'count_map.png']);
    mask = count_map > 0;
    mask([1 end],:) = 0;
    mask(:, [1 end]) = 0;
    mask = imerode(mask, strel('disk', 2));

    gmin = min(nailfold_im(mask));
    gmax = max(nailfold_im(mask));
    figure;
    subplot(2,1,1); 
    imgray(nailfold_im); caxis([gmin gmax]);
    title(selected_sequences{i_seq});
    subplot(2,1,2);
    imgray(vessel_pred);

    plot(ac(selected_distal,1), ac(selected_distal,2), 'rx');
    plot(ac(selected_non_distal,1), ac(selected_non_distal,2), 'b+');
    plot(dxy(:,1), dxy(:,2), 'g--');

end
%%
fid = fopen('C:\isbe\nailfold\data\wellcome_study\sequences_to_correct.txt', 'wt');
for i_seq = 1:num_seqs
    display(['Processing sequence ' num2str(i_seq)]);
    seq_dir = sequence_dirs{i_seq};
    if ~exist([seq_dir 'sequence_frames_data.dat'], 'file') || i_seq==348
        continue;
    end
    [sequence_data] = read_processed_sequence_from(...
        [seq_dir 'sequence_data\sequence_data.dat']);
    write_processed_sequence_to(...
        [seq_dir 'sequence_data\sequence_data_w.dat'], sequence_data);
    
    dividers = seq_dir == '\';
    seq_dir(dividers) = '/';
    fprintf(fid, '%s\n', [seq_dir 'sequence_data/sequence_data_w.dat']);
    
end
fclose(fid);
        

