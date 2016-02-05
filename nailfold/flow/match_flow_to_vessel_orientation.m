function [] = match_flow_to_vessel_orientation(sequence_names)

flow_data_dir = 'N:\Nailfold Capillaroscopy\wellcome\flow_data\';
flow_results_dir = 'C:\isbe\nailfold\data\wellcome_study\flow_results\';
flow_metrics_dir = 'C:\isbe\nailfold\data\wellcome_study\flow_metrics\';
capillary_data_dir = 'C:\isbe\nailfold\data\wellcome_study\capillary_data\';

create_folder(flow_metrics_dir);

connect_thresh = 0.2;
num_ori_bins = 36;
do_plot = 0;
%
for i_seq = 1:numel(sequence_names)
    
    if isempty(sequence_names{i_seq})
        continue;
    end
    
    seq_dir = [sequence_names{i_seq} '\'];%sequence_dirs{i_seq};
    seq_name = seq_dir;
    dividers = seq_name == '\';
    pos = find(dividers, 4, 'last');
    seq_name(dividers) = '_';
    seq_name = seq_name(pos(1)+1:end-1);
    
    display(['processing sequence ' num2str(i_seq) ': ' seq_name]);
    
    vessel_ori = imread([seq_dir '\capillary_data\vessels_o_pred.png']);
    vessel_pred = imread([seq_dir '\capillary_data\vessels_v_pred.png']);
    [sequence_data] = read_processed_sequence_from([seq_dir '\sequence_data\sequence_data.dat']);
    
    if isfield(sequence_data, 'final_mosaic_size')
        vessel_ori = imresize(vessel_ori, fliplr(sequence_data.final_mosaic_size));
        vessel_pred = imresize(vessel_pred, fliplr(sequence_data.final_mosaic_size));
    else
        nailfold_im = imread([seq_dir '\sequence_data\full_mosaic.png']);
        vessel_ori = imresize(vessel_ori, size(nailfold_im));
        vessel_pred = imresize(vessel_pred, size(nailfold_im));
    end
        
    vessel_data = load([capillary_data_dir seq_name '_capillary_data.mat'],...
        'resize_factor', 'apex_measures');
    
    for i_seg = 1:sequence_data.num_segments
    
        flow_list = dir([flow_results_dir seq_name '*s' zerostr(i_seg,2) '*.mat']);
        if isempty(flow_list)
            continue;
        end
        segment_mask = imread([seq_dir '\registered' zerostr(i_seg-1,2) '\segment_mask.png']);
       
        if do_plot && i_seq <= 1
            rf_list = dir([seq_dir '\registered' zerostr(i_seg-1,2) '\frame*.png']);
            frame0 = imread([seq_dir '\registered' zerostr(i_seg-1,2) '\' rf_list(1).name]);
            figure; 
            a1 = subplot(1,2,1); imgray(segment_mask);
            a2 = subplot(1,2,2); imgray(frame0);
        end
        
        for i_flow = 1:length(flow_list)
           
            f = load([flow_data_dir flow_list(i_flow).name]);
            flow_results = u_load([flow_results_dir flow_list(i_flow).name]);
            weighted_width = vessel_data.apex_measures.distal.mean_weighted_width(f.apex_idx);
            
            frame_ori = vessel_ori(...
                (f.y_min:f.y_max)+floor(f.frames_offset(2)),...
                (f.x_min:f.x_max)+floor(f.frames_offset(1)),:);
            frame_ori = rgb2complex(frame_ori, [], 1, [], 0);
            frame_pred = double(vessel_pred(...
                (f.y_min:f.y_max)+floor(f.frames_offset(2)),...
                (f.x_min:f.x_max)+floor(f.frames_offset(1)),:)) / 100;
            
            segment_frame_mask = segment_mask(f.y_min:f.y_max, f.x_min:f.x_max);
            
            if do_plot && i_seq <= 1
                plot(a1, [f.x_min f.x_max f.x_max f.x_min f.x_min], [f.y_min f.y_min f.y_max f.y_max f.y_min]);
                plot(a2, [f.x_min f.x_max f.x_max f.x_min f.x_min], [f.y_min f.y_min f.y_max f.y_max f.y_min]);
            end

            ax = vessel_data.apex_measures.distal.apex_xy(f.apex_idx,1)/vessel_data.resize_factor ...
                - (floor(f.frames_offset(1)) + f.x_min);
            ay = vessel_data.apex_measures.distal.apex_xy(f.apex_idx,2)/vessel_data.resize_factor ...
                - (floor(f.frames_offset(2)) + f.y_min);

            frame_mask = bwselect(frame_pred > connect_thresh, ax, ay, 4) &...
                segment_frame_mask;
            frame_pred_masked = frame_pred;
            frame_pred_masked(~frame_mask) = 0;
            total_vessel_prob = sum(frame_pred_masked(:));

            flow_angle = exp(2i*angle(flow_results.flowPyramidEst{1}));
            flow_angle_diffs = frame_pred_masked .* frame_ori .* flow_angle;
            flow_angle_weights = abs(flow_angle_diffs);
            flow_angle_dirs = angle(flow_angle_diffs);
    
            weighted_flow_rate = sum(frame_pred_masked(:).*abs(flow_results.flowPyramidEst{1}(:))) /...
                total_vessel_prob;
    
            mean_error = mean(abs(flow_angle_dirs(:)));
            mean_weighted_error = sum(flow_angle_weights(:) .* abs(flow_angle_dirs(:))) /...
                sum(flow_angle_weights(:));
    
            g_theta_idx = flow_angle_dirs/(2*pi) + 0.5;
            g_theta_idx = g_theta_idx * num_ori_bins;
    
            %Convert theta vals in into integer indices - there may still be some 0s,
            %which should be moved to the final bin
            g_theta_idx_c = ceil(g_theta_idx);
            g_theta_w = g_theta_idx_c - g_theta_idx;
    
            g_theta_idx_c(~g_theta_idx_c) = num_ori_bins;   
            g_theta_idx_f = g_theta_idx_c - 1;
            g_theta_idx_f(~g_theta_idx_f) = num_ori_bins;
    
            flow_errors_hist = full(...
                sparse(1, g_theta_idx_c, (1-g_theta_w).*flow_angle_weights, 1, num_ori_bins) + ...
                sparse(1, g_theta_idx_f, g_theta_w.*flow_angle_weights, 1, num_ori_bins));
    
            if do_plot %&& i_seq <= 10 && 2
                figure;
                subplot(2,3,1); imgray(complex2rgb(frame_ori, [], [], [], 1));
                subplot(2,3,2); imgray(complex2rgb(conj(flow_results.flowPyramidEst{1}.^2), [], [], [], 1));
                title(['MFR = ' num2str(weighted_flow_rate, 3)]);
                subplot(2,3,3); imgray(complex2rgb(flow_angle_diffs, [], [], [], 1));
                subplot(2,3,4); bar(flow_errors_hist);
                title(['ME = ' num2str(mean_error,3) ' MWE = ' num2str(mean_weighted_error, 3)]);
                subplot(2,3,5); imgray(frame_pred); plot(ax, ay, 'rx');
                subplot(2,3,6); imgray(frame_pred_masked);
            end
    
            save([flow_metrics_dir flow_list(i_flow).name],...
                'mean_error', 'mean_weighted_error', 'weighted_flow_rate',...
                'flow_errors_hist', 'total_vessel_prob', 'weighted_width');
            
        end
    end
end