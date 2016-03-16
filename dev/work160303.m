flow_data_list = dir('N:\Nailfold Capillaroscopy\wellcome\flow_data\*.mat');
flow_results_list = dir('N:\Nailfold Capillaroscopy\wellcome\flow_results\*.mat');

flow_data_names = {flow_data_list(:).name}';
flow_results_names = {flow_results_list(:).name}';

missing = setdiff(flow_data_names, flow_results_names);
%%
num_missing = length(missing);
cropped_present = false(num_missing,1);

warning('off', 'MATLAB:load:variableNotFound');
for i_ve = 1:length(missing)
    fd = load(['N:\Nailfold Capillaroscopy\wellcome\flow_data\' missing{i_ve}], 'cropped_frames');
    cropped_present(i_ve) = isfield(fd, 'cropped_frames'); clear fd;
end

display(sum(cropped_present));
%%
frames_to_process = false(num_missing,1);
for i_ve = 1:num_missing
    if cropped_present(i_ve); continue; end;
    
    display(['processing vessel ' num2str(i_ve)]);
    
    load(['N:\Nailfold Capillaroscopy\wellcome\flow_data\' missing{i_ve}], 'frames_i');
    
    rows = any(frames_i(:,:,1),2);
    cols = any(frames_i(:,:,1),1);
    cropped_frames = frames_i(rows, cols, :);

    if sum(rows) < 16 || sum(cols) < 16
        save(['N:\Nailfold Capillaroscopy\wellcome\flow_data\' missing{i_ve}], ...
            'cropped_frames', '-append');
        continue;
    end
    
    [vessel_transforms] = ...
        register_tiles_features(cropped_frames, ...
                            'ref_type', 'mosaic',...
                            'region', 'all',...
                            'theta_range', 0, ...
                            'offset_lim', 20, ...
                            'mosaic', mean(cropped_frames,3),...
                            'sigma', 6,...
                            'tile_masks', [],...
                            'debug', 0);
    [~, ~, ~, cleaned_frames, edge_mask] = ...
        create_mosaic(cropped_frames, vessel_transforms);

    rows = ~all(edge_mask,2);
    cols = ~all(edge_mask);
    cropped_frames = cleaned_frames(rows, cols, :);

    save(['N:\Nailfold Capillaroscopy\wellcome\flow_data\' missing{i_ve}], ...
        'vessel_transforms', 'cropped_frames', 'edge_mask', '-append');
    
    if all(size(cropped_frames)>=16)
        frames_to_process(i_ve) = 1;
    end
end
%%
flow_results_dir = 'N:\Nailfold Capillaroscopy\wellcome\flow_results\';
processed = 0;
for i_ve = 1:num_missing
    if frames_to_process(i_ve)
        display(['processing vessel ' num2str(processed)]);
        frames = load(['N:\Nailfold Capillaroscopy\wellcome\flow_data\' missing{i_ve}], 'cropped_frames');
        
        g_lims = prctile(double(frames.cropped_frames(:)), [1 99]);
        g_range = g_lims(2) - g_lims(1);

        flow_results = [];
        [flow_results.flowPyramidEst, flow_results.flowConfidence] = ...
            estimate_flow_multilevel(255*(frames.cropped_frames-g_lims(1))/g_range, [], [], 1:3);

        flow_results.g_lims = g_lims;

        save([flow_results_dir flow_list(i_ves).name], 'flow_results');
        processed = processed + 1;
    end
end
    
%%
g_prob = gaussian_filters_1d(2);
g_prob = g_prob / sum(g_prob);

for i_ve = 1:num_vessels
    flow_metrics = u_load([flow_metrics_dir flow_list(i_ve).name]);
    
    
    smooth_ves = conv2(g_prob', g_prob, double(flow_metrics.vessel_pred)/100, 'same');
    smooth_ori = conv2(g_prob', g_prob, rgb2complex(flow_metrics.vessel_ori,[],1), 'same');

    %make a vessel mask connected to the apex
    frame_vessels_mask = smooth_ves > 0.25;
    frame_vessels_mask = bwselect(frame_vessels_mask, flow_metrics.apex_xy(1), flow_metrics.apex_xy(2));
    thin_mask = bwmorph(frame_vessels_mask, 'thin', 'inf');   
    flow_metrics.shape_score = abs(mean(smooth_ori(thin_mask)));
    
    
end
%%

