function [] = make_flow_comparison_video(sequence_name, segment_idx, vessel_idx, speed_factor, varargin)
%MAKE_FLOW_COMPARISON_VIDEO *Insert a one line summary here*
%   [] = make_flow_comparison_video(varargin)
%
% MAKE_FLOW_COMPARISON_VIDEO uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Optional Arguments:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 15-Feb-2016
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'raw_data_dir',         'N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\',... 
    'capillary_data_dir',   'C:\isbe\nailfold\data\wellcome_study\capillary_data\',...
    'flow_data_dir',        'N:\Nailfold Capillaroscopy\wellcome\flow_data\',...
    'flow_results_dir',     'N:\Nailfold Capillaroscopy\wellcome\flow_results\',...
    'video_dir',            'C:\isbe\nailfold\data\wellcome_study\flow_videos\flow_comparisons\',...
    'prob_sigma',           2,...
    'ori_sigma',            2,...
    'width_sigma',          2,...
    'flow_sigma',           2,...
    'num_cells',            8,...
    'scale_cell_width',     0.4,...
    'fixed_cell_width',     8, ...
    'solid_cells',          true,...
    'adjacent_videos',      false,...
    'temp_frames_dir',      [],... if empty, uses the system temporary dir
    'frame_rate',           15,...
    'plot', 0);
clear varargin;

%Load and reshape vessel prediction images
vessel_pred = imread([args.raw_data_dir sequence_name '\capillary_data\vessels_v_pred.png']);
vessel_ori = imread([args.raw_data_dir sequence_name '\capillary_data\vessels_o_pred.png']);
vessel_width = imread([args.raw_data_dir sequence_name '\capillary_data\vessels_w_pred.png']);

seq_name = sequence_name;
seq_name(sequence_name == '\') = '_';
load([args.capillary_data_dir seq_name '_capillary_data.mat']);

[sequence_data] = read_processed_sequence_from([args.raw_data_dir sequence_name '\sequence_data\sequence_data.dat']);
mosaic_w = sequence_data.final_mosaic_size(1);
mosaic_h = sequence_data.final_mosaic_size(2);

vessel_pred = imresize(vessel_pred, [mosaic_h mosaic_w]);
vessel_ori = imresize(vessel_ori, [mosaic_h mosaic_w]);
vessel_width = imresize(vessel_width, [mosaic_h mosaic_w]);

if ~isnumeric(segment_idx) && strcmpi(segment_idx, 'all')
   segment_idx = zeros(0,1); 
   vessel_idx = zeros(0,1);
   for i_cap = 1:length(apex_measures.distal.flow_names)
       for i_f = 1:length(apex_measures.distal.flow_names{i_cap})
            segment_idx(end+1,1) = str2double(apex_measures.distal.flow_names{i_cap}{i_f}(2:3)); %#ok
            vessel_idx(end+1,1) = str2double(apex_measures.distal.flow_names{i_cap}{i_f}(6:7)); %#ok
       end
   end
    
end
%%
for i_vessel = 1:length(segment_idx)
    %Load in flow data and predictions
    vessel_name = [seq_name '_s' zerostr(segment_idx(i_vessel),2) '_v' zerostr(vessel_idx(i_vessel),2)];
    load([args.flow_data_dir vessel_name '.mat'],...
        'vessel_transforms', 'apex_idx', 'cropped_frames',...
        'x_min', 'y_min', 'edge_mask', 'frames_offset');
    load([args.flow_results_dir vessel_name '.mat']);
    
    %Smooth the flow field
    g_flow = gaussian_filters_1d(args.flow_sigma);
    g_flow = g_flow / sum(g_flow);    
    smooth_flow = conv2(g_flow', g_flow, flow_results.flowPyramidEst{1}, 'same');
    
    %Get bounding box for cropped frames in mosaic, and extract patches of
    %vessel predictions
    [~, ~, dx, dy] = ...
        mosaic_limits(size(flow_results.flowPyramidEst{1}), vessel_transforms);

    x_start = floor(x_min + find(any(~edge_mask,1),1) - (dx + 1));
    y_start = floor(y_min + find(any(~edge_mask,2),1) - (dy + 1));

    [flow_h flow_w] = size(flow_results.flowPyramidEst{1});
    cropped_rows = y_start+(1:flow_h);
    cropped_cols = x_start+(1:flow_w);

    %Now compute information about these frames - first get vessel and
    %orientation patches for the correctly translated bundign box
    pred_rows = cropped_rows+floor(frames_offset(2));
    pred_cols = cropped_cols+floor(frames_offset(1));
    valid_pred_rows = pred_rows >= 1 & pred_rows <= size(vessel_pred,1);
    valid_pred_cols = pred_cols >= 1 & pred_cols <= size(vessel_pred,2);

    frame_ori = zeros(flow_h, flow_w, 3);
    frame_ori(valid_pred_rows, valid_pred_cols, :) = vessel_ori(...
        pred_rows(valid_pred_rows),...
        pred_cols(valid_pred_cols),:);
    frame_ori = rgb2complex(frame_ori, [], 1, [], 0);
    frame_pred = zeros(flow_h, flow_w);
    frame_pred(valid_pred_rows, valid_pred_cols) = double(vessel_pred(...
        pred_rows(valid_pred_rows),...
        pred_cols(valid_pred_cols))) / 100;

    frame_width = zeros(flow_h, flow_w);
    frame_width(valid_pred_rows, valid_pred_cols) = double(vessel_width(...
        pred_rows(valid_pred_rows),...
        pred_cols(valid_pred_cols)));

    %Make a frame mask from the vessel predictions, selected pixels
    %connected to the vessel apex
    ax = apex_measures.distal.apex_xy(apex_idx,1)/resize_factor ...
        - (floor(frames_offset(1)) + x_start);
    ay = apex_measures.distal.apex_xy(apex_idx,2)/resize_factor ...
        - (floor(frames_offset(2)) + y_start);

    connect_thresh = 0.5;
    frame_vessels_mask = frame_pred > connect_thresh;
    frame_apex_mask = bwselect(frame_vessels_mask, ax, ay, 4);
    
    if isempty(frame_apex_mask)
        continue;
    end

    % %Find the candidates for the centre of the vessel
    % g_prob = gaussian_filters_1d(args.prob_sigma);
    % g_prob = g_prob / sum(g_prob);
    % smooth_pred = conv2(g_prob', g_prob, frame_pred, 'same');
    % vessel_centre = mb_non_maximal_supp(smooth_pred, angle(frame_ori)/2);
    % vessel_centre(~frame_apex_mask) = 0;
    % vessel_centre_mask = vessel_centre > 0;
    % [centre_y centre_x] = find(vessel_centre);
    % 
    % fxy = smooth_flow(vessel_centre_mask);
    % 
    % figure; show_flow_as('rgb', smooth_flow); hold all;
    % quiver(centre_x, centre_y, 10*real(fxy), 10*imag(fxy));
    %
    % figure;
    % subplot(2,3,1); imgray(frame_pred);
    % plot(ax, ay, 'rx');
    % subplot(2,3,2); imgray(frame_apex_mask);
    % plot(ax, ay, 'rx');
    % subplot(2,3,3); imgray(complex2rgb(frame_ori, [], [], [], 1));
    % plot(ax, ay, 'rx');
    % subplot(2,3,4); imgray(frame_width);
    % plot(ax, ay, 'rx');
    % subplot(2,3,5); imgray(vessel_centre);
    % plot(ax, ay, 'rx');
    % subplot(2,3,6); show_flow_as('rgb', flow_results.flowPyramidEst{1}); hold on;
    % plot(ax, ay, 'rx');
    %
    %
    %Get vessel path by moving away from the apex in the direction of the flow
    %field
    %Set max vessel path to 2 x the diagonal of the patch
    max_dist = 2*sqrt(flow_h^2 + flow_w^2);
    max_steps = ceil(max_dist);%flow_h*flow_w;
    art_path = zeros(max_steps,3);
    ven_path = zeros(max_steps,3);

    %Go forwards from apex (i.e venous path)
    ven_path(1,:) = [ax ay 0];
    fxy_curr = complex(0,0);
    for i_step = 2:max_steps
        ix = ven_path(i_step-1,1);
        iy = ven_path(i_step-1,2);
        dist = ven_path(i_step-1,3);

        if ix > 1 && ix < flow_w && iy > 1 && iy < flow_h && dist < max_dist;
            
            %Get new flow direction, but only take this path if some angle
            %consistency is maintained, otherwise assume we're crossing the
            %path of the other limb (which has taken precendence of this
            %limb in the flow predictions)
            fxy_new = smooth_flow(round(iy), round(ix));
            if i_step == 2 || abs(angle(conj(fxy_new)*fxy_curr)) < pi/3;
                fxy_curr = fxy_new;
            else
                display('Crossing vessel path!');
            end
            fm = abs(fxy_curr) + 1e-6;
            ven_path(i_step,1) = ix + real(fxy_curr) / fm;
            ven_path(i_step,2) = iy + imag(fxy_curr) / fm;
            ven_path(i_step,3) = dist + 1;
        else
            ven_path(i_step:end,:) = [];
            break;
        end
    end

    %Go backwards from apex (i.e arterial path)
    max_dist = -max_dist;
    art_path(1,:) = [ax ay 0];
    fxy_curr = complex(0,0);
    for i_step = 2:max_steps
        ix = art_path(i_step-1,1);
        iy = art_path(i_step-1,2);
        dist = art_path(i_step-1,3);

        if ix > 1 && ix < flow_w && iy > 1 && iy < flow_h && dist > max_dist;
            fxy_new = smooth_flow(round(iy), round(ix));
            if i_step == 2 || abs(angle(conj(fxy_new)*fxy_curr)) < pi/3;
                fxy_curr = fxy_new;
            else
                display('Crossing vessel path!');
            end
            fm = abs(fxy_curr) + 1e-6;
            art_path(i_step,1) = ix - real(fxy_curr) / fm;
            art_path(i_step,2) = iy - imag(fxy_curr) / fm;
            art_path(i_step,3) = dist - 1;
        else
            art_path(i_step:end,:) = [];
            break;
        end
    end

    %Concatenate the two paths (flipping the arterial) to complete the vessel
    %path, then use the complete path to make a distance map for the vessel
    vessel_path = [flipud(art_path); ven_path(2:end,:)];

    if args.plot
        
        [vessel_y vessel_x] = find(frame_apex_mask);
        all_vd = griddata(vessel_path(:,1), vessel_path(:,2), vessel_path(:,3), ...
            vessel_x, vessel_y, 'nearest'); %#ok
    
        dist_map = zeros(flow_h, flow_w);
        dist_map(frame_apex_mask) = all_vd;
    
        figure; 

        subplot(1,2,1); show_flow_as('rgb', flow_results.flowPyramidEst{1}); hold on;
        plot(vessel_path(:,1), vessel_path(:,2), 'k--');

        subplot(1,2,2); imgray(dist_map);
        plot(vessel_path(:,1), vessel_path(:,2), 'r--');
    end

    %--------------------------------------------------------------------------
    %Now we can make our video
    %--------------------------------------------------------------------------

    %Scale and resize the frames to be video compatible (and match the flow
    %field)
    frames = cropped_frames;
    [fh fw num_frames] = size(frames);
    g_lims = prctile(double(frames(:)), [1 99]);
    g_range = g_lims(2) - g_lims(1);
    frames = 255*(double(frames)-g_lims(1))/g_range;
    frames = mb_pad(frames, [flow_h - fh, flow_w - fw, 0], 'replicate', 'post');

    %Set the cell radius - fixed or based on the frame width?
    if args.scale_cell_width
        cell_radius2 = (frame_width(round(ax), round(ay)) * args.scale_cell_width * 0.5).^2;
        cell_radius = sqrt(cell_radius2);
    else
        cell_radius = args.fixed_cell_width;
        cell_radius2 = cell_radius^2;
    end

    %Set the first point, the select equally spaced further points on the
    %vessel
    first_pt = ceil(cell_radius) + ceil(cell_radius*rand);

    cell_idx = round(linspace(first_pt, size(vessel_path,1), args.num_cells+1));
    cell_idx(end) = [];
    cell_pts = repmat(vessel_path(cell_idx,1:2), 1, 1, 2);

    % OLD METHOD: select starting points across the vessel width, has annoying
    % habit of projecting cells outside the vessel path
    % cell_dists = vessel_path(cell_idx,3);
    % cell_pts = zeros(args.num_cells,2,2);
    % for i_cell = 1:args.num_cells
    %     candidate_cells = all_vd == cell_dists(i_cell);
    %     
    %     vxi = vessel_x(candidate_cells);
    %     vyi = vessel_y(candidate_cells);
    %     r_idx = ceil(length(vxi)*rand);
    %     cell_pts(i_cell,1,:) = vxi(r_idx);
    %     cell_pts(i_cell,2,:) = vyi(r_idx);
    % end

    %Set up points matrices to compute cell masks
    xx = repmat(1:flow_w, flow_h, 1);
    yy = repmat((1:flow_h)', 1, flow_w);

    %Make a tempoary directory to store frame in
    if isempty(args.temp_frames_dir)
        temp_frames_dir = tempname;
    else
        temp_frames_dir = args.temp_frames_dir;
    end
    mkdir(temp_frames_dir);

    if rand > 0.5
        scaling = [1.0 speed_factor];
    else
        scaling = [speed_factor 1.0];
    end

    %Make the video frames
    for i_fr = 1:num_frames;

        if args.adjacent_videos
            vid_frame = zeros(flow_h, 2*flow_w + 16, 3, 'uint8'); 
        else
            vid_frame = zeros(flow_h, flow_w, 3, 'uint8'); 
        end

        for i_dim = 1:2
            %Get the raw frame
            frame_i = frames(:,:,i_fr);

            %For each cell
            for i_cell = 1:args.num_cells

                %Make cell mask at current position
                ix = cell_pts(i_cell,1,i_dim);
                iy = cell_pts(i_cell,2,i_dim);
                xy2 = (xx-ix).^2 + (yy-iy).^2;
                cell_mask = xy2 < cell_radius2;
                if ~args.solid_cells
                    cell_mask = cell_mask & xy2 > (cell_radius2/4);
                end
                frame_i(cell_mask) = 255;    

                %If we're still in a valid position, move to next position
                %based on flow estimate
                if ix > 1 && ix < flow_w && iy > 1 && iy < flow_h;
                    fxyi = smooth_flow(round(iy), round(ix));
                    cell_pts(i_cell,1,i_dim) = ix + scaling(i_dim)*real(fxyi);
                    cell_pts(i_cell,2,i_dim) = iy + scaling(i_dim)*imag(fxyi);
                else
                    %Otherwise start at the beginning
                    cell_pts(i_cell,:,i_dim) = vessel_path(first_pt,1:2);
                end
            end

            %Put frames either on top or by one another
            if args.adjacent_videos
                vid_frame(:, (1:flow_w) + (i_dim-1)*(flow_w + 16), 1) = uint8(frame_i);
                vid_frame(:, (1:flow_w) + (i_dim-1)*(flow_w + 16), 2) = uint8(frames(:,:,i_fr));
                vid_frame(:, (1:flow_w) + (i_dim-1)*(flow_w + 16), 3) = uint8(frames(:,:,i_fr));
            else
                vid_frame(:, :, i_dim) = uint8(frame_i);
            end

        end
        if ~args.adjacent_videos
            vid_frame(:, :, 3) = uint8(frames(:,:,i_fr));
        end

        %Write the frame to the temporary directory
        imwrite(vid_frame,...
            [temp_frames_dir '\frame' zerostr(i_fr,4) '.bmp']);
    end

    %Make path to save video
    create_folder(args.video_dir);
    if args.adjacent_videos
        video_path = [args.video_dir vessel_name ...
            '_f1_' num2str(scaling(1),3) '_f2_' num2str(scaling(2),3) '_adjacent.mp4'];
    else
        video_path = [args.video_dir vessel_name ...
            '_f1_' num2str(scaling(1),3) '_f2_' num2str(scaling(2),3) '_overlay.mp4'];
    end
    cmd = ['ffmpeg -y -r ' num2str(args.frame_rate) ' -i "' temp_frames_dir '\frame%04d.bmp" -c:v libx264 -preset slow -crf 18 -an "' video_path '"'];
    system(cmd);

    delete([temp_frames_dir '\*']);
    rmdir(temp_frames_dir);
end