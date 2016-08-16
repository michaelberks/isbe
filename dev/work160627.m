root_dir = 'N:\Nailfold Capillaroscopy\Wellcome\';
seq_name = '035wellcome_2015_03_23_R2_12_24_58';
seg_name = '_s07';

vessels = dir([root_dir 'flow_results\' seq_name seg_name '*.mat']);

frames_dir = 'N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\035wellcome\2015_03_23\R2_12_24_58\registered06\';

frame = imread([frames_dir 'frame03838.png']);
%%
flow_map = zeros(size(frame));
flow_mask = false(size(frame));
for i_ve = 1:5
    flow_data = load([root_dir 'flow_data\' vessels(i_ve).name],...
        'x_*_cropped', 'y_*_cropped');
    
    flow_results = load([root_dir 'flow_results\' vessels(i_ve).name]);
    
    flow_metrics = load([root_dir 'flow_metrics\' vessels(i_ve).name]);
    vessel_mask = flow_metrics.flow_metrics.vessel_pred > 50;
    vessel_mask = bwselect(vessel_mask,...
        flow_metrics.flow_metrics.apex_xy(1), flow_metrics.flow_metrics.apex_xy(2));
    
    vessel_sz = [...
        flow_data.y_max_cropped - flow_data.y_min_cropped + 1
        flow_data.x_max_cropped - flow_data.x_min_cropped + 1];
    
    vessel_r = 1:vessel_sz(1);
    vessel_c = 1:vessel_sz(2);
    
    flow_map(...
        flow_data.y_min_cropped - 1 + vessel_r,...
        flow_data.x_min_cropped - 1 + vessel_c) = max(flow_map(...
        flow_data.y_min_cropped - 1 + vessel_r,...
        flow_data.x_min_cropped - 1 + vessel_c),...
        flow_results.flow_results.flowPyramidEst{1}(...
        vessel_r, vessel_c));
    
    flow_mask(...
        flow_data.y_min_cropped - 1 + vessel_r,...
        flow_data.x_min_cropped - 1 + vessel_c) = flow_mask(...
        flow_data.y_min_cropped - 1 + vessel_r,...
        flow_data.x_min_cropped - 1 + vessel_c) | ...
        vessel_mask(vessel_r, vessel_c);
end

flow_mask_rgb = cat(3, flow_mask, flow_mask, flow_mask);
flow_rgb = complex2rgb(flow_map,[],[],[],1);
flow_rgb2 = complex2rgb(flow_map);
frame_rgb = double(cat(3, frame, frame, frame)) / 255;

figure;
subplot(1,2,1); imgray(frame);
subplot(1,2,2); imgray(flow_rgb);
%%
frames = dir([frames_dir 'frame*.png']);
num_frames = length(frames);

steps = linspace(0.25, 1, num_frames/2);
steps = [fliplr(steps) steps];
%%
for i_fr = 1:20:num_frames;
    prct = steps(i_fr);
    frame_i = imread([frames_dir frames(i_fr).name]);    
    overlay_rgb = double(cat(3,frame_i,frame_i,frame_i))/255;
    overlay_rgb(flow_mask_rgb) = ...
        prct*overlay_rgb(flow_mask_rgb) + (1-prct)*flow_rgb(flow_mask_rgb);
    overlay_rgb = uint8(255*overlay_rgb);
    figure; imgray(overlay_rgb);
end
%%
for i_fr = 1:20:num_frames;
    prct = steps(i_fr);
    frame_i = imread([frames_dir frames(i_fr).name]);    
    overlay_rgb = double(cat(3,frame_i,frame_i,frame_i))/255;
    overlay_rgb = ...
        prct*overlay_rgb + (1-prct)*flow_rgb;
    overlay_rgb = uint8(255*overlay_rgb);
    figure; imgray(overlay_rgb);
end
%%
for i_fr = 1:20:num_frames;
    prct = steps(i_fr);
    frame_i = imread([frames_dir frames(i_fr).name]);    
    overlay_rgb = double(cat(3,frame_i,frame_i,frame_i))/255;
    overlay_rgb = ...
        prct*overlay_rgb + (1-prct)*flow_rgb2;
    overlay_rgb = uint8(255*overlay_rgb);
    figure; imgray(overlay_rgb);
end
%%
flow_gif = 'C:\isbe\flow.gif';

for i_fr = 1:num_frames
    frame_i = imread([frames_dir frames(i_fr).name]);
    frame_i = cat(3,frame_i,frame_i,frame_i);
    
    [im, map] = rgb2ind(frame_i, 256);
    if i_fr == 1
        imwrite(im, map, flow_gif, 'gif', 'WriteMode', 'overwrite', 'Loop', Inf, 'DelayTime', 0.033);
    else
        imwrite(im, map, flow_gif, 'gif', 'WriteMode', 'append', 'DelayTime', 0.033);
    end
end
%%
steps = linspace(0.4, 1, num_frames/2);
steps = [fliplr(steps) steps];
flow_gif = 'C:\isbe\epsrc_talk_June2016\flow_overlay.gif';
for i_fr = 1:2:num_frames
    prct = steps(i_fr);
    frame_i = imread([frames_dir frames(i_fr).name]);    
    overlay_rgb = double(cat(3,frame_i,frame_i,frame_i))/255;
    overlay_rgb = ...
        prct*overlay_rgb + (1-prct)*flow_rgb2;
    overlay_rgb = uint8(255*overlay_rgb);
    
    [im, map] = rgb2ind(overlay_rgb, 256);
    if i_fr == 1
        imwrite(im, map, flow_gif, 'gif', 'WriteMode', 'overwrite', 'Loop', Inf, 'DelayTime', 0.066);
    else
        imwrite(im, map, flow_gif, 'gif', 'WriteMode', 'append', 'DelayTime', 0.066);
    end
end

