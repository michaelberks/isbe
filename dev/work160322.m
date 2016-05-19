%%
data_dir = 'C:\isbe\nailfold\data\rsa_study\set12g\';
contour_dir = [data_dir 'vessel_contours\'];
flow_map_dir = [data_dir 'flow_maps\'];
flow_data_dir = [data_dir 'syn_flow_data\'];
flow_results_dir = [data_dir 'syn_flow_results\'];
flow_metrics_dir = [data_dir 'syn_flow_metrics\'];
flow_params_dir = [data_dir 'syn_flow_params\'];
flow_videos_dir = [data_dir 'syn_flow_videos\'];
%%

flow_list = dir([flow_results_dir 'normalapex0905*.mat']);
flow_names = {flow_list(:).name}';

noise = zeros(10,1);
max_flows = zeros(10,1);
flow_results = cell(10,1);

for i_rpt = 1:10
    
    fr = load([flow_results_dir flow_names{i_rpt}]);
    flow_results{i_rpt} = fr.flow_results.flowPyramidEst{1};
    
    fm = load([flow_metrics_dir flow_names{i_rpt}]);
    
    
    fp = load([flow_params_dir flow_names{i_rpt}]);
    
    noise(i_rpt) = fp.speckle_noise;
    max_flows(i_rpt) = fp.max_flow;
    
    fd = load([flow_data_dir flow_names{i_rpt}], 'frames');
    
    if i_rpt == 1
        flow_frames = zeros(size(fd.frames,1), size(fd.frames,2), size(fd.frames,3), 10);
        flow_metrics = fm.flow_metrics;
    else
        flow_metrics(i_rpt) = fm.flow_metrics;
    end
    flow_frames(:,:,:,i_rpt) = fd.frames;
    
    if ~exist([flow_videos_dir flow_names{i_rpt} '.mp4'], 'file')
        make_flow_video(flow_frames(:,:,:,i_rpt), [flow_videos_dir flow_names{i_rpt} '.mp4'], 30);
    end
    
    if 0 == 1
        figure; 
        subplot(1,3,1); imgray(mean(flow_frames(:,:,:,i_rpt),3));
        subplot(1,3,2); imgray(flow_metrics(i_rpt).vessel_ori);
        subplot(1,3,3); imgray(complex2rgb(flow_results{i_rpt}));
    end
end
%%
for i_rpt = 1:10
    figure; 
    subplot(1,3,1); imgray(mean(flow_frames(:,:,:,i_rpt),3));
    subplot(1,3,2); imgray(flow_metrics(i_rpt).vessel_ori);
    subplot(1,3,3); imgray(complex2rgb(flow_results{i_rpt}, [], 8));
end
%%
figure; 
subplot(1,2,1); axis equal; hold all;
plot([0 7], [0 7], 'k--');
for i_rpt = 1:10
    text(true_weighted_flow_rates(i_rpt), flow_metrics(i_rpt).weighted_flow_rate,...
        num2str(i_rpt));
end

subplot(1,2,2); axis equal; hold all;
plot([0 7], [0 7], 'k--');
for i_rpt = 1:10
    text(true_weighted_flow_rates(i_rpt), prctile_flow_rates(:,i_rpt,6),...
        num2str(i_rpt));
end
%%
new_flows = [4.00 6.00 8.00 10.0 12.0 8.00 8.00 8.00 8.00 4.00 6.00 10.0 12.0 2.00 2.00 14.0 14.0 12.0 12.0 12.0 2.00 4.00 6.00 8.00 10.0 12.0];
new_noise = [0.03 0.03 0.03 0.03 0.03 0.01 0.02 0.04 0.05 0.05 0.05 0.05 0.05 0.03 0.05 0.03 0.05 0.01 0.02 0.04 0.07 0.07 0.07 0.07 0.07 0.07];
vessel_name = 'normalapex0905_vessel';
num_rpts = length(new_flows);
new_results = cell(num_rpts,1);
%%
for i_rpt = 1:num_rpts

    %Randomly sample max flow and noise properties
    vessel_contrast = 1;%2*rand;
    speckle_noise = new_noise(i_rpt);
    max_flow = new_flows(i_rpt);

    %Make synthetic frames and save
    [frames] = make_synthetic_flow_frames(...
        vessel_name, max_flow,...
        'vessel_contrast', vessel_contrast,...
        'speckle_noise', speckle_noise,...
        'flow_map_dir', flow_map_dir);
    
    %Compute flow estimates  
    fr = estimate_flow_multilevel(255*frames, [], [], 1:3);
    new_results{i_rpt} = fr{1};
    
    %Compute flow metrics
    [fm] =  compute_vessel_flow_rf(...
        'frames', frames,...
        'flow_results', new_results{i_rpt},...
        'rescale_factor', 4/3,...
        'patch_contrast_scale', 40,...
        'model_root', 'C:\isbe\nailfold\models\');
    if i_rpt == 1
        new_frames = zeros(size(frames,1), size(frames,2), size(frames,3), num_rpts);
        new_metrics = fm;
    else
        new_metrics(i_rpt) = fm;
    end
    new_frames(:,:,:,i_rpt) = frames;
    
    figure; 
    subplot(1,3,1); imgray(mean(frames,3));
    title(['MFR = ' num2str(new_metrics(i_rpt).weighted_flow_rate, 3)]);
    subplot(1,3,2); imgray(new_metrics(i_rpt).vessel_pred);
    plot(new_metrics(i_rpt).apex_xy(:,1), new_metrics(i_rpt).apex_xy(:,2),'rx');
    title(['MWE = ' num2str(new_metrics(i_rpt).mean_weighted_error, 3)]);
    subplot(1,3,3); imgray(complex2rgb(new_results{i_rpt}, [], [], [], 1));
    title(['Flow ratio = ' num2str(new_metrics(i_rpt).vessel_flow / new_metrics(i_rpt).background_flow, 3)]);
        
end
%%
mask = imerode(new_metrics(6).frame_apex_mask, strel('disk', 5));
figure; imgray(mask);
%%
load([flow_map_dir 'normalapex0905_vessel_fm.mat']);
flow_map = max(flowmap, [], 3);
mean_gt_flow = mean(abs(flow_map(mask_thin)));
%%
mean_gt_flows = zeros(num_rpts,1);
mean_estimated_flows = zeros(num_rpts,1);
prctile_estimated_flows = zeros(num_rpts,6);
for i_rpt = 1:num_rpts
    mean_estimated_flows(i_rpt,1) = mean(abs(new_results{i_rpt}(mask_thin)));
    mean_gt_flows(i_rpt,1) = mean_gt_flow * new_flows(i_rpt);
    prctile_estimated_flows(i_rpt,:) = prctile(abs(new_results{i_rpt}(mask_thin)), 50:10:100);
end

%%
figure; hold all;  axis equal;
plot(mean_gt_flows([14 1:5 16]), mean_estimated_flows([14 1:5 16]), '-o');
plot(mean_gt_flows([15 10 11 9 12 13 17]), mean_estimated_flows([15 10 11 9 12 13 17]), '-o');
plot(mean_gt_flows(21:26), mean_estimated_flows(21:26), '-o');
xlabel('Mean GT flow');
ylabel('Mean estimated flow');
legend({'Noise = 0.03', 'Noise = 0.05', 'Noise = 0.07'});
%%
figure;
for i_pct = 1:6
    
    subplot(2,3,i_pct); hold all; axis equal;
    plot(mean_gt_flows([14 1:5 16]), prctile_estimated_flows([14 1:5 16], i_pct), '-o');
    plot(mean_gt_flows([15 10 11 9 12 13 17]), prctile_estimated_flows([15 10 11 9 12 13 17], i_pct), '-o');
    plot(mean_gt_flows(21:26), prctile_estimated_flows(21:26, i_pct), '-o');
    plot([0 max(mean_gt_flows)], [0 max(mean_gt_flows)], 'k--'); 
    axis([0 max(mean_gt_flows) 0 max(mean_gt_flows)]);
    xlabel('Mean GT flow');
    ylabel('Mean estimated flow');
    legend({'Noise = 0.03', 'Noise = 0.05', 'Noise = 0.07'}, 'location', 'northwest');
end
%%
figure; hold all;
plot(new_noise([6 7 3 8 9 24]), mean_estimated_flows([6 7 3 8 9 24 ]), '-o');
plot(new_noise([18 19 5 20 13 26]), mean_estimated_flows([18 19 5 20 13 26]), '-o');
xlabel('GT Noise');
ylabel('mean estimated flow');
legend({'Max flow = 8.0', 'Max flow = 12.0'});
%%
figure; hold all; 
plot(prctile(abs(new_results{5}(mask_thin)), 0:100), 0:100, '-o');
plot(prctile(abs(new_results{13}(mask_thin)), 0:100), 0:100, '-o');
plot(prctile(abs(new_results{26}(mask_thin)), 0:100), 0:100, '-o');
%%
new_results4 = cell(num_rpts,1);
for i_rpt = 1:num_rpts

    %Compute flow estimates  
    fr = estimate_flow_multilevel(255*new_frames(:,:,:,i_rpt), [], [], 1:4);
    new_results4{i_rpt} = fr{1};
    
    %Compute flow metrics
    [fm] =  compute_vessel_flow_rf(...
        'frames', new_frames(:,:,:,i_rpt),...
        'flow_results', new_results4{i_rpt},...
        'rescale_factor', 4/3,...
        'patch_contrast_scale', 40,...
        'model_root', 'C:\isbe\nailfold\models\');
    if i_rpt == 1
        new_metrics4 = fm;
    else
        new_metrics4(i_rpt) = fm;
    end
    
    figure; 
    subplot(1,3,1); imgray(mean(frames,3));
    title(['MFR = ' num2str(new_metrics4(i_rpt).weighted_flow_rate, 3)]);
    subplot(1,3,2); imgray(new_metrics4(i_rpt).vessel_pred);
    plot(new_metrics4(i_rpt).apex_xy(:,1), new_metrics4(i_rpt).apex_xy(:,2),'rx');
    title(['MWE = ' num2str(new_metrics4(i_rpt).mean_weighted_error, 3)]);
    subplot(1,3,3); imgray(complex2rgb(new_results{i_rpt}, [], [], [], 1));
    title(['Flow ratio = ' num2str(new_metrics4(i_rpt).vessel_flow / new_metrics4(i_rpt).background_flow, 3)]);
        
end
%%
mean_estimated_flows = zeros(num_rpts,1);
prctile_estimated_flows = zeros(num_rpts,6);
for i_rpt = 1:num_rpts
    mean_estimated_flows(i_rpt,1) = mean(abs(new_results4{i_rpt}(mask_thin)));
    prctile_estimated_flows(i_rpt,:) = prctile(abs(new_results4{i_rpt}(mask_thin)), 50:10:100);
end
%%
figure; hold all; axis equal;
plot(mean_gt_flows([14 1:5 16]), mean_estimated_flows([14 1:5 16]), '-o');
plot(mean_gt_flows([15 10 11 9 12 13 17]), mean_estimated_flows([15 10 11 9 12 13 17]), '-o');
plot(mean_gt_flows(21:26), mean_estimated_flows(21:26), '-o');
plot([0 max(mean_gt_flows)], [0 max(mean_gt_flows)], 'k--'); 
axis([0 max(mean_gt_flows) 0 max(mean_gt_flows)]);
xlabel('Mean GT flow');
ylabel('Mean estimated flow');
legend({'Noise = 0.03', 'Noise = 0.05', 'Noise = 0.07'});
%%
figure;
for i_pct = 1:6
    
    subplot(2,3,i_pct); hold all; axis equal;
    plot(mean_gt_flows([14 1:5 16]), prctile_estimated_flows([14 1:5 16], i_pct), '-o');
    plot(mean_gt_flows([15 10 11 9 12 13 17]), prctile_estimated_flows([15 10 11 9 12 13 17], i_pct), '-o');
    plot(mean_gt_flows(21:26), prctile_estimated_flows(21:26, i_pct), '-o');
    plot([0 max(mean_gt_flows)], [0 max(mean_gt_flows)], 'k--'); 
    axis([0 max(mean_gt_flows) 0 max(mean_gt_flows)]);
    xlabel('Mean GT flow');
    ylabel('Mean estimated flow');
    legend({'Noise = 0.03', 'Noise = 0.05', 'Noise = 0.07'}, 'location', 'northwest');
end
%%
figure; hold all;
plot(new_noise([6 7 3 8 9 24]), mean_estimated_flows([6 7 3 8 9 24 ]), '-o');
plot(new_noise([18 19 5 20 13 26]), mean_estimated_flows([18 19 5 20 13 26]), '-o');

xlabel('GT Noise');
ylabel('mean estimated flow');
legend({'Max flow = 8.0', 'Max flow = 12.0'});






    
    