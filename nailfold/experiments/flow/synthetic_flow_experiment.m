%--------------------------
% Synthetic flow experiment
%--------------------------
%%
data_dir = 'C:\isbe\nailfold\data\rsa_study\set12g\';
contour_dir = [data_dir 'vessel_contours\'];
flow_map_dir = [data_dir 'flow_maps\'];
flow_data_dir = [data_dir 'syn_flow_data\'];
flow_results_dir = [data_dir 'syn_flow_results\'];
flow_metrics_dir = [data_dir 'syn_flow_metrics\'];
flow_params_dir = [data_dir 'syn_flow_params\'];

create_folder(flow_map_dir);

contour_list = dir([contour_dir '*_vc.mat']);
num_vessels = length(contour_list);
%%
for i_ve = 1:num_vessels
    display(['Making flow map ' num2str(i_ve)]);
    load([contour_dir contour_list(i_ve).name],...
        'vessel_centre', 'inner_edge', 'outer_edge', 'apex_idx');
    apex_idx = sort(apex_idx);
    
    %Make a flowmap for the vessel shape
    [flowmap, mask, vessel_centre, widths] = ...
        create_flowmap_profile(vessel_centre, inner_edge, outer_edge, apex_idx, []);
    
    save([flow_map_dir contour_list(i_ve).name(1:end-6) 'fm.mat'],...
        'flowmap', 'mask', 'vessel_centre', 'widths');
end
%%
est_flow_to_ori_error_mean = nan(num_vessels, 10);
est_flow_to_ori_error_weighted = nan(num_vessels, 10);
est_flow_to_true_error_mean = nan(num_vessels, 10);
weighted_flow_rates = nan(num_vessels, 10);
prctile_flow_rates = nan(num_vessels,10,6);
prctile_matched_flow_rates = nan(num_vessels,10,6);
weighted_matched_flow_rates = nan(num_vessels, 10);
total_vessel_probs = nan(num_vessels, 10);
mean_widths = nan(num_vessels, 10);
shape_scores = nan(num_vessels, 10);
vessel_flow_rates = nan(num_vessels, 10);
bg_flow_rates = nan(num_vessels, 10);

true_flow_rates = nan(num_vessels, 10);
true_weighted_flow_rates = nan(num_vessels, 10);
true_flow_to_ori_error_mean = nan(num_vessels, 10);
true_flow_to_ori_error_weighted = nan(num_vessels, 10);
true_widths = nan(num_vessels, 10); 

true_noise = nan(num_vessels, 10);
true_mean_flow = nan(num_vessels, 10);

g_prob = gaussian_filters_1d(2);
g_prob = g_prob / sum(g_prob);
%%
for i_ve = 1:num_vessels
    
    display(['Processing vessel ' num2str(i_ve)]);
    vessel_name = contour_list(i_ve).name(1:end-6);
    
    load([flow_map_dir vessel_name 'fm.mat'], 'flowmap', 'mask', 'widths');
    flat_map = max(flowmap, [], 3);
    flat_mask = any(mask,3);
    flowmap = max(flowmap,[],3)  / nanmean(abs(flat_map(flat_mask)));
        
    for i_rpt = 1:10
         
        if ~exist([flow_metrics_dir vessel_name 'rpt' zerostr(i_rpt, 2) '.mat'], 'file');
            continue;
        end
        
        load([flow_results_dir vessel_name 'rpt' zerostr(i_rpt, 2) '.mat']);
        load([flow_metrics_dir vessel_name 'rpt' zerostr(i_rpt, 2) '.mat']);
        load([flow_params_dir vessel_name 'rpt' zerostr(i_rpt, 2) '.mat']);
        
        apex_estimated_flow = flow_results.flowPyramidEst{1}(flow_metrics.frame_apex_mask);
        prctile_flow_rates(i_ve,i_rpt,:) = prctile(abs(apex_estimated_flow), 50:10:100);
       
        flow_map_i = flowmap * mean_flow;
        flow_map_i(size(flow_metrics.frame_apex_mask,1), size(flow_metrics.frame_apex_mask,2)) = 0;
        flat_mask(size(flow_metrics.frame_apex_mask,1), size(flow_metrics.frame_apex_mask,2)) = 0;
        valid_gt_pts = flat_mask(flow_metrics.frame_apex_mask);

        %est_flow_to_ori_error_mean(i_ve,i_rpt) = flow_metrics.mean_error;
        est_flow_to_ori_error_weighted(i_ve,i_rpt) = flow_metrics.mean_weighted_error;
        weighted_flow_rates(i_ve,i_rpt) = flow_metrics.weighted_flow_rate;
        total_vessel_probs(i_ve,i_rpt) = flow_metrics.total_vessel_prob;
        mean_widths(i_ve,i_rpt) = flow_metrics.mean_width;
        shape_scores(i_ve,i_rpt) = flow_metrics.shape_score;
        vessel_flow_rates(i_ve,i_rpt) = flow_metrics.vessel_flow;
        bg_flow_rates(i_ve,i_rpt) = flow_metrics.background_flow;
        true_noise(i_ve,i_rpt) = speckle_noise;
        true_mean_flow(i_ve,i_rpt) = mean_flow;
        true_widths(i_ve,i_rpt) = mean(widths);

        smooth_ves = conv2(g_prob', g_prob, double(flow_metrics.vessel_pred)/100, 'same');
        smooth_ori = conv2(g_prob', g_prob, rgb2complex(flow_metrics.vessel_ori,[],1), 'same');

        apex_gt_flow = flow_map_i(flat_mask & flow_metrics.frame_apex_mask);
        %apex_gt_flow(isnan(apex_gt_flow)) = 0;
        if any(isnan(apex_gt_flow))
            display('Something gone wrong!');
        end
        
        apex_pred = smooth_ves(flow_metrics.frame_apex_mask);
        apex_ori = smooth_ori(flow_metrics.frame_apex_mask);

        true_flow_velocity = abs(apex_gt_flow);
        true_flow_angle = exp(2i*angle(apex_gt_flow));
        flow_angle_diffs = apex_pred(valid_gt_pts) .* apex_ori(valid_gt_pts) .* true_flow_angle;
        flow_angle_weights = abs(flow_angle_diffs);
        flow_angle_dirs = angle(flow_angle_diffs);% / 2

        true_flow_rates(i_ve,i_rpt) = mean(abs(apex_gt_flow));
        true_weighted_flow_rates(i_ve,i_rpt) = sum(apex_pred(valid_gt_pts).*true_flow_velocity) /...
            sum(apex_pred(valid_gt_pts));

        true_flow_to_ori_error_mean(i_ve,i_rpt) = mean(abs(flow_angle_dirs));
        true_flow_to_ori_error_weighted(i_ve,i_rpt) = sum(flow_angle_weights .* abs(flow_angle_dirs)) /...
            sum(flow_angle_weights);
        
        est_flow_velocity = abs(apex_estimated_flow);
        est_flow_angle = exp(2i*angle(apex_estimated_flow));
        flow_angle_diffs = apex_pred .* apex_ori .* est_flow_angle;
        flow_angle_dirs = angle(flow_angle_diffs);%/2
        flow_angle_weights = abs(flow_angle_diffs).*cos(flow_angle_dirs);

        matched_pts = abs(flow_angle_dirs) < pi/3;
        prctile_matched_flow_rates(i_ve,i_rpt,:) =...
            prctile(abs(apex_estimated_flow(matched_pts)), 50:10:100);
        
        weighted_matched_flow_rates(i_ve,i_rpt) =...
            sum(est_flow_velocity .* flow_angle_weights) / sum(flow_angle_weights);
        est_flow_to_ori_error_mean(i_ve,i_rpt) = mean(abs(flow_angle_dirs));
        
        flow_angle_diffs = true_flow_angle .* conj(est_flow_angle(valid_gt_pts));
        flow_angle_dirs = angle(flow_angle_diffs);
        est_flow_to_true_error_mean(i_ve,i_rpt) = mean(abs(flow_angle_dirs));

        if false && i_rpt == 1
            figure;
            subplot(2,3,1); imgray(flow_metrics.vessel_pred); 
            plot(flow_metrics.apex_xy(1), flow_metrics.apex_xy(2), 'rx');
            title(['ME = ' num2str(flow_metrics.mean_error,3) ' MWE = ' num2str(flow_metrics.mean_weighted_error, 3)]);
            subplot(2,3,2); imgray(flow_metrics.vessel_ori); 
            plot(flow_metrics.apex_xy(1), flow_metrics.apex_xy(2), 'rx');
            title(['MFR = ' num2str(flow_metrics.weighted_flow_rate, 3)]);        
            subplot(2,3,3); imgray(flow_metrics.vessel_ori); 
            plot(flow_metrics.apex_xy(1), flow_metrics.apex_xy(2), 'rx');
            title(['Flow ratio = ' num2str(flow_metrics.vessel_flow / flow_metrics.background_flow, 3)]);

            subplot(2,3,4); imgray(flow_metrics.frame_apex_mask);
            plot(flow_metrics.apex_xy(1), flow_metrics.apex_xy(2), 'rx');
            subplot(2,3,5); imgray(complex2rgb(flow_map_i, [], [], [], 1));
            plot(flow_metrics.apex_xy(1), flow_metrics.apex_xy(2), 'rx');
            subplot(2,3,6); imgray(complex2rgb(flow_results.flowPyramidEst{1}, [], [], [], 1));
            plot(flow_metrics.apex_xy(1), flow_metrics.apex_xy(2), 'rx');
        end
    end
    
end
%%
%lets ignore giants, they don't work too well
valid = [true(2130,1); false(450,1); true(1920,1)];
%%
figure; 
subplot(1,2,1); hold on; axis equal;
plot(vessel_flow_rates, weighted_flow_rates, 'rx'); axis equal;
plot([0 7], [0 7], 'k--');
axis([0 7 0 7]);
title('Estimated flow - mean vs weighted mean');
xlabel('Estimated flow (mean)');
ylabel('Estimated flow (weighted mean)');
subplot(1,2,2); hold on; axis equal;
plot(weighted_flow_rates, weighted_matched_flow_rates, 'rx'); axis equal;
plot([0 7], [0 7], 'k--');
axis([0 7 0 7]);
title('Estimated flow, weighted mean - matched vs unmatched');
xlabel('Estimated flow (unmatched)');
ylabel('Estimated flow (matched)');
%%
pts = true_weighted_flow_rates > 0;

figure; 
subplot(1,2,1); hold on; axis equal;
plot(true_flow_rates, weighted_flow_rates, 'rx');
plot([0 12], [0 12], 'k--');
[tf, ef] = ...
    kernel_smoother(true_flow_rates(pts), weighted_flow_rates(pts), 50);
plot(tf, ef, 'b-', 'linewidth', 2);
axis([0 12 0 12]);
title('True vs estimated flow - mean in apex ROI');
xlabel('True flow');
ylabel('Estimated flow (weighted mean)');

subplot(1,2,2); hold on; axis equal;
plot(true_flow_rates, weighted_matched_flow_rates, 'rx')
plot([0 12], [0 12], 'k--');
[tf, ef] = ...
    kernel_smoother(true_flow_rates(pts), weighted_matched_flow_rates(pts), 50);
plot(tf, ef, 'b-', 'linewidth', 2);
axis([0 12 0 12]);
title('True vs estimated flow - weighted mean in apex ROI');
xlabel('True flow');
ylabel('Estimated flow (matched + weighted mean)');
%%
pts = true_weighted_flow_rates > 0;

figure; 
subplot(1,2,1); hold on;
plot(true_flow_rates, est_flow_to_true_error_mean, 'rx');
[tf, ef] = ...
    kernel_smoother(true_flow_rates(pts), est_flow_to_true_error_mean(pts), 50);
plot(tf, ef, 'b-', 'linewidth', 2);
axis([0 12 0 pi/2]);
title('Flow velocity vs direction error');
xlabel('True flow');
ylabel('|Est. flow direction - true flow direction|');

subplot(1,2,2); hold on;
plot(true_flow_rates, est_flow_to_ori_error_mean, 'rx')
[tf, ef] = ...
    kernel_smoother(true_flow_rates(pts), est_flow_to_ori_error_mean(pts), 50);
plot(tf, ef, 'b-', 'linewidth', 2);
axis([0 12 0 pi/2]);
title('Flow velocity vs direction error');
xlabel('True flow');
ylabel('|Est. flow direction - est. orientation|');
%%
figure
subplot(2,2,1); hold on; axis equal;
plot(true_flow_to_ori_error_mean, est_flow_to_ori_error_mean, 'rx')
plot([0 pi/2], [0 pi/2], 'k--');
title('True vs estimated flow direction errors - mean in apex ROI');
xlabel('Direction error (true flow vs est ori)');
ylabel('Direction error (est flow vs est ori)');

subplot(2,2,2); hold on; axis equal;
plot(est_flow_to_true_error_mean, est_flow_to_ori_error_mean, 'rx')
plot([0 pi/2], [0 pi/2], 'k--');
title('True vs estimated flow direction errors - mean in apex ROI');
xlabel('Direction error (true flow vs est flow)');
ylabel('Direction error (est flow vs est ori)');

subplot(2,2,3); hold on; axis equal;
plot(true_flow_to_ori_error_mean, est_flow_to_true_error_mean, 'rx')
plot([0 pi/2], [0 pi/2], 'k--');
title('True vs estimated flow direction errors - mean in apex ROI');
xlabel('Direction error (true flow vs est ori)');
ylabel('Direction error (true flow vs est flow)');
%%
figure
subplot(2,2,1); hold on; axis equal;
plot(true_flow_to_ori_error_mean, est_flow_to_ori_error_mean, 'rx')
plot([0 pi/2], [0 pi/2], 'k--');
title('True vs estimated flow direction errors - mean in apex ROI');
xlabel('Direction error (true flow vs est ori)');
ylabel('Direction error (est flow vs est ori)');

%%
figure; 
for i_pct = 1:6
    subplot(2,3,i_pct); hold on; axis equal;
    plot(true_flow_rates, prctile_flow_rates(:,:,i_pct), 'rx');
    plot([0 12], [0 12], 'k--');
    pf = prctile_flow_rates(:,:,i_pct);
    [tf, ef] = ...
        kernel_smoother(true_flow_rates(pts), pf(pts), 50);
    plot(tf, ef, 'b-', 'linewidth', 2);
    
    title('True vs estimated flow - %tile of apex ROI');
    xlabel('True flow');
    ylabel(['Estimated flow - ' num2str(40 + 10*i_pct) 'th']);
end
%%
figure; 
for i_pct = 1:6
    subplot(2,3,i_pct); hold on; axis equal;
    plot(true_flow_rates, prctile_matched_flow_rates(:,:,i_pct), 'rx');
    plot([0 12], [0 12], 'k--');
    pf = prctile_matched_flow_rates(:,:,i_pct);
    [tf, ef] = ...
        kernel_smoother(true_flow_rates(pts), pf(pts), 50);
    plot(tf, ef, 'b-', 'linewidth', 2);
    
    title('True vs estimated flow - %-ile of apex ROI');
    xlabel('True flow');
    ylabel(['Estimated flow - ' num2str(40 + 10*i_pct) 'th']);
end
%%
figure; 
subplot(1,2,1); hold on; axis equal;
plot(true_flow_rates, prctile_flow_rates(:,1), 'rx');
plot([0 12], [0 12], 'k--');

title('True vs estimated flow - mean in apex ROI');
xlabel('True flow');
ylabel('Estimated flow');

subplot(1,2,2); hold on; axis equal;
plot(true_flow_rates, prctile_flow_rates(:,4), 'rx');
plot([0 12], [0 12], 'k--');

title('True vs estimated flow - mean in apex ROI');
xlabel('True flow');
ylabel('Estimated flow');
%%
error_prctiles = prctile(est_flow_to_ori_error_weighted, linspace(0, 100, 9));
figure;
for ii = 1:8
    pts = est_flow_to_ori_error_weighted > error_prctiles(ii) & est_flow_to_ori_error_weighted <= error_prctiles(ii+1);
    subplot(2,4,ii); axis equal; hold all;
    plot(true_weighted_flow_rates(pts), prctile_matched_flow_rates(pts,5), 'x');
    plot([0 12], [0 12], 'k--');
    plot([0 12], [2 2], 'k--');
    xlabel('True flow');
    ylabel('Estimated flow');
end
%%
interesting = find(est_flow_to_ori_error_weighted < 0.2201 & true_flow_rates > 8);
for idx = 1:length(interesting)
    i_ve = ceil(interesting(idx) / 10);
    i_rpt = interesting(idx) - 10*(i_ve-1);
    
    vessel_name = contour_list(i_ve).name(1:end-6);
    
    load([flow_map_dir vessel_name 'fm.mat'], 'flowmap');  
    load([flow_results_dir vessel_name 'rpt' zerostr(i_rpt, 2) '.mat']);
    load([flow_metrics_dir vessel_name 'rpt' zerostr(i_rpt, 2) '.mat']);
    load([flow_params_dir vessel_name 'rpt' zerostr(i_rpt, 2) '.mat']);
    flow_map_i = max(flowmap,[],3) * mean_flow;
    
    figure;
    subplot(2,3,1); imgray(flow_metrics.vessel_pred); 
    plot(flow_metrics.apex_xy(1), flow_metrics.apex_xy(2), 'rx');
    title(['ME = ' num2str(flow_metrics.mean_error,3) ' MWE = ' num2str(flow_metrics.mean_weighted_error, 3)]);
    subplot(2,3,2); imgray(complex2rgb(flow_map_i, [], [], [], 1)); 
    plot(flow_metrics.apex_xy(1), flow_metrics.apex_xy(2), 'rx');
    title(['MFR = ' num2str(flow_metrics.weighted_flow_rate, 3)]);        
    subplot(2,3,3); imgray(complex2rgb(flow_results.flowPyramidEst{1}, [], [], [], 1)); 
    plot(flow_metrics.apex_xy(1), flow_metrics.apex_xy(2), 'rx');
    title(['True FR = ' num2str(true_weighted_flow_rates(interesting(idx)), 3)]); 

    subplot(2,3,4); imgray(flow_metrics.vessel_ori);
    plot(flow_metrics.apex_xy(1), flow_metrics.apex_xy(2), 'rx');
    subplot(2,3,5); imgray(abs(flow_map_i)); caxis([0 nanmax(abs(flow_map_i(:)))]); colorbar;
    plot(flow_metrics.apex_xy(1), flow_metrics.apex_xy(2), 'rx');
    subplot(2,3,6); imgray(abs(flow_results.flowPyramidEst{1})); colorbar;
    plot(flow_metrics.apex_xy(1), flow_metrics.apex_xy(2), 'rx');
end
%%
interesting = find(est_flow_to_ori_error_weighted < 0.2201 & true_flow_rates < 2 & true_flow_rates > 0);
for idx = 1:length(interesting)
    i_ve = ceil(interesting(idx) / 10);
    i_rpt = interesting(idx) - 10*(i_ve-1);
    
    vessel_name = contour_list(i_ve).name(1:end-6);
    
    load([flow_map_dir vessel_name 'fm.mat'], 'flowmap');  
    load([flow_results_dir vessel_name 'rpt' zerostr(i_rpt, 2) '.mat']);
    load([flow_metrics_dir vessel_name 'rpt' zerostr(i_rpt, 2) '.mat']);
    load([flow_params_dir vessel_name 'rpt' zerostr(i_rpt, 2) '.mat']);
    
    flat_map = max(flowmap, [], 3);
    flat_mask = any(mask,3);
    flowmap = max(flowmap,[],3)  / nanmean(abs(flat_map(flat_mask)));
    flow_map_i = flowmap * mean_flow;
    
    figure;
    subplot(2,3,1); imgray(flow_metrics.vessel_pred); 
    plot(flow_metrics.apex_xy(1), flow_metrics.apex_xy(2), 'rx');
    title(['ME = ' num2str(flow_metrics.mean_error,3) ' MWE = ' num2str(flow_metrics.mean_weighted_error, 3)]);
    subplot(2,3,2); imgray(complex2rgb(flow_map_i, [], [], [], 1)); 
    plot(flow_metrics.apex_xy(1), flow_metrics.apex_xy(2), 'rx');
    title(['MFR = ' num2str(flow_metrics.weighted_flow_rate, 3)]);        
    subplot(2,3,3); imgray(complex2rgb(flow_results.flowPyramidEst{1}, [], [], [], 1)); 
    plot(flow_metrics.apex_xy(1), flow_metrics.apex_xy(2), 'rx');
    title(['True FR = ' num2str(true_weighted_flow_rates(interesting(idx)), 3)]); 

    subplot(2,3,4); imgray(flow_metrics.vessel_ori);
    plot(flow_metrics.apex_xy(1), flow_metrics.apex_xy(2), 'rx');
    subplot(2,3,5); imgray(abs(flow_map_i)); caxis([0 nanmax(abs(flow_map_i(:)))]); colorbar;
    plot(flow_metrics.apex_xy(1), flow_metrics.apex_xy(2), 'rx');
    subplot(2,3,6); imgray(abs(flow_results.flowPyramidEst{1})); colorbar;
    plot(flow_metrics.apex_xy(1), flow_metrics.apex_xy(2), 'rx');
end
%%
interesting = find(est_flow_to_ori_error_weighted > 1.3 & valid);
for idx = 1:length(interesting)
    i_ve = ceil(interesting(idx) / 10);
    i_rpt = interesting(idx) - 10*(i_ve-1);
    
    vessel_name = contour_list(i_ve).name(1:end-6);
    
    load([flow_map_dir vessel_name 'fm.mat'], 'flowmap', 'mask');  
    load([flow_results_dir vessel_name 'rpt' zerostr(i_rpt, 2) '.mat']);
    load([flow_metrics_dir vessel_name 'rpt' zerostr(i_rpt, 2) '.mat']);
    load([flow_params_dir vessel_name 'rpt' zerostr(i_rpt, 2) '.mat']);
    flat_map = max(flowmap, [], 3);
    flat_mask = any(mask,3);
    flowmap = max(flowmap,[],3)  / nanmean(abs(flat_map(flat_mask)));
    flow_map_i = flowmap * mean_flow;
    
    figure;
    subplot(2,3,1); imgray(flow_metrics.vessel_pred); 
    plot(flow_metrics.apex_xy(1), flow_metrics.apex_xy(2), 'rx');
    title(['ME = ' num2str(flow_metrics.mean_error,3) ' MWE = ' num2str(flow_metrics.mean_weighted_error, 3)]);
    subplot(2,3,2); imgray(complex2rgb(flow_map_i, [], [], [], 1)); 
    plot(flow_metrics.apex_xy(1), flow_metrics.apex_xy(2), 'rx');
    title(['MFR = ' num2str(flow_metrics.weighted_flow_rate, 3)]);        
    subplot(2,3,3); imgray(complex2rgb(flow_results.flowPyramidEst{1}, [], [], [], 1)); 
    plot(flow_metrics.apex_xy(1), flow_metrics.apex_xy(2), 'rx');
    title(['True FR = ' num2str(true_weighted_flow_rates(interesting(idx)), 3)]); 

    subplot(2,3,4); imgray(flow_metrics.vessel_ori);
    plot(flow_metrics.apex_xy(1), flow_metrics.apex_xy(2), 'rx');
    subplot(2,3,5); imgray(abs(flow_map_i)); caxis([0 nanmax(abs(flow_map_i(:)))]); colorbar;
    plot(flow_metrics.apex_xy(1), flow_metrics.apex_xy(2), 'rx');
    subplot(2,3,6); imgray(abs(flow_results.flowPyramidEst{1})); colorbar;
    plot(flow_metrics.apex_xy(1), flow_metrics.apex_xy(2), 'rx');
end
%%
pts = valid & true_mean_flow > 1;
figure; 
subplot(2,2,1); hold on;
plot(true_weighted_flow_rates(pts), true_flow_to_ori_error_weighted(pts), 'rx');
[flow_centres, smoothed_err] = ...
    kernel_smoother(true_weighted_flow_rates(pts), true_flow_to_ori_error_weighted(pts), 50);
plot(flow_centres, smoothed_err, 'linewidth', 2);
 
title('Flow direction error vs velocity - mean in apex ROI');
ylabel('Flow direction error (using true flow)');
xlabel('True flow');

% subplot(2,2,2); hold on;
% plot(true_weighted_flow_rates(pts), est_flow_to_ori_error_weighted(pts), 'rx')
% [flow_centres, smoothed_errs] = ...
%     kernel_smoother(true_weighted_flow_rates(pts), est_flow_to_ori_error_weighted(pts), 50);
% plot(flow_centres, smoothed_errs, 'linewidth', 2);
% 
% title('Flow error vs velocity - mean in apex ROI');
% ylabel('Flow direction error (using estimated flow)');
% xlabel('True flow');

subplot(2,2,2); hold on;
plot(true_mean_flow(pts), est_flow_to_ori_error_weighted(pts), 'rx')
[flow_centres, smoothed_errs] = ...
    kernel_smoother(true_mean_flow(pts), est_flow_to_ori_error_weighted(pts), 50);
plot(flow_centres, smoothed_errs, 'linewidth', 2);

title('Flow error vs velocity - mean in apex ROI');
ylabel('Flow direction error (using estimated flow)');
xlabel('True flow');

subplot(2,2,3); hold on;
plot(weighted_flow_rates(pts), est_flow_to_ori_error_weighted(pts), 'rx')
[flow_centres, smoothed_err] = ...
    kernel_smoother(weighted_flow_rates(pts), est_flow_to_ori_error_weighted(pts), 50);
plot(flow_centres, smoothed_err, 'linewidth', 2);

title('Flow error vs velocity - mean in apex ROI');
ylabel('Flow direction error (using estimated flow)');
xlabel('Estimated flow');

% subplot(2,2,4); hold on;
% plot(abs(weighted_flow_rates(pts) - true_weighted_flow_rates(pts)) ./...
%     true_weighted_flow_rates(pts),...
%     est_flow_to_ori_error_weighted(pts), ...
%     'rx')
% [flow_centres, smoothed_err] = ...
%     kernel_smoother(abs(weighted_flow_rates(pts) - true_weighted_flow_rates(pts)) ./...
%     true_weighted_flow_rates(pts),...
%     est_flow_to_ori_error_weighted(pts), 50);
% plot(flow_centres, smoothed_err, 'linewidth', 2);

subplot(2,2,4); hold on;
plot(abs(prctile_flow_rates(pts,4) - true_weighted_flow_rates(pts)) ./...
    true_weighted_flow_rates(pts),...
    est_flow_to_ori_error_weighted(pts), ...
    'rx')
[flow_centres, smoothed_err] = ...
    kernel_smoother(abs(prctile_flow_rates(pts,4) - true_weighted_flow_rates(pts)) ./...
    true_weighted_flow_rates(pts),...
    est_flow_to_ori_error_weighted(pts), 50);
plot(flow_centres, smoothed_err, 'linewidth', 2);

title('Flow direction error vs velocity error - mean in apex ROI');
ylabel('Flow direction error (using estimated flow)');
xlabel('%|True flow - estimated flow|');
%%
pts = true_weighted_flow_rates > 0.5;
pts(1:10:end) = 0;

figure;
subplot(1,2,1); hold on;
plot(mean_widths(pts), ...
    abs(weighted_flow_rates(pts) - true_weighted_flow_rates(pts)) ./...
    true_weighted_flow_rates(pts), 'rx')

[err_centres, smoothed_err_diff] = ...
    kernel_smoother(mean_widths(pts),...
    abs(weighted_flow_rates(pts) - true_weighted_flow_rates(pts)) ./...
    true_weighted_flow_rates(pts),...
    50);
plot(err_centres, smoothed_err_diff, 'linewidth', 2);
title('Vessel width vs velocity error - mean in apex ROI');
ylabel('%|True flow - estimated flow|');
xlabel('Vessel width');

subplot(1,2,2); hold on;
plot(true_noise(pts), ...
    abs(weighted_flow_rates(pts) - true_weighted_flow_rates(pts)) ./...
    true_weighted_flow_rates(pts), 'rx')
[noise_centres, smoothed_err_diff] = ...
    kernel_smoother(true_noise(pts),...
    abs(weighted_flow_rates(pts) - true_weighted_flow_rates(pts)) ./...
    true_weighted_flow_rates(pts),...
    50);
plot(noise_centres, smoothed_err_diff, 'linewidth', 2);

title('Added noise vs velocity error - mean in apex ROI');
ylabel('%|True flow - estimated flow|');
xlabel('Added nois');
%%
figure; hold on;
plot(true_weighted_flow_rates, weighted_flow_rates, 'rx')
title('True vs estimated flow - weighted mean in apex ROI');
xlabel('True flow');
ylabel('Estimated flow');
plot([0 12], [0 12], 'k--');
%%
%Look at relative flow rates for each vessel
sorted_mean_flow = reshape(true_mean_flow, 10, 450);
sorted_flows = reshape(weighted_matched_flow_rates, 10, 450);
sorted_noise = reshape(true_noise, 10, 450);
sorted_errors = reshape(est_flow_to_ori_error_weighted, 10, 450);
sorted_widths = reshape(mean_widths, 10, 450);
sorted_bg_flows = reshape(bg_flow_rates, 10, 450);

[sorted_mean_flow, sort_order] = sort(sorted_mean_flow);
for i_ve = 1:450
    sorted_flows(:,i_ve) = sorted_flows(sort_order(:,i_ve),i_ve);
    sorted_noise(:,i_ve) = sorted_noise(sort_order(:,i_ve),i_ve);
    sorted_errors(:,i_ve) = sorted_errors(sort_order(:,i_ve),i_ve);
    sorted_widths(:,i_ve) = sorted_widths(sort_order(:,i_ve),i_ve);
    sorted_bg_flows(:,i_ve) = sorted_bg_flows(sort_order(:,i_ve),i_ve);
end

figure; hold all;
plot(sorted_mean_flow, sorted_flows, '-');

relative_mean_flow = bsxfun(@rdivide, sorted_mean_flow, sorted_mean_flow(5,:));
relative_flows = bsxfun(@rdivide, sorted_flows, sorted_flows(5,:));
relative_widths = sorted_widths;%reshape(true_widths, 10, 450);
relative_noise = sorted_noise;
relative_errors = sorted_errors;
relative_bg_flows = sorted_bg_flows;
abs_flows = repmat(sorted_flows(5,:), 10, 1);


figure; hold all;
plot(relative_mean_flow, relative_flows, 'x');
%%
correct_order = relative_flows > 1;

X = [...
    %reshape(relative_noise(:,valid_ve), [], 1) ...
    reshape(relative_widths(:,valid_ve), [], 1) ...
    reshape(relative_errors(:,valid_ve), [], 1) ...
    reshape(relative_bg_flows(:,valid_ve), [], 1) ...
    reshape(abs_flows(:,valid_ve), [], 1) ...
    reshape(relative_flows(:,valid_ve), [], 1)];


y = reshape(relative_mean_flow(:,valid_ve), [], 1);

mdl = fitglm(X, y, 'linear')
%%
X = [...
    reshape(relative_noise(:,valid_ve), [], 1) ...
    reshape(relative_widths(:,valid_ve), [], 1) ...
    reshape(relative_errors(:,valid_ve), [], 1) ...
    reshape(relative_bg_flows(:,valid_ve), [], 1) ...
    reshape(abs_flows(:,valid_ve), [], 1) ...
    reshape(relative_mean_flow(:,valid_ve), [], 1)];

%y = reshape(correct_order([1:4 6:10],valid_ve), [], 1);
y = reshape(relative_flows(:,valid_ve), [], 1);

% discard = X > 10;
% X(discard) = [];
% y(discard) = [];

mdl = fitglm(X, y, 'linear')

