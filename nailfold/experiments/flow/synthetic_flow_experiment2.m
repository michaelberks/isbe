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

contour_list = [dir([contour_dir 'enlarged*_vc.mat']); dir([contour_dir 'normal*_vc.mat'])];
num_vessels = length(contour_list);
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
         
        if ~exist([flow_metrics_dir vessel_name 'rpt' zerostr(i_rpt + 10, 2) '.mat'], 'file');
            continue;
        end
        
        load([flow_results_dir vessel_name 'rpt' zerostr(i_rpt + 10, 2) '.mat']);
        load([flow_metrics_dir vessel_name 'rpt' zerostr(i_rpt + 10, 2) '.mat']);
        load([flow_params_dir vessel_name 'rpt' zerostr(i_rpt + 10, 2) '.mat']);
        
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
figure; 
subplot(2,2,1);
hold on;
plot(true_noise(:,1:5), weighted_matched_flow_rates(:,1:5), 'rx');
plot(true_noise(1,1:5), mean(weighted_matched_flow_rates(:,1:5)), 'b-', 'linewidth', 2);
xlabel('Noise');
ylabel('Estimated flow rates');
title('Noise vs estimated flow, true flow = 3 px/frame');
 
subplot(2,2,2);
hold on;
plot(true_noise(:,6:10), weighted_matched_flow_rates(:,6:10), 'rx');
plot(true_noise(1,1:5), mean(weighted_matched_flow_rates(:,6:10)), 'b-', 'linewidth', 2);
xlabel('Noise');
ylabel('Estimated flow rates');
title('Noise vs estimated flow, true flow = 6 px/frame');

subplot(2,2,3);
hold on;
plot(true_noise(:,1:5), prctile_matched_flow_rates(:,1:5,4), 'rx');
plot(true_noise(1,1:5), mean(prctile_matched_flow_rates(:,1:5,4)), 'b-', 'linewidth', 2);
xlabel('Noise');
ylabel('Estimated flow rates - 90%');
title('Noise vs estimated flow, true flow = 3 px/frame');
 
subplot(2,2,4);
hold on;
plot(true_noise(:,6:10), prctile_matched_flow_rates(:,6:10,4), 'rx');
plot(true_noise(1,1:5), mean(prctile_matched_flow_rates(:,6:10,4)), 'b-', 'linewidth', 2);
xlabel('Noise');
ylabel('Estimated flow rates - 90th %');
title('Noise vs estimated flow, true flow = 6 px/frame');
%%
figure; 
subplot(2,2,1);
hold on;
plot(true_noise(:,1:5), est_flow_to_true_error_mean(:,1:5), 'rx');
plot(true_noise(1,1:5), mean(est_flow_to_true_error_mean(:,1:5)), 'b-', 'linewidth', 2);
xlabel('Noise');
ylabel('|Est. flow direction - true flow direction|');
title('Noise vs flow orientation error, true flow = 3 px/frame');
 
subplot(2,2,2);
hold on;
plot(true_noise(:,6:10), est_flow_to_true_error_mean(:,6:10), 'rx');
plot(true_noise(1,1:5), mean(est_flow_to_true_error_mean(:,6:10)), 'b-', 'linewidth', 2);
xlabel('Noise');
ylabel('|Est. flow direction - true flow direction|');
title('Noise vs flow orientation error, true flow = 6 px/frame');

subplot(2,2,3);
hold on;
plot(true_noise(:,1:5), est_flow_to_ori_error_mean(:,1:5), 'rx');
plot(true_noise(1,1:5), mean(est_flow_to_ori_error_mean(:,1:5)), 'b-', 'linewidth', 2);
xlabel('Noise');
ylabel('|Est. flow direction - est. orientation|');
title('Noise vs flow orientation error, true flow = 3 px/frame');
 
subplot(2,2,4);
hold on;
plot(true_noise(:,6:10), est_flow_to_ori_error_mean(:,6:10), 'rx');
plot(true_noise(1,1:5), mean(est_flow_to_ori_error_mean(:,6:10)), 'b-', 'linewidth', 2);
xlabel('Noise');
ylabel('|Est. flow direction - est. orientation|');
title('Noise vs flow orientation error, true flow = 6 px/frame');


