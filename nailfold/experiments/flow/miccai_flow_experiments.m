%--------------------------------------------------------------------------
% Script for producing results for 2016 submission to MICCAI
%--------------------------------------------------------------------------
%%
%Define directory paths we'll use

%Data
flow_results_dir = 'N:\Nailfold Capillaroscopy\Wellcome\flow_results\';
flow_metrics_dir = 'C:\isbe\nailfold\data\wellcome_study\flow_metrics\';
flow_data_dir = 'N:\Nailfold Capillaroscopy\Wellcome\flow_data\';
capillary_data_dir = 'C:\isbe\nailfold\data\wellcome_study\capillary_data\';

%Ouput
paper_dir = 'C:\isbe\matlab_code\isbe\papers\2016miccai_flow\';
fig_dir = [paper_dir 'figs\'];

%%
% Load in results
load('C:\isbe\nailfold\data\wellcome_study\sequence_names.mat');
load('C:\isbe\nailfold\data\wellcome_study\subject_summary.mat');
load('C:\isbe\nailfold\data\wellcome_study\results\auto_stats.mat');

hc_idx = strcmpi(subject_summary(:,3), '0');
pr_idx = strcmpi(subject_summary(:,3), '1');
ss_idx = strcmpi(subject_summary(:,3), '2') | strcmpi(subject_summary(:,3), '3');
%%

flow_list = dir([flow_metrics_dir '*.mat']);
%
num_vessels = length(flow_list);
mean_errors = zeros(num_vessels,1);
mean_weighted_errors = zeros(num_vessels,1);
weighted_flow_rates = zeros(num_vessels,1);
total_vessel_probs = zeros(num_vessels,1);
mean_widths = zeros(num_vessels,1);
shape_scores = zeros(num_vessels,1);
vessel_flow_rates = zeros(num_vessels,1);
bg_flow_rates = zeros(num_vessels,1);
%
for i_ve = 1:num_vessels
    f = u_load([flow_metrics_dir flow_list(i_ve).name]);
    mean_errors(i_ve,1) = f.mean_error;
    mean_weighted_errors(i_ve,1) = f.mean_weighted_error;
    weighted_flow_rates(i_ve,1) = f.weighted_flow_rate;
    total_vessel_probs(i_ve,1) = f.total_vessel_prob;
    mean_widths(i_ve,1) = f.mean_width;
    shape_scores(i_ve,1) = f.shape_score;
    vessel_flow_rates(i_ve,1) = f.vessel_flow;
    bg_flow_rates(i_ve,1) = f.background_flow;
end
flow_ratios = vessel_flow_rates ./ bg_flow_rates;
%%
save C:\isbe\nailfold\data\wellcome_study\results\flow_metrics_data.mat
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
hist(vessel_flow_rates, linspace(0,8,100));
title('Histogram of mean flow_rates');

figure;
hist(bg_flow_rates, linspace(0,8,100));
title('Histogram of mean background flow_rates');

figure;
hist(vessel_flow_rates ./ bg_flow_rates, linspace(0,8,100));
title('Histogram of mean background flow_rates');

figure;
subplot(1,2,1);
plot(mean_weighted_errors, weighted_flow_rates,'rx');
xlabel('Weighted flow error');
ylabel('Weighted flow rate');

subplot(1,2,2);
plot(mean_weighted_errors, shape_scores,'rx');
xlabel('Weighted flow error');
ylabel('Shape score');

figure;
subplot(1,2,1);
plot(mean_weighted_errors, weighted_flow_rates,'rx');
xlabel('Weighted flow error');
ylabel('Weighted flow rate');

subplot(1,2,2);
plot(shape_scores, weighted_flow_rates,'rx');
xlabel('Shape score');
ylabel('Weighted flow rate');

figure;
plot(mean_widths, weighted_flow_rates,'rx');
xlabel('Vessel width');
ylabel('Weighted flow rate');

figure;
plot(mean_weighted_errors, flow_ratios,'rx');
xlabel('Weighted flow error');
ylabel('Flow ratio');

figure;
subplot(1,2,1);
plot(mean_weighted_errors, vessel_flow_rates,'rx');
xlabel('Weighted flow error');
ylabel('Vessel flow rate');

subplot(1,2,2);
plot(mean_weighted_errors, bg_flow_rates,'rx');
xlabel('Weighted flow error');
ylabel('Background flow');
%%
selected_features = {...
    'vessel_density',...%'mean_mean_width',...
    'mean_mean_width',...%adjusted_width
    'max_mean_width',...
    'mean_connected_orientation_dispersion',...
    'dispersion_connected_orientation',...%
    'mean_mean_flow'};%,...
    %'mean_weighted_flow',...
    %'mean_vessel_flow'};

feature_display_names = {...
    'Capillary density',...%'Mean width',..
    'Mean width',...
    'Max width',...
    'Shape score',...
    'Derangement score',...%
    'Mean flow',...
    'Weighted flow',...
    'Vessel flow'};

num_features = length(selected_features);
all_stats = zeros(112, num_features);

all_roc_pts = zeros(100, 2, num_features);
all_roc_auc = zeros(1, num_features);
all_roc_se = zeros(1, num_features);

all_dist_means = zeros(3, num_features);
all_dist_se = zeros(3, num_features);

for i_f = 1:num_features
    dist_f = auto_stats.(selected_features{i_f});
    valid_ims = sum(auto_stats.present & auto_stats.gradeable...
        & ~isnan(dist_f) & (dist_f < inf), 2);       
    dist_f = naNsum(dist_f, 2) ./  valid_ims;
    all_stats(:,i_f) = dist_f;
    
    all_dist_means(1,i_f) = mean(dist_f(hc_idx));
    all_dist_means(2,i_f) = mean(dist_f(pr_idx));
    all_dist_means(3,i_f) = mean(dist_f(ss_idx));
    
    all_dist_se(1,i_f) = 1.96*std(dist_f(hc_idx)) / sqrt(sum(hc_idx));
    all_dist_se(2,i_f) = 1.96*std(dist_f(pr_idx)) / sqrt(sum(pr_idx));
    all_dist_se(3,i_f) = 1.96*std(dist_f(ss_idx)) / sqrt(sum(ss_idx));
    
    if mean(dist_f(ss_idx)) < mean(dist_f(~ss_idx))
        dist_f = -dist_f;
    end   
    [all_roc_pts(:,:,i_f), all_roc_auc(:,i_f), ~, ~, all_roc_se(:,i_f)] =...
        calculate_roc_curve(dist_f, ss_idx);
 
end
%%
markers = 'xosv^*';
figure('windowstyle', 'normal'); axis equal; axis([0 1 0 1]); hold all;
label = cell(num_features,1);
for i_f = 1:num_features
    plot(all_roc_pts(:,1,i_f), all_roc_pts(:,2,i_f), ...
        ['-' markers(i_f)], 'linewidth', 2);
    label{i_f} = [feature_display_names{i_f} ': A_z = ' num2str(all_roc_auc(i_f),2) ' \pm ' num2str(all_roc_se(i_f),2)];
end
legend(label, 'location', 'southeast');

title('ROC curves for individual biomarkers', 'fontsize', 18);
xlabel('FPR', 'fontsize', 18);
ylabel('TPR', 'fontsize', 18);
set(gca, 'fontsize', 14);
exportfig([fig_dir 'individual_rocs.png']);
%%
ypred_orig = zeros(112,1);
ypred_flow = zeros(112,1);
for i_sub = 1:112
    sub_idx = [1:i_sub-1 i_sub+1:112];
    mdl_orig = stepwiseglm(all_stats(sub_idx,1:5), ss_idx(sub_idx), 'linear','Distribution','binomial','link','logit'); %stepwise
    mdl_flow = stepwiseglm(all_stats(sub_idx,:), ss_idx(sub_idx), 'linear','Distribution','binomial','link','logit');
    
    ypred_orig(i_sub) = predict(mdl_orig, all_stats(i_sub,1:5));
    ypred_flow(i_sub) = predict(mdl_flow, all_stats(i_sub,:));
end
%%
figure('windowstyle', 'normal'); axis equal; axis([0 1 0 1]); hold all;

[roc_pts, auc, ~, ~, auc_se] = calculate_roc_curve(ypred_orig, ss_idx);
plot(roc_pts(:,1), roc_pts(:,2), 'r-o', 'linewidth', 2);
label{1} = ['Structure: A_z = ' num2str(auc,'%3.3f') ' \pm ' num2str(auc_se,'%2.2g')];

[roc_pts, auc, ~, ~, auc_se] = calculate_roc_curve(ypred_flow, ss_idx);
plot(roc_pts(:,1), roc_pts(:,2), 'b-s', 'linewidth', 2);
label{2} = ['Structure + flow: A_z = ' num2str(auc,'%3.3f') ' \pm ' num2str(auc_se,'%2.2g')];
legend(label, 'location', 'southeast');

title('ROC curves for combined biomarkers', 'fontsize', 18);
xlabel('FPR', 'fontsize', 18);
ylabel('TPR', 'fontsize', 18);
set(gca, 'fontsize', 14);
exportfig([fig_dir 'combined_rocs.png']);


%%
sequence_name = ['N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\035wellcome\2015_03_23\R2_12_24_58\'...
    'sequence_frames_data.dat'];
seq_name = '035wellcome_2015_03_23_R2_12_24_58';
sequence = load(sequence_name);
[segments_s, segments_ns, motor_x, motor_y, motor_z, sharpness] = ...
    get_stationary_segments(sequence, 120);
plot_segment_traces(1000, segments_s, segments_ns, motor_x, motor_y, motor_z);

exportfig('C:\isbe\matlab_code\isbe\papers\2016miccai_flow\figs\sequence_trace.png');
%%
vessel_list = dir([flow_results_dir seq_name '*s07*.mat']);
for i_ve = 1:length(vessel_list)
    load([flow_data_dir vessel_list(i_ve).name], 'cropped_frames');
    load([flow_results_dir vessel_list(i_ve).name], 'flow_results');
    
    write_im_from_colormap(mean(cropped_frames,3), ...
        [fig_dir 'vessel' zerostr(i_ve,2) '.png'], gray(256));
    imwrite(complex2rgb(flow_results.flowPyramidEst{1}, [], [], [], 1),...
        [fig_dir 'vessel_fm' zerostr(i_ve,2) '.png']);
end
%%
%--------------------------------------------------------------------------
%Make table of results
%--------------------------------------------------------------------------
table_path = [paper_dir 'results_table.txt'];
%Write a table from this
o_txt = {...
    '\specialcell{Capillary density\\$mm^{-1}$}   ';...
    '\specialcell{Mean width\\${\mu}m$}           ';...
    '\specialcell{Max width\\${\mu}m$}            ';...
    '\specialcell{Shape score\\(no units)}        ';...
    '\specialcell{Derangment score\\(no units)}   ';...
    '\specialcell{Mean flow velocity\\$mms^{-1}$} '};

fid = fopen(table_path, 'wt');
fprintf(fid, '%s \n', '\begin{tabular*}{0.95\textwidth}{@{\extracolsep{\fill} } l r r r r}');
fprintf(fid, '%s \n', '\toprule');
fprintf(fid, '%s \n', '%');

fprintf(fid, '%s \n', 'Biomarker & \multicolumn{3}{c}{Subject group means ($95\%$ CI)} 	           & \multicolumn{1}{c}{ROC $A_z$}	\\');
fprintf(fid, '%s \n', '          & \multicolumn{1}{c}{HC} & \multicolumn{1}{c}{PRP} & \multicolumn{1}{c}{SSc} & \multicolumn{1}{c}{HC,PRP $v$ SSc}  	\\');
for i_f = 1:6
    fprintf(fid, '%s ', o_txt{i_f});
    fprintf(fid, '& %3.3g $\\pm$ %3.2f        ', all_dist_means(1,i_f), all_dist_se(1,i_f));
    fprintf(fid, '& %3.3g $\\pm$ %3.2f        ', all_dist_means(2,i_f), all_dist_se(2,i_f));
    fprintf(fid, '& %3.3g $\\pm$ %3.2f        ', all_dist_means(3,i_f), all_dist_se(3,i_f));
    fprintf(fid, '& %3.3g $\\pm$ %3.2f        ', all_roc_auc(i_f), all_roc_se(i_f));
    fprintf(fid, '%s \n', '\\');
end

fprintf(fid, '%s \n', '%');
fprintf(fid, '%s \n', '\bottomrule \noalign{\smallskip}');
fprintf(fid, '%s \n', '\end{tabular*}');
fclose(fid);
