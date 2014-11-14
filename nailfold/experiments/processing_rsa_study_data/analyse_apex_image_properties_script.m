%%
%file_names = {'1-36', '37-270', '271-393', '394-657'};
file_names = {''};

for i_s = 1:31
    %Load in the partial stats file
    v_stats_i = u_load(['C:\isbe\nailfold\data\rsa_study\models\apex_image_statistics_c_' zerostr(i_s, 2) '.mat']);
    
    %Workout which apices haven't been filled yet
    missing_apices = ~v_stats_i.vessel_width.apex_truth;

    %Remove these apices from the relevant 
    for i_v = {'vessel_prob', 'vessel_curv', 'vessel_corr', 'vessel_width'}
        for i_a = {'apex_mean', 'apex_max', 'apex_val'}
            if isfield(v_stats_i.(i_v{1}), i_a{1})
                v_stats_i.(i_v{1}).(i_a{1})(missing_apices) = [];
                display(['Discarding ' i_v{1} '.' i_a{1}]);
            end
        end
    end    
    v_stats_i.vessel_width.apex_truth(missing_apices) = [];
    v_stats_i.vessel_centres.apex_d2c(missing_apices) = [];
    v_stats_i.vessel_properties.shape(missing_apices) = [];
    v_stats_i.vessel_properties.size(missing_apices) = [];
    v_stats_i.vessel_properties.im_num(missing_apices) = [];
    v_stats_i.vessel_properties.apex_pos(missing_apices,:) = [];
    
    
    if i_s == 1;
        v_stats = v_stats_i;
    else
        %Combine the relevant fields
        for i_v = {'vessel_prob', 'vessel_curv', 'vessel_corr', 'vessel_width', 'vessel_centres'}
            for i_a = {'apex_mean', 'apex_max', 'apex_val', 'apex_truth', 'apex_d2c'}
                if isfield(v_stats_i.(i_v{1}), i_a{1})
                    display(['Combining ' i_v{1} '.' i_a{1}]);
                    v_stats.(i_v{1}).(i_a{1}) = ....
                        [v_stats_i.(i_v{1}).(i_a{1}); v_stats.(i_v{1}).(i_a{1})];
                end
            end
            if isfield(v_stats_i.(i_v{1}), 'background')
                v_stats.(i_v{1}).background = ....
                    v_stats.(i_v{1}).background + v_stats_i.(i_v{1}).background;
            end
        end
        v_stats.vessel_properties.shape = ...
            [v_stats_i.vessel_properties.shape; v_stats.vessel_properties.shape];
        v_stats.vessel_properties.size = ...
            [v_stats_i.vessel_properties.size; v_stats.vessel_properties.size];
        v_stats.vessel_properties.im_num = ...
            [v_stats_i.vessel_properties.im_num; v_stats.vessel_properties.im_num];
        v_stats.vessel_properties.apex_pos = ...
            [v_stats_i.vessel_properties.apex_pos; v_stats.vessel_properties.apex_pos];
    end
    clear v_stats_i;
end
%%
if ispc
    rsa_dir = 'rsa_study/';
else
    rsa_dir = [];
end

model_dir = [nailfoldroot 'data/' rsa_dir 'models/apex_templates/'];
image_dir = [nailfoldroot 'data/' rsa_dir 'test/images/'];
prob_dir = [nailfoldroot 'data/' rsa_dir 'test/predictions/detection/rf_classification/182321/'];
ori_dir = [nailfoldroot 'data/' rsa_dir 'test/predictions/orientation/rf_regression/182263/'];
width_dir = [nailfoldroot 'data/' rsa_dir 'test/predictions/width/rf_regression/182367/'];
centre_dir = [nailfoldroot 'data/' rsa_dir 'test/vessel_centres/'];

apex_radius = 5;
max_apex_guess = 5000;
vessel_prob_smoothing_sigma = 2;
curvature_smoothing_sigma = 2;
strong_vessel_thresh = 0.25;
curv_max = 0.5;
d2c_max = 100;
%%
vessel_prob_bins = linspace(0, 1, 100);
vessel_corr_bins = linspace(-1, 1, 100);
vessel_curv_bins = linspace(0, curv_max, 100);
vessel_d2c_bins = linspace(0, d2c_max, 100);

include_apexes = v_stats.vessel_width.apex_truth <= 50;
%include_apexes = true(size(v_stats.vessel_curv.apex_max));

vessel_prob_bg_cdf = cumsum(sum(v_stats.vessel_prob.background));
vessel_prob_bg_cdf = vessel_prob_bg_cdf / vessel_prob_bg_cdf(end);
figure; hold on;
plot(vessel_prob_bins, vessel_prob_bg_cdf);
plot(prctile(v_stats.vessel_prob.apex_max(include_apexes), 0:100), (0:100)/100, 'r');
title('CDF of vessel probabilities for non-apex pixels');
xlabel('Vessel probabilities');
ylabel('CDF %');

vessel_corr_bg_cdf = cumsum(sum(v_stats.vessel_corr.background));
vessel_corr_bg_cdf = vessel_corr_bg_cdf / vessel_corr_bg_cdf(end);
figure; hold on;
plot(vessel_corr_bins, vessel_corr_bg_cdf);
plot(prctile(v_stats.vessel_corr.apex_max(include_apexes), 0:100), (0:100)/100, 'r');
title('CDF of apex template correlation scores for non-apex pixels');
xlabel('Apex template correlation scores');
ylabel('CDF %');

vessel_curv_bg_cdf = cumsum(sum(abs(v_stats.vessel_curv.background)));
vessel_curv_bg_cdf = vessel_curv_bg_cdf / vessel_curv_bg_cdf(end);
figure; hold on;
plot(vessel_curv_bins, vessel_curv_bg_cdf);
plot(prctile(abs(v_stats.vessel_curv.apex_max(include_apexes)), 0:100), (0:100)/100, 'r');
title('CDF of orientation curvature for non-apex pixels');
xlabel('Orientation curvature');
ylabel('CDF %');

% vessel_d2c_bg_cdf = cumsum(sum(v_stats.vessel_centres.background));
% vessel_d2c_bg_cdf = vessel_d2c_bg_cdf / vessel_d2c_bg_cdf(end);
% figure; hold on;
% plot(vessel_d2c_bins, vessel_d2c_bg_cdf);
% plot(prctile(v_stats.vessel_centres.apex_d2c, 0:100), (0:100)/100, 'r');
% title('CDF of distance to centreline for non-apex pixels');
% xlabel('Distance to centreline');
% ylabel('CDF %');

vessel_prob_corr_bg_cdf = cumsum(sum(v_stats.vessel_prob_corr.background));
vessel_prob_corr_bg_cdf = vessel_prob_corr_bg_cdf / vessel_prob_corr_bg_cdf(end);
figure; hold on;
plot(vessel_corr_bins, vessel_prob_corr_bg_cdf);
plot(prctile(v_stats.vessel_prob.apex_max(include_apexes) .* v_stats.vessel_corr.apex_max(include_apexes), 0:100), (0:100)/100, 'r');
title('CDF of vessel probabilities multiplied by correlation scores for non-apex pixels');
xlabel('Vessel probabilities multiplied by correlation scores');
ylabel('CDF %');

vessel_prob_curv_bg_cdf = cumsum(sum(v_stats.vessel_prob_curv.background));
vessel_prob_curv_bg_cdf = vessel_prob_curv_bg_cdf / vessel_prob_curv_bg_cdf(end);
figure; hold on;
plot(vessel_curv_bins, vessel_prob_curv_bg_cdf);
plot(prctile(v_stats.vessel_prob.apex_max(include_apexes) .* v_stats.vessel_curv.apex_max(include_apexes), 0:100), (0:100)/100, 'r');
title('CDF of vessel probabilities multiplied by curvature for non-apex pixels');
xlabel('Vessel probabilities multiplied by curvature');
ylabel('CDF %');
%%
figure; plot(vessel_prob_bins, sum(v_stats.vessel_prob.background))

figure; plot(vessel_corr_bins, sum(v_stats.vessel_corr.background))
figure; plot(vessel_curv_bins, sum(v_stats.vessel_curv.background))
figure; plot(vessel_d2c_bins, sum(v_stats.vessel_centres.background))
%
figure; hist(v_stats.vessel_prob.apex_max, vessel_prob_bins);
figure; hist(v_stats.vessel_corr.apex_max, vessel_corr_bins);
figure; hist(v_stats.vessel_curv.apex_max, vessel_curv_bins);
figure; hist(v_stats.vessel_centres.apex_d2c, vessel_d2c_bins);
%%
figure; plot(v_stats.vessel_width.apex_val, v_stats.vessel_width.apex_truth, 'rx'); hold on;
plot([0 55], [0 55]);
title('Predicted width vs marked width at apexes');
xlabel('Predicted width');
ylabel('Marked width');

figure; plot(v_stats.vessel_prob.apex_val, v_stats.vessel_width.apex_truth, 'rx');
title('Vessel probability vs marked width at apexes');
xlabel('Vessel probability');
ylabel('Marked width');

figure; plot(v_stats.vessel_corr.apex_val, v_stats.vessel_width.apex_truth, 'rx');
title('Apex correlation vs marked width');
xlabel('Apex correlation score');
ylabel('Marked width');

figure; plot(v_stats.vessel_curv.apex_val, v_stats.vessel_width.apex_truth, 'rx');
title('Apex curvature vs marked width');
xlabel('Predicted width');
ylabel('Marked width');

figure; plot(v_stats.vessel_centres.apex_d2c, v_stats.vessel_width.apex_truth, 'rx');
title('Distance from apex to nearest centreline vs marked width');
xlabel('Distance from apex to nearest centreline');
ylabel('Marked width');
%%
corr_vs_prob = hist3([v_stats.vessel_prob.apex_max v_stats.vessel_corr.apex_max], [11 10]);
figure; imgray(corr_vs_prob); colormap(jet(256));

curv_vs_prob = hist3([v_stats.vessel_prob.apex_max v_stats.vessel_curv.apex_max], [11 10]);
figure; imgray(curv_vs_prob); colormap(jet(256));
%%
%Display curvatures
figure; hold all; axis equal;
for kappa = [0.01 0.02:0.02:0.2];
    x = 0;
    y = 0;
    theta = pi/2;
    ii = 1;
    while theta < 3*pi/2
        x(ii+1) = x(ii) + cos(theta);
        y(ii+1) = y(ii) + sin(theta);
        theta = theta + kappa;
        ii = ii + 1;
    end
    plot(x, y);
end
legend(strcat('\kappa = ', cellstr(num2str(([0.01 0.02:0.02:0.2])'))));
%%
pred_list = dir([prob_dir '*.mat']);
include_apexes = find(v_stats.vessel_curv.apex_max < 0.03);
include_im_nums = v_stats.vessel_properties.im_num(include_apexes);
im_nums = unique(include_im_nums);
%%
for i_im = im_nums(9:10)'%length(include_apexes)
    
    vessel_prob = u_load([prob_dir pred_list(i_im).name]);
    vessel_ori = u_load([ori_dir pred_list(i_im).name]);
    figure; imgray(complex2rgb(vessel_prob.*vessel_ori));
 
%         im_name = pred_list(i_im).name(1:end-9);
%         vessel_im = u_load([image_dir im_name '.mat']);
%         [rows cols] = size(vessel_im);
%         rows = round(.4*rows):round(.6*rows);
%         cols = round(.25*cols):round(.75*cols);
%         g_min = min(min(vessel_im(rows,cols)));
%         g_max = max(max(vessel_im(rows,cols)));
%         figure; imgray(vessel_im); caxis([g_min g_max]);

    for i_ap = include_apexes(include_im_nums==i_im)'
        apex_xy = v_stats.vessel_properties.apex_pos(i_ap,:);
        plot(apex_xy(1), apex_xy(2), 'g*');     
    end
        
end
%%
load('C:\isbe\nailfold\data\rsa_study\test\vessel_centres\10147c_vc.mat');
vessel_im = u_load([image_dir '10147c.mat']);
vessel_prob = u_load([prob_dir '10147c_pred.mat']);
vessel_ori = u_load([ori_dir '10147c_pred.mat']);

g = gaussian_filters_1d(vessel_prob_smoothing_sigma);
g = g / sum(g);
vessel_prob = conv2(g', g, vessel_prob, 'same');

include_pts = (vessel_centre_prob>.5) & (vessel_centre_curv>.05);

figure; imgray(vessel_im); caxis([-20 10]);
plot(vessel_centre_x, vessel_centre_y, 'r.', 'markersize', 2);
plot(vessel_centre_x(include_pts), vessel_centre_y(include_pts), 'g.', 'markersize', 4);

[sc_pi pc_pi] = shape_context_prob_track_mult(vessel_ori.*vessel_prob, vessel_prob, vessel_ori,...
    vessel_centre_y(include_pts), vessel_centre_x(include_pts), 'num_streams', 1e3);

%%
cum_probs = zeros(100);
max_curv = max(v_stats.vessel_curv.apex_max);
num_apexes = length(v_stats.vessel_prob.apex_max);
for ii = 1:100
    p_thr = ii / 100;
    p_counts = v_stats.vessel_prob.apex_max > p_thr;
    for jj = 1:100
        c_thr =  max_curv*jj / 100;
        
        cum_probs(ii,jj) = sum(p_counts & v_stats.vessel_curv.apex_max > c_thr) / num_apexes;
    end
end
figure; imgray(cum_probs);
ylabel('Vessel probs');
xlabel('Vessel curvature');
%%
cum_probs = zeros(100);
max_corr = max(v_stats.vessel_corr.apex_max);
num_apexes = length(v_stats.vessel_prob.apex_max);
for ii = 1:100
    p_thr = ii / 100;
    p_counts = v_stats.vessel_prob.apex_max > p_thr;
    for jj = 1:100
        c_thr =  max_corr*jj / 100;
        
        cum_probs(ii,jj) = sum(p_counts & v_stats.vessel_corr.apex_max > c_thr) / num_apexes;
    end
end
figure; imgray(cum_probs);
ylabel('Vessel probs');
xlabel('Vessel correlation');

%%
%Lets set up a workflow...
include_apexes = v_stats.vessel_width.apex_truth <= 50;
total_apexes = sum(include_apexes);

%Now select all points with a corr > 0.4. These are almost certainly
%vessels and saves time from further processing
corr_apexes_idx = include_apexes & (v_stats.vessel_corr.apex_max > 0.4);
corr_apexes = sum(corr_apexes_idx);
remaining_apexes = sum(~corr_apexes_idx);

figure; hold on;
plot(prctile(v_stats.vessel_prob.apex_max(corr_apexes_idx), 0:100), (0:100)/100, 'b');
plot(prctile(v_stats.vessel_prob.apex_max(~corr_apexes_idx), 0:100), (0:100)/100, 'r');
title('CDF of vessel probabilities for non-apex pixels');
xlabel('Vessel probabilities');
ylabel('CDF %');

figure; hold on;
plot(prctile(v_stats.vessel_corr.apex_max(corr_apexes_idx), 0:100), (0:100)/100, 'b');
plot(prctile(v_stats.vessel_corr.apex_max(~corr_apexes_idx), 0:100), (0:100)/100, 'r');
title('CDF of apex template correlation scores for non-apex pixels');
xlabel('Apex template correlation scores');
ylabel('CDF %');


figure; hold on;
plot(prctile(abs(v_stats.vessel_curv.apex_max(corr_apexes_idx)), 0:100), (0:100)/100, 'b');
plot(prctile(abs(v_stats.vessel_curv.apex_max(~corr_apexes_idx)), 0:100), (0:100)/100, 'r');
title('CDF of orientation curvature for non-apex pixels');
xlabel('Orientation curvature');
ylabel('CDF %');


        
        
        
        