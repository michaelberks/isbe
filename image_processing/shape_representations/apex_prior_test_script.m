vessel_prob_smoothing_sigma = 2;

g = gaussian_filters_1d(vessel_prob_smoothing_sigma);
g = g / sum(g);
rsa_dir = 'rsa_study';

ori_dir = [nailfoldroot 'data/' rsa_dir '/set12/predictions/orientation/rf_regression/222835/'];
prob_dir = [nailfoldroot 'data/' rsa_dir '/set12/predictions/detection/rf_classification/222836/'];
vessel_size = 'normal';

%Get list of vessel names
v_list = dir([ori_dir, vessel_size '*.mat']);

load('C:\isbe\nailfold\data\rsa_study\models\apex_curv_prior_model.mat', 'curv_stats');
prior_left = bsxfun(@rdivide, curv_stats.l.marked_curv_by_step.hist, sum(curv_stats.l.marked_curv_by_step.hist,2));
prior_right = bsxfun(@rdivide, curv_stats.r.marked_curv_by_step.hist, sum(curv_stats.r.marked_curv_by_step.hist,2));

ori_hist_r_smoothed = imfilter(curv_stats.r.pred_ori_by_step.hist, conv(g,g), 'circular');
ori_hist_r_smoothed = imfilter(ori_hist_r_smoothed, conv(g,g)', 'replicate');

ori_hist_l_smoothed = imfilter(curv_stats.l.pred_ori_by_step.hist, conv(g,g), 'circular');
ori_hist_l_smoothed = imfilter(ori_hist_l_smoothed, conv(g,g)', 'replicate');

apex_prior = cat(3, ori_hist_l_smoothed, ori_hist_r_smoothed);

num_streams = [1e2 5e2 1e3 5e3];
apex_prior_scores_all = zeros(50,10,3);
v_count = 0;
%
for i_sz = {'normal', 'enlarged'}
    
    vessel_size = i_sz{1};
    apex_dir = [nailfoldroot 'data/' rsa_dir '/apexes/' vessel_size '/'];
    contour_dir = [nailfoldroot 'data/' rsa_dir '/vessel_contours/' vessel_size '/'];

    %Get list of vessel names
    v_list = dir([ori_dir, vessel_size '*.mat']);
    %
    for i_ve = 1:25%1:20%
        v_count = v_count + 1;
    
        %load vessel
        apex_name = v_list(i_ve).name(length(vessel_size)+(1:8));
        contour_struc = load([contour_dir apex_name '_vessel_contour.mat']);
        apex_struc = load([apex_dir apex_name '.mat']);
        vessel_struc = u_load([apex_dir apex_name '_vessel.mat']);
        
        display(['Processing ' apex_name ' from ' vessel_size]);

        prob_patch = u_load([prob_dir v_list(i_ve).name]);
        prob_patch = conv2(g', g, prob_patch, 'same');
        ori_patch = u_load([ori_dir v_list(i_ve).name]);       

        %Resample the vessel pts to be equalled spaced, 1 pixel apart
        v_pts = spline_contour(contour_struc.vessel_centre, [], 1);

        %Correct apex coordinates frame
        apex_xy(:,1) = apex_struc.apex_xy(:,1) +...
            apex_struc.apex_properties.sc - vessel_struc.vessel_properties.sc;
        apex_xy(:,2) = apex_struc.apex_xy(:,2) +...
            apex_struc.apex_properties.sr - vessel_struc.vessel_properties.sr;

        %Find the vessel pt nearets the apex        
        dists = sum(bsxfun(@minus, v_pts, mean(apex_xy)).^2, 2);
        [~, apex_idx] = min(dists);

        for i_num = 1:4
            
            for i_rep = 1:10
                [apex_prior_scores] = apex_prior_prob_track_mult(...
                            prob_patch, ori_patch,...
                            round(v_pts(apex_idx,2)), round(v_pts(apex_idx,1)),...
                            'num_streams', num_streams(i_num), 'stopping_prob', 0.1,...
                            'apex_prior', apex_prior);
                        
                apex_prior_scores_all(v_count,i_rep,i_num) = sum(apex_prior_scores(:));
            end
        end
    end        
end
%%
save('C:\isbe\nailfold\data\rsa_study\models\apex_prior_stream_repeat_test_results.mat', 'apex_prior_scores_all');
mean_scores  = squeeze(mean(apex_prior_scores_all,2));
std_scores = squeeze(std(apex_prior_scores_all,1,2));

figure; hold all;
plot(mean_scores(:,1), std_scores(:,1), 'x');
plot(mean_scores(:,2), std_scores(:,2), '+');
plot(mean_scores(:,3), std_scores(:,3), '*');
plot(mean_scores(:,4), std_scores(:,4), 'o');
