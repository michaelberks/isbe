retroot = [asymmetryroot 'data\retinograms\DRIVE\test\predictions\orientations\'];

mkdir([retroot 'g2d\analytic\orientations']);
mkdir([retroot 'g2d\analytic\responses']);
mkdir([retroot 'g2d\analytic\scales']);


for ii = 1:20
    %load retinogram and merge RGB channels
    ret = u_load([asymmetryroot 'data\retinograms\DRIVE\test\images_extended\' zerostr(ii,2) '_test_ext.mat']);
    ret = rgb2gray(ret);
    

    %Compute repsonses for mono and save
    [response_map ori_map scale_map norm_response_map] = gaussian_2nd_derivative_line(ret, [0.5 1 2 4 8]);
    save_uint8([retroot 'g2d\analytic\orientations\' zerostr(ii,3) '_test_ori.mat'], ori_map);
    save_uint8([retroot 'g2d\analytic\responses\' zerostr(ii,3) '_test_response.mat'], response_map);
    save_uint8([retroot 'g2d\analytic\responses\' zerostr(ii,3) '_test_response_n.mat'], norm_response_map);
    save_uint8([retroot 'g2d\analytic\scales\' zerostr(ii,3) '_test_scale.mat'], scale_map);
    

end
%%
retroot = [asymmetryroot 'data\retinograms\DRIVE\test\'];
load([retroot '\orientations\all_gt_orientations.mat']);

    
predicted_orientations = [];
predicted_orientationsd = [];
orientation_responses1 = [];
orientation_responses2 = [];
orientation_responsesd = [];
for jj = 1:20

    %load foveal mask and vessel mask and indices to vessel centres
    foveal_mask = u_load([retroot 'foveal_masks\' zerostr(jj,2) '_test_f_mask']);
    vessel_mask = u_load([retroot 'vessel_masks\' zerostr(jj,2) '_test_v_mask.mat']);
    vessel_mask(~foveal_mask) = false;

    %load in orientation map
    ori_map = load_uint8([retroot 'predictions\orientations\g2d\analytic\orientations\' zerostr(jj,3) '_test_ori.mat']);
    ori_map(~vessel_mask) = NaN;
    
    ori_mapd = load_uint8([retroot 'predictions\orientation\g2d\rf_3\' zerostr(jj,2) '_test_ext_class.mat']);
    ori_mapd(~vessel_mask) = NaN;


    %load the response map
    response_map1 = load_uint8([retroot 'predictions\orientations\g2d\analytic\responses\' zerostr(jj,3) '_test_response.mat']);
    response_map2 = load_uint8([retroot 'predictions\orientations\g2d\analytic\responses\' zerostr(jj,3) '_test_response_n.mat']);
    response_mapd = abs(ori_mapd);
    
    %convert the orientation to complex form
    ori_map = exp(complex(0, 2*ori_map));

    %Get the orientation predictions and responses
    predicted_orientations = [predicted_orientations; ori_map(vessel_mask)]; %#ok
    predicted_orientationsd = [predicted_orientationsd; ori_mapd(vessel_mask)]; %#ok
    orientation_responses1 = [orientation_responses1; response_map1(vessel_mask)]; %#ok
    orientation_responses2 = [orientation_responses2; response_map2(vessel_mask)]; %#ok
    orientation_responsesd = [orientation_responsesd; response_mapd(vessel_mask)]; %#ok


end
%Now compute the orientation errors for the set
[prediction_errs_rad] = ori_error(predicted_orientations, gt_orientations);
[prediction_errs_radd] = ori_error(predicted_orientationsd, gt_orientations);
%%
orientation_responses1(isnan(gt_orientations)) = [];
orientation_responses2(isnan(gt_orientations)) = [];
orientation_responsesd(isnan(gt_orientations)) = [];
orientation_responsesr = orientation_responses1 ./ orientation_responses2;

orientation_responses1a = abs(orientation_responses1);
orientation_responses2a = abs(orientation_responses2);
orientation_responsesra = abs(orientation_responsesr);
orientation_responses1a = orientation_responses1a / max(orientation_responses1a);
orientation_responses2a = orientation_responses2a / max(orientation_responses2a);
orientation_responsesra = orientation_responsesra / max(orientation_responsesra);

orientation_responses1s = (orientation_responses1 - min(orientation_responses1)) / (max(orientation_responses1) - min(orientation_responses1));
orientation_responses2s = (orientation_responses2 - min(orientation_responses2)) / (max(orientation_responses2) - min(orientation_responses2));
orientation_responsesrs = (orientation_responsesr - min(orientation_responsesr)) / (max(orientation_responsesr) - min(orientation_responsesr));

prediction_errs = 180*abs(prediction_errs_rad)/pi;
prediction_errsd = 180*abs(prediction_errs_radd)/pi;


[responses_x1a, responses_y1a] = ...
    kernel_smoother(orientation_responses1a, prediction_errs, 100);
[responses_x2a, responses_y2a] = ...
    kernel_smoother(orientation_responses2a, prediction_errs, 100);
[responses_xra, responses_yra] = ...
    kernel_smoother(orientation_responsesra, prediction_errs, 100);
[responses_x1s, responses_y1s] = ...
    kernel_smoother(orientation_responses1s, prediction_errs, 100);
[responses_x2s, responses_y2s] = ...
    kernel_smoother(orientation_responses2s, prediction_errs, 100);
[responses_xrs, responses_yrs] = ...
    kernel_smoother(orientation_responsesrs, prediction_errs, 100);
[responses_xd, responses_yd] = ...
    kernel_smoother(orientation_responsesd, prediction_errsd, 100);
%
figure; hold all;
plot(responses_x1a, responses_y1a);
plot(responses_x2a, responses_y2a);
plot(responses_xra, responses_yra);
plot(responses_x1s, responses_y1s);
plot(responses_x2s, responses_y2s);
plot(responses_xrs, responses_yrs);
plot(responses_xd, responses_yd);
legend({...
    '|\lambda_1|',...
    '|\lambda_2|',...
    '|\lambda_1/\lambda_2|',...
    '\lambda_1',...
    '\lambda_2',...
    '\lambda_1/\lambda_2',...
    'D'});
title('Kernel estimate of Mean Absolute Error in orientation estimate for varying combinations of G" filter responses');
xlabel('Filter response, scaled between 0 - 1');
ylabel('MAE (degrees)');
figure; hist([orientation_responses1a orientation_responsesd], 100);
legend({'|\lambda_1|', 'D'});
title('Histogram of filter responses at vessel pixels');
%%
figure; hist(abs(orientation_responses1)-abs(orientation_responses2),100);
figure; hist(abs(orientation_responses2), 100);
%%
figure; plot3(orientation_responses1, orientation_responses2, prediction_errs, 'r.');

