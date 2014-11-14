[g,dg,ddg] = gaussian_filters_1d(4,5*4);

gxy1 = (g'*ddg) - (ddg'*g);
gxy2 = 2*(dg'*dg);

ddgm = ddg - mean(ddg);
dgm = dg - mean(dg);

gn	= g / sum(abs(g));
dgn	= dgm / sum(abs(dgm));
ddgn	= ddgm / sum(abs(ddgm));

gxym1 = (g'*ddgm) - (ddgm'*g);
gxym2 = 2*dg'*dgm;

gxyn1 = (gn'*ddgn) - (ddgn'*gn);
gxyn2 = 2*(dgn'*dgn);

[gs,dgs,ddgs] = gaussian_filters_1d(4,3*4);

gxys1 = (gs'*ddgs) - (ddgs'*gs);
gxys2 = 2*(dgs'*dgs);

figure; hold all; plot(g); plot(dg); plot(ddg); title('Unormalised');
figure; hold all; plot(g); plot(dgm); plot(ddgm); title('Shifted');
figure; hold all; plot(gn); plot(dgn); plot(ddgn); title('Shifted and normalised');

figure; imgray(gxy1);
figure; imgray(gxy2);

figure; hold all; title('Unormalised');
plot(-20:20, gxy1(21,:));
plot(sqrt(2)*(-20:20), diag(gxy2));

figure; hold all; title('Shifted');
plot(-20:20, gxym1(21,:));
plot(sqrt(2)*(-20:20), diag(gxym2));

figure; hold all; title('Shifted and normalised');
plot(-20:20, gxyn1(21,:));
plot(sqrt(2)*(-20:20), diag(gxyn2));

figure; hold all; title('Un-shifted, short support');
plot(-12:12, gxys1(13,:));
plot(sqrt(2)*(-12:12), diag(gxys2)); set(gca, 'xlim', [-30 30]);
%%
retroot = [asymmetryroot('shared'),'data\retinograms\DRIVE\test\'];
load([retroot,'orientations\all_gt_orientations.mat']);
% pts_im_list = [];
% for ii = 1:20
%     pts_im_list = [pts_im_list; ii*ones(pts_per_image(ii),1)]; %#ok;
% end
% pts_im_list(isnan(gt_orientations)) = [];
discard = isnan(gt_orientations);
gt_orientations(discard) = [];
num_pts = size(gt_orientations,1);

plot_orientations = angle(gt_orientations)/2;

errors_all_scales = zeros(num_pts,4);
responses_all_scales = zeros(num_pts,4);

for ss = 1:4;
    predicted_orientations = [];
    responses = [];
    for jj = 1:20

        sigma = 2^(ss-1);
        
        ret = u_load([asymmetryroot('shared'),'data\retinograms\DRIVE\test\images_extended\' zerostr(jj,2) '_test_ext.mat']);
        ret = rgb2gray(ret);

        %Compute repsonses for mono and save
        [response_map, ori_map] = gaussian_2nd_derivative_line(ret, sigma);

        %convert the orientation to complex form
        ori_map = exp(complex(0, 2*ori_map));

        %load foveal mask and vessel mask and indices to vessel centres
        foveal_mask = u_load([retroot 'foveal_masks\' zerostr(jj,2) '_test_f_mask']);
        vessel_mask = u_load([retroot 'vessel_masks\' zerostr(jj,2) '_test_v_mask.mat']);
        vessel_mask(~foveal_mask) = false;

        %Get the orientation predictions and responses
        predicted_orientations = [predicted_orientations; ori_map(vessel_mask)]; %#ok
        responses = [responses; response_map(vessel_mask)]; %#ok

    end
    predicted_orientations(discard) = []; %#ok
    responses(discard) = []; %#ok
    
    %Now compute the orientation errors for the set
    [prediction_errs all_stats] = ori_error(predicted_orientations, gt_orientations);
    errors_all_scales(:,ss) = prediction_errs;
    responses_all_scales(:,ss) = responses;
    
    %display the results
    display(['Median angular error: ' num2str(median(abs(prediction_errs))*180/pi)]);

    
    figure; plot(plot_orientations , prediction_errs, 'b.');
    title(['\sigma = ' num2str(sigma) ', median angular error: ' num2str(median(abs(prediction_errs))*180/pi)]);
end

rows = (1:num_pts)';

[~, optimal_scale] = min(abs(errors_all_scales),[],2);
[~, response_scale] = max(abs(responses_all_scales), [], 2);
[~, sigma_scale] = max(bsxfun(@times, abs(responses_all_scales), [1 2 4 8].^2), [], 2);%

optimal_errors = errors_all_scales( sub2ind([num_pts 4], rows, optimal_scale) );
response_errors = errors_all_scales( sub2ind([num_pts 4], rows, response_scale) );
sigma_errors = errors_all_scales( sub2ind([num_pts 4], rows, sigma_scale) );

figure; plot(plot_orientations , optimal_errors, 'b.');
title(['Optimal scale, median angular error: ' num2str(median(abs(optimal_errors))*180/pi)]);

figure; plot(plot_orientations , response_errors, 'b.');
title(['Response scale - valid, median angular error: ' num2str(median(abs(response_errors))*180/pi)]);

figure; plot(plot_orientations , sigma_errors, 'b.');
title(['Response scale - fudged, median angular error: ' num2str(median(abs(sigma_errors))*180/pi)]);

% %Compute bootstrap errors
% median_error_bt = zeros(2000,1);
% for jj = 1:2000
%     boot_samp = ceil(20*rand(20,1));
%     pts_idx = ismember(pts_im_list, boot_samp);
%     median_error_bt(jj) = median(abs(prediction_errs(pts_idx)));
% end
% median_ci = 180*grpstats(median_error_bt, ones(2000,1), 'meanci')/pi;
% 
% %display the results
% display([pred_decomp{ii} ' ' pred_type{ii} ': ' num2str(median(abs(prediction_errs))*180/pi) ...
%     ' (' num2str(median_ci(1)) ', ' num2str(median_ci(2)) ')' ]);

