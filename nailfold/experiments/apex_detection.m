%--------------------------------------------------------------------------
% Script demonstrating template matching of vessel apexes in nailfold
% images
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%% 1) Create template from the training vessels
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% a) Basic mean template, no aligment of vessel shapes
load('C:\isbe\nailfold\playground\model_experiments\vessel_names.mat', 'train_names');

num_vessels = length(train_names);
gd1_patches = zeros(33,33,num_vessels);
gd2_patches = zeros(33,33,num_vessels);
img_patches = zeros(33,33,num_vessels);

for ii = 1:num_vessels
       
    %Sample patch from 1st deriv.
    mag_1d = u_load(['C:\isbe\nailfold\playground\model_experiments\images\g1d\vessel_top' train_names{ii} '.mat']);

    %Sample patch from 2nd deriv.
    mag_2d = u_load(['C:\isbe\nailfold\playground\model_experiments\images\g2d\vessel_top' train_names{ii} '.mat']);
    
    %Sample patch from 2nd deriv.
    img_patch = imread(['C:\isbe\nailfold\playground\model_experiments\images\orig\vessel_top' train_names{ii} '.png']);

    %Load in pts
    v_pts = u_load(['C:\isbe\nailfold\playground\model_experiments\points_m\vessel_top' train_names{ii} '.mat']);

    gd1_patches(:,:,ii) = abs(sample_window(mag_1d, 33, round(v_pts(16,2)), round(v_pts(16,1))));
    gd2_patches(:,:,ii) = sample_window(mag_2d, 33, round(v_pts(16,2)), round(v_pts(16,1)));
    img_patches(:,:,ii) = sample_window(img_patch, 33, round(v_pts(16,2)), round(v_pts(16,1)));
end

%Compute the mean of the 2 sets of patches to use as templates
apex_template1 = mean(gd1_patches,3);
apex_template2 = mean(gd2_patches,3);
apex_templatei = mean(img_patches,3);

%save the templates for later use
save([nailfoldroot 'data\apex_detection\apex_templates_unaligned.mat'], 'apex_template*');

figure; imagesc(apex_template1); axis image;
title('Template for 1st derivative patches');
figure; imagesc(apex_template2); axis image;
title('Template for 2nd derivative patches');
figure; imagesc(apex_templatei); axis image;
title('Template for image patches');
%%
%--------------------------------------------------------------------------
% b) Now align the shapes, allowing rotation and scaling
load('C:\isbe\nailfold\playground\model_experiments\vessel_names.mat', 'train_names');

num_vessels = length(train_names);
u_shapes = zeros(num_vessels, 62);
areas = zeros(num_vessels, 1);

for ii = 1:num_vessels
    %Load in pts
    v_pts = u_load(['C:\isbe\nailfold\playground\model_experiments\points_m\vessel_top' train_names{ii} '.mat']);
    u_shapes(ii,:) = v_pts(:)';
    areas(ii) = polyarea(v_pts(:,1), v_pts(:,2));
end
[a_shapes, a_scales, mean_target, a_rots, a_trans a_origins] = align_shapes(u_shapes, 'area', mean(areas));
mean_shape = [mean(a_shapes(:,1:31))', mean(a_shapes(:,32:62))'];

% v = mean_shape(1,:) - mean_shape(end,:);
% theta = atan(v(2)/v(1));
theta = -mean(acos(squeeze(a_rots(1,1,:))));
theta_rot = [cos(theta) sin(theta); -sin(theta) cos(theta)];

figure; hold on;
plot(a_shapes(:,1:31)', a_shapes(:,32:62)'); axis ij equal;
plot(mean_shape(:,1), mean_shape(:,2), 'k', 'linewidth', 2)

mean_shape = mean_target * theta_rot';
%mean_shape = [mean_shape(:,1) - (mean_shape(end,1)+mean_shape(1,1))/2 mean_shape(:,2)-mean_shape(16,2)];
mean_shape = [mean_shape(:,1) - mean_shape(16,1) mean_shape(:,2)-mean_shape(16,2)];
figure; plot(mean_shape(:,1), mean_shape(:,2)); axis ij equal;
%
% b1) Extract image patches for each apex given rotation and scaling to
% mean

x = repmat(-24:24, 49, 1);
y = repmat(-24:24, 49, 1)';
xy = [x(:) y(:)];

gd1_patches_a = zeros(49,49,num_vessels);
gd2_patches_a = zeros(49,49,num_vessels);
img_patches_a = zeros(49,49,num_vessels);

for ii = 1:num_vessels
    %Load in pts
    v_pts = u_load(['C:\isbe\nailfold\playground\model_experiments\points_m\vessel_top' train_names{ii} '.mat']);
    
    %Sample patch from 1st deriv.
    mag_1d = abs(u_load(['C:\isbe\nailfold\playground\model_experiments\images\g1d\vessel_top' train_names{ii} '.mat']));

    %Sample patch from 1st deriv.
    mag_2d = u_load(['C:\isbe\nailfold\playground\model_experiments\images\g2d\vessel_top' train_names{ii} '.mat']);
    
    %Sample patch from 2nd deriv.
    img_patch = double(imread(['C:\isbe\nailfold\playground\model_experiments\images\orig\vessel_top' train_names{ii} '.png']));
    
    xya = (xy * a_rots(:,:,ii)' / a_scales(ii))* theta_rot;
    xa = reshape(xya(:,1) + v_pts(16,1), 49, 49);
    ya = reshape(xya(:,2) + v_pts(16,2), 49, 49);
    
    gd1_patches_a(:,:,ii) = interp2(mag_1d, xa, ya, 'bilinear');
    gd2_patches_a(:,:,ii) = interp2(mag_2d, xa, ya, 'bilinear');
    img_patches_a(:,:,ii) = interp2(img_patch, xa, ya, 'bilinear');
    
    if ii < 0
        figure;
        subplot(2,2,1); imgray(mag_1d); 
        plot(xa, ya, 'rx'); plot(v_pts(:,1), v_pts(:,2), 'g'); 
        plot(v_pts(16,1), v_pts(16,2), 'co');
        
        subplot(2,2,2); imgray(mag_2d); 
        plot(xa, ya, 'rx'); plot(v_pts(:,1), v_pts(:,2), 'g'); 
        plot(v_pts(16,1), v_pts(16,2), 'co');
        
        subplot(2,2,3); imgray(gd1_patches_a(:,:,ii));
        plot(mean_shape(:,1)+24, mean_shape(:,2)+12);
        subplot(2,2,4); imgray(gd2_patches_a(:,:,ii));
        plot(mean_shape(:,1)+24, mean_shape(:,2)+12);
    end

end

%Compute the mean of the 2 sets of patches to use as templates
circle_mask = x.^2 + y.^2 > 24^2;
apex_template1_r = zeros(49,49, 11);
apex_template2_r = zeros(49,49, 11);
apex_template1_a = mean(gd1_patches_a,3);
apex_template2_a = mean(gd2_patches_a,3);
apex_templatei_a = mean(img_patches_a,3);

apex_template1_a(circle_mask) = NaN;
apex_template2_a(circle_mask) = NaN;
apex_templatei_a(circle_mask) = NaN;
%
% b2) Save rotated copies of the template
figure; 
subplot(1,2,1); imagesc(apex_templatei_a); axis image; hold on;
plot(mean_shape(:,1)+24, mean_shape(:,2)+24, 'k', 'linewidth', 2);
title('Template for 1st derivative patches');
subplot(1,2,2); imagesc(apex_template2_a); axis image; hold on;
plot(mean_shape(:,1)+24, mean_shape(:,2)+24, 'k', 'linewidth', 2);
title('Template for 2nd derivative patches');
    
    
for ii = 1:11
    theta = (ii - 6)*6;
    at1a = imrotate(apex_template1_a, theta, 'bilinear', 'crop');
    at1a(circle_mask) = NaN;
    at2a = imrotate(apex_template2_a, theta, 'bilinear', 'crop');
    at2a(circle_mask) = NaN;
    
    apex_template1_r(:,:,ii) = at1a;
    apex_template2_r(:,:,ii) = at2a;

end

%save the templates for later use
save([nailfoldroot 'data\apex_detection\mean_shape.mat'], 'mean_shape');
save([nailfoldroot 'data\apex_detection\apex_templates_aligned.mat'], 'apex_template*a');
save([nailfoldroot 'data\apex_detection\apex_templates_rotated.mat'], 'apex_template*r');

%%
%--------------------------------------------------------------------------
%% 2) Apply the templates to test images
%--------------------------------------------------------------------------

%-------------------------------------------------------------------------- 
% a) test unaligned mean template, NOT allowing the template
% to rotate or scale
load([nailfoldroot 'data\apex_detection\apex_templates_unaligned.mat'], 'apex_template*');
nf_files = dir('C:\isbe\nailfold\images\anonymous_oct\annotations_qt\*.txt');

all_candidate_labels12 = [];
all_candidate_labels1 = [];
all_candidate_labels2 = [];
all_candidate_labelsi = [];
all_candidate_vals12 = [];
all_candidate_vals1 = [];
all_candidate_vals2 = [];
all_candidate_valsi = [];
num_gt = 0;

for nn = 1:6%7:12
    %LOad in nailfold
    nailfold = imread(['C:\isbe\nailfold\images\anonymous_oct\bmp\' nf_files(nn).name(1:end-11) '.bmp']);   
    nailfold = nailfold(:,:,1);
    nailfold_mask = imerode(nailfold > 10 & nailfold < 250, strel('disk', 20));
    
    %Compute the gaussian derivatives
    [mag_2d] = gaussian_2nd_derivative_line(nailfold, 4);
    [mag_1d] = abs(gaussian_1st_derivative_gradient(nailfold, 2));

    %Apply normalised cross correlations using the 2 templates
    C1 = mb_normxcorr2(apex_template1, mag_1d);
    C2 = mb_normxcorr2(apex_template2, mag_2d);
    Ci = mb_normxcorr2(apex_templatei, nailfold);

    %Look for local maxima in the 2 independent NCC maps
    [maxima_pos1, maxima_vals1] = local_image_maxima(C1, 20, nailfold_mask, 0);
    [maxima_pos2, maxima_vals2] = local_image_maxima(C2, 20, nailfold_mask, 0);
    [maxima_posi, maxima_valsi] = local_image_maxima(Ci, 20, nailfold_mask, 0);

    %Now multiply the 2 NCC maps to get a single map and again look for local
    %maxima
    C12 = C1 .* C2;
    C12(~nailfold_mask) = 0;
    [maxima_pos12, maxima_vals12] = local_image_maxima(C12, 20, [], 0);
         
%     %Plot the maxima in the nailfold image, with marker size scaled according
%     %to NCC value
%     n_min = min(nailfold(nailfold_mask));
%     n_max = max(nailfold(nailfold_mask)); 
%     figure; imagesc(nailfold); axis image off; hold on; colormap(gray(256)); caxis([n_min n_max]);
%     for ii = 1:length(maxima_vals12)
%         ms = 1 + 20*(maxima_vals12(ii));
%         plot(maxima_pos12(ii,1), maxima_pos12(ii,2), 'r+', 'MarkerSize', ms);
%     end
    
    %Load in ground truth apex locations
    apexes = u_load(['C:\isbe\nailfold\data\apexes\' nf_files(nn).name(1:end-4) '_apexes.mat']);
    apex_gt = zeros(length(apexes),2);
    for ii = 1:length(apexes)
        apex_gt(ii,:) = apexes{ii}(16,:);
    end
    num_gt = num_gt + length(apexes);
    
    %Workout whether candidate are hits or not
    [candidate_labels12, gt_hits12] = evaluate_apex_candidates(apex_gt, maxima_pos12, 15);
    [candidate_labels1] = evaluate_apex_candidates(apex_gt, maxima_pos1, 15);
    [candidate_labels2] = evaluate_apex_candidates(apex_gt, maxima_pos2, 15);
    [candidate_labelsi] = evaluate_apex_candidates(apex_gt, maxima_posi, 15);
    
    all_candidate_labels12 = [all_candidate_labels12; candidate_labels12]; %#ok
    all_candidate_labels1 = [all_candidate_labels1; candidate_labels1]; %#ok
    all_candidate_labels2 = [all_candidate_labels2; candidate_labels2]; %#ok
    all_candidate_labelsi = [all_candidate_labelsi; candidate_labelsi]; %#ok
    
    all_candidate_vals12 = [all_candidate_vals12; maxima_vals12]; %#ok
    all_candidate_vals1 = [all_candidate_vals1; maxima_vals1]; %#ok
    all_candidate_vals2 = [all_candidate_vals2; maxima_vals2]; %#ok
    all_candidate_valsi = [all_candidate_valsi; maxima_valsi]; %#ok
    
%     plot(apex_gt(gt_hits12,1), apex_gt(gt_hits12,2), 'gx');
%     plot(apex_gt(~gt_hits12,1), apex_gt(~gt_hits12,2), 'cx');
%     plot(maxima_pos12(candidate_labels12,1), maxima_pos12(candidate_labels12,2), 'g+');
    
    display([sum(gt_hits12), sum(~gt_hits12)]);
    
    save([nailfoldroot 'data\apex_detection\unaligned\' nf_files(nn).name(1:end-4) '_corr.mat'], 'C1', 'C2');
end
save([nailfoldroot 'data\apex_detection\unaligned\all_candidates.mat'], 'all_candidate_*');
%%
[~, idx12] = sort(all_candidate_vals12, 'descend');
[~, idx1] = sort(all_candidate_vals1, 'descend');
[~, idx2] = sort(all_candidate_vals2, 'descend');
[~, idxi] = sort(all_candidate_valsi, 'descend');

tp_count12 = cumsum(all_candidate_labels12(idx12));
tp_count1 = cumsum(all_candidate_labels1(idx1));
tp_count2 = cumsum(all_candidate_labels2(idx2));
tp_counti = cumsum(all_candidate_labelsi(idxi));

fp_count12 = (1:length(all_candidate_labels12))' - tp_count12;
fp_count1 = (1:length(all_candidate_labels1))' - tp_count1;
fp_count2 = (1:length(all_candidate_labels2))' - tp_count2;
fp_counti = (1:length(all_candidate_labelsi))' - tp_counti;

sensitivity12 = tp_count12 / num_gt;
sensitivity1 = tp_count1 / num_gt;
sensitivity2 = tp_count2 / num_gt;
sensitivityi = tp_counti / num_gt;

%Find out what threshold gives at least 95% sensitivity
thresh12 = all_candidate_vals12(idx12( sum(sensitivity12 > 0.95) + 1 ));
thresh1 = all_candidate_vals1(idx1( sum(sensitivity1 > 0.95) + 1 ));
thresh2 = all_candidate_vals2(idx2( sum(sensitivity2 > 0.95) + 1 ));
save([nailfoldroot 'data\apex_detection\unaligned\thresh.mat'], 'thresh*');

figure; hold all; title('Template matching - unaligned shapes');
plot(log10(fp_count12), sensitivity12);
plot(log10(fp_count1), sensitivity1);
plot(log10(fp_count2), sensitivity2);
plot(log10(fp_counti), sensitivityi);
legend({'G'' * G"', 'G''', 'G"', 'Image'}, 'location', 'southeast');
set(gca, 'ylim', [0 1]);

% figure; hold all; title('Template matching - unaligned shapes');
% plot(fp_count12 / num_gt, sensitivity12);
% plot(fp_count1 / num_gt, sensitivity1);
% plot(fp_count2 / num_gt, sensitivity2);
% legend({'G'' * G"', 'G''', 'G"'}, 'location', 'southeast');
%% ------------------------------------------------------------------------
% b) test aligned mean template, NOT allowing the template
% to rotate or scale
load([nailfoldroot 'data\apex_detection\apex_templates_aligned.mat'], 'apex_template*');
nf_files = dir('C:\isbe\nailfold\images\anonymous_oct\annotations_qt\*.txt');

all_candidate_labels12 = [];
all_candidate_labels1 = [];
all_candidate_labels2 = [];
all_candidate_labelsi = [];
all_candidate_vals12 = [];
all_candidate_vals1 = [];
all_candidate_vals2 = [];
all_candidate_valsi = [];
num_gt = 0;

for nn = 1:6%7:12
    %LOad in nailfold
    nailfold = imread(['C:\isbe\nailfold\images\anonymous_oct\bmp\' nf_files(nn).name(1:end-11) '.bmp']);   
    nailfold = nailfold(:,:,1);
    nailfold_mask = imerode(nailfold > 10 & nailfold < 250, strel('disk', 20));
    
    %Compute the gaussian derivatives
    [mag_2d] = gaussian_2nd_derivative_line(nailfold, 4);
    [mag_1d] = abs(gaussian_1st_derivative_gradient(nailfold, 2));

    %Apply normalised cross correlations using the 2 templates
    C1 = mb_normxcorr2(apex_template1_a, mag_1d);
    C2 = mb_normxcorr2(apex_template2_a, mag_2d);
    Ci = mb_normxcorr2(apex_templatei_a, nailfold);

    %Look for local maxima in the 2 independent NCC maps
    [maxima_pos1, maxima_vals1] = local_image_maxima(C1, 20, nailfold_mask, 0);
    [maxima_pos2, maxima_vals2] = local_image_maxima(C2, 20, nailfold_mask, 0);
    [maxima_posi, maxima_valsi] = local_image_maxima(Ci, 20, nailfold_mask, 0);

    %Now multiply the 2 NCC maps to get a single map and again look for local
    %maxima
    C12 = C1 .* C2;
    C12(~nailfold_mask) = 0;
    [maxima_pos12, maxima_vals12] = local_image_maxima(C12, 20, [], 0);
    
%     maxima_pos1(:,2) = maxima_pos1(:,2) - 12;
%     maxima_pos2(:,2) = maxima_pos2(:,2) - 12;
%     maxima_pos12(:,2) = maxima_pos12(:,2) - 12;
         
    %Plot the maxima in the nailfold image, with marker size scaled according
    %to NCC value
    n_min = min(nailfold(nailfold_mask));
    n_max = max(nailfold(nailfold_mask)); 
    figure; imagesc(nailfold); axis image off; hold on; colormap(gray(256)); caxis([n_min n_max]);
    for ii = 1:length(maxima_vals12)
        ms = 1 + 20*(maxima_vals12(ii));
        plot(maxima_pos12(ii,1), maxima_pos12(ii,2), 'r+', 'MarkerSize', ms);
    end
    
    %Load in ground truth apex locations
    apexes = u_load(['C:\isbe\nailfold\data\apexes\' nf_files(nn).name(1:end-4) '_apexes.mat']);
    apex_gt = zeros(length(apexes),2);
    for ii = 1:length(apexes)
        apex_gt(ii,:) = apexes{ii}(16,:);
    end
    num_gt = num_gt + length(apexes);
    
    %Workout whether candidate are hits or not
    [candidate_labels12, gt_hits12] = evaluate_apex_candidates(apex_gt, maxima_pos12, 15);
    [candidate_labels1] = evaluate_apex_candidates(apex_gt, maxima_pos1, 15);
    [candidate_labels2] = evaluate_apex_candidates(apex_gt, maxima_pos2, 15);
    [candidate_labelsi] = evaluate_apex_candidates(apex_gt, maxima_posi, 15);
    
    all_candidate_labels12 = [all_candidate_labels12; candidate_labels12]; %#ok
    all_candidate_labels1 = [all_candidate_labels1; candidate_labels1]; %#ok
    all_candidate_labels2 = [all_candidate_labels2; candidate_labels2]; %#ok
    all_candidate_labelsi = [all_candidate_labelsi; candidate_labelsi]; %#ok
    
    all_candidate_vals12 = [all_candidate_vals12; maxima_vals12]; %#ok
    all_candidate_vals1 = [all_candidate_vals1; maxima_vals1]; %#ok
    all_candidate_vals2 = [all_candidate_vals2; maxima_vals2]; %#ok
    all_candidate_valsi = [all_candidate_valsi; maxima_valsi]; %#ok
    
    plot(apex_gt(gt_hits12,1), apex_gt(gt_hits12,2), 'gx');
    plot(apex_gt(~gt_hits12,1), apex_gt(~gt_hits12,2), 'cx');
    plot(maxima_pos12(candidate_labels12,1), maxima_pos12(candidate_labels12,2), 'g+');
    
    display([sum(gt_hits12), sum(~gt_hits12)]);
    
    save([nailfoldroot 'data\apex_detection\aligned\' nf_files(nn).name(1:end-4) '_corr.mat'], 'C1', 'C2');
end
save([nailfoldroot 'data\apex_detection\aligned\all_candidates.mat'], 'all_candidate_*');
%%
%num_gt = 111;
load([nailfoldroot 'data\apex_detection\aligned\all_candidates.mat'], 'all_candidate_*');
[~, idx12] = sort(all_candidate_vals12, 'descend');
[~, idx1] = sort(all_candidate_vals1, 'descend');
[~, idx2] = sort(all_candidate_vals2, 'descend');
[~, idxi] = sort(all_candidate_valsi, 'descend');

tp_count12 = cumsum(all_candidate_labels12(idx12));
tp_count1 = cumsum(all_candidate_labels1(idx1));
tp_count2 = cumsum(all_candidate_labels2(idx2));
tp_counti = cumsum(all_candidate_labelsi(idxi));

fp_count12 = (1:length(all_candidate_labels12))' - tp_count12;
fp_count1 = (1:length(all_candidate_labels1))' - tp_count1;
fp_count2 = (1:length(all_candidate_labels2))' - tp_count2;
fp_counti = (1:length(all_candidate_labelsi))' - tp_counti;

sensitivity12 = tp_count12 / num_gt;
sensitivity1 = tp_count1 / num_gt;
sensitivity2 = tp_count2 / num_gt;
sensitivityi = tp_counti / num_gt;

%Find out what threshold gives at least 95% sensitivity
thresh12 = all_candidate_vals12(idx12( sum(sensitivity12 > 0.95) + 1 ));
thresh1 = all_candidate_vals1(idx1( sum(sensitivity1 > 0.95) + 1 ));
thresh2 = all_candidate_vals2(idx2( sum(sensitivity2 > 0.95) + 1 ));
%save([nailfoldroot 'data\apex_detection\aligned\thresh.mat'], 'thresh*');

figure; hold all; title('Template matching - unaligned shapes');
plot(log10(fp_count12), sensitivity12);
plot(log10(fp_count1), sensitivity1);
plot(log10(fp_count2), sensitivity2);
plot(log10(fp_counti), sensitivityi);
legend({'G'' * G"', 'G''', 'G"', 'Image'}, 'location', 'southeast');
set(gca, 'ylim', [0 1]);

%%
% c) test aligned template, allowing template to rotate but NOT scale
load([nailfoldroot 'data\apex_detection\apex_templates_rotated.mat'], 'apex_template*');
nf_files = dir('C:\isbe\nailfold\images\anonymous_oct\annotations_qt\*.txt');

all_candidate_labels12 = [];
all_candidate_labels1 = [];
all_candidate_labels2 = [];
all_candidate_vals12 = [];
all_candidate_vals1 = [];
all_candidate_vals2 = [];
num_gt = 0;

for nn = 7:12
    %LOad in nailfold
    nailfold = imread(['C:\isbe\nailfold\images\anonymous_oct\bmp\' nf_files(nn).name(1:end-11) '.bmp']);   
    nailfold = nailfold(:,:,1);
    nailfold_mask = imerode(nailfold > 10 & nailfold < 250, strel('disk', 20));
    
%     %Compute the gaussian derivatives
%     [mag_2d] = gaussian_2nd_derivative_line(nailfold, 4);
%     [mag_1d] = abs(gaussian_1st_derivative_gradient(nailfold, 2));
% 
%     %Apply normalised cross correlations using the 2 templates
%     C1 = -inf(size(nailfold));
%     C2 = -inf(size(nailfold));
%     for ii = 1:11
%     
%         c1 = mb_normxcorr2(apex_template1_r(:,:,ii), mag_1d);
%         c2 = mb_normxcorr2(apex_template2_r(:,:,ii), mag_2d);
%         
%         C1 = max(C1, c1);
%         C2 = max(C1, c2);
%     end
    
    load([nailfoldroot 'data\apex_detection\aligned_rotations\' nf_files(nn).name(1:end-4) '_corr.mat'], 'C1', 'C2');
    
    %Look for local maxima in the 2 independent NCC maps
    [maxima_pos1, maxima_vals1] = local_image_maxima(C1, 33, nailfold_mask, 0.1);
    [maxima_pos2, maxima_vals2] = local_image_maxima(C2, 33, nailfold_mask, 0.1);

    %Now multiply the 2 NCC maps to get a single map and again look for local
    %maxima
    C12 = C1 .* C2;
    C12(~nailfold_mask) = 0;
    [maxima_pos12, maxima_vals12] = local_image_maxima(C12, 33, [], 0.1);
         
    %Plot the maxima in the nailfold image, with marker size scaled according
    %to NCC value
    n_min = min(nailfold(nailfold_mask));
    n_max = max(nailfold(nailfold_mask)); 
    figure; imagesc(nailfold); axis image off; hold on; colormap(gray(256)); caxis([n_min n_max]);
    for ii = 1:length(maxima_vals12)
        ms = 1 + 20*(maxima_vals12(ii)-0.1);
        plot(maxima_pos12(ii,1), maxima_pos12(ii,2), 'r+', 'MarkerSize', ms);
    end
    
    %Load in ground truth apex locations
    apexes = u_load(['C:\isbe\nailfold\data\apexes\' nf_files(nn).name(1:end-4) '_apexes.mat']);
    apex_gt = zeros(length(apexes),2);
    for ii = 1:length(apexes)
        apex_gt(ii,:) = apexes{ii}(16,:);
    end
    num_gt = num_gt + length(apexes);
    
    %Workout whether candidate are hits or not
    [candidate_labels12] = evaluate_apex_candidates(apex_gt, maxima_pos12, 15);
    [candidate_labels1] = evaluate_apex_candidates(apex_gt, maxima_pos1, 15);
    [candidate_labels2] = evaluate_apex_candidates(apex_gt, maxima_pos2, 15);
    
    all_candidate_labels12 = [all_candidate_labels12; candidate_labels12]; %#ok
    all_candidate_labels1 = [all_candidate_labels1; candidate_labels1]; %#ok
    all_candidate_labels2 = [all_candidate_labels2; candidate_labels2]; %#ok
    
    all_candidate_vals12 = [all_candidate_vals12; maxima_vals12]; %#ok
    all_candidate_vals1 = [all_candidate_vals1; maxima_vals1]; %#ok
    all_candidate_vals2 = [all_candidate_vals2; maxima_vals2]; %#ok
    
    %save([nailfoldroot 'data\apex_detection\aligned_rotations\' nf_files(nn).name(1:end-4) '_corr.mat'], 'C1', 'C2');
end
save([nailfoldroot 'data\apex_detection\aligned_rotations\all_candidates.mat'], 'all_candidate_*');
%
[~, idx12] = sort(all_candidate_vals12, 'descend');
[~, idx1] = sort(all_candidate_vals1, 'descend');
[~, idx2] = sort(all_candidate_vals2, 'descend');

tp_count12 = cumsum(all_candidate_labels12(idx12));
tp_count1 = cumsum(all_candidate_labels1(idx1));
tp_count2 = cumsum(all_candidate_labels2(idx2));

fp_count12 = (1:length(all_candidate_labels12))' - tp_count12;
fp_count1 = (1:length(all_candidate_labels1))' - tp_count1;
fp_count2 = (1:length(all_candidate_labels2))' - tp_count2;

sensitivity12 = tp_count12 / num_gt;
sensitivity1 = tp_count1 / num_gt;
sensitivity2 = tp_count2 / num_gt;

figure; hold all; title('Template matching - unaligned shapes');
plot(log10(fp_count12), sensitivity12);
plot(log10(fp_count1), sensitivity1);
plot(log10(fp_count2), sensitivity2);
legend({'G'' * G"', 'G''', 'G"'}, 'location', 'southeast');
xlabel('Total FPs - log_{10} scale');
ylabel('Sensitivity');
set(gca, 'ylim', [0 1]);

%%
%d) Allow scale variation as well?
load([nailfoldroot 'data\apex_detection\apex_gd_templates.mat'], 'apex_template*');
nf_files = dir('C:\isbe\nailfold\images\anonymous_oct\annotations_qt\*.txt');

nn = 8;
image_path = ['C:\isbe\nailfold\images\anonymous_oct\bmp\' nf_files(nn).name(1:end-11) '.bmp'];
nailfold = imread(image_path);   
nailfold = nailfold(:,:,1);
[rows cols] = size(nailfold);
[mag_1d] = abs(gaussian_1st_derivative_gradient(nailfold, 2));
clear nailfold;
    
C1 = -inf(rows, cols);
for jj = 1:11
    for ii = -5:5
        scale = 1 + ii/10;
        mag_1ds = imresize(mag_1d, scale, 'bilinear');
        c1 = imresize(mb_normxcorr2(apex_template1_a(:,:,jj), mag_1ds), [rows cols], 'bilinear');
        C1 = max(C1, c1);
        %figure; imgray(c1); colormap jet; caxis([-1 1]);
    end
end
figure; imgray(C1); colormap jet; caxis([-1 1]);
%%
%--------------------------------------------------------------------------
%% 3) Create regions to test AAM fitting after initial candidate point selection
%--------------------------------------------------------------------------
%Reload threshold for aligned method
load([nailfoldroot 'data\apex_detection\apex_templates_aligned.mat'], 'apex_template*');
load([nailfoldroot 'data\apex_detection\aligned\thresh.mat'], 'thresh2');
load([nailfoldroot 'data\apex_detection\mean_shape.mat'], 'mean_shape');
nf_files = dir('C:\isbe\nailfold\images\anonymous_oct\annotations_qt\*.txt');

x = repmat(-24:24, 49, 1);
y = repmat(-24:24, 49, 1)';
xy = [x(:) y(:)];

circle_mask = ~(x.^2 + y.^2 > 24^2);
N = sum(circle_mask(:));
template = apex_template2_a(circle_mask);

T = sum(template) / N;
T2 = sum(template.^2) / N;
denom_T = sqrt( max(T2 - T^2,0) );

apex_num = 1;
%             
for nn = 7:12
    %load image
    image_path = ['C:\isbe\nailfold\images\anonymous_oct\bmp\' nf_files(nn).name(1:end-11) '.bmp'];
    nailfold = imread(image_path);   
    nailfold = nailfold(:,:,1);
    nailfold_mask = imerode(nailfold > 10 & nailfold < 250, strel('disk', 20));
    
    %load correlation map for G"
    load([nailfoldroot 'data\apex_detection\aligned\' nf_files(nn).name(1:end-4) '_corr.mat'], 'C2');
    
    %Compute the gaussian derivatives
    [mag_2d] = gaussian_2nd_derivative_line(nailfold, 4);
    
    %Get local maxima using precomputed threshold
    [maxima_pos, maxima_vals] = local_image_maxima(C2, 20, nailfold_mask, thresh2);
    
    %Load in ground truth apex locations
    apexes = u_load(['C:\isbe\nailfold\data\apexes\' nf_files(nn).name(1:end-4) '_apexes.mat']);
    apex_gt = zeros(length(apexes),2);
    for ii = 1:length(apexes)
        apex_gt(ii,:) = apexes{ii}(16,:);
    end
    
    %Workout whether candidates are hits or not
    [candidate_label gt_hits candidate_match] = evaluate_apex_candidates(apex_gt, maxima_pos, 15);
    
    num_maxima = size(maxima_pos,1);
    maxima_scales = zeros(num_maxima,1);
    maxima_thetas = zeros(num_maxima,1);
    
    for ii = 1:num_maxima
        c_max = 0;
        for theta = -15:3:15
            rot = [cosd(theta) sind(theta); -sind(theta) cosd(theta)];

            for scale = 0.8:0.1:1.5

                xya = xy * rot * scale;
                xa = reshape(xya(:,1) + maxima_pos(ii,1), 49, 49);
                ya = reshape(xya(:,2) + maxima_pos(ii,2), 49, 49);
                
                gd2_patch = interp2(mag_2d, xa, ya, 'bilinear');
                gd2_vec = gd2_patch(circle_mask);
                
                A = sum(gd2_vec) / N;
                A2 = sum(gd2_vec.^2) / N;
                
                denom = denom_T * sqrt( max(A2 - A^2,0) );
                numerator = (template' * gd2_vec) / N - A*T;
                
                c = numerator / denom;
                
                if c > c_max
                    maxima_scales(ii) = scale;
                    maxima_thetas(ii) = theta;
                    vessel_xy = mean_shape * rot * scale;
                    c_max = c;
                end
                
                %figure; imgray(gd2_patch);
                %title(['\theta = ' num2str(theta) ', scale = ' num2str(scale) ', C = ' num2str(c)]);
                
            end
        end
        
        sr = max(1, floor(min(vessel_xy(:,2)) - 50 + maxima_pos(ii,2)));
        er = min(size(nailfold,1), ceil(max(vessel_xy(:,2)) + 50 + maxima_pos(ii,2)));
        sc = max(1, floor(min(vessel_xy(:,1)) - 50 + maxima_pos(ii,1)));
        ec = min(size(nailfold,2), floor(max(vessel_xy(:,1)) + 50 + maxima_pos(ii,1)));
        
        vessel_xy(:,1) = vessel_xy(:,1) - sc + maxima_pos(ii,1);
        vessel_xy(:,2) = vessel_xy(:,2) - sr + maxima_pos(ii,2);
        
        %Sample patch from image
        image_patch = nailfold(sr:er, sc:ec);
        
        true_vessel_xy = [];
        if candidate_match(ii)
            true_vessel_xy = apexes{candidate_match(ii)};
            true_vessel_xy(:,1) = true_vessel_xy(:,1) - sc;
            true_vessel_xy(:,2) = true_vessel_xy(:,2) - sr;
            
            if maxima_scales(ii) > 1.1
                figure; imgray(image_patch);
                plot(vessel_xy(:,1), vessel_xy(:,2), 'g');
                plot(true_vessel_xy(:,1), true_vessel_xy(:,2), 'r');
                title(['\theta = ' num2str(maxima_thetas(ii)) ', scale = ' num2str(maxima_scales(ii)) ', C = ' num2str(c_max)]);
            end            
        end
        
        apex_candidate.image_patch = image_patch;
        apex_candidate.vessel_xy = vessel_xy;
        apex_candidate.true_vessel_xy = true_vessel_xy;
        apex_candidate.scale = maxima_scales(ii);
        apex_candidate.theta = maxima_thetas(ii);
        apex_candidate.sr = sr;
        apex_candidate.sc = sc;
        apex_candidate.method = 'template_matching_aligned_g2d';
        apex_candidate.image_path = image_path;
        
        save(['C:\isbe\nailfold\data\aam\candidates\template_matching\aligned\apex' zerostr(apex_num, 4) '.mat'], 'apex_candidate');
        apex_num = apex_num + 1;
        
    end
end

%%
candidate_labels = false(338,1);
num_pts = 31;

fid(1) = fopen('C:\isbe\nailfold\data\aam\candidates\template_matching\aligned\all_candidates_test.smd', 'wt');
fid(2) = fopen('C:\isbe\nailfold\data\aam\candidates\template_matching\aligned\true_candidates_test.smd', 'wt');
fid(3) = fopen('C:\isbe\nailfold\data\aam\candidates\template_matching\aligned\false_candidates_test.smd', 'wt');

for ff = 1:3
    fprintf(fid(ff), '%s \n', '// Text file describing test data for vessel apex models');
    fprintf(fid(ff), '\n');
    fprintf(fid(ff), '%s \n', '// Directory containing images');
    fprintf(fid(ff), '%s \n', 'image_dir: C:/isbe/nailfold/data/aam/candidates/template_matching/aligned/images/orig/');
    fprintf(fid(ff), '%s \n', '// Directory containing points');
    fprintf(fid(ff), '%s \n', 'points_dir: C:/isbe/nailfold/data/aam/candidates/template_matching/aligned/points/');
    fprintf(fid(ff), '\n');
    fprintf(fid(ff), '%s \n', '// Details of points : images');
    fprintf(fid(ff), '\n');
    fprintf(fid(ff), '%s \n', 'training_set:');
    fprintf(fid(ff), '%s \n', '{');
end

for ii = 1:338
    
    %Load in saved candidate details
    load(['C:\isbe\nailfold\data\aam\candidates\template_matching\aligned\apex' zerostr(ii, 4) '.mat'], 'apex_candidate');
    vessel_xy = apex_candidate.vessel_xy;
    image_patch = apex_candidate.image_patch;
    candidate_labels(ii) = ~isempty(apex_candidate.true_vessel_xy);
    clear apex_candidate;
        
    %Write out image patch
    imwrite(image_patch, ['C:\isbe\nailfold\data\aam\candidates\template_matching\aligned\images\orig\candidate_apex' zerostr(ii,4) '.png']);
    
    %Write out a pts file we can read in to VXL
    fid1 = fopen(['C:\isbe\nailfold\data\aam\candidates\template_matching\aligned\points\candidate_apex' zerostr(ii,4) '.pts'], 'wt');
    fprintf(fid1, '%s \n', 'version: 1');
    fprintf(fid1, '%s %d \n', 'n_points:', num_pts);
    fprintf(fid1, '%s \n', '{'); 
    for jj = 1:num_pts
        fprintf(fid1,'%.2f %.2f \n', vessel_xy(jj,1), vessel_xy(jj,2));
    end
    fprintf(fid1, '%s \n', '}');
    fprintf(fid1, 'nailfold: %s \n', image_path);
    fprintf(fid1, '%s %d \n', 'start_row:', sr);
    fprintf(fid1, '%s %d \n', 'start_col: ', sc);
    fclose(fid1);
    
    %Write entry for this image/pts pair in model .smd files
    str = ['candidate_apex' zerostr(ii,4) '.pts : candidate_apex' zerostr(ii,4) '.png'];
    fprintf(fid(1), '%s \n', str);
    if candidate_labels(ii)
        fprintf(fid(2), '%s \n', str);
    else
        fprintf(fid(3), '%s \n', str);
    end
end
for ff = 1:3
    fprintf(fid(ff), '%s \n', '}');
    fclose(fid(ff));
end
save('C:\isbe\nailfold\data\aam\candidates\template_matching\aligned\candidate_labels.mat', 'candidate_labels');

%%
load('C:\isbe\nailfold\data\aam\candidates\template_matching\aligned\candidate_labels.mat', 'candidate_labels');
fid = fopen('C:/isbe/nailfold/data/aam/candidates/template_matching/aligned/model_qualities.txt');
q_txt = textscan(fid, '%s %f', 'delimiter', ':');
fclose(fid);
model_q = q_txt{2}; clear q_txt;
[sorted_model_qualities qidx] = sort(model_q, 'descend');

feature = 'orig';

num_rows = 3;
num_cols = 4;
ii = 1;

while ii <= 338
    figure;
    for row = 1:num_rows, 
        for col = 1:num_cols
            if ii > 338
                break;
            end
            jj = qidx(ii);
            
            load(['C:\isbe\nailfold\data\aam\candidates\template_matching\aligned\apex' zerostr(jj, 4) '.mat'], 'apex_candidate');
            f1 = fopen(['C:\isbe\nailfold\data\aam\candidates\template_matching\aligned\out_points\orig\candidate_apex' zerostr(jj,4) '.pts']);
            textscan(f1, '%[^{]');
            fgetl(f1);
            vessel_str = textscan(f1,'%[^}]');
            fclose(f1);
            test_pts = str2num(vessel_str{1}{1});

            im = imread(['C:\isbe\nailfold\data\aam\candidates\template_matching\aligned\images\orig\candidate_apex' zerostr(jj,4) '.png']);

            axes('units', 'normalized', 'position', [(col-1)/num_cols (num_rows-row)/num_rows 1/num_cols 1/num_rows]);
            
            [r c] = size(im);
            roi = im(round(r/4):round(3*r/4), round(c/4):round(3*c/4));
            
            imagesc(im); axis image off; colormap(gray(256)); caxis([min(roi(:)) max(roi(:))]);
            hold on;
            if ~candidate_labels(jj)
                plot(test_pts(:,1), test_pts(:,2), 'bx');
            else
                plot(apex_candidate.true_vessel_xy(:,1),apex_candidate.true_vessel_xy(:,2), 'rx'); 
                plot(test_pts(:,1), test_pts(:,2), 'gx');
            end
            text(10, 10, sprintf('QoF: %5.2f', sorted_model_qualities(ii)), 'color', 'r');
            ii = ii + 1;
            
        end
    end
    
end
%%

%%
if ~isdir([nailfoldroot 'data\apex_detection\']); mkdir([nailfoldroot 'data\apex_detection\']); end
if ~isdir([nailfoldroot 'images\anonymous_oct\apexes\']); mkdir([nailfoldroot 'images\anonymous_oct\apexes\']); end
%%-------------------------------------------------------------------------
%%Old Stuff....
% %% 1) Locate a set of apexes in 1 image given manual annotations
% % We'll use MB's markup of the KK demo image (because it just happens to be
% % an image I've marked up lots of vessels on)
% 
% %load in the nailfold - note we need to resize this by a factor 3/4 to
% %match the other nailfolds in our set
% nailfold = imresize(imread([nailfoldroot 'images\kk_demo\nailfolds\n3_mb.bmp']), 0.75, 'bilinear');
% 
% %Load in the annotated vessels
% v = u_load([nailfoldroot 'images\kK_demo\nailfolds\n3_mb_vessels_all.mat']);
% vessels = [];
% for ii = 1:length(v)
%     if size(v{ii},1) > 1;       
%         %Discard duplicate points
%         keep = [true; any(diff(v{ii}),2)];
%         vessel = 0.75*[v{ii}(keep,1) v{ii}(keep,2)];
% 
%         %Sample points evenly along the vessel
%         dists = cumsum([0; sqrt(sum(diff(vessel).^2, 2))]);
%         vessels{end+1,1} = interp1(dists, vessel, linspace(0, dists(end), floor(dists(end))), 'linear'); %#ok
%     end
% end
% num_vessels = length(vessels);
% 
% %Now choose each vessel apex as the point with highest y-val (I know this
% %may not strictly be accurate given local deviations, rotations etc - but
% %it appears to work more robustly than more complex methods such as
% %searching for points of maximum curvature etc)
% figure; hold on; axis equal ij;
% vessel_apex = zeros(num_vessels, 2);
% for ii = 1:num_vessels
%     [y_min min_idx] = min(vessels{ii}(:,2));
%     vessel_apex(ii,:) = [vessels{ii}(min_idx,1) y_min];
%     plot(vessels{ii}(:,1), vessels{ii}(:,2)); 
%     plot(vessel_apex(ii,1), vessel_apex(ii,2), 'rx');          
% end
% 
% %% 2) Now compute Gaussian 1st + 2nd derivatives of the nailfold and sample
% % a patch of each about each apex
% 
% %Compute the Gaussian derivatives
% [mag_2d, ori_2d] = gaussian_clover_line(nailfold, 4);
% [mag_1d] = gaussian_1st_derivative_gradient2(nailfold, 2);
% 
% f1 = figure;
% f2 = figure;
% 
% gd1_patches = zeros(33,33,num_vessels);
% gd2_patches = zeros(33,33,num_vessels);
% for ii = 1:num_vessels
%     
%     %Sample patch from 1st deriv.
%     mag_patch = sample_window(mag_1d, 33, round(vessel_apex(ii,2)), round(vessel_apex(ii,1)));
%     mag_patch = (mag_patch - mean(mag_patch(:))) / std(mag_patch(:)); %Normalise patch
%     gd1_patches(:,:,ii) = mag_patch;
%     
%     %Sample patch from 1st deriv.
%     mag_patch = sample_window(mag_2d, 33, round(vessel_apex(ii,2)), round(vessel_apex(ii,1)));
%     mag_patch = (mag_patch - mean(mag_patch(:))) / std(mag_patch(:));  %Normalise patch
%     gd2_patches(:,:,ii) = mag_patch;
% 
%     if ii <= 12
%         figure(f1);
%         subplot(3,4,ii); imagesc(gd1_patches(:,:,ii)); axis image; colormap(gray(256));
%         figure(f2);
%         subplot(3,4,ii); imagesc(gd2_patches(:,:,ii)); axis image; colormap(gray(256));
%     end
% end
% 
% %Compute the mean of the 2 sets of patches to use as templates
% apex_template1 = mean(gd1_patches,3);
% apex_template2 = mean(gd2_patches,3);
% 
% %save the templates for later use
% save([nailfoldroot 'data\apex_detection\apex_gd_templates.mat'], 'apex_template*');
% 
% figure; imagesc(apex_template1); axis image;
% title('Template for 1st derivative patches');
% figure; imagesc(apex_template2); axis image;
% title('Template for 2nd derivative patches');
% 
% %% 3) Use the template to try and find vessel apexes in a new nailfold
% 
% %load in nailfold
% nailfold = imread([nailfoldroot 'images\ncm\andrea Murray.bmp']);
% nailfold = nailfold(:,:,1);
% nailfold_mask = imerode(nailfold > 10 & nailfold < 250, strel('disk', 20));
% 
% %Compute the gaussian derivatives
% [mag_2d] = gaussian_clover_line(nailfold, 4);
% [mag_1d] = gaussian_1st_derivative_gradient2(nailfold, 2);
% 
% %Apply normalised cross correlations using the 2 templates
% C1 = normxcorr2(apex_template1, mag_1d);
% C1 = C1(17:end-16, 17:end-16); %Need to discard floor(patch_size)/2 from start and end
% 
% C2 = normxcorr2(apex_template2, mag_2d);
% C2 = C2(17:end-16, 17:end-16); %Need to discard floor(patch_size)/2 from start and end
% 
% %Look for local maxima in the 2 independent NCC maps
% [maxima_pos1, maxima_vals1] = local_image_maxima(C1, 33, nailfold_mask, 0.5);
% [maxima_pos2, maxima_vals2] = local_image_maxima(C2, 33, nailfold_mask, 0.5);
% 
% %Display the maps and plot the local maxima
% figure; 
% subplot(3,1,1); imagesc(nailfold); axis image; hold on; colormap(gray(256));
% plot(maxima_pos1(:,1), maxima_pos1(:,2), 'r+');
% plot(maxima_pos2(:,1), maxima_pos2(:,2), 'gx');
% subplot(3,1,2); imagesc(C1); axis image; hold on;
% subplot(3,1,3); imagesc(C2); axis image; hold on;
% 
% %Now multiply the 2 NCC maps to get a single map and again look for local
% %maxima
% C12 = C1 .* C2;
% C12(~nailfold_mask) = 0;
% [maxima_pos12, maxima_vals12] = local_image_maxima(C12, 33, [], 0.2);
% 
% %Display a heat map of the combined tempate matching
% figure; imagesc(C12); axis image off; colormap(hot(256)); caxis([0 max(maxima_vals12)]);
% 
% %Plot the maxima in the nailfold image, with marker size scaled according
% %to NCC value
% n_min = min(nailfold(nailfold_mask));
% n_max = max(nailfold(nailfold_mask)); 
% figure; imagesc(nailfold); axis image off; hold on; colormap(gray(256)); caxis([n_min n_max]);
% for ii = 1:length(maxima_vals12)
%     ms = 1 + 20*(maxima_vals12(ii)-0.2);
%     plot(maxima_pos12(ii,1), maxima_pos12(ii,2), 'r+', 'MarkerSize', ms);
% end
% %%
% %Repeat the above test for each nailfold in the anonymous OCT set
% 
% load([nailfoldroot 'data\apex_detection\apex_gd_templates.mat'], 'apex_template*');
% 
% n_files = dir([nailfoldroot 'images\anonymous_oct\bmp\*.bmp']);
% do_plot = 0;
% for nn = 15
%     %load in nailfold
%     nailfold = imread([nailfoldroot 'images\anonymous_oct\bmp\' n_files(nn).name]);
%     nailfold = nailfold(:,:,1);
%     nailfold_mask = imerode(nailfold > 10 & nailfold < 250, strel('disk', 20));
% 
%     %Compute the gaussian derivatives
%     [mag_2d] = gaussian_clover_line(nailfold, 4);
%     [mag_1d] = gaussian_1st_derivative_gradient2(nailfold, 2);
% 
%     %Apply normalised cross correlations using the 2 templates
%     C1 = normxcorr2(apex_template1, mag_1d);
%     C1 = C1(17:end-16, 17:end-16); %Need to discard floor(patch_size)/2 from start and end
% 
%     C2 = normxcorr2(apex_template2, mag_2d);
%     C2 = C2(17:end-16, 17:end-16); %Need to discard floor(patch_size)/2 from start and end
% 
%     %Look for local maxima in the 2 independent NCC maps
%     [maxima_pos1, maxima_vals1] = local_image_maxima(C1, 33, nailfold_mask, 0.5);
%     [maxima_pos2, maxima_vals2] = local_image_maxima(C2, 33, nailfold_mask, 0.5);
% 
%     %Now multiply the 2 NCC maps to get a single map and again look for local
%     %maxima
%     C12 = C1 .* C2;
%     C12(~nailfold_mask) = 0;
%     [maxima_pos12, maxima_vals12] = local_image_maxima(C12, 33, [], 0.2);
%     
%     if do_plot
%         figure; 
%         %Display a heat map of the combined tempate matching
%         subplot(2,1,1); imagesc(C12); axis image off; colormap(hot(256)); caxis([0 max(maxima_vals12)]);
% 
%         %Plot the maxima in the nailfold image, with marker size scaled according
%         %to NCC value
%         n_min = min(nailfold(nailfold_mask));
%         n_max = max(nailfold(nailfold_mask)); 
%         subplot(2,1,2); imagesc(nailfold); axis image off; hold on; colormap(gray(256)); caxis([n_min n_max]);
%         for ii = 1:length(maxima_vals12)
%             ms = 1 + 20*(maxima_vals12(ii)-0.2);
%             plot(maxima_pos12(ii,1), maxima_pos12(ii,2), 'r+', 'MarkerSize', ms);
%         end
%     end
%     
%     %%
%     save([nailfoldroot 'images\anonymous_oct\apexes\' n_files(nn).name(1:end-4) '_apexes.mat'], 'maxima_*');
% end
% %%
% nf_files = dir('C:\isbe\nailfold\images\anonymous_oct\annotations_qt\*.txt');
% 
% for nn = 8:12
%     image_path = ['C:\isbe\nailfold\images\anonymous_oct\bmp\' nf_files(nn).name(1:end-11) '.bmp'];
%     nailfold = imread(image_path);   
%     nailfold = nailfold(:,:,1);
%     nailfold_mask = imerode(nailfold > 10 & nailfold < 250, strel('disk', 20));
% 
%     load(['C:\isbe\nailfold\data\apexes\' nf_files(nn).name(1:end-4) '_apexes.mat'], 'apexes');
%     
%     %Compute the gaussian derivatives
%     [mag_2d ori_2d] = gaussian_2nd_derivative_line(nailfold, 4);
%     vessel_nms = mb_non_maximal_supp(mag_2d, ori_2d);
%     vessel_nms(~nailfold_mask) = 0;
%     [yv xv] = find(vessel_nms);
%     
%     n_min = min(nailfold(nailfold_mask));
%     n_max = max(nailfold(nailfold_mask)); 
%     figure; imgray(nailfold); caxis([n_min n_max]);
%     plot(xv, yv, 'g.', 'markersize', 2);
%     
%     for ii = 1:length(apexes)
%         plot(apexes{ii}(:,1), apexes{ii}(:,2), 'r');
%     end
% end