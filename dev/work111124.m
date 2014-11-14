load Z:\asym\data\retinograms\dRIVE\training\vessel_masks\21_training_v_mask.mat
gts = bwmorph(vessel_mask, 'skel', 'inf');
[yv xv] = find(gts);
num_pts = length(xv);
%%
% offs = -2:2;
% count_mask = ones(5); count_mask(2:4,2:4) = 0;

offs = -1:1;
count_mask = ones(3); count_mask(2,2) = 0;
is_branch = false(num_pts,1);
%
for ii = 1:num_pts
    
    patch = gts(yv(ii) + offs, xv(ii) + offs) & count_mask;
    label = bwlabel(patch, 4);
    is_branch(ii) = max(label(:)) >= 3;
end

figure; imagesc(gts); axis image; colormap gray; hold on;
plot(xv(is_branch), yv(is_branch), 'rx');
%%
d_root = 'Z:\asym\data\retinograms\dRIVE\test\';
%d_root = 'C:\isbe\asymmetry_project\data\retinograms\dRIVE\test\';

tp_counts = zeros(20, 102);
fp_counts = zeros(20, 102);
t_counts = zeros(20, 1);
f_counts = zeros(20, 1);

forest = '51101';

for ii = 1:20
    vessel_prob = load_uint8([d_root 'images_extended\results\' forest '\' zerostr(ii,2) '_test_ext_class.mat']);
    %vessel_prob = abs(load_uint8(['Z:\asym\data\retinograms\DRIVE\test\predictions2\dt\rf_3\' zerostr(ii,2) '_ori.mat']));
    %figure; imagesc(vessel_prob); axis image; colormap(gray(256));
    
    v_mask = u_load([d_root 'vessel_masks\' zerostr(ii,2) '_test_v_mask.mat']);
    f_mask = u_load([d_root 'foveal_masks\' zerostr(ii,2) '_test_f_mask.mat']);

    %v_mask = bwmorph(v_mask, 'skel', 'inf');
    
    %Compute ROC counts for image
    [roc_pts auc tp_count fp_count] = calculate_roc_curve(vessel_prob(f_mask), v_mask(f_mask),(-1:100)/100);
    %figure; plot(roc_pts(:,1), roc_pts(:,2)); axis([0 1 0 1]); title(['AUC: ' num2str(auc)]);
    %Increment total counts
    tp_counts(ii,:) = tp_count;
    fp_counts(ii,:) = fp_count;
    t_counts(ii) = sum(v_mask(f_mask));
    f_counts(ii) = sum(~v_mask(f_mask));

end

%Compute ROC points for complete set of data
roc_pts = [sum(fp_counts)' / sum(f_counts) sum(tp_counts)' / sum(t_counts)];

%Compute AUC for ROC curve
auc = sum( (roc_pts(2:end,1)-roc_pts(1:end-1,1)) .* roc_pts(1:end-1,2)) + ...
        0.5*sum( (roc_pts(2:end,1)-roc_pts(1:end-1,1)) .* (roc_pts(2:end,2)-roc_pts(1:end-1,2)) );
    
figure; plot(roc_pts(:,1), roc_pts(:,2)); axis([0 1 0 1]);
title([ forest ' AUC: ' num2str(auc)]);
%%
emap = [0 0 0; jet(255)];
forest1 = '49756';
forest2 = '51101';
for ii = 1:20
    vessel_prob1 = load_uint8(['Z:\asym\data\retinograms\dRIVE\test\images_extended\results\' forest1 '\' zerostr(ii,2) '_test_ext_class.mat']);
    vessel_prob2 = load_uint8(['Z:\asym\data\retinograms\dRIVE\test\images_extended\results\' forest2 '\' zerostr(ii,2) '_test_ext_class.mat']);
    
    %figure; imagesc(vessel_prob); axis image; colormap(gray(256));
    
    v_mask = u_load(['Z:\asym\data\retinograms\drIVE\test\vessel_masks\' zerostr(ii,2) '_test_v_mask.mat']);
    f_mask = u_load(['Z:\asym\data\retinograms\drIVE\test\foveal_masks\' zerostr(ii,2) '_test_f_mask.mat']);

    v_mask = bwmorph(v_mask, 'skel', 'inf');
    
    error_im1 = -vessel_prob1;
    error_im1(~v_mask) = nan;
    
    error_im2 = -vessel_prob2;
    error_im2(~v_mask) = nan;
    figure; 
    subplot(1,2,1); imagesc(error_im1); axis image; colormap(emap);
    subplot(1,2,2); imagesc(error_im2); axis image; colormap(emap);

end