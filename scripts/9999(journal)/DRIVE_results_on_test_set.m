retroot = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\test\';
fov_mask_dir = [retroot 'fov_masks\'];
fg_mask_dir = [retroot 'vessel_masks\'];
%%
pred_dir = [retroot 'predictions\detection\rf_classification\'];
label_dir = [retroot 'vessel_masks\'];

rf_codes = cell(20, 4);
% rf_codes( 1,:) = {'8892', 'g2', '1', 'orig'};
% rf_codes( 2,:) = {'8893', 'g2', '3', 'orig'};
% rf_codes( 3,:) = {'8894', 'dt', '1', 'orig'};
% rf_codes( 4,:) = {'8895', 'dt', '3', 'orig'};
% rf_codes( 5,:) = {'9151', 'g2', '1', '0.25'};
% rf_codes( 6,:) = {'9152', 'g2', '3', '0.25'};
% rf_codes( 7,:) = {'9153', 'dt', '1', '0.25'};
% rf_codes( 8,:) = {'9154', 'dt', '3', '0.25'};
% rf_codes( 9,:) = {'9155', 'g2', '1', '0.50'};
% rf_codes(10,:) = {'9156', 'g2', '3', '0.50'};
% rf_codes(11,:) = {'9157', 'dt', '1', '0.50'};
% rf_codes(12,:) = {'9158', 'dt', '3', '0.50'};

rf_codes( 1,:) = {'9162',  'dt', '3', 'orig'};
rf_codes( 2,:) = {'88122', 'gabor_ri', '3', 'orig'};
rf_codes( 3,:) = {'24303', 'gabor_a', '3', 'orig'}; %Old version of filters '11480'
rf_codes( 4,:) = {'88126', 'gh2da', '3', 'orig'};
rf_codes( 5,:) = {'12394', 'mono', '3', 'orig'};
rf_codes( 6,:) = {'42382', 'gabor_i', '3', 'orig'};
rf_codes( 7,:) = {'50981', 'gabor_big_scale', '3', 'orig'};
rf_codes( 8,:) = {'51417', 'gabor_big_angle', '3', 'orig'};

% rf_codes(11,:) = {'9159', 'dt', '1', 'orig'};
% rf_codes(12,:) = {'9161', 'g2', '1', 'orig'};
% rf_codes(13,:) = {'', 'g12d', '1', 'orig'};
% rf_codes(14,:) = {'', 'mono', '1', 'orig'};
% rf_codes(15,:) = {'', 'linop', '1', 'orig'};
% rf_codes(16,:) = {'', 'haar', '1', 'orig'};
% rf_codes(17,:) = {'', 'gabor_a', '1', 'orig'}; %Old version of filters '11480'
% rf_codes(18,:) = {'', 'g2di', '1', 'orig'};
% rf_codes(19,:) = {'', 'g1d', '1', 'orig'};
% rf_codes(20,:) = {'', 'gabor_i', '1', 'orig'};
%
auc_all = [];
roc_pts_all = [];
for ii = [1 5]
    if 0%exist([pred_dir rf_codes{ii,1} '\roc\roc_stats.mat'], 'file')
        load([pred_dir rf_codes{ii,1} '\roc\roc_stats.mat'])
    else
        [roc_pts, auc, auc_ind] = compute_roc_image_set([pred_dir rf_codes{ii,1} '\'], label_dir, 'fov_mask_dir', fov_mask_dir);
        create_folder([pred_dir rf_codes{ii,1} '\roc\']);
        save([pred_dir rf_codes{ii,1} '\roc\roc_stats.mat'], 'roc_pts', 'auc', 'auc_ind');
    end
    roc_pts_all(:,:,ii) = roc_pts; %#ok
    auc_all(ii,1) = auc; %#ok
    auc_all(ii,2:21) = auc_ind; %#ok
end

if exist([retroot 'predictions\detection\staal\roc\roc_stats.mat'], 'file')
    load([retroot 'predictions\detection\staal\roc\roc_stats.mat'])
else
    [roc_pts, auc, auc_ind] = compute_roc_image_set([retroot 'predictions\detection\staal\'], label_dir, 'fov_mask_dir', fov_mask_dir);
    create_folder([pred_dir rf_codes{ii,1} '\roc\']);
    save([pred_dir rf_codes{ii,1} '\roc\roc_stats.mat'], 'roc_pts', 'auc', 'auc_ind');
end
%
rfs = find(auc_all(:,1));
for i_rf = rfs'
    display(['Decomp: ' rf_codes{i_rf, 2} ', w = ' rf_codes{i_rf,3} ', sampling mode: ' rf_codes{i_rf, 4}...
        ', Az = ' num2str(auc_all(i_rf,1),3) ' +/- ' num2str(std(auc_all(i_rf,2:21)),3)]);
end
%%
%Compute Staal's ROC curve
pred_dir = [retroot 'predictions\detection\staal\'];
if exist([pred_dir 'roc\roc_stats.mat'], 'file')
    load([pred_dir '\roc\roc_stats.mat'])
else
    [roc_pts, auc, auc_ind] = compute_roc_image_set(pred_dir, label_dir, 'fov_mask_dir', fov_mask_dir);
    create_folder([pred_dir '\roc\']);
    save([pred_dir '\roc\roc_stats.mat'], 'roc_pts', 'auc', 'auc_ind');
end
%%
warning('off', 'ori_error:nans');
pred_dir = [retroot 'predictions\orientation\rf_regression\'];
label_dir = [retroot 'orientations\'];

load([label_dir 'gt\all_gt_orientations.mat'], 'gt_orientations', 'vessel_centres');
rf_codes = cell(4, 4);
rf_codes( 1,:) = {'42813', 'dt', '3', '0.75'};
rf_codes( 2,:) = {'13295', 'mono', '3', '0.75'};
rf_codes( 3,:) = {'24308', 'gabor_a', '3', '0.75'}; %Old version of filters '11480'
rf_codes( 4,:) = {'88171', 'gh2da', '3', '0.75'};

error_medians = zeros(4,2);
error_means = zeros(4,2);

for ii = 1:4

    if exist([pred_dir rf_codes{ii,1} '\errors\error_stats.mat'], 'file')
        load([pred_dir rf_codes{ii,1} '\errors\error_stats.mat']);
    else
        [orientation_errors, ~, error_stats] =...
            compute_image_orientation_errors([pred_dir rf_codes{ii,1} '\'], fg_mask_dir,...
            'gt_orientations', gt_orientations, 'fov_mask_dir', fov_mask_dir,...
            'save_dir', [pred_dir rf_codes{ii,1} '\errors\']);
    end
    centre_stats = ori_error_stats(orientation_errors(vessel_centres));
    
    error_medians(ii,1) = error_stats.abs_median;
    error_medians(ii,2) = centre_stats.abs_median;
    
    error_means(ii,1) = error_stats.abs_mean;
    error_means(ii,2) = centre_stats.abs_mean;
    
    display(['Errors for ' rf_codes{ii,2} ', w = ' rf_codes{ii,3} ', sampling = ' rf_codes{ii,4}...
        ', abs mean: ' num2str(error_stats.abs_mean,3) ', abs_median: ' num2str(error_stats.abs_median,3)]); 
       
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%*************************** STARE ****************************************
%%
%
retroot = 'C:\isbe\asymmetry_project\data\retinograms\STARE\training\';
fov_mask_dir = [retroot 'fov_masks_eroded\'];
fg_mask_dir = [retroot 'vessel_masks\'];
warning('off', 'ori_error:nans');
%%
pred_dir = [retroot 'predictions\orientation\rf_regression\'];
label_dir = [retroot 'orientations\'];

load([label_dir 'gt\all_gt_orientations.mat'], 'gt_orientations', 'vessel_centres');

rf_codes = cell(4,4);
rf_codes( 1,:) = {'42813', 'dt', '3', '0.75'};
rf_codes( 2,:) = {'13295', 'mono', '3', '0.75'};
rf_codes( 3,:) = {'24308', 'gabor_a', '3', '0.75'}; %Old version of filters '11480'
rf_codes( 4,:) = {'88171', 'gh2da', '3', '0.75'};

error_medians = zeros(1,2);
error_means = zeros(1,2);

for ii = 1:4

    if exist([pred_dir rf_codes{ii,1} '\errors\error_stats.mat'], 'file')
        load([pred_dir rf_codes{ii,1} '\errors\error_stats.mat']);
    else
        [orientation_errors, ~, error_stats] =...
            compute_image_orientation_errors([pred_dir rf_codes{ii,1} '\'], fg_mask_dir,...
            'label_dir', label_dir, 'fov_mask_dir', fov_mask_dir,...
            'save_dir', [pred_dir rf_codes{ii,1} '\errors\']);
    end
    
    error_medians(ii,1) = error_stats.abs_median;
    error_means(ii,1) = error_stats.abs_mean;
    
    display(['Errors for ' rf_codes{ii,2} ', w = ' rf_codes{ii,3} ', sampling = ' rf_codes{ii,4}...
        ', abs mean: ' num2str(error_stats.abs_mean,3) ', abs_median: ' num2str(error_stats.abs_median,3)]); 
       
end
%%
pred_dir = [retroot 'predictions\detection\rf_classification\old\'];
label_dir = [retroot 'vessel_masks\'];

rf_codes = cell( 4, 4);

rf_codes( 1,:) = {'gh2da', 'dt', '3', '0.75'};
rf_codes( 2,:) = {'gabor', 'mono', '3', '0.75'};
rf_codes( 3,:) = {'dt', 'gabor_a', '3', '0.75'}; %Old version of filters '11480'
rf_codes( 4,:) = {'mono', 'gh2da', '3', '0.75'};

roc_pts_all = [];
auc_all = [];
%
for ii = 1:4
    if 0%exist([pred_dir rf_codes{ii,1} '\roc\roc_stats.mat'], 'file')
        load([pred_dir rf_codes{ii,1} '\roc\roc_stats.mat'])
    else
        [roc_pts, auc, auc_ind] = ...
            compute_roc_image_set([pred_dir rf_codes{ii,1} '\'],...
            label_dir, 'fov_mask_dir', fov_mask_dir);
        create_folder([pred_dir rf_codes{ii,1} '\roc\']);
        save([pred_dir rf_codes{ii,1} '\roc\roc_stats.mat'], 'roc_pts', 'auc', 'auc_ind');
    end
    roc_pts_all(:,:,ii) = roc_pts; %#ok
    auc_all(ii,1) = auc; %#ok
    auc_all(ii,2:21) = auc_ind; %#ok
end
%
rfs = find(auc_all(:,1));
for i_rf = rfs'
    display(['Decomp: ' rf_codes{i_rf, 2} ', w = ' rf_codes{i_rf,3} ', sampling mode: ' rf_codes{i_rf, 4}...
        ', Az = ' num2str(auc_all(i_rf,1),3) ' +/- ' num2str(auc_all(i_rf,2),3)]);
end
%%
rf_codes = cell( 4, 4);

rf_codes( 1,:) = {'gh2da', 'dt', '3', '0.75'};
rf_codes( 2,:) = {'gabor', 'mono', '3', '0.75'};
rf_codes( 3,:) = {'dt', 'gabor_a', '3', '0.75'}; %Old version of filters '11480'
rf_codes( 4,:) = {'mono', 'gh2da', '3', '0.75'};

% rf_codes( 1,:) = {'88126', 'gaussian', '3', '0.75'};
% rf_codes( 2,:) = {'24303', 'gabor', '3', '0.75'};
% rf_codes( 3,:) = {'9162', 'dt', '3', '0.75'}; %Old version of filters '11480'
% rf_codes( 4,:) = {'12394', 'mono', '3', '0.75'};

p_place = [2 3 5 6];
for ii = 1:20
    
    r_im = u_load([retroot 'images\' zerostr(ii,2) '_training.mat']);
    m_im = u_load([retroot 'vessel_masks\' zerostr(ii,2) '_training_v_mask.mat']);
    
    figure;
    subplot(2,3,1); imgray(r_im);
    subplot(2,3,4); imgray(m_im);
    
    for i_p = 1:4
        p_im = u_load([pred_dir  'old\' rf_codes{i_p,1} '\' zerostr(ii,2) '_training_pred.mat']);
        %p_im2 = u_load([pred_dir 'old\' rf_codes{i_p,1} '\' zerostr(ii,2) '_training_pred.mat']);
        %subplot(2,3,p_place(i_p)); imgray(p_im - p_im2); colormap jet;
        subplot(2,3,p_place(i_p)); imgray(p_im); colormap jet;
    end
        
end