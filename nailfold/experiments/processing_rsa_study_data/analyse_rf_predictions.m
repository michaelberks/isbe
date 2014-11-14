warning('off', 'ori_error:nans');
pred_dir = 'C:\isbe\nailfold\data\rsa_study\set2\predictions\orientation\rf_regression\';
label_dir = 'C:\isbe\nailfold\data\rsa_study\set2\orientations\';
fov_mask_dir = 'C:\isbe\nailfold\data\rsa_study\set2\fov_masks\';
fg_mask_dir = 'C:\isbe\nailfold\data\rsa_study\set2\vessel_masks\';
fgc_mask_dir = 'C:\isbe\nailfold\data\rsa_study\set2\vessel_centre_masks\';
%%
rf_codes = cell(4, 2);
rf_codes( 1,:) = {'219121', 'OrigIms'};
rf_codes( 2,:) = {'219122', 'NormIms'};
rf_codes( 3,:) = {'219133', 'OrigImsAllVess'};
rf_codes( 4,:) = {'222210', 'RandBG'};
rf_codes( 5,:) = {'222835', 'RandBG?'};

error_medians = zeros(4,2);
error_means = zeros(4,2);
%%
for ii = 5%1:5

    if 0%exist([pred_dir rf_codes{ii,1} '\errors\error_stats.mat'], 'file')
        load([pred_dir rf_codes{ii,1} '\errors\error_stats.mat']);
    else
        [orientation_errors, ~, error_stats] =...
            compute_image_orientation_errors([pred_dir rf_codes{ii,1} '\'], fgc_mask_dir,...
            'label_dir', label_dir, 'fov_mask_dir', fov_mask_dir,...
            'save_dir', [pred_dir rf_codes{ii,1} '\errors\']);
    end
    %centre_stats = ori_error_stats(orientation_errors(vessel_centres));
    
    error_medians(ii,1) = error_stats.abs_median;
    %error_medians(ii,2) = centre_stats.abs_median;
    
    error_means(ii,1) = error_stats.abs_mean;
    %error_means(ii,2) = centre_stats.abs_mean;
    
    display(['Errors for ' rf_codes{ii,2} ...
        ', abs mean: ' num2str(error_stats.abs_mean,3) ', abs_median: ' num2str(error_stats.abs_median,3)]); 
       
end
%%
pred_dir = 'C:\isbe\nailfold\data\rsa_study\set2\predictions\detection\rf_classification\';
label_dir = 'C:\isbe\nailfold\data\rsa_study\set2\vessel_centre_masks\';
rf_codes = cell(2, 2);
rf_codes( 1,:) = {'219128', 'NormIms'};
rf_codes( 2,:) = {'219131', 'OrigIms'};
rf_codes( 3,:) = {'182321', 'OrigIms'};
rf_codes( 4,:) = {'182321', 'OrigIms'};

auc_allc = [];
roc_pts_allc = [];
%
for ii = 1:3
    if 0%exist([pred_dir rf_codes{ii,1} '\roc\roc_stats.mat'], 'file')
        load([pred_dir rf_codes{ii,1} '\roc\roc_stats.mat'])
    else
        [roc_pts, auc, auc_ind] = compute_roc_image_set([pred_dir rf_codes{ii,1} '\'], label_dir, 'fov_mask_dir', fov_mask_dir);
        %create_folder([pred_dir rf_codes{ii,1} '\roc\']);
        %save([pred_dir rf_codes{ii,1} '\roc\roc_stats.mat'], 'roc_pts', 'auc', 'auc_ind');
    end
    roc_pts_allc(:,:,ii) = roc_pts; %#ok
    auc_allc(ii,:) = [auc auc_ind']; %#ok
end
%%
warning('off', 'ori_error:nans');
fov_mask_dir = 'C:\isbe\nailfold\data\rsa_study\set12\fov_masks\';
fg_mask_dir = 'C:\isbe\nailfold\data\rsa_study\set12\vessel_masks\';
fgc_mask_dir = 'C:\isbe\nailfold\data\rsa_study\set12\vessel_centre_masks\';

pred_dir = 'C:\isbe\nailfold\data\rsa_study\set12\predictions\detection\rf_classification\';
label_dir = 'C:\isbe\nailfold\data\rsa_study\set12\vessel_centre_masks\';
rf_codes = cell(2, 2);
rf_codes( 1,:) = {'222836', 'NormIms'};
rf_codes( 2,:) = {'247131', 'OrigIms'};
rf_codes( 3,:) = {'256372', 'OrigIms'};
rf_codes( 4,:) = {'257273', 'OrigIms'};

auc_allc = [];
roc_pts_allc = [];
%
for ii = 1:4
    if exist([pred_dir rf_codes{ii,1} '\roc\roc_stats.mat'], 'file')
        load([pred_dir rf_codes{ii,1} '\roc\roc_stats.mat'])
    else
        [roc_pts, auc, auc_ind] = compute_roc_image_set([pred_dir rf_codes{ii,1} '\'], label_dir, 'fov_mask_dir', fov_mask_dir);
        %create_folder([pred_dir rf_codes{ii,1} '\roc\']);
        %save([pred_dir rf_codes{ii,1} '\roc\roc_stats.mat'], 'roc_pts', 'auc', 'auc_ind');
    end
    roc_pts_allc(:,:,ii) = roc_pts; %#ok
    auc_allc(ii,:) = [auc auc_ind']; %#ok
end
%%
fov_mask_dir = 'C:\isbe\nailfold\data\rsa_study\set1g\fov_masks\';
fg_mask_dir = 'C:\isbe\nailfold\data\rsa_study\set1g\vessel_masks\';
fgc_mask_dir = 'C:\isbe\nailfold\data\rsa_study\set1g\vessel_centre_masks\';

pred_dir = 'C:\isbe\nailfold\data\rsa_study\set1g\predictions\detection\rf_classification\';
label_dir = 'C:\isbe\nailfold\data\rsa_study\set1g\vessel_centre_masks\';
rf_codes = cell(2, 2);
rf_codes( 1,:) = {'222836', 'NormIms'};
rf_codes( 2,:) = {'256372', 'OrigIms'};
rf_codes( 3,:) = {'257273', 'OrigIms'};

auc_allc = [];
roc_pts_allc = [];
%
for ii = 1:3
    if exist([pred_dir rf_codes{ii,1} '\roc\roc_stats.mat'], 'file')
        load([pred_dir rf_codes{ii,1} '\roc\roc_stats.mat'])
    else
        [roc_pts, auc, auc_ind] = compute_roc_image_set([pred_dir rf_codes{ii,1} '\'], label_dir, 'fov_mask_dir', fov_mask_dir);
        %create_folder([pred_dir rf_codes{ii,1} '\roc\']);
        %save([pred_dir rf_codes{ii,1} '\roc\roc_stats.mat'], 'roc_pts', 'auc', 'auc_ind');
    end
    roc_pts_allc(:,:,ii) = roc_pts; %#ok
    auc_allc(ii,:) = [auc auc_ind']; %#ok
end
%%


