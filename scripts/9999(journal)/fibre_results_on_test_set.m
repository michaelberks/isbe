fibre_root = 'C:\isbe\asymmetry_project\data\fibre\test\';
fov_mask_dir = [fibre_root 'fov_masks\'];
fg_mask_dir = [fibre_root 'vessel_masks\'];
%%
pred_dir = [fibre_root 'predictions\detection\rf_classification\'];
label_dir = [fibre_root 'fibre_masks\'];
tol_mask_dir = [fibre_root 'fibre_masks_dilated\'];

rf_codes = cell(4, 4);

rf_codes( 1,:) = {'88311', 'gabor', '3', '0.75'};
rf_codes( 2,:) = {'88443', 'dt', '3', '0.75'};
rf_codes( 3,:) = {'88444', 'gh2d', '3', '0.75'};
rf_codes( 4,:) = {'88445', 'mono', '3', '0.75'};

roc_pts = [];
auc = zeros(0,1);
auc_all = zeros(0,101);
%
c_args.fov_mask_dir = fov_mask_dir;
%c_args.tol_mask_dir = tol_mask_dir;
c_args.centre_only = 0;
c_args.thin = 0;
%%
for ii = 1:4
    if 0%exist([pred_dir rf_codes{ii,1} '\roc\roc_stats.mat'], 'file')
        load([pred_dir rf_codes{ii,1} '\roc\roc_stats.mat'])
    else
        [roc_pts, auc, auc_ind] = compute_roc_image_set(...
            [pred_dir rf_codes{ii,1} '\'], label_dir, c_args);
        
        create_folder([pred_dir rf_codes{ii,1} '\roc\']);
        save([pred_dir rf_codes{ii,1} '\roc\roc_stats.mat'], 'roc_pts', 'auc', 'auc_ind');
    end
    roc_pts_all(:,:,ii) = roc_pts; %#ok
    auc_all(ii,1) = auc; %#ok
    auc_all(ii,2:end) = auc_ind; %#ok
end
%%
figure; hold all;
rfs = find(auc_all(:,1));
for i_rf = rfs'
    display(['Decomp: ' rf_codes{i_rf, 2} ', w = ' rf_codes{i_rf,3} ', sampling mode: ' rf_codes{i_rf, 4}...
        ', Az = ' num2str(auc_all(i_rf,1),3) ' +/- ' num2str(std(auc_all(i_rf,2:21)),3)]);
    plot(roc_pts_all(:,1,i_rf), roc_pts_all(:,2,i_rf), '-x');
    
    [~, idx] = min(abs(sum(roc_pts_all(:,:,i_rf),2)-1));
    display([roc_pts_all(idx,:,i_rf) mean([roc_pts_all(idx,1,i_rf) 1-roc_pts_all(idx,2,i_rf)])]);
end
%plot([0 1], [1 0], 'k--');
%%
%--------------------------------------------------------------------------
% Orientation
%--------------------------------------------------------------------------
fibre_root = 'C:\isbe\asymmetry_project\data\fibre\test\';
fov_mask_dir = [fibre_root 'fov_masks\'];
fg_mask_dir = [fibre_root 'fibre_masks\'];
pred_dir = [fibre_root 'predictions\orientation\rf_regression\'];
label_dir = [fibre_root 'orientation_maps\'];

warning('off', 'ori_error:nans');
%

rf_codes = cell(4,4);
rf_codes( 1,:) = {'99189', 'gabor', '3', '0.75'};
rf_codes( 2,:) = {'99190', 'dt', '3', '0.75'};
rf_codes( 3,:) = {'99191', 'gh2da', '3', '0.75'}; %Old version of filters '11480'
rf_codes( 4,:) = {'99192', 'mono', '3', '0.75'};

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
f_list = dir([fibre_root 'images\*.mat']);
m_list = dir([fibre_root 'fibre_masks\*.mat']);
p_list = dir([pred_dir rf_codes{1,1} '\*.mat']);
for ii = 1:10
    
    f_im = u_load([fibre_root 'images\' f_list(ii).name]);
    m_im = u_load([fibre_root 'fibre_masks\' m_list(ii).name]);
    p_im1 = u_load([pred_dir rf_codes{1,1} '\' p_list(ii).name]);
    p_im2 = u_load([pred_dir rf_codes{2,1} '\' p_list(ii).name]);
    p_im3 = u_load([pred_dir rf_codes{3,1} '\' p_list(ii).name]);
    p_im4 = u_load([pred_dir rf_codes{4,1} '\' p_list(ii).name]);

    figure;
    subplot(2,3,1); imgray(f_im);
    subplot(2,3,4); imgray(m_im);
    subplot(2,3,2); imgray(p_im1);
    subplot(2,3,3); imgray(p_im2);
    subplot(2,3,5); imgray(p_im3);
    subplot(2,3,6); imgray(p_im4);
end
%%
f_list = dir([fibre_root 'images\*.mat']);
m_list = dir([fibre_root 'fibre_masks\*.mat']);
p_list = dir([pred_dir rf_codes{1,1} '\*.mat']);

p_place = [2 3 5 6];
for ii = 1:10
    
    f_im = u_load([fibre_root 'images\' f_list(ii).name]);
    m_im = u_load([fibre_root 'fibre_masks\' m_list(ii).name]);
    
    figure;
        subplot(2,3,1); imgray(f_im);
        subplot(2,3,4); imgray(m_im);
    
    for i_p = 1:4
        p_im = u_load([pred_dir rf_codes{1,1} '\' p_list(ii).name]);
        %p_im = bwmorph(p_im > 0.5, 'thin', 'inf');
        subplot(2,3,p_place(i_p)); imgray(p_im.^4); colormap hot;
    end
        
end
%%
sens = roc_pts(:,1);
spec = 1 - roc_pts(:,2);

    
