%Script, examining the results of various tests applied to the training set 
% of DRIVE retinograms
%--------------------------------------------------------------------------
retroot = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\training\';
fov_mask_dir = [retroot 'fov_masks\'];
fg_mask_dir = [retroot 'vessel_masks\'];
%%
gt_orientations = [];
vessel_centres = true(0,1);
for ii = 21:40
    gt_ori = load_uint8([retroot '/orientations/' zerostr(ii,2) '_ori1.mat']);
    fg_mask = u_load([fg_mask_dir zerostr(ii,2) '_training_v_mask.mat']);
    fov_mask = u_load([fov_mask_dir zerostr(ii,2) '_training_f_mask.mat']);
    fg_mask = fg_mask & fov_mask;
    centre_mask = bwmorph(fg_mask, 'thin', 'inf');
    
    gt_orientations = [gt_orientations; gt_ori(fg_mask)]; %#ok
    vessel_centres = [vessel_centres; centre_mask(fg_mask)]; %#ok
end
mkdir([retroot 'orientations\gt\']);
save([retroot 'orientations\gt\all_gt_orientations.mat'], 'gt_orientations', 'vessel_centres');
%%
%--------------------------------------------------------------------------
% Detection...
pred_dir = [retroot 'predictions\detection\rf_classification\'];
label_dir = [retroot 'vessel_masks\'];

rf_codes = cell(20, 4);

rf_codes( 1,:) = {'8895', 'dt', '3', 'orig'};
rf_codes( 2,:) = {'8893', 'g2', '3', 'orig'};
rf_codes( 3,:) = {'11476', 'g12d', '3', 'orig'};
rf_codes( 4,:) = {'11477', 'mono', '3', 'orig'};
rf_codes( 5,:) = {'11478', 'linop', '3', 'orig'};
rf_codes( 6,:) = {'11479', 'haar', '3', 'orig'};
rf_codes( 7,:) = {'24066', 'gabor_a', '3', 'orig'}; %Old version of filters '11480'
rf_codes( 8,:) = {'38797', 'g2di', '3', 'orig'};
rf_codes( 9,:) = {'24222', 'g1d', '3', 'orig'};
rf_codes(10,:) = {'38798', 'gabor_i', '3', 'orig'};

rf_codes(11,:) = {'8894', 'dt', '1', 'orig'};
rf_codes(12,:) = {'8892', 'g2', '1', 'orig'};
rf_codes(13,:) = {'11476', 'g12d', '1', 'orig'};
rf_codes(14,:) = {'11477', 'mono', '1', 'orig'};
rf_codes(15,:) = {'11478', 'linop', '1', 'orig'};
rf_codes(16,:) = {'11479', 'haar', '1', 'orig'};
rf_codes(17,:) = {'24066', 'gabor_a', '1', 'orig'}; %Old version of filters '11480'
rf_codes(18,:) = {'24221', 'g2di', '1', 'orig'};
rf_codes(19,:) = {'24222', 'g1d', '1', 'orig'};
rf_codes(20,:) = {'24336', 'gabor_i', '1', 'orig'};


rf_codes(21,:) = {'24072', 'gabor_r', '3', 'orig'};
rf_codes(22,:) = {'33429', 'dt_r', '3', 'orig'};

rf_codes(23,:) = {'34237', 'g2d_5', '5', 'orig'};
rf_codes(24,:) = {'34015', 'dt_r', '3', 'orig'};

%rf_codes(22,:) = {'29920', 'gabor_m', '3', 'orig'};

% rf_codes(32,:) = {'30345', 'gabor_1', '3', 'orig'};
% rf_codes(33,:) = {'30346', 'gabor_2', '3', 'orig'};
% rf_codes(34,:) = {'30347', 'gabor_3', '3', 'orig'};
% rf_codes(35,:) = {'30348', 'gabor_4', '3', 'orig'};
% rf_codes(36,:) = {'30349', 'gabor_5', '3', 'orig'};
% rf_codes(37,:) = {'33429', 'dt_r', '3', 'orig'};
% rf_codes(38,:) = {'33430', 'dt_m', '3', 'orig'};
% 
% rf_codes(11,:) = {'24336', 'g2da', '3', 'orig'};
% rf_codes( 5,:) = {'9151', 'g2', '1', '0.25'};
% rf_codes( 6,:) = {'9152', 'g2', '3', '0.25'};
% rf_codes( 7,:) = {'9153', 'dt', '1', '0.25'};
% rf_codes( 8,:) = {'9154', 'dt', '3', '0.25'};
% rf_codes( 9,:) = {'9155', 'g2', '1', '0.50'};
% rf_codes(10,:) = {'9156', 'g2', '3', '0.50'};
% rf_codes(11,:) = {'9157', 'dt', '1', '0.50'};
% rf_codes(12,:) = {'9158', 'dt', '3', '0.50'};
% rf_codes(13,:) = {'9159', 'g2', '1', '0.75'};
% rf_codes(14,:) = {'9160', 'g2', '3', '0.75'};
% rf_codes(15,:) = {'9161', 'dt', '1', '0.75'};
% rf_codes(16,:) = {'9162', 'dt', '3', '0.75'};
% rf_codes(17,:) = {'9473', 'g2', '1', '1.00'};
% rf_codes(18,:) = {'9474', 'g2', '3', '1.00'};
% rf_codes(19,:) = {'9475', 'dt', '1', '1.00'};
% rf_codes(20,:) = {'9476', 'dt', '3', '1.00'};

% rf_codes(26,:) = {'12393', 'g12d', '3', '0.75'};
% rf_codes(27,:) = {'12394', 'mono', '3', '0.75'};
% rf_codes(28,:) = {'12395', 'linop', '3', '0.75'};
% rf_codes(29,:) = {'12396', 'haar', '3', '0.75'};
% rf_codes(30,:) = {'24303', 'gabor_a', '3', '0.75'}; %Old version of filters '12397'
% rf_codes(31,:) = {'12595', 'g2dg', '3', 'orig'};
% rf_codes(32,:) = {'13508', 'g2da', '1', 'orig'};
% rf_codes(33,:) = {'13509', 'g2da', '3', 'orig'};
% rf_codes(34,:) = {'13569', 'g2da_r', '1', 'orig'};
% rf_codes(35,:) = {'13570', 'g2da_m', '1', 'orig'};
% rf_codes(36,:) = {'13571', 'g2da_r', '3', 'orig'};
% rf_codes(37,:) = {'13572', 'g2da_m', '3', 'orig'};

roc_pts_all = [];
auc_all = [];
%%
for ii = [3 13]%30:36%[1:4 21:30]
    if exist([pred_dir rf_codes{ii,1} '\roc\roc_stats.mat'], 'file')
        load([pred_dir rf_codes{ii,1} '\roc\roc_stats.mat'])
    else
        [roc_pts, auc, auc_ind] = compute_roc_image_set([pred_dir rf_codes{ii,1} '\'], label_dir, 'fov_mask_dir', fov_mask_dir);
        create_folder([pred_dir rf_codes{ii,1} '\roc\']);
        save([pred_dir rf_codes{ii,1} '\roc\roc_stats.mat'], 'roc_pts', 'auc', 'auc_ind');
    end
    roc_pts_all(:,:,ii) = roc_pts; %#ok
    auc_all(ii,1) = auc; %#ok
    auc_all(ii,2) = std(auc_ind); %#ok
end
%
rfs = find(auc_all(:,1));
for i_rf = rfs'
    display(['Decomp: ' rf_codes{i_rf, 2} ', w = ' rf_codes{i_rf,3} ', sampling mode: ' rf_codes{i_rf, 4}...
        ', Az = ' num2str(auc_all(i_rf,1),3) ' +/- ' num2str(auc_all(i_rf,2),3)]);
end

%%
%--------------------------------------------------------------------------
%Centre orientation
warning('off', 'ori_error:nans');
pred_dir = [retroot 'predictions\centre_orientation\rf_regression\'];
label_dir = [retroot 'orientations\'];

load([label_dir 'all_gt_orientations.mat'], 'gt_orientations', 'vessel_centres');
    
rf_codes = cell(20, 4);
rf_codes( 1,:) = {'8298', 'g2', '1', 'orig'};
rf_codes( 2,:) = {'8321', 'g2', '3', 'orig'};
rf_codes( 3,:) = {'8328', 'dt', '1', 'orig'};
rf_codes( 4,:) = {'8331', 'dt', '3', 'orig'};
rf_codes( 5,:) = {'9832', 'dt', '3', '0.25'};
rf_codes( 6,:) = {'9976', 'dt', '3', '0.50'};
rf_codes( 7,:) = {'9977', 'dt', '3', '0.75'};
rf_codes( 8,:) = {'9978', 'dt', '3', '1.00'};
rf_codes( 9,:) = {'10035', 'g2a', '1', 'orig'};
rf_codes(10,:) = {'10204', 'g2a', '3', 'orig'};

%
for ii = [2 4 10]

    display(['Errors for ' rf_codes{ii,2} ', w = ' rf_codes{ii,3} ', sampling = ' rf_codes{ii,4}]); 
    [ori_errors, ~, error_stats] =...
        compute_image_orientation_errors([pred_dir rf_codes{ii,1} '\'], fg_mask_dir,...
        'gt_orientations', gt_orientations, 'fov_mask_dir', fov_mask_dir);
    display(error_stats);
    
    centre_stats = ori_error_stats(ori_errors(vessel_centres));
    display(centre_stats);
end
%%
%--------------------------------------------------------------------------
%All vessel orientation
warning('off', 'ori_error:nans');
pred_dir = [retroot 'predictions\orientation\rf_regression\'];
label_dir = [retroot 'orientations\'];

load([label_dir 'gt\all_gt_orientations.mat'], 'gt_orientations', 'vessel_centres');
    
rf_codes = cell(7, 4);
rf_codes( 1,:) = {'12242', 'dt', '3', 'orig'};
rf_codes( 2,:) = {'12243', 'g2', '3', 'orig'};
rf_codes( 3,:) = {'12244', 'g12', '3', 'orig'};
rf_codes( 4,:) = {'12245', 'mono', '3', 'orig'};
rf_codes( 5,:) = {'12246', 'linop', '3', 'orig'};
rf_codes( 6,:) = {'12247', 'haar', '3', 'orig'};
rf_codes( 7,:) = {'24069', 'gabor', '3', 'orig'};
rf_codes( 8,:) = {'24293', 'g1d', '3', 'orig'};
rf_codes( 9,:) = {'38788', 'g2di', '3', 'orig'};
rf_codes(10,:) = {'38789', 'gabori', '3', 'orig'};

rf_codes(11,:) = {'38809', 'dt', '1', 'orig'};
rf_codes(12,:) = {'38810', 'g2', '1', 'orig'};
rf_codes(13,:) = {'38811', 'g12', '1', 'orig'};
rf_codes(14,:) = {'38812', 'mono', '1', 'orig'};
rf_codes(15,:) = {'38813', 'linop', '1', 'orig'};
rf_codes(16,:) = {'38814', 'haar', '1', 'orig'};
rf_codes(17,:) = {'38815', 'gabor', '1', 'orig'};
rf_codes(18,:) = {'38816', 'g1d', '1', 'orig'};
rf_codes(19,:) = {'38817', 'g2di', '1', 'orig'};
rf_codes(20,:) = {'38818', 'gabori', '1', 'orig'};

% rf_codes(11,:) = {'30469', 'gabor_1', '3', 'orig'};
% rf_codes(12,:) = {'30470', 'gabor_2', '3', 'orig'};
% rf_codes(13,:) = {'30471', 'gabor_4', '3', 'orig'};
% rf_codes(14,:) = {'30472', 'gabor_8', '3', 'orig'};
% rf_codes(15,:) = {'30473', 'gabor_16', '3', 'orig'};


%
error_medians = zeros(9,2);
error_means = zeros(9,2);
%%
for ii = [3 13]

    if exist([pred_dir rf_codes{ii,1} '\errors\error_stats.mat'], 'file')
        load([pred_dir rf_codes{ii,1} '\errors\error_stats.mat']);
    else
        [orientation_errors, ~, error_stats] =...
            compute_image_orientation_errors([pred_dir rf_codes{ii,1} '\'], fg_mask_dir,...
            'gt_orientations', gt_orientations, 'fov_mask_dir', fov_mask_dir);
    end
    centre_stats = ori_error_stats(orientation_errors(vessel_centres));
    
    error_medians(ii,1) = error_stats.abs_median;
    error_medians(ii,2) = centre_stats.abs_median;
    
    error_means(ii,1) = error_stats.abs_mean;
    error_means(ii,2) = centre_stats.abs_mean;
    
    display(['Errors for ' rf_codes{ii,2} ', w = ' rf_codes{ii,3} ', sampling = ' rf_codes{ii,4}...
        ', abs mean: ' num2str(error_stats.abs_mean,4) ', abs_median: ' num2str(error_stats.abs_median,4)]); 
       
end

%%
% Old (non-normalised) version of gabor filters
% rf_codes(32,:) = {'12787', 'gabor_ag', '3', 'orig'};%Old version of filters
% rf_codes(25,:) = {'11480', 'gabor_a', '3', 'orig'}; %Old version of filters
% rf_codes(30,:) = {'12397', 'gabor_a', '3', '0.75'}; %Old version of filters
% rf_codes(39,:) = {'14285', 'gabor_r', '1', 'orig'};%Old version of filters
% rf_codes(40,:) = {'14286', 'gabor_r', '3', 'orig'};%Old version of filters
% rf_codes(41,:) = {'14287', 'gabor_rr', '1', 'orig'};%Old version of filters
% rf_codes(42,:) = {'14288', 'gabor_rm', '1', 'orig'};%Old version of filters
% rf_codes(43,:) = {'14289', 'gabor_rr', '3', 'orig'};%Old version of filters
% rf_codes(44,:) = {'14290', 'gabor_rm', '3', 'orig'};%Old version of filters
% rf_codes(45,:) = {'14831', 'gabor_a', '1', 'orig'};%Old version of filters
% rf_codes(46,:) = {'14832', 'gabor_ar', '1', 'orig'};%Old version of filters
% rf_codes(47,:) = {'14833', 'gabor_am', '1', 'orig'};%Old version of filters
% rf_codes(48,:) = {'14834', 'gabor_ar', '3', 'orig'};%Old version of filters
% rf_codes(49,:) = {'14835', 'gabor_am', '3', 'orig'};%Old version of filters