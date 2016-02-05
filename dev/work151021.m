base_dir = 'A:\D_2.2_Datasets\125_375_Single_MLO_Dataset_VOLPARA\';
control_list = dir([base_dir 'CONTROLS\ASSURE_CONTROLS*']);
cancer_list = dir([base_dir 'CANCERS\ASSURE_CANCERS*']);

num_controls = length(control_list);
num_cancers = length(cancer_list);
%
control_im_names = cell(num_controls,1);
control_mask_names = cell(num_controls,1);
for i_case = 1:num_controls
    im_name = dir([base_dir 'CONTROLS\' control_list(i_case).name '\*ML*hint_densityMap*.pgm']);
    mask_name = dir([base_dir 'CONTROLS\' control_list(i_case).name '\*ML*hint_Segmentation.pgm']);
    
    control_im_names{i_case} = [base_dir 'CONTROLS\' control_list(i_case).name '\' im_name(end).name];
    control_mask_names{i_case} = [base_dir 'CONTROLS\' control_list(i_case).name '\' mask_name(end).name];
end

cancer_im_names = cell(num_cancers,1);
cancer_mask_names = cell(num_cancers,1);
for i_case = 1:num_cancers
    im_name = dir([base_dir 'CANCERS\' cancer_list(i_case).name '\*ML*hint_densityMap*.pgm']);
    mask_name = dir([base_dir 'CANCERS\' cancer_list(i_case).name '\*ML*hint_Segmentation.pgm']);
    
    cancer_im_names{i_case} = [base_dir 'CANCERS\' cancer_list(i_case).name '\' im_name(end).name];
    cancer_mask_names{i_case} = [base_dir 'CANCERS\' cancer_list(i_case).name '\' mask_name(end).name];
end
    
case_im_names = [control_im_names; cancer_im_names];
case_mask_names = [control_mask_names; cancer_mask_names];

num_cases = length(case_im_names);
rand_i = randperm(num_cases);

num_train = num_cases*0.8;
num_test = num_cases - num_train;
train_case_idx = rand_i(1:num_train);
test_case_idx = rand_i(num_train+1:end);
%%
pts_per_img = 100;
total_samples = num_train*pts_per_img;

all_responses = zeros(total_samples, 12);
xy = zeros(total_samples, 2);
curr_sample = 0;

decomposition_args.decomp_type = 'dt';      %Use the dual-tree
decomposition_args.win_size = 1;            %Window size about pixel to sample features
decomposition_args.levels = 6;            %Levels to include
decomposition_args.feature_shape = 'rect';  %Keep this as rect, although there are other things you can play with (see sample_dt_data)
decomposition_args.feature_type = 'conj';   %Again, 'conj' works best, but can use 'mag' for magnitude only, or 'phase' (see convert_complex_representation)
decomposition_args.do_max = 0;              %Keep all 6 sub-bands (0) or only keep the maximum (1)
decomposition_args.rotate = 0;              %Try and make rotation invariant (best left switched off)
decomposition_args.use_nag = 0;             %Use the NAG toolbox to interpolate (best left swicthed off)
decomposition_args.normalise = 0;

for i_idx = 1:num_train
    i_case = train_case_idx(i_idx);
    
    mam_vdm = imread(case_im_names{i_case});
    mask = imread(case_mask_names{i_case});
    dual_tree = compute_dual_tree(mam_vdm, 6);       

    [rows cols] = find(mask==4);
    %[rows cols] = find(mask>-1);
    rand_i = randperm(length(rows),pts_per_img);
    rows = rows(rand_i);
    cols = cols(rand_i);
    sample_idx = curr_sample+(1:pts_per_img);

    [responses] = compute_filter_responses(mam_vdm, decomposition_args);
    all_responses(sample_idx,:) = sample_image_features(responses, rows, cols, decomposition_args); 
    xy(sample_idx,1) = cols;
    xy(sample_idx,2) = rows;      
    curr_sample = curr_sample + pts_per_img;
end
%
opts = statset('MaxIter', 1000);
num_k = 8;
[k_idx, k_centres, sumd] = kmeans(all_responses, num_k, 'EmptyAction', 'drop', 'Replicates', 10, 'Options', opts);
save C:\isbe\density\Assure\k_compressed.mat k_idx k_centres sumd all_responses 
%save C:\isbe\density\Assure\k_all.mat k_idx k_centres sumd all_responses 
%%
colors = jet(num_k);
%
for i_idx = 1:20
    i_case = test_case_idx(i_idx);
    
    mam_vdm = imread(case_im_names{i_case});
    mask = imread(case_mask_names{i_case});
    
    dual_tree = compute_dual_tree(mam_vdm, 6);

    [rows cols] = find(mask > -1);
    [responses] = compute_filter_responses(mam_vdm, decomposition_args);
    sampled_features = sample_image_features(responses, rows, cols, decomposition_args);
    %
    k_dists = zeros(size(sampled_features,1), num_k);
    for i_k = 1:num_k
        k_dists(:,i_k) = sum(bsxfun(@minus, sampled_features, k_centres(i_k,:)).^2,2);
    end
    [~,assigned_k] = min(k_dists,[],2);
    
    b_inner = bwboundaries(mask == 4);
    b_outer = bwboundaries(mask > 0);

    %
    figure; 
    subplot(1,2,1); imgray(mam_vdm);
    subplot(1,2,2); axis equal ij; hold all;
    for i_k = 1:num_k
        idx = assigned_k == i_k;
        plot(cols(idx), rows(idx), '.', 'markeredgecolor', colors(i_k,:));
    end
    plot(b_inner{1}(:,2), b_inner{1}(:,1), 'k', 'linewidth', 2);
    plot(b_outer{1}(:,2), b_outer{1}(:,1), 'k', 'linewidth', 2);
    
    set(gca, 'xlim', [0 size(mam_vdm,2)], 'ylim', [0 size(mam_vdm,1)]);
    exportfig(['C:\isbe\density\assure\figs\k_compressed\test' zerostr(i_idx,2) '.png']);
end
%%
decomposition_args.decomp_type = 'dt';      %Use the dual-tree
decomposition_args.win_size = 1;            %Window size about pixel to sample features
decomposition_args.levels = 1:6;            %Levels to include
decomposition_args.feature_shape = 'rect';  %Keep this as rect, although there are other things you can play with (see sample_dt_data)
decomposition_args.feature_type = 'conj';   %Again, 'conj' works best, but can use 'mag' for magnitude only, or 'phase' (see convert_complex_representation)
decomposition_args.do_max = 0;              %Keep all 6 sub-bands (0) or only keep the maximum (1)
decomposition_args.rotate = 0;              %Try and make rotation invariant (best left switched off)
decomposition_args.use_nag = 0;             %Use the NAG toolbox to interpolate (best left swicthed off)
decomposition_args.normalise = 0;

mkdir c:\isbe\density\assure\375_125\images
mkdir c:\isbe\density\assure\375_125\full_masks
mkdir c:\isbe\density\assure\375_125\compressed_masks
mkdir c:\isbe\density\assure\375_125\label_masks

for i_case = 1:num_cases
    
    mam_vdm = imread(case_im_names{i_case});
    mask = imread(case_mask_names{i_case});   
    
    [~, fname, ext] = fileparts(case_im_names{i_case});
    
%     dual_tree = compute_dual_tree(mam_vdm, 6);
% 
%     [rows cols] = find(mask == 4);
%     [responses] = compute_filter_responses(mam_vdm, decomposition_args);
%     samples_X = sample_image_features(responses, rows, cols, decomposition_args);
%     if i_case <= num_controls
%         samples_y = false(size(rows));
%     else
%         samples_y = true(size(rows));
%     end
%         
%     sample_name = [];
%     save(sample_name, 'samples_X', 'samples_y');
end
%%
im1 = imread(case_im_names{1});
mask1 = imread(case_mask_names{1});
mask1(mask1==4) = 2;

im2 = imread(case_im_names{2});
mask2 = imread(case_mask_names{2});
mask2(mask2==4) = 2;

points_mask = false(size(mask1));
points_mask(32:32:end,32:32:end) = 1;
points_mask(~(mask1==2)) = 0;
figure; imgray(points_mask);
%points_mask(381, 193) = 1;
%%
[bcs1, bcs2] = BreastCoodSystem.imageCorrespondentPoints(im1, im2, mask1, mask2, find(points_mask), 1, 0, 2, 3);