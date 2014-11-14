function [] = train_hog_apex_detectors(varargin)

args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'data_path',            [],...
    'data_dir',             [nailfoldroot 'data/rsa_study/set12g/'],...
    'model_id',             [],...
    'model_root',           [unixenv('DATA_ROOT',[]) unixenv('MODEL_ROOT',[nailfoldroot,'models/apex'])], ...
    'model_name',           'rf',...
    'num_pos_samples',      1e4,...
    'neg_ratio',            1,...
    'num_trees',            unixenv('NUM_TREES',100), ...
    'selected_ims',         [],...
    'labels_dir',           'apex_class_data',...
    'hog_dir',              'vessel_hogs',...
    'num_cells',            8,...
    'cell_sz',              8,... %Size of HoG cells in blocks
    'block_sz',             [2 2],...%Size of blocks in cells
    'num_ori_bins',         9,... %Number of bins in orientation histograms
    'norm_method',          'l1-sqrt',... %Method for local normalisation
    'block_spacing',        8,... %Separation between blocks in pixels - controls level of overlap
    'gradient_operator',    [-1 0 1],...
    'spatial_sigma',        0, ...
    'angle_wrap',           1,...
    'base_width',           20,...
    'dist_thresh',          24,...
    'save_data',            0);
    
%Set up file lists
labels_dir = [args.data_dir '/' args.labels_dir '/'];
labels_list = dir([labels_dir '*.mat']);
hog_dir = [args.data_dir '/' args.hog_dir '/'];
hog_list = dir([hog_dir '*.mat']);
num_ims = length(labels_list);

if isempty(args.model_id)
    model_id = datestr(now, 30);
else
    model_id = args.model_id;
end
rf_class_dir = [args.model_root '/classification/' model_id '/'];
rf_offset_x_dir = [args.model_root '/offset_x/' model_id '/'];
rf_offset_y_dir = [args.model_root '/offset_y/' model_id '/'];

if ~isempty(args.data_path)
    load(args.data_path, 'train_X', 'train_c', 'train_o', 'samples_per_image');
    n_pos_so_far = sum(samples_per_image(:,1)); %#ok
else

    %Don't worry, if we've asked for too much data, the full size storage will
    %be created, but we'll only sample what we can from each image, then discard the excess
    %space at the end
    num_pos_samples = args.num_pos_samples;
    num_neg_samples = round(num_pos_samples*args.neg_ratio);
    n_pos_so_far = 0;
    n_neg_so_far = 0;
    idx_p = 0; %Staring index for pos samples
    idx_n = num_pos_samples; %Starting index for neg samples

    samples_per_image = zeros(num_ims, 2);
    
    pos_samples_count = zeros(num_ims,2);
    for i_im = 1:num_ims

        display(['Processing image ' num2str(i_im) ' of ' num2str(num_ims)]);

        %Load the HoG data and class labels/apex offsets
        load([labels_dir labels_list(i_im).name], 'apex_class', 'apex_offsets');
        apex_hog = u_load([hog_dir hog_list(i_im).name]);

        %For the first image, create storage containers for all the data
        if i_im == 1
            hog_size = size(apex_hog, 2);
            train_X = zeros(num_pos_samples+num_neg_samples, hog_size);
            train_c = [true(num_pos_samples,1); false(num_neg_samples, 1)];
            train_o = zeros(num_pos_samples, 2);
            train_idx = zeros(num_pos_samples+num_neg_samples, 2);
        end   

        %Randomyl sample the number of positive and negative samples to take
        %from this image
        pos_samples_i = compute_samples_per_image(num_pos_samples, n_pos_so_far, num_ims, i_im);
        %neg_samples_i = compute_samples_per_image(num_neg_samples, n_neg_so_far, num_ims, i_im);

        %Make sure these don't exceed the #samples we actually have
        total_pos_i = sum(apex_class);
        total_neg_i = sum(~apex_class);
        pos_samples_i = min(pos_samples_i, total_pos_i);
        %neg_samples_i = min(neg_samples_i, total_neg_i);
        neg_samples_i = min(pos_samples_i, total_neg_i);
        n_pos_so_far = n_pos_so_far + pos_samples_i;
        n_neg_so_far = n_neg_so_far + neg_samples_i;
        samples_per_image(i_im,:) = [pos_samples_i neg_samples_i];
        pos_samples_count(i_im,:) = [pos_samples_i, total_pos_i];
        
        %Update the indices both in the storage and where to sample from the
        %image
        [idx_p rand_idx_p] = ...
            update_sample_indices(idx_p, pos_samples_i, apex_class, total_pos_i);
        [idx_n rand_idx_n] = ...
            update_sample_indices(idx_n, neg_samples_i, ~apex_class, total_neg_i);

        %sample the data/labels into the storage containers
        train_X(idx_p,:) = apex_hog(rand_idx_p,:);
        train_X(idx_n,:) = apex_hog(rand_idx_n,:);
        train_o(idx_p,:) = apex_offsets(rand_idx_p,:); %#ok

        %Also keep a record of which points we've sample
        train_idx(idx_p,1) = i_im;
        train_idx(idx_n,1) = i_im;
        train_idx(idx_p,2) = rand_idx_p;
        train_idx(idx_n,2) = rand_idx_n;
    end

    %Discard any points we've not filled
    missing_pos = (n_pos_so_far+1):num_pos_samples;
    missing_neg = (num_pos_samples+n_neg_so_far+1):(num_pos_samples+num_neg_samples);
    train_X([missing_pos missing_neg],:) = [];
    train_c([missing_pos missing_neg],:) = [];
    train_idx([missing_pos missing_neg],:) = []; %#ok
    train_o(missing_pos,:) = [];

    %Save the data if we've been given a save directory
    if args.save_data
        rf_data_dir = [args.model_root '/training_data/' model_id '/'];
        create_folder(rf_data_dir);
        save([rf_data_dir 'data.mat'], 'train_X', 'train_c', 'train_o', 'train_idx', 'samples_per_image');
        clear train_idx; %we don't need this again
    end
end

%Now build some forests with this data!
warning('off', 'ASYM:unexpectedArgument');
rf_class_args.prediction_type = 'rf_classification';
rf_class_args.n_trees = args.num_trees;
rf_class_args.d = [];
rf_class_args.w_prior = 0;
rf_class_args.impure_thresh = 1.0000e-004;
rf_class_args.split_min = 100;
rf_class_args.end_cut_min = 25;
rf_class_args.do_ubound = 0;
rf_class_args.quiet = 1;
rf_class_args.overwrite = 0;
rf_class_args.minimise_size = 0;
rf_class_args.split_criterion = 'gdi';
rf_class_args.var_criterion = 'mabs';

rf_class_args.sampling_args.sampling_method = 'sample_saved_training_data';
rf_class_args.decomposition_args = [];

%Save the RF classification args (before we add in the data)
create_folder(rf_class_dir);
save([rf_class_dir 'rf_args.mat'], 'rf_class_args');
rf_class_args.tree_dir = [rf_class_dir 'trees/']; 

%Add the data/labels and train the forest
rf_class_args.sampling_args.y = train_c;
rf_class_args.sampling_args.X = train_X;
apex_class_rf = random_forest_class_train(rf_class_args);

%Now we can clear the class args from memory, along with the negative data
clear rf_class_args train_c;
train_X(n_pos_so_far+1:end,:) = [];

%Create args for the regressor to predict apex offsets - we'll do this separately for
% x and y offsets
rf_off_args.prediction_type = 'rf_regression';
rf_off_args.n_trees = args.num_trees;
rf_off_args.d = [];
rf_off_args.w_prior = 0;
rf_off_args.impure_thresh = 1.0000e-008;
rf_off_args.split_min = 100;
rf_off_args.end_cut_min = 0;
rf_off_args.do_ubound = 0;
rf_off_args.quiet = 1;
rf_off_args.do_circular = [];
rf_off_args.overwrite = 1;
rf_off_args.minimise_size = 0;
rf_off_args.split_criterion = 'ssq';
rf_off_args.var_criterion = 'ssq';

rf_off_args.sampling_args.sampling_method = 'sample_saved_training_data';
rf_off_args.decomposition_args = [];

%Save the RF regression args (before we add in the data)
create_folder(rf_offset_x_dir);
create_folder(rf_offset_y_dir);
save([rf_offset_x_dir 'rf_args.mat'], 'rf_off_args');
save([rf_offset_y_dir 'rf_args.mat'], 'rf_off_args');

%Add the features to args
rf_off_args.sampling_args.X = train_X;

%Do the x (horizontal) offsets
rf_off_args.tree_dir = [rf_offset_x_dir 'trees/'];            
rf_off_args.sampling_args.y = train_o(:,1);
apex_offset_x_rf = random_forest_reg_train(rf_off_args);

%Do the y (vertical) offsets
rf_off_args.tree_dir = [rf_offset_y_dir 'trees/'];            
rf_off_args.sampling_args.y = train_o(:,2);
apex_offset_y_rf = random_forest_reg_train(rf_off_args);

%Load the trees into each RF structure
for i_tree = 1:args.num_trees
    apex_class_rf.trees{i_tree} = ...
        u_load([apex_class_rf.tree_root apex_class_rf.tree_dir apex_class_rf.trees{i_tree}]);
    apex_offset_x_rf.trees{i_tree} = ...
        u_load([apex_offset_x_rf.tree_root apex_offset_x_rf.tree_dir apex_offset_x_rf.trees{i_tree}]);
    apex_offset_y_rf.trees{i_tree} = ...
        u_load([apex_offset_y_rf.tree_root apex_offset_y_rf.tree_dir apex_offset_y_rf.trees{i_tree}]);
end
    
%Save the RFs
save([rf_class_dir args.model_name '.mat'], 'apex_class_rf');
save([rf_offset_x_dir args.model_name '.mat'], 'apex_offset_x_rf');
save([rf_offset_y_dir args.model_name '.mat'], 'apex_offset_y_rf');

function samples_per_image = compute_samples_per_image(n_samples, n_sampled_so_far, n_images, curr_image)

N = n_samples - n_sampled_so_far;
p = 1 / (n_images - curr_image + 1);
samples_per_image = sample_from_binomial(N, p, 1);

function [idx rand_idx] = update_sample_indices(idx, n_samples, class_labels, n_data)

%Update the storage indices
idx = idx(end) + (1:n_samples);
    
%Get a random sample of indices from the data
rand_idx = find(class_labels);
rp = randperm(n_data);
rand_idx = rand_idx(rp(1:n_samples));