function [] = train_hog_apex_rescorer(varargin)

args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    {'image_names'},        ...
    'data_path',            [],...
    'data_dir',             [nailfoldroot 'data/rsa_study/master_set/'],...
    'model_id',             [],...
    'model_root',           [unixenv('DATA_ROOT',[]) unixenv('MODEL_ROOT',[nailfoldroot,'models/apex'])], ...
    'model_name',           'rf',...
    'num_trees',            unixenv('NUM_TREES',100), ...
    'labels_dir',           'labels',...
    'hog_dir',              'apex_hogs',...
    'save_data',            0);
    
%Set up file lists
labels_dir = [args.data_dir '/' args.labels_dir '/'];
hog_dir = [args.data_dir '/' args.hog_dir '/'];

if isempty(args.model_id)
    model_id = datestr(now, 30);
else
    model_id = args.model_id;
end
rf_class_dir = [args.model_root '/rescoring/' model_id '/'];


if ~isempty(args.data_path)
    load(args.data_path, 'train_X', 'train_c');
else
    
    num_images = length(args.image_names);
    num_pos_samples = 0;
    num_neg_samples = 0;
    for i_im = 1:num_images
        load([labels_dir args.image_names{i_im} '_label.mat'], 'candidates_class');
        num_pos_samples = num_pos_samples + sum(candidates_class);
        num_neg_samples = num_neg_samples + sum(~candidates_class);
    end
    num_samples = num_pos_samples + num_neg_samples;
    
    idx = 0;
    for i_im = 1:num_images

        display(['Processing image ' num2str(i_im) ' of ' num2str(num_images)]);

        %Load the HoG data and class labels/apex offsets
        load([labels_dir args.image_names{i_im} '_label.mat'], 'candidates_class');
        candidate_hog = u_load([hog_dir args.image_names{i_im} '_hog.mat']);       

        %For the first image, create storage containers for all the data
        if i_im == 1
            hog_size = size(candidate_hog, 2);
            train_X = zeros(num_samples, hog_size);
            train_c = false(num_samples,1);
            train_idx = zeros(num_samples,2);
        end 
        
        %sample the data/labels into the storage containers
        num_candidates = length(candidates_class);
        im_idx = 1:num_candidates;
        idx = idx(end) + im_idx;
                
        train_X(idx,:) = candidate_hog;
        train_c(idx,:) = candidates_class;
        train_idx(idx,1) = i_im;
        train_idx(idx,2) = im_idx;
    end

    %Save the data if we've been given a save directory
    if args.save_data
        rf_data_dir = [args.model_root '/training_data/' model_id '/'];
        create_folder(rf_data_dir);
        save([rf_data_dir 'data.mat'], 'train_X', 'train_c', 'train_idx');
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

%Load the trees into each RF structure
for i_tree = 1:args.num_trees
    apex_class_rf.trees{i_tree} = ...
        u_load([apex_class_rf.tree_root apex_class_rf.tree_dir apex_class_rf.trees{i_tree}]);
end
    
%Save the RFs
save([rf_class_dir args.model_name '.mat'], 'apex_class_rf');