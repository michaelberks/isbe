function [] = build_vessel_apex_predictor(varargin)
%BUILD_VESSEL_APEX_PREDICTOR *Insert a one line summary here*
%   [] = extract_vessel_centres()
%
% Inputs:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 18-Jun-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'task_id',              unixenv('SGE_TASK_ID',1), ...
    'num_jobs',             unixenv('NUM_JOBS',100), ...
    'data_dir',             [nailfoldroot 'data/rsa_study/test/'],...
    'feature_im_dir',       'predictions/detection/rf_classification/257273/',...
    'centre_dir',           'vessel_centres/',...
    'hog_dir',              'vessel_hogs/',...
    'max_size',             unixenv('MAX_SIZE', 1000),...
    'smoothing_sigma',      2,...
    'num_cells',            8,...
    'cell_sz',              8,... %Size of HoG cells in blocks
    'block_sz',             [2 2],...%Size of blocks in cells
    'num_ori_bins',         9,... %Number of bins in orientation histograms
    'norm_method',          'l1-sqrt',... %Method for local normalisation
    'block_spacing',        8,... %Separation between blocks in pixels - controls level of overlap
    'gradient_operator',    [-1 0 1],...
    'spatial_sigma',        0, ...
    'angle_wrap',           1,...
    'base_width',           20, ...
    'make_parts_folders',   0);

%Form full directory paths and create folder for HoGs
hog_dir = [args.data_dir args.hog_dir ];

x_list = [hog_dir '*hog_X*'];
y_list = [hog_dir '*hog_y*'];

%
if length(x_list) ~= length(y_list)
    error('Number of output and feature files don''t match');
end

num_ims = length(x_list);

if isempty(args.selected_idx)
    selected_ims = 1:num_ims;
else
    selected_ims = args.selected_ims;
    length(selected_ims);
end
    

%Loop through the images to count how many positive and negative samples we
%have
data_counts = zeros(num_ims,2);
for i_im = selected_ims
    load([hog_dir y_list(i_im).name], 'apex_class');
    data_counts(im_num,1) = sum(apex_class);
    data_counts(im_num,2) = sum(~apex_class);
end
train_n = sum(data_counts(:,1));

%Pre-allocate the outputs
train_c = [true(train_n,1); false(train_n, 1)];
train_y = zeros(train_n, 2);

curr_idx = 0;
for i_im = selected_ims
    
    load([hog_dir y_list(i_im).name], 'apex_class', 'apex_offsets');
    apex_hog = u_load([hog_dir x_list(i_im).name], 'apex_hog');
    
    if i_im == 1;
        train_X = zeros(2*train_n, size(apex_hog,2));
    end
    
    %Get indices for where we'll place the positive and negative data
    idx = curr_idx + (1:data_counts(i_im,1));
    curr_idx = idx(end);

    n_idx = find(~apex_class);
    r_idx = randperm(data_counts(i_im,2));
    n_idx = n_idx(r_idx(1:data_counts(i_im,1)));

    %Sample the hog features for both pos and neg samples, and the apex
    %offsets for the pos samples
    train_X(idx,:) = apex_hog(apex_class,:);
    train_X(idx + train_n,:) = apex_hog(n_idx,:);
    train_y(idx,:) = apex_offsets(apex_class,:); %#ok
    
end

%Now build the predictors   
warning('off', 'ASYM:unexpectedArgument');

%The pos/neg classifier
rf_args.prediction_type = 'rf_classification';
rf_args.n_trees = 100;
rf_args.d = [];
rf_args.w_prior = 0;
rf_args.impure_thresh = 1.0000e-004;
rf_args.split_min = 100;
rf_args.end_cut_min = 25;
rf_args.do_ubound = 0;
rf_args.quiet = 1;
rf_args.overwrite = 0;
rf_args.minimise_size = 0;
rf_args.split_criterion = 'gdi';
rf_args.var_criterion = 'mabs';

rf_args.sampling_args.sampling_method = 'sample_saved_training_data';
rf_args.decomposition_args = [];

rf_dir = ['C:\isbe\nailfold\models\vessel\apex_location\class' datestr(now, 30) '\'];
rf_args.tree_dir = [rf_dir 'trees/'];            

rf_args.sampling_args.y = train_c;
rf_args.sampling_args.X = train_X;
apex_class_rf = random_forest_class_train(rf_args);   
    
%The regressors for the apex offset
rf_args.prediction_type = 'rf_regression';
rf_args.n_trees = 100;
rf_args.d = [];
rf_args.w_prior = 0;
rf_args.impure_thresh = 1.0000e-008;
rf_args.split_min = 100;
rf_args.end_cut_min = 0;
rf_args.do_ubound = 0;
rf_args.quiet = 1;
rf_args.do_circular = [];
rf_args.overwrite = 1;
rf_args.minimise_size = 0;
rf_args.split_criterion = 'ssq';
rf_args.var_criterion = 'ssq';

rf_args.sampling_args.sampling_method = 'sample_saved_training_data';
rf_args.decomposition_args = [];

rf_args.sampling_args.X = train_X(1:train_n,:);

rf_dir = ['C:\isbe\nailfold\models\vessel\apex_location\offset_x\' datestr(now, 30) '\'];
rf_args.tree_dir = [rf_dir 'trees/'];            
rf_args.sampling_args.y = train_y(:,1);
apex_offset_x_rf = random_forest_reg_train(rf_args); 

rf_dir = ['C:\isbe\nailfold\models\vessel\apex_location\offset_y\' datestr(now, 30) '\'];
rf_args.tree_dir = [rf_dir 'trees/'];            
rf_args.sampling_args.y = train_y(:,2);
apex_offset_y_rf = random_forest_reg_train(rf_args);

%Load the trees into the forests and save

    
    
