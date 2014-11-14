function [random_forest abnormal_roi_abnormalvotes normal_roi_abnormalvotes abnormal_roi_totalvotes normal_roi_totalvotes] = mb_random_forest_image_train(varargin)

%MB_RANDOM_FOREST_CLASS_TRAIN given input data and a set of class labels,
%builds a random forest for classification
%   [random_forest] = mb_train_random_forest_class(varargin)
%
% MB_RANDOM_FOREST_CLASS_TRAIN uses the U_PACKARGS interface function
% to allow optional arguments to be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%   sampling_method: Function name used to extract each set of training
%       data from the global population. It is assumed this produces an N x D
%       matrix of training data (where each row is a datapoint consisting of D
%       input variables) and an N x 1 vector of class labels for each data
%       point
%
%   sampling_method_args: structure of arguments used by the sampling
%       method
%
%   d:  The size of the random subset of variables to test at each
%           node. Note if d = 0, the algorithm becomes standard
%           bagging
%
%   tree_dir: path to the directory in which trees should be saved to
%
% Optional Arguments:
%   n_trees: Number of trees in the forest to build (default 100)
%
%   prior:  Prior probabilities for each class in the data
%           (default proportional representation from input)
%
%   cost:   The relative cost of misclassications (default equal
%           cost across all classes)
%
%   split_criterion: either {'gdi'}, 'twoing' or 'deviance'. The function
%           used to determine the optimal split for each variable
%
%   split_min: Minimum number of datapoints at a node, below which no
%           further splitting is performed (default 1 - i.e. build full trees)
%
%   prune:  Turn tree pruning on/{off}
%
%   names:  Add text labels for each class
%
%
% Outputs:
%   random_forest: structure containing the following fields
%       trees: Nx1 cell, each containing a single classification tree
%       D = number of input variables;
%       nclasses = number of classes in output labels;
%       classname = value/name for each class (as opposed to ordinal
%       index);
%       prior = copy of input prior
%       cost = copy of input cost;
%       names = copy of input names;
%
% Example: [random_forest] = mb_random_forest_class_train('X', X, 'y', y);
%
% Notes: This function implements a version of random forests in which each treat
% is trained from data randomly sampled (without replacement) from a population of
% data stored on the hard-disk. This adds a computational overhead before
% computing each tree, but allows arbitrarily large populations of data.
%
% Because it is assumed the input datasets will be large, trees in the
% forest are not stored in memory but written to disk, with the filepath
% saved in the forest structure
%
% For a more standard implementation of random forests (as described by
% Breiman) see MB_RANDOM_FOREST_CLASS_TRAIN_BOOT. This version also
% includes the option of calculating out-of-bag error rates and creating a
% proximity matrix
%
% The code currently works for classification trees,
% but should also be implemented for regression trees.
%
% The function only works for continuous variables however functionality
% for categorical should also be added.
%
% For more information on the random forests algorithm see:
%   Random Forests, Leo Breiman, Machine Learning, 45, 5-32, 2001
%
% See also: MB_TREE_CLASS_TRAIN MB_RANDOM_FOREST_CLASS_PREDICT
%
% Created: 13-Oct-2009
% Author: Michael Berks
% Email : michael.berks@postgrad.man.ac.uk
% Phone : +44 (0)161 275 1241
% Copyright: (C) University of Manchester

% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    {'sampling_method',... % the mandatory arguments
    'sampling_method_args',...
    'd',...
    'tree_dir'}, ...
    'X_test', [],...
    'y_test', [],...
    'do_err', 0,...
    'tree_root', [],...
    'n_trees', 100,...
    'prior', [],... % the optional arguments
    'cost', [],...
    'split_criterion', 'gdi',...
    'split_min', 1, ...
    'prune', 0,...
    'minimise_size', 1,...
    'names', [],...
    'quiet', 0,...
    'save_path', []);

%Get input arguments for whole forest from fields of args
sampling_method = args.sampling_method;
sampling_method_args = args.sampling_method_args;
X_test = args.X_test;
y_test = args.y_test;
do_err = args.do_err;
minimise_size = args.minimise_size;
n_trees = args.n_trees;
tree_dir = [args.tree_root args.tree_dir];

if ~isempty(tree_dir)
    if ~strcmp(tree_dir(end), '/') && ~strcmp(tree_dir(end), '\');
        tree_dir = [tree_dir filesep];
    end
    if ~exist(tree_dir, 'dir')
        mkdir(tree_dir);
    end
end

%Now copy input args to tree_args to pass to individual tree building
tree_args.random_m = args.d;
tree_args.prior = args.prior;
tree_args.cost= args.cost;
tree_args.split_criterion = args.split_criterion;
tree_args.split_min = args.split_min;
tree_args.prune = args.prune;
tree_args.names = args.names;

num_train = sampling_method_args.num_train;
num_test = sampling_method_args.num_test;

%Create random forest structure and pre-allocate space for each tree
random_forest.trees = cell(n_trees, 1);

%If we're testing as we go, pre-allocate for error trace
if ~isempty(X_test) || do_err
    random_forest.error_trace = zeros(n_trees, 1);
    random_forest.tree_error = zeros(n_trees, 1);
end

%Compute individual classification trees for boostrap samples of X
for ii = 1:n_trees
    
    if ~args.quiet
        display(['Building tree ', num2str(ii), ' of ', num2str(n_trees)]);
    end
    
    tic;
    
    % Get training data (and test data if necessary)
    if do_err
        [X y X_test y_test row col image_idx_abnormal image_idx_normal] = feval(sampling_method, sampling_method_args);
        
        abnormal_test_idx = image_idx_abnormal(num_train+1:end);
        normal_test_idx = image_idx_normal(num_train+1:end);
        
        if ~exist('abnormal_roi_abnormalvotes', 'var')
            abnormal_roi_abnormalvotes = zeros(row, col);
            normal_roi_abnormalvotes = zeros(row, col);
            abnormal_roi_totalvotes = zeros(row, col);
            normal_roi_totalvotes = zeros(row, col);
        end
    else
        [X y] = feval(sampling_method, sampling_method_args);
    end
    
    t = toc;
    display(['Time sampling data = ', num2str(t)]);
    
    if ii == 1
        % Get number of dimensions in samples
        D = size(X,2);
        
        %Convert input classes to class labels as returned by tree classifiers
        [y_idx, y_vals] = mb_grp2idx(y);
        y = y_vals(y_idx);
        nclasses = max(y_idx);
        
        if ~isempty(X_test)
            y_pred = zeros(size(X_test,1), nclasses);
            y_test_idx = mb_grp2idx(y_test);
            y_test = y_vals(y_test_idx);
        end
        
        %Could check here that D and nclasses match for each dataset - but for
        % assume they always do
    end
    
    tic;
    %build a tree for the sampled data set
    tree = mb_tree_class_train(X, y, tree_args); clear X y;
    t = toc;
    display(['Time building tree = ', num2str(t)]);
    
    %If we've been given test data, compute test error rate given the new
    %tree
    if ~isempty(X_test)
        
        %Get predictions from new tree
        [y_tree] = mb_tree_predict(tree, X_test);
        y_double = str2double(y_tree);
        
        temp_abnormaltotalvotes = zeros(row, col);
        temp_abnormaltotalvotes(abnormal_test_idx) = 1;
        abnormal_roi_totalvotes = abnormal_roi_totalvotes + temp_abnormaltotalvotes;
        
        temp_normaltotalvotes = zeros(row, col);
        temp_normaltotalvotes(normal_test_idx) = 1;
        normal_roi_totalvotes = normal_roi_totalvotes + temp_normaltotalvotes;
      
        abnormal_test_idx(y_double(1:num_test)==0) = [];
        temp_abnormalvotes = zeros(row, col);
        temp_abnormalvotes(abnormal_test_idx) = 1;
        abnormal_roi_abnormalvotes = abnormal_roi_abnormalvotes + temp_abnormalvotes;
        
        normal_test_idx(y_double(num_test+1:end)==0) = [];
        temp_normalvotes = zeros(row, col);
        temp_normalvotes(normal_test_idx) = 1;
        normal_roi_abnormalvotes = normal_roi_abnormalvotes + temp_normalvotes;
        
        
        %For each class increment the class votes
        for jj = 1:nclasses
            y_pred(:, jj) = y_pred(:, jj) + ...
                strcmp(y_tree, y_vals(jj));
        end
        
        %Get index to class with maximum votes
        [dummy idx_y_pred] = max(y_pred, [], 2); clear dummy
        
        %Convert class index to class labels and compute error rate
        random_forest.error_trace(ii) = 1 - mean(strcmp(y_vals(idx_y_pred), y_test));
        random_forest.tree_error(ii) = 1 - mean(strcmp(y_tree, y_test));
        
    end
    
    %To save space, throwaway fields from the tree not required for each
    %tree (e.g. some fields such as nclasses need only be saved once per
    %forest, not in every tree)
    tree = rmfield(tree, {'prior', 'cost', 'names', 'nclasses'});
    if minimise_size
        tree = rmfield(tree,...
            {'classcount', 'classprob','nodeprob', 'nodeerr', 'risk', 'nodesize'}); %#ok
    end
    
    %Set filename at which to save tree
    tree_name = ['rf_tree', zerostr(ii, 4), '.mat'];
    save([tree_dir tree_name], 'tree'); clear tree;
    
    % Save the filepath to the tree
    random_forest.trees{ii} = tree_name;
    
end

%Save fields to random_forest output structure
random_forest.D = D;
random_forest.d = args.d;
random_forest.tree_dir = args.tree_dir;
random_forest.tree_root = args.tree_root;
random_forest.nclasses = nclasses;
random_forest.classname = y_vals;
random_forest.prior = tree_args.prior;
random_forest.cost = tree_args.cost;
random_forest.names = tree_args.names;

if ~isempty(args.save_path)
    save(args.save_path, 'random_forest');
end

