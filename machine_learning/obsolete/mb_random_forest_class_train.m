function [predictor] = mb_random_forest_class_train(varargin)
%MB_RANDOM_FOREST_CLASS_TRAIN given input data and a set of class labels,
%builds a random forest for classification
%   [predictor] = mb_train_random_forest_class(varargin)
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
%   predictor: structure containing the following fields
%       trees: Nx1 cell, each containing a single classification tree
%       D = number of input variables;
%       nclasses = number of classes in output labels;
%       classname = value/name for each class (as opposed to ordinal
%       index);
%       prior = copy of input prior
%       cost = copy of input cost;
%       names = copy of input names;
%
% Example: [predictor] = mb_random_forest_class_train('X', X, 'y', y);
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
   {'sampling_args',... % for sampling training data
    'decomposition_args',... % for sampling training features
    'prediction_type',...   
    'tree_dir'}, ...
    'd', [],...
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
    'minimise_size', 0,...
    'random_m', [], ...
    'names', [],...
    'quiet', 0,...
    'overwrite', 0,...
    'predictor_path', []);
clear varargin;

%Get input arguments for whole forest from fields of args
sampling_args = args.sampling_args;
sampling_method = sampling_args.sampling_method;
minimise_size = args.minimise_size;
n_trees = args.n_trees;

create_folder(args.tree_dir);

%Create random forest structure and pre-allocate space for each tree
predictor.trees = cell(n_trees, 1);

if args.overwrite
    start_tree = 1;
else
    %see if there any trees already in the tree dir (e.g. from an aborted
    %job)
    old_trees = length(dir([args.tree_dir '/rf_tree*.mat']));
    for ii = 1:old_trees
        predictor.trees{ii} = sprintf('rf_tree%04d.mat', ii);
    end
    start_tree = old_trees + 1;
end

%Now copy input args to tree_args to pass to individual tree building
tree_args = rmfield(args, {'sampling_args', 'decomposition_args'});

X_test = args.X_test;
y_test = args.y_test;
do_err = args.do_err;

%If we're testing as we go, pre-allocate for error trace
if ~isempty(X_test) || do_err
    predictor.error_trace = zeros(n_trees, 1);
    predictor.tree_error = zeros(n_trees, 1);
end

%Compute individual classification trees for boostrap samples of X
tb = timebar('limit', n_trees);
for ii = start_tree:n_trees
    if ~args.quiet
        display(['Building tree ', num2str(ii), ' of ', num2str(n_trees)]);
    end
    % Get training data (and test data if necessary)
    tic;
    sampling_method_args.task_id = ii;
    if do_err
        [X y X_test y_test] = feval(sampling_method, args);
    else
        [X y] = feval(sampling_method, args);
    end
    display(['Time sampling data = ', num2str(toc)]);
    
    if ii == start_tree
        % Get number of dimensions in samples
        D = size(X,2);

        %Convert input classes to class labels as returned by tree classifiers
        [y_idx, y_vals] = mb_grp2idx(y);
        %y = y_vals(y_idx);
        nclasses = max(y_idx);
        
        if ~isempty(X_test)
            y_pred = zeros(size(X_test,1), nclasses);
            %y_test_idx = mb_grp2idx(y_test);
            %y_test = y_vals(y_test_idx);
        end

        if isempty(tree_args.random_m)
            tree_args.random_m = round(sqrt(D));
            args.d = round(sqrt(D));
        end
        %Could check here that D and nclasses match for each dataset - but for
        % assume they always do
    end
    
    fprintf('Building tree...');
    tic;
    % Build a tree for the sampled data set
    tree = mb_tree_class_train(X, y, tree_args); clear X y;
    clear X y; %#ok
    fprintf('done (took %0.2f sec)\n', toc);

    %If we've been given test data, compute test error rate given the new
    %tree
    if ~isempty(X_test)
        %Get predictions from new tree
        [y_tree] = mb_tree_predict(tree, X_test);
       
        %For each class increment the class votes
        for jj = 1:nclasses
            y_pred(:, jj) = y_pred(:, jj) + strcmp(y_tree, y_vals(jj));
        end
        
        %Get index to class with maximum votes
        [dummy idx_y_pred] = max(y_pred, [], 2); clear dummy

        %Convert class index to class labels and compute error rate
        predictor.error_trace(ii) = 1 - mean(strcmp(y_vals(idx_y_pred), y_test));
        predictor.tree_error(ii) = 1 - mean(strcmp(y_tree, y_test));
    end
    
    %To save space, throwaway fields from the tree not required for each
    %tree (e.g. some fields such as nclasses need only be saved once per
    %forest, not in every tree)
    tree = rmfield(tree, {'prior', 'cost', 'names', 'nclasses'});
    if minimise_size
        tree = rmfield(tree, {'classcount', 'classprob','nodeprob', ...
                              'nodeerr', 'risk', 'nodesize'} ); %#ok
    end
    
	% Set filename (e.g. 'tree0003.mat') at which to save tree
    tree_name = sprintf('rf_tree%04d.mat',ii);
    save([tree_args.tree_dir tree_name], 'tree'); 
    clear tree;
    
    % Save the filepath to the tree
    predictor.trees{ii} = tree_name;

    timebar(tb,'advance');
end
timebar(tb,'close'); clear tb;

%Save fields to predictor output structure
predictor.D = D;
predictor.d = args.d;
predictor.tree_dir = args.tree_dir;
predictor.tree_root = args.tree_root;
predictor.nclasses = nclasses;
predictor.classname = y_vals;
predictor.prior = tree_args.prior;
predictor.cost = tree_args.cost;
predictor.names = tree_args.names;

if ~isempty(args.predictor_path)
    save(args.predictor_path, 'predictor');
end

