function [predictor] = mb_random_forest_reg_train(varargin)
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
%   tree_dir: path to the directory in which trees should be saved to
%
% Optional Arguments:
%   d:  The size of the random subset of variables to test at each
%           node. Note if d = 0, the algorithm becomes standard
%           bagging
%
%   n_trees: Number of trees in the forest to build (default 100)
%
%   prior:  Prior probabilities for each class in the data
%           (default proportional representation from input)
%
%   cost:   The relative cost of misclassications (default equal
%           cost across all classes)
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
%       n_classes = number of classes in output labels;
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
    ... % Mandatory arguments
   {'sampling_args',... % for sampling training data
    'decomposition_args',... % for sampling training features
    'prediction_type',...   
    'tree_dir'}, ...
    ... % Optional arguments
    'd', [], ...
    'tree_root', [],...
    'n_trees', 100,...
    'split_criterion', 'dabs',...
    'var_criterion', 'dabs',...
    'split_min', 100, ...
    'end_cut_min', 25, ...
    'do_ubound', 1,...
    'do_circular', [],...
    'w_prior', 0, ...
    'impure_thresh', 1e4,...
    'prune', 0,...
    'mod', 0, ...
    'random_m', [], ...
    'minimise_size', 1,...
    'names', [],...
    'overwrite', 0,...
    'predictor_path', []);

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

%Compute individual classification trees for boostrap samples of X
tb = timebar('limit', n_trees);
for ii = start_tree:n_trees
	timebar(tb,'title', sprintf('Sampling tree #%02i',ii));
    
    %if args.quiet
    fprintf('Building tree %d of %d\n', ii, n_trees);
    tic;
    % Get training data
    sampling_args.task_id = ii;
    [X y] = feval(sampling_method, args);
    t = toc;
    display(['Time sampling data = ', num2str(t)]);

    % Get number of dimensions in samples
    D = size(X,2);
    if isempty(tree_args.random_m)
        tree_args.random_m = round(sqrt(D));
    elseif strcmpi(tree_args.random_m, 'f')
        tree_args.random_m = [round(sqrt(D)) 1];
    end

    fprintf('Building tree...');
    tic;
    %build a tree for the sampled data set
    tree = mb_tree_reg_train(X, y, tree_args); 
    clear X y; %#ok
    fprintf('done (took %0.2f sec)\n', toc);
    
    %To save space, throwaway fields from the tree not required for each
    %tree (e.g. some fields such as n_classes need only be saved once per
    %forest, not in every tree)
    if minimise_size
        %tree = rmfield(tree,...
        %    {'nodeprob', 'nodeerr', 'risk', 'nodesize'}); %#ok
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

%Save fields to random_forest output structure
predictor.D = D;
predictor.d = tree_args.random_m;
predictor.tree_dir = tree_args.tree_dir;
predictor.tree_root = tree_args.tree_root;
predictor.regression_method = tree_args.prediction_type;

% We probably don't need to do this here - it can go in the main
% build_predictor function/script.
if ~isempty(args.predictor_path)
    save(args.predictor_path, 'predictor');
end

