function [random_forest] = mb_random_forest_reg_train_boot(varargin)
%MB_RANDOM_FOREST_REG_TRAIN_BOOTSTRAP given input data and a set of class labels,
%builds a random forest for regression using the original bootstrap method
%suggested by Breiman
%   [random_forest] = mb_random_forest_class_train_boot(varargin)
%
% MB_RANDOM_FOREST_CLASS_TRAIN_BOOT uses the U_PACKARGS interface function
% to allow optional arguments to be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%   X:      N x D matrix of input data, where each row is a
%           datapoint consisting of D input variables
%   y:      N x 1 vector of regression values for each data point
%
% Optional Arguments:
%   n_trees: Number of trees in the forest to build (default 100)
%
%   d:  The size of the random subset of variables to test at each
%           node. Note if d = 0, the algorithm becomes standard bagging
%           (default ceil(sqrt(D)) )
%
%   do_proximities: Flag whether or not to compute proximities between all
%           data points (note this creates an NxN matrix in the output so
%           may not be a good idea if working memory is tight - default
%           off)
%
%   do_oob_trace: Flag whether or not to compute oob errors as each new
%           tree is added (default off);
%
%   split_criterion: ...
%
%   split_min: Minimum number of datapoints at a node, below which no
%           further splitting is performed (default 1 - i.e. build full trees)
%
%   prune:  Turn tree pruning on/{off}
%
%
% Outputs:
%   random_forest: structure containing the following fields
%       trees: Nx1 cell, each containing a single classification tree
%       D = number of input variables;
%       oob_error = out-of-bag error estimate
%       oob_votes = out-of-bag class votes for each point N;
%       oob_trace = if set, the trace of oob_errors as each new tree is
%           added to the forest
%
% Example: [random_forest] = mb_random_forest_reg_train('X', X, 'y', y);
%
% Notes: This function is set up to train a random forest based on
% the basic Breiman method. That is, for each tree a bootsrap set of equal
% size is created and used to train the tree. At each node in the tree a
% random subset of input variables (of pre-selected size) are tested to
% find the best split. However it is intended that additional options for 
% tree building (e.g. taking random linear combinations of variables as features)
% will be added
%
% For a method in which subsets are selected without replacement from a
% hard disk repository of data (thus allowing arbitrarily large sets data)
% see MB_RANDOM_FOREST_REG_TRAIN
%
% The function only works for continuous variables however functionality
% for categorical should also be added.
%
%
% For more information on the random forests algorithm see:
%   Random Forests, Leo Breiman, Machine Learning, 45, 5-32, 2001
%
% See also: MB_TREE_REG_TRAIN MB_RANDOM_FOREST_REG_PREDICT_BOOT
%
% Created: 23-Feb-2011
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
			 '0', ... % non-strict mode
			 {'X',... % the mandatory arguments
             'y'}, ...
             'n_trees', 100,...
             'd', [],...
             'do_oob', 1,...
			 'do_oob_trace', 0,...
             'do_proximities', 0,...
             'oob_proximities', 0,...
             'split_min', 1, ...
             'prune', 0,...
             'minimise_size', 0,...
             'names', [],...
             'quiet', 0);

%Store input data in variable X and remove X from fields of args
X = args.X; args = rmfield(args, 'X');
y = args.y; args = rmfield(args, 'y');
minimise_size = args.minimise_size; args = rmfield(args, 'minimise_size');
n_trees = args.n_trees; args = rmfield(args, 'n_trees');
do_oob = args.do_oob; args = rmfield(args, 'do_oob');
do_oob_trace = args.do_oob_trace; args = rmfield(args, 'do_oob_trace');
do_proximities = args.do_proximities; args = rmfield(args, 'do_proximities');
oob_proximities = args.oob_proximities; args = rmfield(args, 'oob_proximities');
quiet = args.quiet; args = rmfield(args, 'quiet');

%workout number of data points and input variables
[N D] = size(X);

%Now copy input args to tree_args to pass to individual tree building
tree_args = args; clear args;

%work out number of variables to select randomly at each node in tree
if isempty(tree_args.d)
    tree_args.random_m = ceil(sqrt(D));
    %args.random_m = ceil(sqrt(D));
else
    tree_args.random_m = tree_args.d;
    %args.random_m = args.d;
end
random_forest.d = tree_args.d;
tree_args = rmfield(tree_args, 'd');

%Create random forest structure and pre-allocate space for each tree
random_forest.trees = cell(n_trees, 1);
    
%pre-allocate space for out-of-bag (oob) error estimates
%oob_counts = zeros(N, 1);
y_oob = [];
oob_trace = [];
if do_oob
    oob_sum = zeros(N, 1);
    oob_counts = zeros(N, 1);
    if do_oob_trace
        oob_trace = zeros(n_trees, 1);
    end
end

%if we've been asked to compute proximities pre-allocate now
random_forest.proximities = [];
if do_proximities
    random_forest.proximities = zeros(N);
    
    %if we're only using oob proximities we need to count co-occurrences of
    %pts being oob, pre-allocate space for this now
    if oob_proximities
        oob_prox_counts = zeros(N);
    end
end
%Compute individual classification trees for boostrap samples of X
for ii = 1:n_trees
    
    if ~quiet
        display(['Building tree ', num2str(ii), ' of ', num2str(n_trees)]);
    end
    
    %compute bootstrap idx (sample integers from 1 to N with replacement)
    boot_idx = ceil(N*rand(N,1));
    
    %build a tree for the bootstrap set
    tree = mb_tree_reg_train(X(boot_idx,:), y(boot_idx,:), tree_args);
    
    %To save space, throwaway fields from the tree not required for each
    %tree - generally this is a bad idea though - only use in ermengency!
    if minimise_size
        tree = rmfield(tree,...
            {'nodeprob', 'nodeerr', 'risk', 'nodesize'});
    end
    
    %May want to manipulate the fields of each tree into a better
    %structure than below, but for now save in n_trees x 1 cell
    random_forest.trees{ii} = tree;
    
    if do_oob
        %Do out-of-bag error stuff:
        
        %First work out index to out-of-bag samples (i.e. those not in the bootstrap
        %set)
        oob_idx = setdiff((1:N)', boot_idx);
    
        %Get predictions for OOB samples
        %oob_y_tree = mb_tree_predict(tree, X(oob_idx,:));
        [y_tree nodes] = mb_tree_predict(tree, X);
    
        %Update sum of oob predictions
        oob_sum(oob_idx) = oob_sum(oob_idx) + y_tree(oob_idx);

        %Update count of oob predictions
        oob_counts(oob_idx) = oob_counts(oob_idx) + 1;
        
        %If asked for, compute current oob error
        if do_oob_trace
            %Work out which data points we have OOB estimates for
            oob_exists = oob_counts > 0;
            
            %Get y predictions for the points that have existing oob
            %predictions
            y_oob = oob_sum(oob_exists) ./ oob_counts(oob_exists);
            
            %Get RMS error of current y predictions and store in OOB trace
            oob_trace(ii) = sqrt(sum((y_oob - y(oob_exists)).^2)) / sum(oob_exists);
        end
    end
    
    %Compute proximity measure - i.e. count the number of times data points
    %end up at the same leaf node
    if do_proximities
        leaf_nodes = unique(nodes);
             
        for jj = 1:length(leaf_nodes)
            leaf = leaf_nodes(jj);
            leaf_pts = double(nodes == leaf);
            %Check if we're only calculating proximities for out-of-bag
            %data pts
            if oob_proximities
                %if so, set the bootstrap selected pts to zero
                leaf_pts(boot_idx) = 0;
            end
            random_forest.proximities = random_forest.proximities + ...
                leaf_pts * leaf_pts';
        end
        if oob_proximities
            %count co-occurences of oob samples
            prox_pts = zeros(N,1);
            prox_pts(oob_idx) = 1;
            oob_prox_counts = oob_prox_counts + prox_pts * prox_pts';
        end
    end
    
    %Do variable importance stuff - not yet implemented   
end

if do_proximities
    %Divide proximities by n_trees to get measure between 0 and 1
    if oob_proximities
        random_forest.proximities = random_forest.proximities ./ oob_prox_counts;
    else
        random_forest.proximities = random_forest.proximities / n_trees;
    end
    
    %we could subtract from 1 to get a distance measure, but to be consistent with
    %breiman code we don't for now
    %random_forest.proximities = 1 - random_forest.proximities;
end

oob_error = [];
if do_oob
    %Compute final OOB error
   oob_exists = oob_counts > 0;
            
    %Get final y predictions for the points that have existing oob
    %predictions
    y_oob = oob_sum(oob_exists) ./ oob_counts(oob_exists);

    %Get RMS error of final y predictions and store in OOB error
    oob_error = sqrt(sum((y_oob - y(oob_exists)).^2)) / sum(oob_exists);
end

%Save fields to random_forest output structure
random_forest.D = D;
random_forest.oob_error = oob_error;
random_forest.y_oob = y_oob;
random_forest.oob_trace = oob_trace;

