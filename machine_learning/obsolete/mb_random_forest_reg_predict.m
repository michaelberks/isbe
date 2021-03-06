function [y_fit y_trees proximities] = mb_random_forest_reg_predict(forest, X)
%MB_RANDOM_FOREST_CLASS_PREDICT given a random forest and input data,
%predict the set of class labels for the data
%   [random_forest] = mb_random_forest_class_predict(forest, X)
%
%   X:      N x d matrix of input data, where each row is a
%           datapoint consisting of d input variables
%
%   forest: Random forest structure as generated by
%   MB_RANDOM_FOREST_CLASS_TRAIN
%
% Outputs:
%   y_fit:  Set of class labels predicted for each data point (note these
%           are currently returned as text labels in an Nx1 cell string
%           array. To convert to numeric values call str2double(y_fit),
%           which will produce an Nx1 vector of equivalent numeric values
%
%   votes:  Votes for each class for each data point (note by dividing by
%           the number of trees in the forest these can be used to
%           constructed output probabilities (as opposed to hard labels)
%
% Example:
%
% Notes:
%
% See also: MB_RANDOM_FOREST_CLASS_TRAIN MB_TREE_PREDICT
%
% Created: 13-Oct-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

%workout number of data points and input variables
N = size(X, 1);

%Work out number of trees in forest and number of output classes
ntrees = length(forest.trees);
    
%pre-allocate space tree regression values
%y_fit = zeros(N, 1);
y_trees = zeros(N, ntrees);

if nargout > 2
    proximities = zeros(N);
end
for ii = 1:ntrees

    %load tree
    tree = u_load([forest.tree_root forest.tree_dir forest.trees{ii}]);
    
    %Get predictions for data from this tree
    [y_tree nodes] = mb_tree_predict(tree, X);
    y_trees(:,ii) = y_tree;
    
    %Add to running sum;
    %y_fit = y_fit + y_tree;
    %y_fit = mod(y_fit - mb_mod(y_fit - y_tree,180)/ii,180);
    
    if nargout > 2
    %Compute proximity measure - i.e. count the number of times data points
    %end up at the same leaf node
        leaf_nodes = unique(nodes);
        for jj = 1:length(leaf_nodes)
            leaf = leaf_nodes(jj);
            leaf_pts = double(nodes == leaf);
            proximities = proximities + (leaf_pts * leaf_pts');
        end
    end   
end

if isreal(y_trees)
	% if angle predicted then wraparound
	y_trees = pi*y_trees/90; %[0, 2pi]
	x_mean = mean(cos(y_trees),2);
	y_mean = mean(sin(y_trees),2);
	y_fit = mod(90*atan2(y_mean, x_mean)/pi,180);
else
	% otherwise compute mean orientation vector
	y_fit = mean(y_trees,2);
end

%If we need to calculate proximities, divide by ntrees to get measure
%between 0 and 1
if nargout > 2
    proximities = proximities / ntrees;
end