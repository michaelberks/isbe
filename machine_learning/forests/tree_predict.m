function [y_fit, nodes] = tree_predict(tree, X, return_probs)
%TREE_PREDICT Non-recursive tree prediction
%   [y_fit] = pt_tree_predict(tree,X)
%
% Inputs:
%	tree:   Classification or regression tree, as created by
%	CLASS_TREE_TRAIN or REG_TREE_TRAIN
%
%	X:      N x d matrix of input data, where each row is a
%           datapoint consisting of d input variables
%
% Outputs:
%	y_fit: predicted class labels for each data point
%
%   nodes: the destination leaf node in the tree for each data point
%
%
% Example:
%
% Notes: Adapted from the EVAL function for the matlab CLASSREGTREE object.
% Note categorical variables are not yet supported
%
% See also: TREE_CLASS_TRAIN TREE_REG_TRAIN CLASSREGTREE
%
% Created: 15-Oct-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester
if nargin < 3
    return_probs = false;
end

missing_rows = any(isnan(X), 2);
if any(missing_rows)
    warning('TreePredict:missing_data', ['Removing ' num2str(sum(missing_rows))...
        ' rows containing NaNs from data']);
    X(missing_rows,:) = [];
end

[N, nvars] = size(X);
if nvars ~= tree.npred
   error('The number of columns in X must match the number of variables in the tree');
end

if ~isempty(tree.prunelist)
   prunelist = tree.prunelist;
else
   prunelist = repmat(Inf,size(tree.node));
end

% initialise queue: to begin, assign all probes to the root node
node_queue(1,:) = { 1:N, 1 };
nodes = zeros(N,1);

while ~isempty(node_queue)
	% get current node number and probes assigned to it
	thisnode = node_queue{end,2};
	rows = node_queue{end,1};

	% pop from the queue
	% we will either not need it or replace it with its two children
	node_queue(end,:) = [];
	
	% get parameters for this node
	splitvar      = tree.var(thisnode);
	cutoff        = tree.cut(thisnode);
	kids          = tree.children(thisnode,:);
	%catsplit      = tree.catsplit; Ignore categorical variables for now
	prunelevel    = prunelist(thisnode);

	if splitvar==0
		% terminal node
		nodes(rows) = thisnode;
		continue;
	else
		% non-terminal
		
		% Determine if this point goes left, goes right, or stays here
		if splitvar > 0                % continuous variable
			x = X(rows,splitvar);   
			ispresent = ~isnan(x);
			isleft = (x(ispresent) < cutoff);
			isright = ~isleft;
% 		else                         % categorical variable
% 			isleft = ismember(x,catsplit{cutoff,1});
% 			isright = ismember(x,catsplit{cutoff,2});
% 			ismissing = ~(isleft | isright);
		end

		subrows = rows(ispresent);  

		% left child node
		if any(isleft)
			node_queue(end+1,:) = { subrows(isleft) kids(1) };
		end

		% right child node
		if any(isright)
			node_queue(end+1,:) = { subrows(isright) kids(2) };
		end

		% some missing, treat as leaf
		if ~all(ispresent)
			node_queue(rows(~ispresent)) = thisnode;
		end	
		
	end % if splitvar==0
end % while

%If a classification tree assign class labels to class indices
if isequal(tree.method, 'classification')
    if return_probs
        y_fit = tree.classprob(nodes,:);
    else
        %Get class indices/regression values from tree for this set of nodes
        id = tree.class(nodes);
        y_fit = tree.classname(id);
    end
else
   y_fit = tree.class(nodes);
end

if any(missing_rows)
    temp = y_fit;
    y_fit = nan(length(missing_rows), size(y_fit,2));
    y_fit(~missing_rows,:) = temp;
end