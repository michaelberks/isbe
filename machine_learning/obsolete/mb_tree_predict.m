function [y_fit, nodes] = mb_tree_predict(tree, X, return_probs)
%MB_TREE_PREDICT *Insert a one line summary here*
%   [y_fit] = mb_tree_predict(tree,X)
%
% Inputs:
%	tree:   Classification or regression tree, as created by
%	MB_CLASS_TREE_TRAIN or (training regression trees not yet implemented)
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
% See also: MB_CLASS_TREE_TREE CLASSREGTREE
%
% Created: 15-Oct-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester
if nargin < 3
    return_probs = false;
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

%Use recursive function to work out which leaf node each row of X is
%assigned to
nodes = doapply(tree, X, 1:N, 1, zeros(N, 1), 0, prunelist, 1, 1);

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


%------------------------------------------------
function nodes = doapply(tree,X,rows,thisnode,nodes,subtrees,prunelist,endcol,count)
%DOAPPLY Apply classification rule to specified rows starting at a node.
%   This is a recursive function.  Starts at top node, then recurses over
%   child nodes.  THISNODE is the current node at each step.
%
%   NODES has one row per observation and one column per subtree.
%
%   X, NODES, PRUNELIST, and SUBTREES are the same in each recursive call
%   as they were in the top-level call.  ROWS describes the subset of X and
%   NODES to consider.  1:ENDCOL are colums of NODES and the elements of
%   SUBTREES to consider.

splitvar      = tree.var(thisnode);
cutoff        = tree.cut(thisnode);
kids          = tree.children(thisnode,:);
%catsplit      = tree.catsplit; Ignore categorical variables for now
prunelevel    = prunelist(thisnode);

% For how many of the remaining trees is this a terminal node?
if splitvar==0      % all, if it's terminal on the unpruned tree
   ncols = endcol;
else                % some, if it's terminal only after pruning
   ncols = sum(subtrees(1:endcol) >= prunelevel);
end
if ncols>0          % for those trees, assign the node level now
   nodes(rows,(endcol-ncols+1:endcol)) = thisnode;
   endcol = endcol - ncols;
end

% Now deal with non-terminal nodes
if endcol > 0
   % Determine if this point goes left, goes right, or stays here
   x = X(rows,abs(splitvar));   
   if splitvar > 0                % continuous variable
      isleft = (x < cutoff);
      isright = ~isleft;
      ismissing = isnan(x);
%    else                         % categorical variable
%       isleft = ismember(x,catsplit{cutoff,1});
%       isright = ismember(x,catsplit{cutoff,2});
%       ismissing = ~(isleft | isright);
   end

   subrows = rows(isleft & ~ismissing);  % left child node
   if ~isempty(subrows)
      nodes = doapply(tree,X,subrows,kids(1),nodes,subtrees,prunelist,endcol);
   end
   
   subrows = rows(isright & ~ismissing); % right child node
   if ~isempty(subrows)
      nodes = doapply(tree,X,subrows,kids(2),nodes,subtrees,prunelist,endcol);
   end

   subrows = rows(ismissing);            % missing, treat as leaf
   if ~isempty(subrows)
      nodes(subrows,1:endcol) = thisnode;
   end
end