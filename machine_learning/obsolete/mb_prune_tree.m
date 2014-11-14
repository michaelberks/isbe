function [tree] = mb_prune_tree(varargin)
%MB_PRUNE_TREE prunes a classification/regression tree as specified
%   [] = mb_prune_tree(varargin)
%
% MB_PRUNE_TREE uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%   tree:   The tree to be pruned
%
% Optional Arguments:
%   nodes:  The list of nodes to prune (default [])
%
%   level: the level below which is pruned (default [])
%
% If neither 'level' or 'nodes' are set, the tree is 'optimally' pruned
%
% Outputs:
%
%   tree:   the pruned tree
%
% Example:
%
% Notes:
%   T2 = MB_PRUNE_TREE(T1,'level',LEVEL) takes a decision tree T1 and a pruning
%   level LEVEL, and returns the decision tree T2 pruned to that level.
%   The value LEVEL=0 means no pruning.  Trees are pruned based on an
%   optimal pruning scheme that first prunes branches giving less
%   improvement in error cost.
%
%   T2 = MB_PRUNE_TREE(T1,'nodes',NODES) prunes the nodes listed in the NODES
%   vector from the tree.  Any T1 branch nodes listed in NODES become
%   leaf nodes in T2, unless their parent nodes are also pruned.  The
%   VIEW method can display the node numbers for any node you select.
%
%   T2 = MB_PRUNE_TREE(T1) returns the decision tree T2 that is the full,
%   unpruned T1, but with optimal pruning information added.  This is
%   useful only if you created T1 by pruning another tree, or by using
%   the CLASSREGTREE function with pruning set 'off'.  If you plan to prune
%   a tree multiple times along the optimal pruning sequence, it is more
%   efficient to create the optimal pruning sequence first.
%
% See also: PRUNE
%
% Created: 13-Oct-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
			 '0', ... % non-strict mode
			 {'tree'},... % the mandatory arguments
             'nodes', [],... % the optional arguments
             'level', []);

% Process inputs
tree = args.tree;
level = args.level;
nodes = args.nodes;
   
% Now do the pruning
if isempty(tree.alpha) && isempty(nodes)
    tree = getpruneinfo(tree);                % may need optimal prune sequence
end
pruned = false;
if ~isempty(level)
    tree = subtree(tree,level);               % remove stuff using optimal sequence
    pruned = true;
elseif ~isempty(nodes)
    tree = prunenodes(tree,nodes);            % remove children of specified nodes
    pruned = true;
end
if ~isempty(tree.prunelist) && pruned
    tree.prunelist = [];
    tree = getpruneinfo(tree);                % recompute prune info from scratch
end

% ------------------------------------------------------------------
function tree = getpruneinfo(tree)
%GETPRUNEINFO Get optimal pruning information and store into decision tree.

% Start from the smallest tree that minimizes R(alpha,T) for alpha=0
N = length(tree.node);
parent     = tree.parent;
children   = tree.children;

isleaf = tree.var(:)==0;
nleaves = sum(isleaf);
adjfactor = 1 + 100*eps;

% Work up from the bottom of the tree to compute, for each branch node,
% the number of leaves under it and sum of their costs
treatasleaf = isleaf';
nodecost = tree.risk;
costsum = nodecost;
nodecount = double(isleaf);
while(true)
    % Find "twigs" which I define as branches with two leaf children
    branches = find(~treatasleaf);
    twig = branches(sum(treatasleaf(children(branches,:)),2) == 2);
    if isempty(twig), break; end;    % worked our way up to the root node

    % Add the costs and sizes of the two children, give to twig
    kids = children(twig,:);
    costsum(twig)   = sum(costsum(kids'),1)';
    nodecount(twig) = sum(nodecount(kids'),1)';
    treatasleaf(twig) = 1;
end

% Now start pruning to generate a sequence of smaller trees
whenpruned = zeros(N,1);
branches = find(~isleaf);
prunestep = 0;
allalpha = zeros(N,1);
ntermnodes = zeros(N,1);
ntermnodes(1) = nleaves;
while(~isempty(branches))
   prunestep = prunestep + 1;
   
   % Compute complexity parameter -- best tree minimizes cost+alpha*treesize
   alpha = max(0, nodecost(branches) - costsum(branches)) ./ ...
           max(eps,nodecount(branches) - 1);
   bestalpha = min(alpha);
   toprune = branches(alpha <= bestalpha*adjfactor);

   % Mark nodes below here as no longer on tree
   wasleaf = isleaf;
   kids = toprune;
   while ~isempty(kids)
      kids = children(kids,:);
      kids = kids(kids>0);
      kids(isleaf(kids)) = [];
      isleaf(kids) = 1;
   end
   newleaves = toprune(~isleaf(toprune));
   isleaf(toprune) = 1;

   % Remember when branch was pruned, also perhaps leaves under it
   whenpruned(isleaf~=wasleaf & whenpruned==0) = prunestep;
   whenpruned(toprune) = prunestep;   % this branch was pruned

   % Update costs and node counts
   for j=1:length(newleaves)          % loop over branches that are now leaves
      node = newleaves(j);
      diffcost  = nodecost(node) - costsum(node);
      diffcount = nodecount(node) - 1;
      while(node~=0)                  % work from leaf up to root
         nodecount(node) = nodecount(node) - diffcount;
         costsum(node)   = costsum(node) + diffcost;
         node = parent(node);         % move to parent node
      end
   end

   allalpha(prunestep+1) = bestalpha;
   ntermnodes(prunestep+1) = nodecount(1);
   
   % Get list of branches on newly pruned tree
   branches = find(~isleaf);
end

tree.prunelist  = whenpruned;
tree.alpha      = allalpha(1:prunestep+1);
tree.ntermnodes = ntermnodes(1:prunestep+1);

% ------------------------------------------------------------
function tree = subtree(tree,p)
%SUBTREE Get subtree from tree indexed by pruning point.

whenpruned = tree.prunelist;
v = find(whenpruned>0 & whenpruned<=p);
if ~isempty(v)
   tree = prunenodes(tree,v);
end

% ------------------------------------------------------------
function tree = prunenodes(tree,branches)
%PRUNENODES Prune selected branch nodes from tree.

N = length(tree.node);

% Find children of these branches and remove them
parents = branches;
tokeep = true(N,1);
kids = [];
while(true)
   newkids = tree.children(parents,:);
   newkids = newkids(:);
   newkids = newkids(newkids>0 & ~ismember(newkids,kids));
   if isempty(newkids), break; end
   kids = [kids; newkids];
   tokeep(newkids) = 0;
   parents = newkids;
end

% Convert branches to leaves by removing split rule and children
tree.var(branches) = 0;
tree.cut(branches) = 0;
tree.children(branches,:) = 0;

% Get new node numbers from old node numbers
ntokeep = sum(tokeep);
nodenums = zeros(N,1);
nodenums(tokeep) = (1:ntokeep)';

% Reduce tree, update node numbers, update child/parent numbers
%reduce any fields with the same number of rows as tree.node - doing this
%automatically means we don't have to keep changing this function if we add
%in optional fields during tree construction
tree_fields = fieldnames(tree);
for ii = 1:length(tree_fields)
    if ~strcmp(tree_fields{ii}, 'node')
        tf = tree.(tree_fields{ii});
        if size(tf,1) == N
            tree.(tree_fields{ii}) = tf(tokeep,:);
        end
    end
end

%Now update node numbers in tree.node
tree.node = (1:ntokeep)';

%Finally update parent and children to match the new node numbers
mask = tree.parent>0;
tree.parent(mask) = nodenums(tree.parent(mask));
mask = tree.children>0;
tree.children(mask) = nodenums(tree.children(mask));
