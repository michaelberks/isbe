function [tree] = mb_tree_class_train(X, y, varargin)
%MB_TREE_CLASS_TRAIN given input data and set of class labels, construct a
%single classification tree
%   [tree] = mb_class_tree_train(X, y, varargin)
%
% Inputs:
%   X:      N x d matrix of input data, where each row is a
%           datapoint consisting of d input variables
%   y:      N x 1 vector of class labels for each data points
%
% In addition MB_TREE_CLASS_TRAIN uses the U_PACKARGS interface function
% and so optional arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Optional Arguments:
%  prior:  Prior probabilities for each class in the data
%           (default proportional representation from input)
%
%   cost:   The relative cost of misclassications (default equal
%           cost across all classes)
%
%   random_m: The number of variables to randomly select at each node, from
%           which the optimal split will be chosen
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
% Outputs:
%   tree:   structure containing the decision tree
%
% Example: [tree] = mb_class_tree_train(X, y)
%
% Notes: This function is a implementation adopted from the Matlab function
%   CLASSREGTREE. However MB_TREE_CLASS_TRAIN allows more flexibility in
%   how splits are made at each node. In particular we can choose to select
%   a random subset of input variables from which to calculate the optimal
%   split at each node. It is intended that further functionality (e.g. 
%   choosing linear combinations of variables as features) will be added.
%   This allows "random forests" to built from the trees
%
% See also: MB_TREE_PREDICT CLASSREGTREE MB_RANDOM_FOREST_CLASS_TRAIN
%
% Created: 13-Oct-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
			 '0', ... % non-strict mode
			 'prior', [],... % the optional arguments
             'cost', [],...
             'random_m', 0,...
             'split_criterion', 'gdi',...
             'split_min', 10, ...
             'prune', 0,...
             'names', []);

%Store input data in variable X and remove X from fields of args
%X = args.X; args = rmfield(args, 'X');

%Store input labels in variable y and remove y from fields of args
%y = args.y; args = rmfield(args, 'y');

%Check we have the correct number of labels for data points
if numel(y) ~= size(X,1);
    error('Y (labels) must have the same number of elements as rows in X (data).');
end

%Find data points containing NaNs and discard from data and labels
t = any(isnan(X),2);
if any(t)
   X(t,:) = [];
   y(t) = [];
end

%Compute number of remaining data points and number of input varaiables
[N,nvars] = size(X);

%Get names variables
if ~isempty(args.names)
    names = args.names; args = rmfield(args, 'names');
    if ischar(names)
        names = cellstr(names);
    end
    if ~iscellstr(names) || numel(names)~=nvars
        error('stats:treefit:BadNames',...
              'NAMES must be a character array or cell array with %d strings.',nvars);
    end
else
    names = [];
end

switch(args.split_criterion)
%                Criterion function   Is it an impurity measure?
%                ------------------   --------------------------
    case 'gdi',      critfun = @gdi;      impurity = 1;
    case 'twoing',   critfun = @twoing;   impurity = 0;
    case 'deviance', critfun = @deviance; impurity = 1;
    otherwise,     error('Bad value for ''splitcriterion'' parameter.')
end
   
%Group labels
if islogical(y)
    y = double(y);
end
[y, cnames] = mb_grp2idx(y);   %find groups only after NaNs removed from X

%Check if we have any missing labels (stored as NaNs)
if any(isnan(y))
    t = isnan(y);
    y(t) = [];
    X(t,:) = [];
    N = size(X,1);
end

%Workout number of classes
nclasses = max(y);

%Create binary matrix, C(i,j)==1 means data point i is in class j
C = false(N,nclasses);
C(sub2ind([N nclasses],(1:N)',y)) = 1;

%Work out the numbers of points in each class
Nj = sum(C,1);

% Tree structure fields:
%  .method     method
%  .node       node number
%  .parent     parent node number
%  .class      class assignment for points in this node if treated as a leaf
%  .var        column j of X matrix to be split, or 0 for a leaf node,
%              or -j to treat column j as categorical
%  .cut        cutoff value for split (Xj<cutoff goes to left child node),
%              or index into catsplit if var is negative
%  .children   matrix of child nodes (2 cols, 1st is left child)
%  .nodeprob   probability p(t) for this node
%  .nodeerr    resubstitution error estimate r(t) for this node
%  .risk       R(t) = p(t)*r(t)
%  .nodesize   number of points at this node
%  .prunelist  list of indices that define pruned subtrees.  One entry per
%              node.  If prunelist(j)=k then, at the kth level of pruning,
%              the jth node becomes a leaf (or drops off the tree if its
%              parent also gets pruned).
%  .alpha      vector of complexity parameters for each pruning cut
%  .ntermnodes vector of terminal node counts for each pruning cut
%  .catsplit   call array for categorical splits,
%              left categories in column 1 and right categories in column 2
%  .classprob  vector of class probabilities
%  .classname  names of each class
%  .classcount count of members of each class
%  .nclasses   number of classes
%  .cost       misclassification cost

nodenumber = zeros(N,1);
parent = zeros(N,1);
yfitnode = zeros(N,1);
cutvar = zeros(N,1);
cutpoint = zeros(N,1);
children = zeros(N,2);
nodeprob = zeros(N,1);
resuberr = zeros(N,1);
risk = zeros(N,1);
nodesize = zeros(N,1);
classprob = zeros(N,nclasses);
classcount = zeros(N,nclasses);
nodenumber(1) = 1;
assignednode = ones(N,1);
nextunusednode = 2;

% Get default or specified prior class probabilities
if ~isempty(args.prior)
    prior = args.prior;
    prior = prior(:)';
    prior = prior / sum(prior);
else
    prior = ones(1,nclasses) / nclasses;
end

% Get default or specified misclassification costs
if ~isempty(args.cost)
    cost = args.cost;
    havecosts = true;
else
    cost = ones(nclasses) - eye(nclasses);
    havecosts = false;
end

% Adjust priors if required to take misclassification costs into account
adjprior = prior;
if havecosts
  Cj = sum(cost,2)';
  pc = Cj .* prior;
  adjprior = pc / sum(pc);
end

% Keep processing nodes until done
tnode = 1;
while(tnode < nextunusednode)
    
    % Record information about this node
    %noderows = find(assignednode==tnode);
    noderows = assignednode==tnode;
    
    %Nnode = length(noderows);
    Nnode = sum(noderows);
    Cnode = C(noderows,:);
       
    % Compute class probabilities and related statistics for this node
    Njt = sum(Cnode,1);    % number in class j at node t
    Pjandt = prior .* Njt ./ Nj;
    Pjgivent = Pjandt / sum(Pjandt);
    misclasscost = Pjgivent * cost;
    [mincost,nodeclass] = min(misclasscost);
    yfitnode(tnode) = nodeclass;
    Pt = sum(Pjandt);
    nodeprob(tnode) = Pt;
    classprob(tnode,:) = Pjgivent;
    classcount(tnode,:) = Njt;
    pratio = adjprior ./ Nj;
    impure = sum(Pjgivent>0)>1;
  
    bestcrit          = -Inf;
    nodesize(tnode)   = Nnode;
    resuberr(tnode)   = mincost;
    risk(tnode)       = nodeprob(tnode) * resuberr(tnode);
    cutvar(tnode)     = 0;
    cutpoint(tnode)   = 0;
    children(tnode,:) = 0;
   
   % Consider splitting this node
    if (Nnode >= args.split_min) && impure      % split only large impure nodes
        %Xnode = X(noderows,:);
        bestvar = 0;
        bestcut = 0;

        % Find the best of all possible splits - if args.random_m is
        % non-zero, randomly select M = random_m integers first
        if args.random_m
            %Note: avoid randsample because of STATS license, use randperm
            %instead
            poss_vars = randperm(nvars);
            poss_vars = poss_vars(1:args.random_m);
        else
            args.random_m = nvars;
            poss_vars = 1:nvars;
        end
        for ii = 1:args.random_m;
            jvar = poss_vars(ii);
            [x,idx] = sort(X(noderows,jvar)); % get sorted jth x variable
            
             % Determine if there's anything to split along this variable
             maxeps = max(eps(x(1)), eps(x(end)));
             if x(1)+maxeps > x(end)
                continue;
             end
             %rows = find(x(1:end-1)+maxeps < x(2:end));
             rows = x(1:end-1)+maxeps < x(2:end);
             if ~any(rows)
                continue;
             end

            Ccum = cumsum(Cnode(idx,:)); % cum. class counts
            
            %Now do the actual splitting
            [critval,cutval] = ...
                compute_critval(x, Ccum, rows, pratio, Pt, impurity, critfun, bestcrit);

            % Change best split if this one is best so far
            if critval>bestcrit
                bestcrit = critval;
                bestvar = jvar;
                bestcut = cutval;
            end
        end

        % Split this node using the best rule found
        if bestvar~=0
            %x = Xnode(:,bestvar);

            cutvar(tnode) = bestvar;
            cutpoint(tnode) = bestcut;
            %leftside = x<=bestcut;
            %rightside = ~leftside;
            leftside = X(:,bestvar)<=bestcut;
           
            children(tnode,:) = nextunusednode + (0:1);
            %assignednode(noderows(leftside)) = nextunusednode;
            %assignednode(noderows(rightside)) = nextunusednode+1;
            assignednode(noderows & leftside) = nextunusednode;
            assignednode(noderows & ~leftside) = nextunusednode+1;
            nodenumber(nextunusednode+(0:1)) = nextunusednode+(0:1)';
            parent(nextunusednode+(0:1)) = tnode;
            nextunusednode = nextunusednode+2;
        end
    end
    tnode = tnode + 1;
end

%Save all output to tree structure
topnode        = nextunusednode - 1;
tree.alpha     = []; %for compatibility with other classregtree fcns
tree.prunelist = []; %for compatibility with other classregtree fcns
tree.catcols   = []; %for compatibility with other classregtree fcns
tree.method    = 'classification';
tree.node      = nodenumber(1:topnode);
tree.parent    = parent(1:topnode);
tree.class     = yfitnode(1:topnode);
tree.var       = cutvar(1:topnode);
tree.cut       = cutpoint(1:topnode);
tree.children  = children(1:topnode,:);
tree.nodeprob  = nodeprob(1:topnode);
tree.nodeerr   = resuberr(1:topnode);
tree.risk      = risk(1:topnode);
tree.nodesize  = nodesize(1:topnode);
tree.npred     = nvars;
tree.names     = names;
tree.prior     = prior;
tree.nclasses  = nclasses;
tree.cost      = cost;
tree.classprob = classprob(1:topnode,:);
tree.classcount= classcount(1:topnode,:);
tree.classname = cnames;
%Can remove bad splits or prune now (commented out for now)
tree = removebadsplits(tree);

if args.prune
   tree = mb_prune_tree('tree', tree);
end

end

%----------------------------------------------------
function v=gdi(p)
%GDI Gini diversity index

v=1-sum(p.^2,2);
end

%----------------------------------------------------
function v=twoing(Pleft, P1, Pright, P2)
%TWOING Twoing index

v = 0.25 * Pleft .* Pright .* sum(abs(P1-P2),2).^2;
end

%----------------------------------------------------
function v=deviance(p)
%DEVIANCE Deviance

v = -2 * sum(p .* log(max(p,eps(class(p)))), 2);
end

%----------------------------------------------------
function [critval,cutval] = compute_critval(x,Ccum,rows,pratio,Pt,impurity,critfun,bestcrit)
%COMPUTE_CRITVAL Get critical value for splitting node in classification tree.
   
% First get all possible split points

% Get arrays showing left/right class membership at each split
%nsplits = length(rows);
nsplits = sum(rows);

% Split between each pair of distinct ordered values
Csplit1 = Ccum(rows,:);
Csplit2 = Ccum(size(Ccum,1)*ones(nsplits,1),:) - Csplit1; % repmat(Ccum(end,:),nsplits,1) - Csplit1;

% Get left/right class probabilities at each split
temp = pratio(ones(nsplits,1),:); %repmat(pratio,nsplits,1);
P1 = temp .* Csplit1;
P2 = temp .* Csplit2;
Ptleft  = sum(P1,2);
Ptright = sum(P2,2);
nclasses = size(P1,2);
wuns = ones(1,nclasses);
P1 = P1 ./ Ptleft(:,wuns);   %repmat(Ptleft,1,nclasses);
P2 = P2 ./ Ptright(:,wuns);  %repmat(Ptright,1,nclasses);

% Get left/right node probabilities
Pleft = Ptleft ./ Pt;
Pright = 1 - Pleft;

% Evaluate criterion as impurity or otherwise
if impurity
   %crit = - Pleft.*feval(critfun,P1);
   crit = - Pleft.*(1-sum(P1.^2,2));
   t = (crit>bestcrit);   % compute 2nd term only if it would make a difference
   if any(t)
      %crit(t) = crit(t) - Pright(t).*feval(critfun,P2(t,:));
      crit(t) = crit(t) - Pright(t).*(1-sum(P2(t,:).^2,2));
   end
else
   crit = feval(critfun, Pleft, P1, Pright, P2);
end

% Return best split point, but bail out early if no improvement
critval = max(crit);
if critval<bestcrit
   cutval = 0;
   return;
end

maxloc = find(crit==critval);
if length(maxloc)>1
   maxloc = maxloc(1+floor(length(maxloc)*rand));
end

cutloc = find(rows, maxloc);
cutval = (x(cutloc(end)) + x(cutloc(end)+1))/2;
end

% --------------------------------------
function tree = removebadsplits(tree)
%REMOVEBADSPLITS Remove splits that contribute nothing to the tree.

N = length(tree.node);
isleaf = (tree.var==0)';   % no split variable implies leaf node
isntpruned = true(1,N);
doprune = false(1,N);
risk = tree.risk';
adjfactor = (1 - 100*eps(class(risk)));

% Work up from the bottom of the tree
while(true)
    % Find "twigs" with two leaf children
    branches = find(~isleaf & isntpruned);
    twig = branches(sum(isleaf(tree.children(branches,:)),2) == 2);
    if isempty(twig)
        break; % must have just the root node left
    end

    % Find twigs to "unsplit" if the error of the twig is no larger
    % than the sum of the errors of the children
    Rtwig = risk(twig);
    kids = tree.children(twig,:);
    Rsplit = sum(risk(kids),2);
    unsplit = Rsplit >= Rtwig'*adjfactor;
    if any(unsplit)
        % Mark children as pruned, and mark twig as now a leaf
        isntpruned(kids(unsplit,:)) = 0;
        twig = twig(unsplit);   % only these to be marked on next 2 lines
        isleaf(twig) = 1;
        doprune(twig) = 1;
    else
        break;
    end
end

% Remove splits that are useless
if any(doprune)
    tree = mb_prune_tree('tree', tree, 'nodes', find(doprune));
end
end