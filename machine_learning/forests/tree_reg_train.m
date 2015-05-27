function [tree assignednode] = tree_reg_train(X, y, varargin)
%TREE_REG_TRAIN given input data and an output variable, construct a
%single regression tree
%   [tree] = tree_reg_train(X, y, varargin)
%
% Inputs:
%   X:      N x d matrix of input data, where each row is a
%           datapoint consisting of d input variables
%   y:      N x 1 vector of class labels for each data points
%
% In addition TREE_REG_TRAIN uses the U_PACKARGS interface function
% and so optional arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Optional Arguments:
%   random_m: The number of variables to randomly select at each node, from
%           which the optimal split will be chosen
%
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
% Example: [tree] = tree_reg_train(X, y)
%
% Notes: This function is a implementation adopted from the Matlab function
%   CLASSREGTREE. However TREE_REG_TRAIN allows more flexibility in
%   how splits are made at each node. In particular we can choose to select
%   a random subset of input variables from which to calculate the optimal
%   split at each node. It is intended that further functionality (e.g. 
%   choosing linear combinations of variables as features) will be added.
%   This allows "random forests" to built from the trees
%
% See also: TREE_PREDICT CLASSREGTREE RANDOM_FOREST_CLASS_TRAIN RANDOM_FOREST_REG_TRAIN
%
% Created: 13-Oct-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
			 '0', ... % non-strict mode
             'random_m', 0,...
             'split_criterion', 'dabs',...
			 'var_criterion', 'mabs',...
             'split_min', 100, ...
             'end_cut_min', 25, ...
			 'w_prior', 0, ...
             'prune', 0,...
             'do_circular', [],...
             'do_ubound', 1,...
             'impure_thresh', 1e-4,...
             'names', [],...
             'debug', 0);

%Check X is numeric
if ~isnumeric(X)
   error('stats:treefit:BadData','X must be a numeric matrix.');
end
%Check y is numeric
if ~isnumeric(y)
   error('stats:treefit:BadData','y must be a numeric matrix. For classification use TREE_CLASS_TRAIN');
end

%Remove any NaN from the data
t = any(isnan(X),2);
t = t | any(isnan(y),2);
if any(t)
   X(t,:) = [];
   y(t,:) = [];
end

%Workout whether to do circular dispersion (i.e. for orientations) or linear
% variance
if isempty(args.do_circular)
    do_circular = ~isreal(y);
else
    do_circular = args.do_circular;
end

%Get number of data samples and variables from X
[N,nvars] = size(X);

%Ensure that args.random_m is within sensible limits and
% workout if we're using the forced phase/magnitude choice on splitting
%variables
if args.random_m(1) || (args.random_m(1) < nvars)
    if length(args.random_m)==1 || ~args.random_m(2)
        split_vars = nvars;
        force_split = 0;
    else
        split_vars = nvars/2;
        force_split = 1;
    end
else
    args.random_m = nvars;
    poss_vars = 1:nvars;
end

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
    case 'midp',        crit_fun = @midp;
    case 'dabs',        crit_fun = @dabs;
    case 'dsq',         crit_fun = @dsq;
    case 'mabs',        crit_fun = @mabs;
    case 'ssq',         crit_fun = @ssq;
    case 'bob',         crit_fun = @bob;
    case 'acos',		crit_fun = @acosv;
    otherwise,     error([args.split_criterion 'is not a recognised ''splitcriterion'' parameter.'])
end

switch(args.var_criterion)
%                Criterion function   Is it an impurity measure?
%                ------------------   --------------------------
    case 'midp',        var_fun = @midp;
    case 'dabs',        var_fun = @dabs;
    case 'dsq',         var_fun = @dsq;
    case 'mabs',        var_fun = @mabs;
    case 'ssq',         var_fun = @ssq;
    case 'acos',		var_fun = @acosv;
    case 'bob',         var_fun = @bob;
    otherwise,     error(['Bad value for ''var_criterion'' parameter: ',args.var_criterion])
end

% tree structure fields:
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

nodenumber = zeros(N,1);
parent = zeros(N,1);
yfitnode = zeros(N,size(y,2));
cutvar = zeros(N,1);
cutpoint = zeros(N,1);
goodness = zeros(N,2); % goodness of split
children = zeros(N,2);
nodeprob = zeros(N,1);
resuberr = zeros(N,1);
risk = zeros(N,1);
nodesize = zeros(N,1);

if ispc && strcmp(get_username,'ptresadern')
	nodesamples = cell(N,1);
	nodeopts = zeros(N,2*args.random_m);
end

nodenumber(1) = 1;

assignednode = ones(N,1);
nextunusednode = 2;

% Keep processing nodes until done
tnode = 1;
while(tnode < nextunusednode)
    
    % Record information about this node
    %noderows = find(assignednode==tnode);
    noderows = assignednode==tnode;
    
    %Nnode = length(noderows);
    Nnode = sum(noderows);
    y_node = y(noderows,:);
   

    % Compute variance and related statistics for this node
	if ispc && strcmp(get_username,'ptresadern') && length(y_node)==1
		% only one example in this leaf - set variance to NaN
		y_mean = y_node;
		node_var = NaN;
		impure = false;
	else
		if do_circular
			%For orientations compute mean of orientation vectors
			y_mean = mean(y_node, 1);

			%The circular variance is then 1 minus the length of this vector
			node_var = 1 - abs(y_mean);

			%Check if impure - is this ratio appropriate for orientations?
			impure = (node_var > args.impure_thresh);
		else
			%For euclidean vectors compute standard mean and variance
			y_mean = mean(y_node); %Mean
			sst = norm(y_node-y_mean)^2;

			node_var = sst / Nnode; %Variance

			%Check if impure
			impure = (node_var > args.impure_thresh * resuberr(1));
		end
	end
   
    yfitnode(tnode,:)	= y_mean;
    nodeprob(tnode)		= Nnode/N;
    bestcrit			= -Inf;
    nodesize(tnode)		= Nnode;
    resuberr(tnode)		= node_var;
    risk(tnode)			= nodeprob(tnode) * resuberr(tnode);
    cutvar(tnode)		= 0;
    cutpoint(tnode)		= 0;
	goodness(tnode,:)	= NaN;
    children(tnode,:)	= 0;
		
	if ispc && strcmp(get_username,'ptresadern')
		nodeopts(tnode,:)	= 0;
		nodesamples{tnode}	= y_node;
	end
	
    % Consider splitting this node
    if (Nnode >= args.split_min) && impure % split only large impure nodes
        bestvar = 0;
        bestcut = 0;

        if args.do_ubound
            % compute an upper bound on the mean vector length (that we are trying 
            % to maximize)
            rows = 1:Nnode-1;
            if do_circular
                %Sort y in order of angles on cicrle
                [sorted_angles sort_idx] = sort( angle(y_node) );
                y_sort = y_node(sort_idx);

                %Get indices of unique values of sorted y
                all_idx = (1:Nnode)';
                uni_idx = [0; all_idx(diff(y_sort)~=0)];
                uni_len = length(uni_idx);

                %Use each unique vlaue as the starting point for 1 cut
                %on the circle then search for the optimal 2nd cut
                ubound = zeros(uni_len,1);
                for i_ref = 1:uni_len
                    %Circularly permute the y values - we should really
                    %shift the sorted angles too, but as we don't use the
                    %cut_val output there's no need...
                    y_cum	= cumsum(circshift(y_sort, -uni_idx(i_ref)));

                    %Get critical value
                    [ubound(i_ref)] = Ocritval(sorted_angles,y_cum,rows,crit_fun,var_fun,0.0);
                end

                % find maximum critical value over all splits
                ubound = max(ubound);
            else
                %Sort y and find optimal split
                [x,idx] = sort(y_node);
                y_cum = cumsum(y_node(idx,:)-y_mean);
                [ubound] = Rcritval(x,y_cum,rows);
            end
        else
            ubound = 0;
        end
        
		% now find the splitting variable and value that is closest to this
		% bound
        
        % Find the best of all possible splits - if args.random_m is
        % non-zero, randomly select M = random_m integers first
        if args.random_m < nvars
            %Note: avoid randsample because of STATS license, use randperm
            %instead
            poss_vars = randperm(split_vars);
            
            %compute offset: will be zero if we're not force splitting, if
            %we are force splitting, a 50/50 chance of 0 or split_vars
            var_offset = force_split*(rand > .5)*split_vars;
            poss_vars = poss_vars(1:args.random_m) + var_offset;
        else
            poss_vars = 1:nvars;
        end
        
		% now find the splitting variable and value that is closest to this
		% bound
        end_cut_min = args.end_cut_min;
        if end_cut_min < 1
            end_cut_min = round(end_cut_min * Nnode);
        end
        %% DEBUGGING
        if args.debug
            figure;
            ax = zeros(12,1);
        end
        
        %%
        for ii = 1:args.random_m
            jvar = poss_vars(ii);
            [x,idx] = sort(X(noderows,jvar)); % get sorted jth x variable

            % Determine if there's anything to split along this variable
            maxeps = max(eps(x(1)), eps(x(end)));
            if x(1)+maxeps > x(end)
                continue;
            end
            
            %Find rows where x changes value
            rows = find(x(1:end-1)+maxeps < x(2:end));
            
            %Modify rows by the minimum end_cut
            rows = rows(rows >= end_cut_min);
            rows = rows(rows < Nnode-end_cut_min);
            
            %Check if we've got any rows left
            if isempty(rows)
                continue;
            end
            
            %% DEBUGGING
            if args.debug
                subplot(4,6,rem(2*jvar-2,24)+1); hold on;
            end
            
            %%
            if do_circular
                % Compute cumulative sum of the orientation vectors sorted
                % by the j-th dimension of x
                y_cum = cumsum(y_node(idx,:));
                
                %Now do the actual splitting
                [critval,cutval]=Ocritval(x,y_cum,rows,crit_fun,var_fun,args.w_prior,args.debug);
            else
                % Compute cumulative sum of the y sorted
                % by the j-th dimension of x
                y_cum = cumsum(y_node(idx,:) - y_mean);
                
                %Now do the actual splitting
                [critval,cutval]=Rcritval(x,y_cum,rows);
            end
    
			opts(:,ii) = [jvar; critval]; % temporary
            
            %% DEBUGGING
            if args.debug
                y_l = y_node(X(noderows,jvar) < cutval);
                y_r = y_node(X(noderows,jvar) >= cutval);

                ori_l = atan2(imag(y_l), real(y_l));
                ori_r = atan2(imag(y_r), real(y_r));

                counts_l = hist(ori_l, linspace(-pi,pi,24));
                counts_r = hist(ori_r, linspace(-pi,pi,24));

                ax(jvar) = subplot(4,6,rem(2*jvar-1,24)+1); hold on;
                bar(linspace(-pi,pi,24), counts_l, 0.5, 'b');
                bar(linspace(-pi,pi,24) + pi/24, counts_r, 0.5, 'r');
                title(['Max value = ' num2str(critval)]);
            end
            %%
						
            % Change best split if this one is best so far
            if critval>bestcrit
                bestcrit = critval;
                bestvar = jvar;
                bestcut = cutval;
            end
        end
        if args.debug
            title(ax(bestvar), ['\bf Best value = ' num2str(bestcrit)]);
        end
        
		if ispc && strcmp(get_username,'ptresadern')
			opts = sortrows(opts',2)';
			opts = opts(:,end:-1:1);
		end
		
        % Split this node using the best rule found
        if bestvar~=0
            cutvar(tnode) = bestvar;
            cutpoint(tnode) = bestcut;
			goodness(tnode,:) = [bestcrit ubound]; % goodness of cut
            
			leftside = X(:,bestvar)<=bestcut;
			children(tnode,:) = nextunusednode + (0:1);
            assignednode(noderows & leftside) = nextunusednode;
            assignednode(noderows & ~leftside) = nextunusednode+1;
            nodenumber(nextunusednode+(0:1)) = nextunusednode+(0:1)';
            parent(nextunusednode+(0:1)) = tnode;
            nextunusednode = nextunusednode+2;
						
			if ispc && strcmp(get_username,'ptresadern')
				nodeopts(tnode,:) = opts(:)';
			end
		end
    end
    tnode = tnode + 1;       
end

topnode        = nextunusednode - 1;
tree.alpha     = []; %for compatibility with other classregtree fcns
tree.prunelist = []; %for compatibility with other classregtree fcns
tree.catcols   = []; %for compatibility with other classregtree fcns
tree.method    = 'regression';
tree.node      = nodenumber(1:topnode);
tree.parent    = parent(1:topnode);
tree.class     = yfitnode(1:topnode,:);
tree.var       = cutvar(1:topnode);
tree.cut       = cutpoint(1:topnode);
tree.children  = children(1:topnode,:);
tree.nodeprob  = nodeprob(1:topnode);
tree.nodeerr   = resuberr(1:topnode);
tree.risk      = risk(1:topnode);
tree.nodesize  = nodesize(1:topnode);
tree.npred     = nvars;
tree.names     = names;
tree.goodness  = goodness(1:topnode,:);

if ispc && strcmp(get_username,'ptresadern')
	tree.outputs   = nodesamples(1:topnode);
	tree.nodeopts  = nodeopts(1:topnode,:);
end

tree = removebadsplits(tree);

if args.prune
   tree = prune_tree('tree', tree);
end
end

%----------------------------------------------------
function [critval,cutval]=Rcritval(x,Ycum,rows)
%RCRITVAL Get critical value for splitting node in regression tree.
   
% First get all possible split points
% Split between each pair of distinct ordered values
YsplitL = Ycum(rows,:);
YsplitR = Ycum(end) - YsplitL;
nL = rows(:); 
nR = numel(x) - nL(:);

% Get left/right sum of squares at each potential split
%     <--- left ssx --->   <--- right ssx -->
ssx = (YsplitL.^2) ./ nL + ...
      (YsplitR.^2) ./ nR;

%Criti
critval = max(ssx);
maxloc = find(ssx==critval);
if length(maxloc)>1
    maxloc = maxloc(1+floor(length(maxloc)*rand));
end

cutloc = rows(maxloc);
cutval = (x(cutloc) + x(cutloc+1))/2;
end

%----------------------------------------------------
function [crit_val,cut_val]=Ocritval(x,uv_cum,rows,crit_fun,var_fun,w_prior,debug)
%OCRITVAL Get critical value for splitting node in regression tree where
%the regression targets are orientations stored as [u v] vectors
   
if ~exist('debug', 'var')
    debug = false;
end

% First get all possible split points
% Split between each pair of distinct ordered values
N = numel(x);

n_left = rows(:);
n_right = N - n_left;

%uv = [uv_cum(1); uv_cum(2:end)-uv_cum(1:end-1)]; - can this be deleted?

uv_sum_left = uv_cum(rows,:);
uv_sum_right = uv_cum(end,:)-uv_sum_left;

%******************Debugging**********************************************
% if 0
%     abs_left = (real(uv_sum_left).^2 + imag(uv_sum_left).^2) ./ (n_left.*n_left);
%     abs_right = (real(uv_sum_right).^2 + imag(uv_sum_right).^2) ./ (n_right.*n_right);
%     mean_abs = (abs_left + abs_right)/2;
%     [max_abs max_absi] = max(mean_abs);
%     [max_abs_l max_abs_li] = max(abs_left);
%     [max_abs_r max_abs_ri] = max(abs_right);
%     ssq_left = (real(uv_sum_left).^2 + imag(uv_sum_left).^2) ./ n_left;
%     ssq_right = (real(uv_sum_right).^2 + imag(uv_sum_right).^2) ./ n_right;
%     mean_ssq = (ssq_left + ssq_right)/2;
%     [max_ssq max_ssqi] = max(mean_ssq);
%     [max_ssq_l max_ssq_li] = max(ssq_left);
%     [max_ssq_r max_ssq_ri] = max(ssq_right);
%     
%     figure;
%     subplot(1,2,1); hold on; 
%     plot(n_left,abs_left, 'g'); plot([1 N], [max_abs_l max_abs_l], 'g:'); plot(n_left([max_abs_li max_abs_li]), [0 max_abs_l], 'g:');
%     plot(n_left,abs_right, 'r'); plot([1 N], [max_abs_r max_abs_r], 'r:'); plot(n_left([max_abs_ri max_abs_ri]), [0 max_abs_r], 'r:');
%     plot(n_left,mean_abs, 'b'); plot([1 N], [max_abs max_abs], 'b:'); plot(n_left([max_absi max_absi]), [0 max_abs], 'b:');
%     subplot(1,2,2); hold on;
%     plot(n_left,ssq_left, 'g'); plot([1 N], [max_ssq_l max_ssq_l], 'g:'); plot(n_left([max_ssq_li max_ssq_li]), [0 max_ssq_l], 'g:');
%     plot(n_left,ssq_right, 'r'); plot([1 N], [max_ssq_r max_ssq_r], 'r:'); plot(n_left([max_ssq_ri max_ssq_ri]), [0 max_ssq_r], 'r:');
%     plot(n_left,mean_ssq, 'b'); plot([1 N], [max_ssq max_ssq], 'b:'); plot(n_left([max_ssqi max_ssqi]), [0 max_ssq], 'b:');
% end
%**************************************************************************

%Compute critical values based on left/right sums
crit_vals = feval(crit_fun, uv_sum_left, uv_sum_right, ...
	                        n_left, n_right, w_prior, debug);

%Get maximum of the summed dispersions and find all indices at the max
%value - if more than one randomly sample from them
max_idx = find(crit_vals == max(crit_vals));
if length(max_idx) > 1
    max_idx = max_idx(ceil(length(max_idx)*rand));
end

% the actual critical value we want is something like the average length of
% the two mean vectors at the cut point
crit_val = feval(var_fun, uv_sum_left(max_idx), uv_sum_right(max_idx),...
                          n_left(max_idx), n_right(max_idx), w_prior, 0);

%Get the row correpsonding to the maximum index
cut_row = rows(max_idx);

%Compute the cut value as the average of the x value at this row and the
%next
cut_val = (x(cut_row) + x(cut_row+1))/2;
end

%% Functions for choosing split point
function crit_vals = midp(sum_left, sum_right, n_left, n_right, w_prior, debug) %#ok
% ignore all data values and tell it to pick the middle one

if ~exist('debug', 'var')
    debug = false;
end

crit_vals = zeros(length(sum_left),1);
crit_vals(floor(length(crit_vals)/2)) = 1;

end

function crit_vals = dabs(sum_left, sum_right, n_left, n_right, w_prior, debug)
% Compute squared length of the left/right mean vectors for any split
if ~exist('debug', 'var')
    debug = false;
end

d_left = (real(sum_left).^2 + imag(sum_left).^2) ./ (n_left.*n_left);
d_right = (real(sum_right).^2 + imag(sum_right).^2) ./ (n_right.*n_right);

% Return the negative difference of the mean squared lengths for any split
crit_vals = -(d_left - d_right).^2;

% add a prior term that favours balanced trees
if (w_prior>0)
	crit_vals = crit_vals - w_prior*((n_left-n_right)/n_left(end)).^2;
end

end

function crit_vals = dsq(sum_left, sum_right, n_left, n_right, w_prior, debug)
% Compute length of the left/right mean vectors for any split
if ~exist('debug', 'var')
    debug = false;
end

d_left = (real(sum_left).^2 + imag(sum_left).^2) ./ n_left;
d_right = (real(sum_right).^2 + imag(sum_right).^2) ./ n_right;

% Return the negative difference of the mean squared lengths for any split
crit_vals = -(d_left - d_right).^2;

% add a prior term that favours balanced trees
if (w_prior>0)
	crit_vals = crit_vals - w_prior*((n_left-n_right)/n_left(end)).^2;
end

end

function crit_vals = bob(sum_left, sum_right, n_left, n_right, w_prior, debug)
% Compute length of the left/right mean vectors for any split
if ~exist('debug', 'var')
    debug = false;
end

mean_left = sum_left ./ n_left;
mean_right = sum_right ./ n_right;

abs_l = abs(mean_left);
abs_r = abs(mean_right);
crit_vals = min(abs_l, abs_r);% -...
    %0.5*abs(mean_left./abs_l + mean_right./abs_r);

if debug
    plot(abs(mean_left), 'r');
    plot(abs(mean_right), 'b');
    plot(0.5*abs(mean_left + mean_right), 'g');
    plot(0.5*(abs(mean_left) + abs(mean_right)), 'm');
    plot(crit_vals, 'k--', 'linewidth', 2);
    set(gca, 'ylim', [0 1]);
end

% add a prior term that favours balanced trees
if (w_prior>0)
	crit_vals = crit_vals - w_prior*((n_left-n_right)/n_left(end)).^2;
end

end



function crit_vals = mabs(sum_left, sum_right, n_left, n_right, w_prior, debug) %#ok
if ~exist('debug', 'var')
    debug = false;
end

% Compute length of the left/right mean vectors for any split
d_left  = abs(sum_left)  ./ n_left;
d_right = abs(sum_right) ./ n_right;

% Return the average length of the two mean vectors for any split
crit_vals = (d_left + d_right)/2;
end

function crit_vals = ssq(sum_left, sum_right, n_left, n_right, w_prior, debug) %#ok
if ~exist('debug', 'var')
    debug = false;
end

% Compute length of the left/right mean vectors for any split
d_left = (real(sum_left).^2 + imag(sum_left).^2) ./ n_left;
d_right = (real(sum_right).^2 + imag(sum_right).^2) ./ n_right;

% Return the sum of the mean squared lengths for any split
crit_vals = d_left + d_right;

if debug
    plot(d_left, 'r');
    plot(d_right, 'b');
    plot(crit_vals/2, 'k--', 'linewidth', 2);
end

end

function crit_vals = acosv(sum_left, sum_right, n_left, n_right, w_prior, debug) %#ok
if ~exist('debug', 'var')
    debug = false;
end

% Compute length of the left/right mean vectors for any split
d_left = abs(sum_left) ./ n_left;
d_right = abs(sum_right) ./ n_right;

% Return the average length of the two mean vectors for any split
crit_vals = 1-(acos((d_left + d_right)/2)/(pi/2));
end

%% Remove bad splits from the tree

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
      break;            % must have just the root node left
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
   tree = prune_tree('tree', tree,'nodes',find(doprune));
end
end
