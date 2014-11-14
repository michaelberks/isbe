% function [y_fit, nodes] = tree_predict_image(tree, image_in)
function [fit_image, nodes_image] = tree_predict_image(varargin)
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

% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    {'image_in',... % the mandatory arguments
    'tree'}, ...
    'win_size', 3, ...
    'num_levels', 6,...
    'do_max', 0,...
    'feature_type', 'all',...
    'max_size', 128,...
    'mask', []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tree = args.tree;

%Constants computed from arguments
pad_w = floor(args.win_size/2); %half size of window size

win_idx = -pad_w:pad_w;

[ROW COL] = size(args.image_in);

% Create DT-CWT of image
if strcmpi(args.feature_type, 'ilp')    
    dt = dtwavexfm2b(args.image_in, args.num_levels+1);
else
    dt = dtwavexfm2b(args.image_in, args.num_levels);
end

r_parts = ceil(ROW / args.max_size);
c_parts = ceil(COL / args.max_size);

fit_image = zeros(size(args.image_in));
nodes_image = zeros(size(args.image_in));

%Go through each segment
for rp = 1:r_parts
    for cp = 1:c_parts
        
        sr = 1 + (rp-1)*args.max_size;
        er = min(rp*args.max_size, ROW);
        sc = 1 + (cp-1)*args.max_size;
        ec = min(cp*args.max_size, COL);
        
        %Get rows/cols subscripts for this part
        [part_cols part_rows] = meshgrid(sc:ec,sr:er);
        part_idx = sub2ind([ROW COL], part_rows, part_cols);
        
        %Check whether we've been given a mask to select specific pixels
        if ~isempty(args.mask)
            
            %Throw away pixels not belonging to the mask
            part_rows(~args.mask(part_idx)) = [];
            part_cols(~args.mask(part_idx)) = [];
            part_idx(~args.mask(part_idx)) = [];
            
            %check whether there's any pixels left to process
            if isempty(part_rows)
                continue;
            end
        end
        
        num_samples_part = numel(part_cols);  
        
        if strcmpi(args.feature_type, 'ilp')
            dt_samples = squeeze(dt_to_pixel_subset(dt, part_rows(:), part_cols(:)));
            
            mags = reshape(abs(dt_samples(:,:,1:args.num_levels)), num_samples_part, []);
            phases = dt_samples(:,:,1:args.num_levels) .* conj(dt_samples(:,:,2:args.num_levels+1).^2);
            phases = reshape(atan2(imag(phases), abs(real(phases))), num_samples_part, []);
            clear dt_samples;
            
            comp = mags.*exp(i*phases);

            [max_mags max_idx] = max(mags, [], 2);
            max_lev = ceil(max_idx/6);

            new_mags = zeros(size(mags));
            new_comp = zeros(size(comp));
            %----------------------------------------------------------------------
            for lev = 1:args.num_levels
                shift_idx = max_lev == lev;
                if any(shift_idx)
                    cols = (1:6)+6*(lev-1);
                    weights = bsxfun(@rdivide, mags(shift_idx, cols), sum(mags(shift_idx, cols),2));

                    for lev2 = 1:args.num_levels
                        for ori = 1:6
                            col = 6*(lev2-1)+ori;
                            new_mags(shift_idx, col) = diag(weights * mags(shift_idx, cols).');
                            new_comp(shift_idx, col) = diag(weights * comp(shift_idx, cols).');
                            weights = circshift(weights,[0 1]);
                        end
                    end
                end
            end

            X = [new_mags angle(new_comp)];
            
        else
            %Make copies of sample rows and cols at positions of local window patch
            win_rows = repmat(part_rows(:)*ones(1,args.win_size) + ones(num_samples_part,1)*win_idx, 1, args.win_size);
            win_cols = kron(part_cols(:)*ones(1,args.win_size) + ones(num_samples_part,1)*win_idx, ones(1,args.win_size));

            %Get interpolated dual-tree coefficients
            dt_samples = dt_to_pixel_subset(dt, win_rows, win_cols); clear win_rows win_cols;
        
            if args.do_max
                %get the maximum response across orientations
                dt_samples = squeeze(max(dt_samples, [], 3));
            end
            
            %Reshape
            dt_samples = reshape(dt_samples, num_samples_part, []);

            %Store as test data
            switch args.feature_type
                case 'all'
                    X = [abs(dt_samples) angle(dt_samples)];

                case 'real'
                    X = real(dt_samples);

                case 'mag'
                    X = abs(dt_samples);

                case 'phase'
                    X = angle(dt_samples);

                otherwise
                    warning(['Feature type: ', args.feature_type, ' not recognised']); %#ok
                    X = [abs(dt_samples) angle(dt_samples)];
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        nodes = doapply(tree, X, 1:N, 1, zeros(N, 1), 0, prunelist, 1);
        
        %Get class indices/regression values from tree for this set of nodes
        id = tree.class(nodes);
        
        %If a classification tree assign class labels to class indices
        if isequal(tree.method, 'classification')
            y_fit = tree.classname(id);
        else
            y_fit = id;
        end
        fit_image(part_idx) = str2double(y_fit);
        nodes_image(part_idx) = nodes;
       
    end
end


%------------------------------------------------
function nodes = doapply(tree,X,rows,thisnode,nodes,subtrees,prunelist,endcol)
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