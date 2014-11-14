function [random_forest] = mb_random_forest_two_images(varargin)

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
%   random_forest: structure containing the following fields
%       trees: Nx1 cell, each containing a single classification tree
%       D = number of input variables;
%       nclasses = number of classes in output labels;
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
    {... % the mandatory arguments
    'image1',...
    'image2',...
    'num_train',...
    'd',...
    'tree_dir'}, ...
    'mask1', [],...
    'mask2', [],...
    'do_test1', 1,...
    'do_test2', 0,...
    'use_probs', 0,...
    'num_levels', 5,...
    'win_size', 3,...
    'do_max', 0,...
    'feature_type', 'all',...
    'max_test_size', 256,...
    'tree_root', [],...
    'n_trees', 100,...
    'prior', [],... % the optional arguments
    'cost', [],...
    'split_criterion', 'gdi',...
    'split_min', 100, ...
    'prune', 0,...
    'minimise_size', 1,...
    'names', [],...
    'quiet', 0,...
    'save_path', [],...
    'incremental_save', 1);

%Get input arguments for whole forest from fields of args
n_trees = args.n_trees;
tree_dir = [args.tree_root args.tree_dir];

if ~isempty(tree_dir)
    if ~strcmp(tree_dir(end), '/') && ~strcmp(tree_dir(end), '\');
        tree_dir = [tree_dir filesep];
    end
    if ~exist(tree_dir, 'dir')
        mkdir(tree_dir);
    end
end

if isempty(args.mask1)
    args.mask1 = true(size(args.image1));
end
if isempty(args.mask2)
    args.mask2 = true(size(args.image2));
end

%Now copy input args to tree_args to pass to individual tree building
tree_args.random_m = args.d;
tree_args.prior = args.prior;
tree_args.cost= args.cost;
tree_args.split_criterion = args.split_criterion;
tree_args.split_min = args.split_min;
tree_args.prune = args.prune;
tree_args.names = args.names;

%Also copy across sampling args (for convenience)
sample_args.num_levels = args.num_levels;
sample_args.win_size = args.win_size;
sample_args.do_max = args.do_max;
sample_args.feature_type = args.feature_type;

%Create random forest structure and pre-allocate space for each tree
random_forest.trees = cell(n_trees, 1);

%Compute individual classification trees for boostrap samples of X
for ii = 1:n_trees
    
    if ~args.quiet
        display(['Building tree ', num2str(ii), ' of ', num2str(n_trees)]);
    end
    
    tic;
    
    % Get training data (and test data if necessary)
    [training_data1 sample_idx1] = ...
        sample_training_data(args.image1, args.mask1, args.num_train, sample_args);
    [training_data2 sample_idx2] = ...
        sample_training_data(args.image2, args.mask2, args.num_train, sample_args);
    y = [true(size(training_data1,1),1); false(size(training_data2,1),1)];
    
    t = toc;
    display(['Time sampling data = ', num2str(t)]);
    
    if ii == 1
        % Get number of dimensions in samples
        D = size(training_data1,2);
        
        %If we're going to classify the images as well as constructing the
        %forest, pre-alloocate voting maps
        if args.do_test1   
            random_forest.image1_votes1 = zeros(size(args.image1));
            random_forest.image1_total_votes = zeros(size(args.image1));
        end
        if args.do_test2   
            random_forest.image2_votes1 = zeros(size(args.image2));
            random_forest.image2_total_votes = zeros(size(args.image2));
        end
        
        %Could check here that D and nclasses match for each dataset - but for
        % assume they always do
    end
    
    tic;
    %build a tree for the sampled data set
    tree = mb_tree_class_train([training_data1; training_data2], y, tree_args);
    t = toc;
    display(['Time building tree = ', num2str(t)]);
    
    tic;
    %Classify all points in the image execpt the training sample
    if args.do_test1
        
        %Make mask of everything in the image except the training sample
        tree_mask = args.mask1;
        tree_mask(sample_idx1) = false;
        
        %Increment the total votes in image1
        random_forest.image1_total_votes(tree_mask) = random_forest.image1_total_votes(tree_mask) + 1;

        %Classify the image and get the new vote counts for class1
        [random_forest.image1_votes1] =...
            classify_image(args.image1, tree_mask, tree, random_forest.image1_votes1, sample_args, args.max_test_size, args.use_probs);
        
    end
    if args.do_test2
        
        %Make mask of everything in the image except the training sample
        tree_mask = args.mask2;
        tree_mask(sample_idx2) = false;
        
        %Increment the total votes in image2
        random_forest.image2_total_votes(tree_mask) = random_forest.image2_total_votes(tree_mask) + 1;

        %Classify the image and get the new vote counts for class1
        [random_forest.image2_votes1] =...
            classify_image(args.image2, tree_mask, tree, random_forest.image2_votes1, sample_args, args.max_test_size, args.use_probs);
        
    end
    t = toc;
    display(['Time classifying images = ', num2str(t)]);
    
    %To save space, throwaway fields from the tree not required for each
    %tree (e.g. some fields such as nclasses need only be saved once per
    %forest, not in every tree)
    tree = rmfield(tree, {'prior', 'cost', 'names', 'nclasses'});
    if args.minimise_size
        tree = rmfield(tree,...
            {'classcount', 'classprob','nodeprob', 'nodeerr', 'risk', 'nodesize'}); %#ok
    end
    
    %Set filename at which to save tree
    tree_name = ['rf_tree', zerostr(ii, 4), '.mat'];
    save([tree_dir tree_name], 'tree'); clear tree;
    
    % Save the filepath to the tree
    random_forest.trees{ii} = tree_name;
    
    %If incremental save set, save forest structure after every tree
    if ~isempty(args.save_path) && (args.incremental_save)
        save(args.save_path, 'random_forest');
    end
    
end

%Save fields to random_forest output structure
random_forest.D = D;
random_forest.d = args.d;
random_forest.tree_dir = args.tree_dir;
random_forest.tree_root = args.tree_root;
random_forest.nclasses = 2;
random_forest.classname = [0 1];
random_forest.prior = tree_args.prior;
random_forest.cost = tree_args.cost;
random_forest.names = tree_args.names;
if ~isempty(args.save_path)
    save(args.save_path, 'random_forest');
end
%----------------- END of MAIN FUNCTION -----------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [training_data sample_idx] = sample_training_data(image_in, mask, num_train, sample_args)

win_size = sample_args.win_size; % the size of scanning window
num_levels = sample_args.num_levels; %number of levels of DT-CWT from which to extract coefficients
do_max = sample_args.do_max; %number of feature vectors to sample for training from each region
feature_type = sample_args.feature_type; %number of feature vectors to sample for training from each region

%Constants computed from arguments
pad_w = floor(win_size/2); %half size of window size
win_idx = -pad_w:pad_w;

%get training and testing idx from abnormal image
[row col] = size(image_in);

sample_idx = randperm(sum(mask(:)))';
sample_idx(num_train+1:end) = [];

image_idx = find(mask);
sample_idx = image_idx(sample_idx);
[rows cols] = ind2sub([row col], sample_idx);

if strcmpi(feature_type, 'ilp')
    
    dt = dtwavexfm2b(image_in, num_levels+1); clear image_in;
    
    dt_samples = squeeze(dt_to_pixel_subset(dt, rows(:), cols(:))); clear dt rr cc;
    
    mags = reshape(abs(dt_samples(:,:,1:num_levels)), num_train, []);
    phases = dt_samples(:,:,1:num_levels) .* conj(dt_samples(:,:,2:num_levels+1).^2);
    phases = reshape(atan2(imag(phases), abs(real(phases))), num_train, []);
    clear dt_samples;

    comp = mags.*exp(i*phases);

    [max_mags max_idx] = max(mags, [], 2);
    max_lev = ceil(max_idx/6);

    new_mags = zeros(size(mags));
    new_comp = zeros(size(comp));
    %----------------------------------------------------------------------
    for lev = 1:num_levels
        shift_idx = max_lev == lev;
        if any(shift_idx)
            lev_cols = (1:6)+6*(lev-1);
            weights = bsxfun(@rdivide, mags(shift_idx, lev_cols), sum(mags(shift_idx, lev_cols),2));

            for lev2 = 1:num_levels
                for ori = 1:6
                    lev_col = 6*(lev2-1)+ori;
                    new_mags(shift_idx, lev_col) = diag(weights * mags(shift_idx, lev_cols).');
                    new_comp(shift_idx, lev_col) = diag(weights * comp(shift_idx, lev_cols).');
                    weights = circshift(weights,[0 1]);
                end
            end
        end
    end

    training_data = [new_mags angle(new_comp)];
else
    
    % Create DT-CWT of image
    dt = dtwavexfm2b(image_in, num_levels); clear image_in;

    %Make copies of sample rows and cols at positions of local window patch
    rr = repmat(rows*ones(1,win_size) + ones(num_train,1)*win_idx, 1, win_size);
    cc = kron(cols*ones(1,win_size) + ones(num_train,1)*win_idx, ones(1,win_size));

    %Get interpolated dual-tree coefficients
    dt_samples = dt_to_pixel_subset(dt, rr, cc); clear dt rr cc;

    if do_max
        %get the maximum response across orientations
        dt_samples = squeeze(max(dt_samples, [], 3)); clear dt;
    end

    temp_samples=reshape(dt_samples, num_train, []);

    switch feature_type
        case 'all'
            training_data = [abs(temp_samples) angle(temp_samples)];

        case 'real'
            training_data = real(temp_samples);

        case 'mag'
            training_data = abs(temp_samples);

        case 'phase'
            training_data = angle(temp_samples);

        otherwise
            warning(['Feature type: ', args.feature_type, ' not recognised']); %#ok
            training_data = [abs(temp_samples) angle(temp_samples)];
    end
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [votes] = classify_image(image_in, mask, tree, votes, sample_args, max_size, use_probs)

win_size = sample_args.win_size; % the size of scanning window
num_levels = sample_args.num_levels; %number of levels of DT-CWT from which to extract coefficients
do_max = sample_args.do_max; %number of feature vectors to sample for training from each region
feature_type = sample_args.feature_type; %number of feature vectors to sample for training from each region

%Constants computed from arguments
pad_w = floor(win_size/2); %half size of window size

win_idx = -pad_w:pad_w;

[ROW COL] = size(image_in);

% Create DT-CWT of image
if strcmpi(feature_type, 'ilp')    
    dt = dtwavexfm2b(image_in, num_levels+1);
else
    dt = dtwavexfm2b(image_in, num_levels);
end

%compute number of parts we need to break image into
r_parts = ceil(ROW / max_size);
c_parts = ceil(COL / max_size);

%Go through each segment
for rp = 1:r_parts
    for cp = 1:c_parts
        
        sr = 1 + (rp-1)*max_size;
        er = min(rp*max_size, ROW);
        sc = 1 + (cp-1)*max_size;
        ec = min(cp*max_size, COL);
        
        %Get rows/cols subscripts for this part
        [part_cols part_rows] = meshgrid(sc:ec,sr:er);
        part_cols = part_cols(:);
        part_rows = part_rows(:);        
        part_idx = sub2ind([ROW COL], part_rows, part_cols);
        
        %Throw away pixels not belonging to the mask
        part_rows(~mask(part_idx)) = [];
        part_cols(~mask(part_idx)) = [];
        part_idx(~mask(part_idx)) = [];
            
        %check whether there's any pixels left to process
        if isempty(part_rows)
            continue;
        end
            
        num_samples_part = numel(part_cols);
        
        if strcmpi(feature_type, 'ilp')
            dt_samples = squeeze(dt_to_pixel_subset(dt, part_rows(:), part_cols(:)));
            
            mags = reshape(abs(dt_samples(:,:,1:num_levels)), num_samples_part, []);
            phases = dt_samples(:,:,1:num_levels) .* conj(dt_samples(:,:,2:num_levels+1).^2);
            phases = reshape(atan2(imag(phases), abs(real(phases))), num_samples_part, []);
            clear dt_samples;
            
            comp = mags.*exp(i*phases);

            [max_mags max_idx] = max(mags, [], 2);
            max_lev = ceil(max_idx/6);

            new_mags = zeros(size(mags));
            new_comp = zeros(size(comp));
            %----------------------------------------------------------------------
            for lev = 1:num_levels
                shift_idx = max_lev == lev;
                if any(shift_idx)
                    cols = (1:6)+6*(lev-1);
                    weights = bsxfun(@rdivide, mags(shift_idx, cols), sum(mags(shift_idx, cols),2));

                    for lev2 = 1:num_levels
                        for ori = 1:6
                            col = 6*(lev2-1)+ori;
                            new_mags(shift_idx, col) = diag(weights * mags(shift_idx, cols).');
                            new_comp(shift_idx, col) = diag(weights * comp(shift_idx, cols).');
                            weights = circshift(weights,[0 1]);
                        end
                    end
                end
            end

            test_data = [new_mags angle(new_comp)];
            
        else
            
            %Make copies of sample rows and cols at positions of local window patch
            win_rows = repmat(part_rows(:)*ones(1,win_size) + ones(num_samples_part,1)*win_idx, 1, win_size);
            win_cols = kron(part_cols(:)*ones(1,win_size) + ones(num_samples_part,1)*win_idx, ones(1,win_size));

            %Get interpolated dual-tree coefficients
            dt_samples = dt_to_pixel_subset(dt, win_rows, win_cols); clear win_rows win_cols;

            if do_max
                %get the maximum response across orientations
                dt_samples = squeeze(max(dt_samples, [], 3));
            end

            %Reshape
            dt_samples = reshape(dt_samples, num_samples_part, []);

            %Store as test data
            switch feature_type
                case 'all'
                    test_data = [abs(dt_samples) angle(dt_samples)];

                case 'real'
                    test_data = real(dt_samples);

                case 'mag'
                    test_data = abs(dt_samples);

                case 'phase'
                    test_data = angle(dt_samples);

                otherwise
                    warning(['Feature type: ', feature_type, ' not recognised']); %#ok
                    test_data = [abs(dt_samples) angle(dt_samples)];
            end
        end

        %Get predictions from tree
        [y_tree] = mb_tree_predict(tree, test_data, use_probs);
        
        %Update probabilities/vote counts
        if use_probs
            votes(part_idx) =  votes(part_idx) + y_tree(:,1);
        else
            votes(part_idx) =  votes(part_idx) + y_tree;
        end
        
        
    end
end
