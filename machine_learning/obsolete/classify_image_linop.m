function [probability_image] = classify_image_linop(varargin)
%CLASSIFY_IMAGE *Insert a one line summary here*
%   [probability_image] = classify_image(image_in,forest)
%
% Inputs:
%      image_in - Image to classify
%
%      forest - Random forest classifier (model)
%
%
% Outputs:
%      probability_image - Pixel-wise probability of belonging to a bar
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 27-Jan-2010
% Author: Michael Berks
% Email : michael.berks@postgrad.man.ac.uk
% Phone : +44 (0)161 275 1241
% Copyright: (C) University of Manchester

% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    {'image_in',... % the mandatory arguments
    'forest'}, ...
    'forest_type', 'isbe',...
    'win_size', 3, ...
    'num_levels', 5,...
    'do_max', 0,...
    'num_angles', 8,...
    'num_trees', [], ...
    'max_size', 128,...
    'use_probs', 0);

%Constants computed from arguments
pad_w = floor(args.win_size/2); %half size of window size

win_idx = -pad_w:pad_w;

[ROW COL] = size(args.image_in);

%compute number of parts we need to break image into
r_parts = ceil(ROW / args.max_size);
c_parts = ceil(COL / args.max_size);

probability_image = zeros(size(args.image_in));
args.image_in = padarray(args.image_in, [pad_w pad_w], 'replicate');

%Go through each segment
for rp = 1:r_parts
    for cp = 1:c_parts
        
        sr = 1 + (rp-1)*args.max_size;
        er = min(rp*args.max_size, ROW);
        sc = 1 + (cp-1)*args.max_size;
        ec = min(cp*args.max_size, COL);
        
        [part_cols part_rows] = meshgrid(sc:ec,sr:er);
        num_samples_part = numel(part_cols);
        
        %Make copies of sample rows and cols at positions of local window patch
        win_rows = repmat(part_rows(:)*ones(1,args.win_size) + ones(num_samples_part,1)*win_idx, 1, args.win_size);
        win_cols = kron(part_cols(:)*ones(1,args.win_size) + ones(num_samples_part,1)*win_idx, ones(1,args.win_size));

        %Now get linop samples from the image
        linop_samples = line_operator_octave_subset(...
            args.image_in, args.num_angles,  args.num_levels, win_rows+pad_w, win_cols+pad_w); clear win_rows win_cols;
    
        if args.do_max
            %get the maximum response across orientations
            linop_samples = squeeze(max(linop_samples, [], 3));
        end
    
        %Reshape
        test_data = reshape(linop_samples, num_samples_part, []);
        
        %Now classify/predict the data using chosen RF method
        switch args.forest_type
            case 'breiman'
                [labels votes] = classRF_predict(test_data, args.forest);
                prob_image_part = reshape(votes(:,1) / args.forest.ntree, size(part_cols));
                
            case 'isbe'
                if ~isempty(args.num_trees)
                    if args.num_trees<=length(args.forest.trees)
                        args.forest.trees(args.num_trees+1:end)=[];
                    else
                        error('The number of trees bigger than the total number of trees in the random forests! ');
                    end
                end
                true_idx = find(str2double(args.forest.classname));
                if args.use_probs
                    [labels votes] = mb_random_forest_prob_predict(args.forest, test_data);
                else
                    [labels votes] = mb_random_forest_class_predict(args.forest, test_data);
                end
                prob_image_part = reshape(votes(:,true_idx) / length(args.forest.trees), size(part_cols)); %#ok
                
            case 'isbe_boot'
                if ~isempty(args.num_trees)
                    if args.num_trees<=length(args.forest.trees)
                        args.forest.trees(args.num_trees+1:end)=[];
                    else
                        error('The number of trees bigger than the total number of trees in the random forests! ');
                    end
                end
                true_idx = find(str2double(args.forest.classname));
                [labels votes] = mb_random_forest_class_predict(args.forest, test_data);
                prob_image_part = reshape(votes(:,true_idx) / length(args.forest.trees), size(part_cols)); %#ok
                true_idx = find(str2double(args.forest.classname));
                [labels votes] = mb_random_forest_class_predict_boot(args.forest, test_data);
                prob_image_part = reshape(votes(:,true_idx) / length(args.forest.trees), size(part_cols)); %#ok
                
            case 'regression'
                [orientations] = mb_random_forest_reg_predict(args.forest, test_data);
                prob_image_part = reshape(orientations, size(part_cols));
                
            otherwise
                error('Type of random forest not recognised');
        end
        
        probability_image(sr:er, sc:ec) = prob_image_part;
        
    end
end

% %For
% curr_sample = 1;
% %loop through image pixels selecting points
% for cc = (1:COL) + pad_w
%     for rr = (1:ROW) + pad_w
%
%
%         %Extract and sample vector
%         sample_vector = pad_dt(rr+win_idx, cc+win_idx, :);
%
%         %Convert sample vector to mag/phase form and save in test data
%         test_data(curr_sample,:)= [abs(sample_vector(:))' angle(sample_vector(:))'];
%         curr_sample = curr_sample + 1;
%     end
% end
% clear pad_dt;
%
% %Now classify the data using either ours or Breiman's code
% if strcmpi(args.forest_type, 'breiman')
%     [labels votes] = classRF_predict(test_data, forest);
%     probability_image = reshape(votes(:,1) / forest.ntree, size(image_in));
% elseif strcmpi(args.forest_type, 'isbe')
%     [labels votes] = mb_random_forest_class_predict(forest, test_data);
%     probability_image = reshape(votes(:,1) / length(forest.trees), size(image_in));
% elseif strcmpi(args.forest_type, 'isbe_boot')
%     [labels votes] = mb_random_forest_class_predict_boot(forest, test_data);
%     probability_image = reshape(votes(:,1) / length(forest.trees), size(image_in));
% elseif strcmpi(args.forest_type, 'regression')
%     [orientations] = mb_random_forest_reg_predict(forest, test_data);
%     probability_image = reshape(orientations, size(image_in));
% else
%     error('Type of random forest not recognised');
% end
%
%
% %probability_image = test_data;
