function [probability_image] = classify_image_monogenic(varargin)
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
    'num_levels', 3,...
    'min_wavelength', 4,...
    'onf', 0.65,...
    'num_trees', [], ...
    'max_size', 128);

%Constants computed from arguments
pad_w = floor(args.win_size/2); %half size of window size

win_idx = -pad_w:pad_w;

[ROW COL] = size(args.image_in);

%compute number of parts we need to break image into
r_parts = ceil(ROW / args.max_size);
c_parts = ceil(COL / args.max_size);

probability_image = zeros(size(args.image_in));
args.image_in = padarray(args.image_in, [pad_w pad_w], 'replicate');

%Now get monogenic samples from the image
[local_amp local_phase] = monogenic(args.image_in, args.num_levels, args.min_wavelength, 2, args.onf);

%save sample into main training data
local_amp = reshape(local_amp(:,:,2:end), [], args.num_levels);
local_phase = reshape(local_phase(:,:,2:end), [], args.num_levels);

%Go through each segment
for rp = 1:r_parts
    for cp = 1:c_parts
        
        sr = 1 + (rp-1)*args.max_size;
        er = min(rp*args.max_size, ROW);
        sc = 1 + (cp-1)*args.max_size;
        ec = min(cp*args.max_size, COL);
        
        [part_cols part_rows] = meshgrid((sc:ec)+pad_w,(sr:er)+pad_w);
        
        image_idx = sub2ind(size(args.image_in), part_rows(:), part_cols(:));
    	num_samples_part = length(image_idx);
    
        %Make copies of indices at positions of local window patch
        patch_idx = image_idx * ones(1, args.win_size^2) + ...
            ones(num_samples_part,1)*...
            (kron((win_idx*(ROW+2*pad_w)), ones(1,args.win_size)) + repmat(win_idx, 1, args.win_size));
    
        test_data = zeros(num_samples_part, args.num_levels*args.win_size^2);
        for level = 1:args.num_levels
            amp_cols = (1:args.win_size^2) + (level-1)*args.win_size^2;
            test_data(:, amp_cols) = ...
                reshape(local_amp(patch_idx,level), num_samples_part, []);

            phase_cols = (1:args.win_size^2) + (args.num_levels + level-1)*args.win_size^2;
            test_data(:, phase_cols) = ...
                reshape(local_phase(patch_idx,level), num_samples_part, []);
        end
        
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
                [labels votes] = mb_random_forest_class_predict(args.forest, test_data);
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
