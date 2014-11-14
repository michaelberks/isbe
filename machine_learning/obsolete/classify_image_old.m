function [probability_image] = classify_image_old(varargin)
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
    'feature_type', 'all',...
    'num_trees', [], ...
    'max_size', 128,...
    'use_probs', 0,...
    'mask', []);

%Constants computed from arguments
pad_w = floor(args.win_size/2); %half size of window size

win_idx = -pad_w:pad_w;

[ROW COL] = size(args.image_in);

% %Create storage for test data
% test_data = zeros(ROW*COL, size_sample_vector);

% Create DT-CWT of image
dt = dtwavexfm2b(args.image_in, args.num_levels);

% %interpolation DT-CWT to full pixel grid
% dt_inter = dt_to_full_image(dt); clear dt;
%
% if args.do_max
%     %get the maximum response across orientations
%     dt_inter = squeeze(max(dt_inter, [], 3)); clear dt;
% end
%
% %pad the edges of subbands and label
% pad_dt = padarray(dt_inter, [pad_w pad_w], 'replicate'); clear max_dt;

%compute number of parts we need to break image into
r_parts = ceil(ROW / args.max_size);
c_parts = ceil(COL / args.max_size);

probability_image = zeros(size(args.image_in));

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
        
        %Now sample the training data - if using clock method we split here, all other methods
        %use the else clause
        if ~isempty(strfind(args.feature_type, 'clock'))

            %Get samples in clock representation
            [dt_samples] = sample_dt_polar13(dt, part_rows(:), part_cols(:), strcmpi(args.feature_type, 'clock_rotate')); 

            %Convert into magnitude and phase
            mags = abs(dt_samples);
            phases = angle(dt_samples); clear dt_samples

            %Discard surplus magnitude columns
            discard_idx = repmat(logical([ones(1,6) zeros(1,78)]), 1, args.num_levels); 
            mags(:,discard_idx) = [];

            %Save into the output data
            test_data = [mags phases];
        else
            num_samples_part = numel(part_cols);

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
                    test_data = [abs(dt_samples) angle(dt_samples)];

                case 'real'
                    test_data = real(dt_samples);

                case 'mag'
                    test_data = abs(dt_samples);

                case 'phase'
                    test_data = angle(dt_samples);

                case 'rotate'

                    %Convert samples into phase and magnitude
                    mags = abs(dt_samples);
                    phases = angle(dt_samples);

                    %Find the sub-band that corresponds to the maximum magnitude in
                    %each row
                    [dummy max_ori] = max(mags, [], 2);
                    max_ori = rem(max_ori-1,6)+1;

                    %Circular shift the data according to the maximum subband
                    for ori = 1:6
                        shift_idx = max_ori == ori;
                        for lev = 1:args.num_levels
                            cols = 6*(lev-1)+(1:6);
                            mags(shift_idx, cols) = circshift(mags(shift_idx, cols), [0 1-ori]);
                            phases(shift_idx, cols) = circshift(phases(shift_idx, cols), [0 1-ori]);
                        end
                    end

                    %Save into the output data
                    test_data = ...
                        [mags phases];

                otherwise
                    warning(['Feature type: ', args.feature_type, ' not recognised']); %#ok
                    test_data = [abs(dt_samples) angle(dt_samples)];
            end
        end
        
        %Now classify/predict the data using chosen RF method
        switch args.forest_type
            case 'breiman'
                [labels votes] = classRF_predict(test_data, args.forest);
                %prob_image_part = reshape(votes(:,1) / args.forest.ntree, size(part_cols));
                prob_image_part = votes(:,1) / args.forest.ntree;
                
            case 'isbe'
                if ~isempty(args.num_trees)
                    if args.num_trees<=length(args.forest.trees)
                        args.forest.trees(args.num_trees+1:end)=[];
                    else
                        error('The number of trees bigger than the total number of trees in the random forests! ');
                    end
                end
                if isnumeric(args.forest.classname)
                    true_idx = find(args.forest.classname);
                else
                    true_idx = find(str2double(args.forest.classname));
                end
                if args.use_probs
                    [labels votes] = mb_random_forest_prob_predict(args.forest, test_data);
                else
                    [labels votes] = mb_random_forest_class_predict(args.forest, test_data);
                end
                %prob_image_part = reshape(votes(:,true_idx) / length(args.forest.trees), size(part_cols)); %#ok
                prob_image_part = votes(:,true_idx) / length(args.forest.trees); %#ok
                
            case 'isbe_boot'
                if ~isempty(args.num_trees)
                    if args.num_trees<=length(args.forest.trees)
                        args.forest.trees(args.num_trees+1:end)=[];
                    else
                        error('The number of trees bigger than the total number of trees in the random forests! ');
                    end
                end
                if isnumeric(args.forest.classname)
                    true_idx = find(args.forest.classname);
                else
                    true_idx = find(str2double(args.forest.classname));
                end
                [labels votes] = mb_random_forest_class_predict_boot(args.forest, test_data);
                %prob_image_part = reshape(votes(:,true_idx) / length(args.forest.trees), size(part_cols)); %#ok
                prob_image_part = votes(:,true_idx) / length(args.forest.trees); %#ok
                
                
            case 'regression'
                %[orientations] = mb_random_forest_reg_predict(args.forest, test_data);
                %prob_image_part = reshape(orientations, size(part_cols));
                prob_image_part = mb_random_forest_reg_predict(args.forest, test_data);
            otherwise
                error('Type of random forest not recognised');
        end
        
        %probability_image(sr:er, sc:ec) = prob_image_part;
        probability_image(part_idx) = prob_image_part;
        
    end
end