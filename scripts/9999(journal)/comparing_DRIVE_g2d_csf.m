function comparing_DRIVE_g2d_csf(repeat, varargin)

args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'exp_name',     unixenv('EXP_NAME', []), ...
    'image_dir',      [asymmetryroot,'data/retinograms/DRIVE/training/images/'], ...
    'vessel_mask_dir',      [asymmetryroot,'data/retinograms/DRIVE/training/vessel_masks/'], ...
    'fov_mask_dir',      [asymmetryroot,'data/retinograms/DRIVE/training/fov_masks/'], ...
    'ori_dir',      [asymmetryroot,'data/retinograms/DRIVE/training/orientations/'], ...
    'width_dir',      [asymmetryroot,'data/retinograms/DRIVE/training/width_maps/'], ...
    'selected_images', [],...
    'num_pts',      unixenv('NUM_SAMPLES', 2000), ...
    'num_trees',    unixenv('NUM_TREES', 2), ...
    'make_data',    unixenv('MAKE_DATA',1), ...
    'num_angles',   unixenv('NUM_ANGLES', 6),...
    'levels',       unixenv('LEVELS', [1 2 4 8 16]),...
    'do_orientation', unixenv('DO_ORIENTATION',1), ...
    'do_detection', unixenv('DO_DETECTION',1), ...
    'do_width', unixenv('DO_WIDTH',1), ...
    'win_sizes', unixenv('WIN_SIZES',[1 3]) ...
);
warning('off', 'load_uint8:missing_variables');

exp_name = args.exp_name;
num_pts = args.num_pts;
win_sizes = args.win_sizes;
levels = args.levels;
num_angles = args.num_angles;
exp_dir = [asymmetryroot 'experiments/DRIVE/comparing_representations/'];

warning('off', 'ASYM:unexpectedArgument');
%
%Make a load of data for bothing training and testing purposes
% Set up arguments for each decomposition type
d_args = cell(2,1);
d_args{1}.decomp_type = {'g2da'};
d_args{2}.decomp_type = {'h2da'};
for i_decomp = 1:2
    d_args{i_decomp}.sigma_range = 1;
    d_args{i_decomp}.num_angles = 18;
    d_args{i_decomp}.do_max = 0;
    d_args{i_decomp}.rotate = 0;
    d_args{i_decomp}.win_size = 3;
    d_args{i_decomp}.normalise = 0;
    d_args{i_decomp}.pca = [];
end

D = get_samples_per_channel(d_args{1});
D1 = D / 18;

%--------------------------------------------------------------------------
% ***************** GENERATE DATA *****************************************
%--------------------------------------------------------------------------

if args.make_data
    for i_level = levels
        
        if 0%exist([exp_dir 'training/' num2str(repeat) '/' num2str(i_level) '/bg_responses_gh2d.mat'], 'file');
            continue;
        end
        
        rand('twister', 1000 * repeat);
        randn('state', 1000 * repeat);
        display(['Generating data for ' exp_name ', level ' num2str(i_level)]);

        d_args{1}.sigma_range = i_level;
        d_args{2}.sigma_range = i_level;

        image_list = dir([args.image_dir '/*.mat']);
        fov_list = dir([args.fov_mask_dir '/*.mat']);
        vessel_list = dir([args.vessel_mask_dir '/*.mat']);
        ori_list = dir([args.ori_dir '/*.mat']);
        width_list = dir([args.width_dir '/*.mat']);

        %Pre-allocate space for the true line parameters
        true_oris_test = zeros(num_pts,1);
        true_widths_test = zeros(num_pts,1);
        true_centre_test = false(num_pts,1);    
        true_oris_train = zeros(num_pts,1);
        true_widths_train = zeros(num_pts,1);
        true_centre_train = false(num_pts,1);

        %Pre-allocate space for the responses
        responses_test = zeros(num_pts, D);
        bg_responses_test = zeros(num_pts, D);
        responses_train = zeros(num_pts, D);
        bg_responses_train = zeros(num_pts, D);

        %Check which images are selected - we'll assume if images are selected the
        %user has managed to index images within the corrcet range
        if isempty(args.selected_images)
            selected_images = 1:length(image_list);
        else
            selected_images = args.selected_images;
        end
        num_images = length(selected_images);

        %loop through each image sampling data
        curr_sample = 1;
        for i_image = 1:num_images
            this_im = selected_images(i_image);

            %Work out the number of samples to take
            num_samples_image = ...
                sample_from_binomial((num_pts + 1 - curr_sample), 1/(num_images+1-i_image), 1);

            if ~num_samples_image
                continue;
            end
            display(['Sampling ' num2str(num_samples_image) ' pixels from image ' num2str(this_im)]);

            %Load in image and masks
            vessel_mask = load_uint8([args.vessel_mask_dir vessel_list(this_im).name]);
            fov_mask = u_load([args.fov_mask_dir fov_list(this_im).name]);
            ori_map = load_uint8([args.ori_dir ori_list(this_im).name]);
            width_map = load_uint8([args.width_dir width_list(this_im).name]);
            ret = rgb2gray(u_load([args.image_dir image_list(this_im).name]));
            centre_mask = bwmorph(vessel_mask, 'thin', 'inf');

            %Check we have enough samples in data
            total_v_pts = sum(vessel_mask(:) & fov_mask(:));
            total_b_pts = sum(~vessel_mask(:) & fov_mask(:));

            %Get random sample of vessel pixels
            v_idx = find(vessel_mask & fov_mask);
            r_idx = randperm(total_v_pts);
            v_idx = v_idx(r_idx(1:2*num_samples_image));
            [v_rows v_cols] = ind2sub(size(ret), v_idx);

            %Get random sample of background pixels
            b_idx = find(~vessel_mask & fov_mask);
            r_idx = randperm(total_b_pts);
            b_idx = b_idx(r_idx(1:2*num_samples_image));
            [b_rows b_cols] = ind2sub(size(ret), b_idx);        

            %Save the sample labels for this image
            sample_idx = curr_sample:num_samples_image+curr_sample-1;

            true_oris_test(sample_idx) = ori_map(v_idx(1:end/2));
            true_widths_test(sample_idx) = width_map(v_idx(1:end/2));
            true_centre_test(sample_idx) = centre_mask(v_idx(1:end/2));

            true_oris_train(sample_idx) = ori_map(v_idx(end/2+1:end));
            true_widths_train(sample_idx) = width_map(v_idx(end/2+1:end));
            true_centre_train(sample_idx) = centre_mask(v_idx(end/2+1:end));

            %For each decomposition compute the responses
            v_responses = 0;
            b_responses = 0;
            for i_decomp = 1:2
                im_responses = compute_filter_responses(ret, d_args{i_decomp});
                v_responses = v_responses + ...
                    sqrt(3 - 2*i_decomp)*sample_image_features(im_responses, v_rows, v_cols, d_args{i_decomp});
                b_responses = b_responses + ...
                    sqrt(3 - 2*i_decomp)*sample_image_features(im_responses, b_rows, b_cols, d_args{i_decomp});
            end
            responses_test(sample_idx, :) = v_responses(1:end/2, :);                
            bg_responses_test(sample_idx, :) = b_responses(1:end/2, :);
            responses_train(sample_idx, :) =  v_responses(end/2+1:end, :);
            bg_responses_train(sample_idx, :) = b_responses(end/2+1:end, :);

            %Update the current sample count
            curr_sample = curr_sample + num_samples_image;

        end
        mkdir([exp_dir 'training/' num2str(repeat) '/' num2str(i_level)]);
        mkdir([exp_dir 'test/' num2str(repeat) '/' num2str(i_level)]);

        save([exp_dir 'test/' num2str(repeat) '/' num2str(i_level) '/responses_gh2d.mat'],...
            'responses_test');
        save([exp_dir 'test/' num2str(repeat) '/' num2str(i_level) '/bg_responses_gh2d.mat'],...
            'bg_responses_test');

        save([exp_dir 'training/' num2str(repeat) '/' num2str(i_level) '/responses_gh2d.mat'],...
            'responses_train');
        save([exp_dir 'training/' num2str(repeat) '/' num2str(i_level) '/bg_responses_gh2d.mat'],...
            'bg_responses_train');        

        if ~exist([exp_dir 'test/' num2str(repeat) '/true_labels.mat'], 'file')
            save([exp_dir 'test/' num2str(repeat) '/true_labels.mat'], 'true_*test');
        end
        if ~exist([exp_dir 'training/' num2str(repeat) '/true_labels.mat'], 'file')
            save([exp_dir 'training/' num2str(repeat) '/true_labels.mat'], 'true_*train');
        end
    end
end
%%
% --------------------------------------------------------------------------
% ********************** DETECTION ****************************************
% --------------------------------------------------------------------------
if args.do_detection
    rf_args.prediction_type = 'rf_regression';
    rf_args.n_trees = args.num_trees;
    rf_args.d = [];
    rf_args.w_prior = 0;
    rf_args.impure_thresh = 1.0000e-004;
    rf_args.split_min = 100;
    rf_args.end_cut_min = 25;
    rf_args.do_ubound = 0;
    rf_args.quiet = 1;
    rf_args.overwrite = 0;
    rf_args.minimise_size = 0;
    rf_args.split_criterion = 'gdi';
    rf_args.var_criterion = 'mabs';

    rf_args.sampling_args.sampling_method = 'sample_saved_training_data';
    rf_args.decomposition_args = [];
    rf_args.sampling_args.y = [true(num_pts,1); false(num_pts,1)];
    
    d_args{1}.feature_type = 'complex';
    for win_size = win_sizes
        bands = 1:(18/num_angles):18;
        
        %set common parameters for all decomp types and then compute vector sizes
        D_level = D1 * num_angles * 2;
        
        new_decomp_args = [];
        if win_size == 1
            new_decomp_args.win_size = 1; 
            D_level = D_level / 9;
        end
        new_decomp_args.feature_type = 'conj';
        new_decomp_args.bands = bands;
        D_levels =  length(levels) * D_level;

        rf_args.sampling_args.X = zeros(2*num_pts, D_levels);
        test_data = zeros(2*num_pts, D_levels);
        cols = 1:D_level;
        for i_level = levels;                

            display(['Sampling from level ' num2str(i_level)]);
            
            %Load data for this decomp type
            fg_data = u_load([exp_dir 'training/' num2str(repeat) '/' num2str(i_level) '/responses_gh2d.mat']);
            bg_data = u_load([exp_dir 'training/' num2str(repeat) '/' num2str(i_level) '/bg_responses_gh2d.mat']);
            training_data_i = [fg_data; bg_data];

            clear fg_data bg_data;
            fg_data = u_load([exp_dir 'test/' num2str(repeat) '/' num2str(i_level) '/responses_gh2d.mat']);
            bg_data = u_load([exp_dir 'test/' num2str(repeat) '/' num2str(i_level) '/bg_responses_gh2d.mat']);
            test_data_i = [fg_data; bg_data]; 
            clear fg_data bg_data;

            new_decomp_args = [];
            if win_size == 1
                new_decomp_args.win_size = 1; 
            end
            new_decomp_args.feature_type = 'conj';
            new_decomp_args.bands = bands; 
            rf_args.sampling_args.X(:, cols) =...
                convert_decomp_form(training_data_i, d_args{1}, new_decomp_args);
            test_data(:, cols) =...
                convert_decomp_form(test_data_i, d_args{1}, new_decomp_args);
            cols = cols + D_level;
        end

        %Set up directory to save RF to
        rf_dir = [exp_dir '/rfs/' num2str(repeat) '/detection/gh2d_all/' exp_name '/' num2str(win_size) '/'];
        rf_args.tree_dir = [rf_dir 'trees/'];
        results_dir = [exp_dir '/results/' num2str(repeat) '/detection/gh2d_all/' exp_name '/' num2str(win_size) '/'];
        mkdir(results_dir);

        %Train
        predictor = random_forest_class_train(rf_args);
        if isfield(predictor, 'D')             
            save([rf_dir 'predictor.mat'], 'predictor');
            display(['************ FOREST for ' exp_name ' complete!! **************']);

            %Test
            [~, votes] = random_forest_class_predict(predictor, test_data);
            predicted_lines = votes(:,2) / rf_args.n_trees;
            save([results_dir 'results1.mat'], 'predicted_lines');
        end
    end
end

%--------------------------------------------------------------------------
% *********************** ORIENTATION + WIDTH *****************************
%--------------------------------------------------------------------------
%
output_types = cell(1,0);

if args.do_orientation
    output_types{1,end+1} = 'orientation';
    load([exp_dir 'training/' num2str(repeat) '/true_labels.mat'], 'true_oris_train');
end
if args.do_width
    output_types{1,end+1} = 'width';
    load([exp_dir 'training/' num2str(repeat) '/true_labels.mat'], 'true_widths_train');
end

if ~isempty(output_types)

    %First build a forest for each original type
    rf_args.prediction_type = 'rf_regression';
    rf_args.n_trees = args.num_trees;
    rf_args.d = [];
    rf_args.w_prior = 0;
    rf_args.impure_thresh = 1.0000e-004;
    rf_args.split_min = 100;
    rf_args.end_cut_min = 25;
    rf_args.do_ubound = 0;
    rf_args.quiet = 1;
    rf_args.do_circular = [];
    rf_args.overwrite = 0;
    rf_args.minimise_size = 0;
    rf_args.split_criterion = 'ssq';
    rf_args.var_criterion = 'ssq';

    rf_args.sampling_args.sampling_method = 'sample_saved_training_data';
    rf_args.decomposition_args = [];
    d_args{1}.feature_type = 'complex';
    
    for output_type = output_types
        for win_size = win_sizes
            switch output_type{1}                       
                case 'width'
                    rf_args.sampling_args.y = true_widths_train;    
                case 'orientation'
                    rf_args.sampling_args.y = true_oris_train;
            end

            %set common parameters for all decomp types and then compute vector sizes
            bands = 1:(18/num_angles):18;
            D_level = D1 * num_angles * 2;          

            new_decomp_args = [];
            if win_size == 1
                new_decomp_args.win_size = 1; 
                D_level = D_level / 9;
            end
            new_decomp_args.feature_type = 'conj';
            new_decomp_args.bands = bands;
            
            D_levels =  length(levels) * D_level; 
            
            rf_args.sampling_args.X = zeros(num_pts, D_levels);
            test_data = zeros(num_pts, D_levels);
            
            cols = 1:D_level;
            for i_level = levels;                
                display(['Sampling from level ' num2str(i_level)]);
                
                %Load data for this decomp type
                training_data_i = u_load([exp_dir 'training/' num2str(repeat) '/' num2str(i_level) '/responses_gh2d.mat']);

                clear fg_data bg_data;
                test_data_i = u_load([exp_dir 'test/' num2str(repeat) '/' num2str(i_level) '/responses_gh2d.mat']);

             
                rf_args.sampling_args.X(:, cols) =...
                    convert_decomp_form(training_data_i, d_args{1}, new_decomp_args);
                test_data(:, cols) =...
                    convert_decomp_form(test_data_i, d_args{1}, new_decomp_args);
                cols = cols + D_level;
            end

            %Set up directory to save RF to
            rf_dir = [exp_dir '/rfs/' num2str(repeat) '/' output_type{1} '/gh2d_all/' exp_name '/' num2str(win_size) '/'];
            rf_args.tree_dir = [rf_dir 'trees/'];
            results_dir = [exp_dir '/results/' num2str(repeat) '/' output_type{1} '/gh2d_all/' exp_name '/' num2str(win_size) '/'];
            mkdir(results_dir);

            %Train
            predictor = random_forest_reg_train(rf_args);    
            if isfield(predictor, 'D')             
                save([rf_dir 'predictor.mat'], 'predictor');         

                %Test
                predictions = random_forest_reg_predict(predictor, test_data);
                save([results_dir 'results.mat'], 'predictions');
            end
        end
    end
end
%%
%--------------------------------------------------------------------------
%**************************************************************************
%--------------------------------------------------------------------------
