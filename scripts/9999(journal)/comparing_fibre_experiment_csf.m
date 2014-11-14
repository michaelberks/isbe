function comparing_fibre_experiment_csf(idx, varargin)

args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'image_dir',      [asymmetryroot,'data/fibre/training/images/'], ...
    'fibre_mask_dir',      [asymmetryroot,'data/fibre/training/fibre_masks/'], ...
    'fov_mask_dir',      [asymmetryroot,'data/fibre/training/fov_masks/'], ...
    'ori_dir',      [asymmetryroot,'data/fibre/training/orientation_maps/'], ...
    'num_images', [],...
    'num_pts',      unixenv('NUM_SAMPLES', 2000), ...
    'num_trees',    unixenv('NUM_TREES', 2), ...
    'make_data',    unixenv('MAKE_DATA',1), ...
    'do_orientation', unixenv('DO_ORIENTATION',1), ...
    'do_detection', unixenv('DO_DETECTION',1), ...
    'do_tests', unixenv('DO_TESTS',1), ...
    'win_sizes', unixenv('WIN_SIZES',[1 3]) ...
);
warning('off', 'load_uint8:missing_variables');

[decomps, repeats] = meshgrid(1:11, 1:10);
i_decomp = decomps(idx);
repeat = repeats(idx);

num_pts = args.num_pts;
win_sizes = args.win_sizes;
exp_dir = [asymmetryroot 'experiments/fibre/comparing_representations/'];

warning('off', 'ASYM:unexpectedArgument');
%
%Make a load of data for bothing training and testing purposes
% Set up arguments for each decomposition type
d_args{1}.decomp_type = {'dt'};
d_args{1}.levels = 1:5;
d_args{1}.feature_shape = 'rect';
d_args{1}.feature_type = 'complex';
d_args{1}.do_max = 0;
d_args{1}.rotate = 0;
d_args{1}.use_nag = 0;     

d_args{3}.decomp_type = {'gabor'};
d_args{3}.num_angles = 6;
d_args{3}.sigma_range = [1 2 4 8 16];	
d_args{3}.do_max = 0;
d_args{3}.rotate = 0;
d_args{3}.feature_type = 'complex';

d_args{4}.decomp_type = {'mono'};
d_args{4}.num_levels = 5;
d_args{4}.min_wavelength = 4;
d_args{4}.onf = 0.65;

d_args{5}.decomp_type = {'g1d'};
d_args{5}.sigma_range = [1 2 4 8 16];
            
d_args{6}.decomp_type = {'g2d'};
d_args{6}.sigma_range = [1 2 4 8 16];
            
d_args{7}.decomp_type = {'g2di'};
d_args{7}.sigma_range = [1 5];
            
d_args{8}.decomp_type = {'h2d'};
d_args{8}.sigma_range = [1 2 4 8 16];

d_args{9}.decomp_type = {'g2da'};
d_args{9}.sigma_range = [1 2 4 8 16];
d_args{9}.num_angles = 6;
d_args{9}.do_max = 0;
d_args{9}.rotate = 0;

d_args{10}.decomp_type = {'gabori'};
d_args{10}.sigma_range = [1 5];
d_args{10}.num_angles = 6;
d_args{10}.do_max = 0;
d_args{10}.feature_type = 'complex';

d_args{2}.decomp_type = {'h2da'};
d_args{2}.sigma_range = [1 2 4 8 16];
d_args{2}.num_angles = 6;
d_args{2}.do_max = 0;
d_args{2}.rotate = 0;

%set common parameters for all decomp types and then compute vector sizes
for i_d = 1:10
    d_args{i_d}.win_size = 3;
    d_args{i_d}.normalise = 0;
    d_args{i_d}.pca = [];
end

D = get_samples_per_channel(d_args{i_decomp});

%--------------------------------------------------------------------------
% ***************** GENERATE DATA *****************************************
%--------------------------------------------------------------------------

if args.make_data
    rand('twister', 1000 * repeat);
    randn('state', 1000 * repeat);
    
    display(['Generating data for ' d_args{i_decomp}.decomp_type{1}]);
    
    image_list = dir([args.image_dir '/*.mat']);
    fov_list = dir([args.fov_mask_dir '/*.mat']);
    fibre_list = dir([args.fibre_mask_dir '/*.mat']);
    ori_list = dir([args.ori_dir '/*.mat']);

    %Check which images are selected - we'll assume if images are selected the
    %user has managed to index images within the corrcet range
    if isempty(args.num_images)
        num_images = length(image_list);
    else
        num_images = args.num_images;
    end
    rp = randperm(num_images);
    selected_images = cell(2,1);
    selected_images{1} = rp(1:floor(num_images/2));
    selected_images{2} = rp(ceil(num_images/2):end);
    data_type = {'training', 'test'};
    for i_data = 1:2
        %Update the number of images for this training example
        num_images = length(selected_images{i_data});
        
        %Pre-allocate space for the true line parameters
        true_oris = zeros(num_pts,1);
        true_centre = false(num_pts,1);    

        %Pre-allocate space for the responses
        responses = zeros(num_pts, D);
        bg_responses = zeros(num_pts, D);

        %loop through each image sampling data
        curr_sample = 1;
        for i_image = 1:num_images
            this_im = selected_images{i_data}(i_image);

            %Work out the number of samples to take
            num_samples_image = ...
                sample_from_binomial((num_pts + 1 - curr_sample), 1/(num_images+1-i_image), 1);

            if ~num_samples_image
                continue;
            end
            
            %Load in image and masks
            fibre_mask = load_uint8([args.fibre_mask_dir fibre_list(this_im).name]);
            fov_mask = u_load([args.fov_mask_dir fov_list(this_im).name]);
            ori_map = load_uint8([args.ori_dir ori_list(this_im).name]);
            ret = u_load([args.image_dir image_list(this_im).name]);
            centre_mask = bwmorph(fibre_mask, 'thin', 'inf');

            %Check we have enough samples in data
            total_v_pts = sum(fibre_mask(:) & fov_mask(:));
            total_b_pts = sum(~fibre_mask(:) & fov_mask(:));
            
            num_samples_image = min(num_samples_image, total_v_pts);
            display(['Sampling ' num2str(num_samples_image) ' pixels from image ' num2str(this_im)]);

            if ~num_samples_image
                continue;
            end
            
            %Get random sample of fibre pixels
            v_idx = find(fibre_mask & fov_mask);
            r_idx = randperm(total_v_pts);
            v_idx = v_idx(r_idx(1:num_samples_image));
            [v_rows v_cols] = ind2sub(size(ret), v_idx);

            %Get random sample of background pixels
            b_idx = find(~fibre_mask & fov_mask);
            r_idx = randperm(total_b_pts);
            b_idx = b_idx(r_idx(1:num_samples_image));
            [b_rows b_cols] = ind2sub(size(ret), b_idx);        

            %Save the sample labels for this image
            sample_idx = curr_sample:num_samples_image+curr_sample-1;

            true_oris(sample_idx) = ori_map(v_idx);
            true_centre(sample_idx) = centre_mask(v_idx);

            %For each decomposition compute the responses
            im_responses = compute_filter_responses(ret, d_args{i_decomp}); 
            responses(sample_idx, :) = ...
                sample_image_features(im_responses, v_rows, v_cols, d_args{i_decomp});
            bg_responses(sample_idx, :) = ...
                sample_image_features(im_responses, b_rows, b_cols, d_args{i_decomp});

            %Update the current sample count
            curr_sample = curr_sample + num_samples_image;

        end
        mkdir([exp_dir data_type{i_data} '/' num2str(repeat)]);
        mkdir([exp_dir data_type{i_data} '/' num2str(repeat)]);

        save([exp_dir data_type{i_data} '/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat'],...
            'responses');
        save([exp_dir data_type{i_data} '/' num2str(repeat) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat'],...
            'bg_responses');        
        save([exp_dir data_type{i_data} '/' num2str(repeat) '/true_labels.mat'], 'true_*');
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
    
    %
    % 1) Original data
    %
    if ismember(0, args.do_tests) && ismember(i_decomp, 1:11);

        for win_size = win_sizes

            fg_data = u_load([exp_dir 'test/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
            bg_data = u_load([exp_dir 'test/' num2str(repeat) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
            test_data = [fg_data; bg_data]; 
            clear fg_data bg_data;
            orig_data = zeros(20*num_pts, size(test_data,2));
            
            for i_repeat = 0:9 
                %Load data for this decomp type
                fg_data = u_load([exp_dir 'training/'...
                    num2str(repeat+i_repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
                bg_data = u_load([exp_dir 'training/'...
                    num2str(repeat+i_repeat) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
                fg_rows = (1:num_pts) + i_repeat*num_pts;
                bg_rows = fg_rows + 4*num_pts;
                orig_data(fg_rows,:) = fg_data;
                orig_data(bg_rows,:) = bg_data; 
                clear fg_data bg_data;
            end
            
            new_decomp_args = [];
            if win_size == 1
                new_decomp_args.win_size = win_size;
            end
            if isfield(d_args{i_decomp}, 'feature_type');
                new_decomp_args.feature_type = 'conj';
            end 
            if ~isempty(new_decomp_args)
                orig_data = convert_decomp_form(orig_data, d_args{i_decomp}, new_decomp_args);
                test_data = convert_decomp_form(test_data, d_args{i_decomp}, new_decomp_args);
            end
            
            for i_pts = 5*num_pts:num_pts:10*num_pts
                %Set up directory to save RF to
                rf_dir = [exp_dir '/rfs/' num2str(repeat)...
                    '/detection/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/' num2str(i_pts) '/'];
                rf_args.tree_dir = [rf_dir 'trees/'];
                results_dir = [exp_dir '/results/' num2str(repeat)...
                    '/detection/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/' num2str(i_pts) '/'];
                mkdir(results_dir);

                %Train
                rf_args.sampling_args.y = [true(i_pts,1); false(i_pts,1)];
                rf_args.sampling_args.X = orig_data([1:i_pts 10*num_pts+(1:i_pts)],:);
                predictor = random_forest_class_train(rf_args);
                if isfield(predictor, 'D')             
                    save([rf_dir 'predictor.mat'], 'predictor');         
                    display(['************ FOREST for ' d_args{i_decomp}.decomp_type{1} ' complete!! **************']);

                    %Test
                    %load([rf_dir 'predictor.mat'], 'predictor');
                    [~, votes, all_votes] = random_forest_class_predict(predictor, test_data);
                    predicted_lines = votes(:,2) / rf_args.n_trees;
                    save([results_dir 'all_votes.mat'], 'all_votes');
                    save([results_dir 'results1.mat'], 'predicted_lines');
                end
            end
        end
    end
    rf_args.sampling_args.y = [true(num_pts,1); false(num_pts,1)];
    %
    if ismember(1, args.do_tests) && ismember(i_decomp, 1:11);

        for win_size = win_sizes

            %Load data for this decomp type
            fg_data = u_load([exp_dir 'training/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
            bg_data = u_load([exp_dir 'training/' num2str(repeat) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
            orig_data = [fg_data; bg_data]; 
            clear fg_data bg_data;
            fg_data = u_load([exp_dir 'test/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
            bg_data = u_load([exp_dir 'test/' num2str(repeat) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
            test_data = [fg_data; bg_data]; 
            clear fg_data bg_data;

            new_decomp_args = [];
            if win_size == 1
                new_decomp_args.win_size = win_size;
            end
            if isfield(d_args{i_decomp}, 'num_angles') && (d_args{i_decomp}.num_angles == 18)
                new_decomp_args.bands = 1:3:16;
            end
            if isfield(d_args{i_decomp}, 'feature_type');
                new_decomp_args.feature_type = 'conj';
            end 
            if ~isempty(new_decomp_args)
                orig_data = convert_decomp_form(orig_data, d_args{i_decomp}, new_decomp_args);
                test_data = convert_decomp_form(test_data, d_args{i_decomp}, new_decomp_args);
            end
            
            %Set up directory to save RF to
            rf_dir = [exp_dir '/rfs/' num2str(repeat)...
                '/detection/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/'];
            rf_args.tree_dir = [rf_dir 'trees/'];
            results_dir = [exp_dir '/results/' num2str(repeat)...
                '/detection/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/'];
            mkdir(results_dir);

            %Train
            rf_args.sampling_args.X = orig_data;
            predictor = random_forest_class_train(rf_args);
            if isfield(predictor, 'D')             
                save([rf_dir 'predictor.mat'], 'predictor');         
                display(['************ FOREST for ' d_args{i_decomp}.decomp_type{1} ' complete!! **************']);

                %Test
                [~, votes, all_votes] = random_forest_class_predict(predictor, test_data);
                predicted_lines = votes(:,2) / rf_args.n_trees;
                save([results_dir 'all_votes.mat'], 'all_votes');
                save([results_dir 'results1.mat'], 'predicted_lines');
            end
        end
    end
    
    %
    % 2) Complex representations
    %
    if ismember(2, args.do_tests) && ismember(i_decomp, [1 3]);

        fg_data = u_load([exp_dir 'test/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
        bg_data = u_load([exp_dir 'test/' num2str(repeat) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
        test_data = [fg_data; bg_data]; 
        clear fg_data bg_data;
        fg_data = u_load([exp_dir 'training/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
        bg_data = u_load([exp_dir 'training/' num2str(repeat) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
        orig_data = [fg_data; bg_data]; 
        clear fg_data bg_data;

        for win_size = win_sizes
            % 1) Complex forms: all (mag/phase), real/imag, conj, mag, phase
            for feature_type = {'real_imag', 'real', 'imag'}%, 'conj', 'mag', 'phase', 'all'}

                rf_dir = [exp_dir '/rfs/' num2str(repeat) '/detection/' d_args{i_decomp}.decomp_type{1} '/feature_types/' feature_type{1} '/' num2str(win_size) '/'];
                results_dir = [exp_dir '/results/' num2str(repeat) '/detection/' d_args{i_decomp}.decomp_type{1} '/feature_types/' feature_type{1} '/' num2str(win_size) '/'];
                mkdir(results_dir);

                %Args to reform data
                new_decomp_args = [];
                new_decomp_args.feature_type = feature_type{1};
                new_decomp_args.win_size = win_size;
                
                %Train
                rf_args.sampling_args.X = convert_decomp_form(orig_data, d_args{i_decomp}, new_decomp_args);
                rf_args.tree_dir = [rf_dir 'trees/'];
                predictor = random_forest_class_train(rf_args);    
                if isfield(predictor, 'D')             
                    save([rf_dir 'predictor.mat'], 'predictor');

                    %Test
                    test_data_i = convert_decomp_form(test_data, d_args{i_decomp}, new_decomp_args);
                    [dummy, votes] = random_forest_class_predict(predictor, test_data_i);
                    predicted_lines = votes(:,2) / rf_args.n_trees;
                    save([results_dir 'results.mat'], 'predicted_lines');
                end
            end
        end
    end

    % 
    % 3) Different number of levels
    %
    if ismember(3, args.do_tests) && ismember(i_decomp, [1 3 5 6]);

        %Load data for this decomp type
        fg_data = u_load([exp_dir 'training/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
        bg_data = u_load([exp_dir 'training/' num2str(repeat) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
        orig_data = [fg_data; bg_data]; 
        clear fg_data bg_data;
        fg_data = u_load([exp_dir 'test/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
        bg_data = u_load([exp_dir 'test/' num2str(repeat) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
        test_data = [fg_data; bg_data]; 
        clear fg_data bg_data;
                
        for i_level = 4:5

            for win_size = win_sizes

                rf_dir = [exp_dir '/rfs/' num2str(repeat) '/detection/' d_args{i_decomp}.decomp_type{1} '/levels/' num2str(i_level) '/' num2str(win_size) '/'];
                results_dir = [exp_dir '/results/' num2str(repeat) '/detection/' d_args{i_decomp}.decomp_type{1} '/levels/' num2str(i_level) '/' num2str(win_size) '/'];
                mkdir(results_dir);

                %Args to reform data
                new_decomp_args = [];
                new_decomp_args.feature_type = 'conj';
                new_decomp_args.levels = i_level;
                new_decomp_args.win_size = win_size;

                %Train       
                rf_args.sampling_args.X = convert_decomp_form(orig_data, d_args{i_decomp}, new_decomp_args);
                rf_args.tree_dir = [rf_dir 'trees/'];
                predictor = random_forest_class_train(rf_args);    
                if isfield(predictor, 'D')             
                    save([rf_dir 'predictor.mat'], 'predictor');

                    %Test
                    test_data_i = convert_decomp_form(test_data, d_args{i_decomp}, new_decomp_args);
                    [dummy, votes] = random_forest_class_predict(predictor, test_data_i);
                    predicted_lines = votes(:,2) / rf_args.n_trees;
                    save([results_dir 'results.mat'], 'predicted_lines');
                end
            end
        end
    end
    %
    % 4) Rotate and do max
    %
    if ismember(4, args.do_tests) && ismember(i_decomp, [1 3 9 11])

        fg_data = u_load([exp_dir 'test/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
        bg_data = u_load([exp_dir 'test/' num2str(repeat) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
        test_data = [fg_data; bg_data]; 
        clear fg_data bg_data;
        
        %Load data for this decomp type
        fg_data = u_load([exp_dir 'training/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
        bg_data = u_load([exp_dir 'training/' num2str(repeat) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
        orig_data = [fg_data; bg_data]; 
        clear fg_data bg_data;

        for win_size = win_sizes

            for reform_type = 1:2

                if win_size==1 && reform_type==2
                    continue;
                end

                new_decomp_args = [];
                if reform_type == 1
                    rf_dir = [exp_dir '/rfs/' num2str(repeat) '/detection/' d_args{i_decomp}.decomp_type{1} '/rotate/' num2str(win_size)  '/'];
                    results_dir = [exp_dir '/results/' num2str(repeat) '/detection/' d_args{i_decomp}.decomp_type{1} '/rotate/' num2str(win_size) '/'];
                    new_decomp_args.rotate = 1;
                else
                    rf_dir = [exp_dir '/rfs/' num2str(repeat) '/detection/' d_args{i_decomp}.decomp_type{1} '/do_max/' num2str(win_size)  '/'];
                    results_dir = [exp_dir '/results/' num2str(repeat) '/detection/' d_args{i_decomp}.decomp_type{1} '/do_max/' num2str(win_size) '/'];
                    new_decomp_args.do_max = 1;
                end
                mkdir(results_dir);

                %Args to reform data
                new_decomp_args.feature_type = 'conj';

                %Train       
                rf_args.sampling_args.X = convert_decomp_form(orig_data, d_args{i_decomp}, new_decomp_args);
                rf_args.tree_dir = [rf_dir 'trees/'];
                predictor = random_forest_class_train(rf_args);    
                if isfield(predictor, 'D')             
                    save([rf_dir 'predictor.mat'], 'predictor'); 

                    %Test
                    test_data_i = convert_decomp_form(test_data, d_args{i_decomp}, new_decomp_args);
                    [dummy, votes] = random_forest_class_predict(predictor, test_data_i);
                    predicted_lines = votes(:,2) / rf_args.n_trees;
                    save([results_dir 'results.mat'], 'predicted_lines');
                end
            end
        end
    end
    %
    % 5) Different ways of stitching Gaussian filters together
    %
    % G" + G'
    if ismember(5, args.do_tests) && i_decomp == 5

%         g1_data = u_load([exp_dir 'training/' num2str(repeat) '/responses_g1d.mat']);
%         g2_data = u_load([exp_dir 'training/' num2str(repeat) '/responses_g2d.mat']);
%         g1_data_bg = u_load([exp_dir 'training/' num2str(repeat) '/bg_responses_g1d.mat']);
%         g2_data_bg = u_load([exp_dir 'training/' num2str(repeat) '/bg_responses_g2d.mat']);
%         rf_args.sampling_args.X = [g1_data g2_data; g1_data_bg g2_data_bg]; 
%         clear g1_data g2_data;
% 
%         g1_data = u_load([exp_dir 'test/' num2str(repeat) '/responses_g1d.mat']);
%         g2_data = u_load([exp_dir 'test/' num2str(repeat) '/responses_g2d.mat']);
%         g1_data_bg = u_load([exp_dir 'test/' num2str(repeat) '/bg_responses_g1d.mat']);
%         g2_data_bg = u_load([exp_dir 'test/' num2str(repeat) '/bg_responses_g2d.mat']);
%         test_data = [g1_data g2_data; g1_data_bg g2_data_bg]; 
%         clear g1_data g2_data;
% 
%         rf_dir = [exp_dir '/rfs/' num2str(repeat) '/detection/g12d/'];
%         results_dir = [exp_dir '/results/' num2str(repeat) '/detection/g12d/'];
%         mkdir(results_dir);
% 
%         rf_args.tree_dir = [rf_dir 'trees/'];
%         predictor = random_forest_class_train(rf_args);    
%         if isfield(predictor, 'D')             
%             save([rf_dir 'predictor.mat'], 'predictor');                             
% 
%             [dummy, votes] = random_forest_class_predict(predictor, test_data);
%             predicted_lines = votes(:,2) / rf_args.n_trees;
%             save([results_dir 'results.mat'], 'predicted_lines');
%         end

        %G" + H"
        for win_size = win_sizes
            new_decomp_args = [];
            new_decomp_args.win_size = win_size;
            h2_data = convert_decomp_form(...
                u_load([exp_dir 'training/' num2str(repeat) '/responses_h2da.mat']),...
                d_args{2}, new_decomp_args);
            g2_data = convert_decomp_form(...
                u_load([exp_dir 'training/' num2str(repeat) '/responses_g2da.mat']),...
                d_args{9}, new_decomp_args);
            h2_data_bg = convert_decomp_form(...
                u_load([exp_dir 'training/' num2str(repeat) '/bg_responses_h2da.mat']),...
                d_args{2}, new_decomp_args);
            g2_data_bg = convert_decomp_form(...
                u_load([exp_dir 'training/' num2str(repeat) '/bg_responses_g2da.mat']),...
                d_args{9}, new_decomp_args);
            rf_args.sampling_args.X = [h2_data g2_data; h2_data_bg g2_data_bg];
            clear h2_data g2_data;

            h2_data = convert_decomp_form(...
                u_load([exp_dir 'test/' num2str(repeat) '/responses_h2da.mat']),...
                d_args{2}, new_decomp_args);
            g2_data = convert_decomp_form(...
                u_load([exp_dir 'test/' num2str(repeat) '/responses_g2da.mat']),...
                d_args{9}, new_decomp_args);
            h2_data_bg = convert_decomp_form(...
                u_load([exp_dir 'test/' num2str(repeat) '/bg_responses_h2da.mat']),...
                d_args{2}, new_decomp_args);
            g2_data_bg = convert_decomp_form(...
                u_load([exp_dir 'test/' num2str(repeat) '/bg_responses_g2da.mat']),...
                d_args{9}, new_decomp_args);
            test_data = [h2_data g2_data; h2_data_bg g2_data_bg];
            clear h2_data g2_data h2_data_bg g2_data_bg;

            rf_dir = [exp_dir '/rfs/' num2str(repeat) '/detection/gh2da/' num2str(win_size)  '/'];
            results_dir = [exp_dir '/results/' num2str(repeat) '/detection/gh2da/' num2str(win_size)  '/'];
            mkdir(results_dir);

            rf_args.tree_dir = [rf_dir 'trees/'];
            predictor = random_forest_class_train(rf_args);    
            if isfield(predictor, 'D')             
                save([rf_dir 'predictor.mat'], 'predictor');

                [dummy, votes] = random_forest_class_predict(predictor, test_data);
                predicted_lines = votes(:,2) / rf_args.n_trees;
                save([results_dir 'results.mat'], 'predicted_lines');
            end
        end
            
    end
    
    if ismember(6, args.do_tests)
        for win_size = win_sizes
            %set common parameters for all decomp types and then compute vector sizes
            all_args = d_args;
            all_decomps = [1 3 4 5 6 8];
            D_all = zeros(max(all_decomps),1);
            for i_d = all_decomps
                all_args{i_d}.win_size = win_size;
                all_args{i_d}.normalise = 0;
                all_args{i_d}.pca = [];
                if isfield(d_args{i_d}, 'feature_type');
                    all_args{i_d}.feature_type = 'conj';
                end

                d_args{i_d}.win_size = 3;
                d_args{i_d}.normalise = 0;
                d_args{i_d}.pca = [];

                D_all(i_d) = get_samples_per_channel(all_args{i_d});
            end

            rf_args.sampling_args.X = zeros(2*num_pts, sum(D_all));
            test_data = zeros(2*num_pts, sum(D_all));
            curr_col = 0;
            for i_d = all_decomps
                cols = curr_col + (1:D_all(i_d));
                curr_col = cols(end);

                %Load data for this decomp type
                fg_data = u_load([exp_dir 'training/' num2str(repeat) '/responses_' d_args{i_d}.decomp_type{1} '.mat']);
                bg_data = u_load([exp_dir 'training/' num2str(repeat) '/bg_responses_' d_args{i_d}.decomp_type{1} '.mat']);
                training_data_i = [fg_data; bg_data];

                clear fg_data bg_data;
                fg_data = u_load([exp_dir 'test/' num2str(repeat) '/responses_' d_args{i_d}.decomp_type{1} '.mat']);
                bg_data = u_load([exp_dir 'test/' num2str(repeat) '/bg_responses_' d_args{i_d}.decomp_type{1} '.mat']);
                test_data_i = [fg_data; bg_data]; 
                clear fg_data bg_data;

                new_decomp_args = [];
                if win_size == 1
                    new_decomp_args.win_size = 1; 
                end
                if isfield(d_args{i_d}, 'feature_type');
                    new_decomp_args.feature_type = 'conj';
                end 
                rf_args.sampling_args.X(:, cols) =...
                    convert_decomp_form(training_data_i, d_args{i_d}, new_decomp_args);
                test_data(:, cols) =...
                    convert_decomp_form(test_data_i, d_args{i_d}, new_decomp_args);
            end

            %Set up directory to save RF to
            rf_dir = [exp_dir '/rfs/' num2str(repeat) '/detection/all/' num2str(win_size) '/'];
            rf_args.tree_dir = [rf_dir 'trees/'];
            results_dir = [exp_dir '/results/' num2str(repeat) '/detection/all/' num2str(win_size) '/'];
            mkdir(results_dir);

            %Train
            predictor = random_forest_class_train(rf_args);
            if isfield(predictor, 'D')             
                save([rf_dir 'predictor.mat'], 'predictor');
                display(['************ FOREST for ' d_args{i_decomp}.decomp_type{1} ' complete!! **************']);

                %Test
                [~, votes] = random_forest_class_predict(predictor, test_data);
                predicted_lines = votes(:,2) / rf_args.n_trees;
                save([results_dir 'results1.mat'], 'predicted_lines');
            end
        end
    end
    
    if ismember(7, args.do_tests) && i_decomp == 2

        %G" + H"
        h2_data = u_load([exp_dir 'training/' num2str(repeat) '/responses_h2da.mat']);
        g2_data = u_load([exp_dir 'training/' num2str(repeat) '/responses_g2da.mat']);
        h2_data_bg = u_load([exp_dir 'training/' num2str(repeat) '/bg_responses_h2da.mat']);
        g2_data_bg = u_load([exp_dir 'training/' num2str(repeat) '/bg_responses_g2da.mat']);
        training_data = ...
            convert_complex_representation_back([g2_data h2_data; g2_data_bg h2_data_bg], 'real_imag'); 
        clear h2_data g2_data h2_data_bg g2_data_bg;

        h2_data = u_load([exp_dir 'test/' num2str(repeat) '/responses_h2da.mat']);
        g2_data = u_load([exp_dir 'test/' num2str(repeat) '/responses_g2da.mat']);
        h2_data_bg = u_load([exp_dir 'test/' num2str(repeat) '/bg_responses_h2da.mat']);
        g2_data_bg = u_load([exp_dir 'test/' num2str(repeat) '/bg_responses_g2da.mat']);
        test_data = ...
            convert_complex_representation_back([g2_data h2_data; g2_data_bg h2_data_bg], 'real_imag'); 
        clear h2_data g2_data g2_data_bg h2_data_bg;

        for feature_type = {'all'}%{'phase', 'mag', 'conj', 'real_imag'}

            rf_dir = [exp_dir '/rfs/' num2str(repeat) '/detection/gh2d/feature_types/' feature_type{1} '/'];
            results_dir = [exp_dir '/results/' num2str(repeat) '/detection/gh2d/feature_types/' feature_type{1} '/'];
            mkdir(results_dir);

            %Args to reform data
            rf_args.sampling_args.X = ...
                convert_complex_representation(training_data, feature_type{1}, 3);
            test_data_i = ...
                convert_complex_representation(test_data, feature_type{1}, 3);

            rf_args.tree_dir = [rf_dir 'trees/'];
            predictor = random_forest_class_train(rf_args);
            if isfield(predictor, 'D')             
                save([rf_dir 'predictor.mat'], 'predictor');

                %Test
                [~, votes] = random_forest_class_predict(predictor, test_data_i);
                predicted_lines = votes(:,2) / rf_args.n_trees;
                save([results_dir 'results1.mat'], 'predicted_lines');
            end
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
    %
    % 0) Original data with increasing num points
    %
    if ismember(0, args.do_tests) && ismember(i_decomp, 1:11);
        for output_type = output_types
            for win_size = win_sizes
                
                %Load data for this decomp type
                test_data = u_load([exp_dir 'test/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);

                orig_data = zeros(4*num_pts, size(test_data,2));
                orig_y = zeros(4*num_pts, 1);
                for i_repeat = 0:3 
                    %Load data for this decomp type
                    fg_data = u_load([exp_dir 'training/'...
                        num2str(repeat+i_repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
                    fg_rows = (1:num_pts) + i_repeat*num_pts;
                    orig_data(fg_rows,:) = fg_data;
                    clear fg_data;
                    
                    switch output_type{1}                       
                        case 'width'
                            load([exp_dir 'training/' num2str(repeat) '/true_labels.mat'], 'true_widths');
                            orig_y(fg_rows,:) = true_widths;    
                        case 'orientation'
                            load([exp_dir 'training/' num2str(repeat) '/true_labels.mat'], 'true_oris');
                            orig_y(fg_rows,:) = true_oris;
                    end
                end
                
                new_decomp_args = [];
                if win_size == 1
                    new_decomp_args.win_size = win_size;
                end
                if isfield(d_args{i_decomp}, 'num_angles') && (d_args{i_decomp}.num_angles == 18)
                    new_decomp_args.bands = 1:3:16;
                end
                if isfield(d_args{i_decomp}, 'feature_type');
                    new_decomp_args.feature_type = 'conj';
                end 
                if ~isempty(new_decomp_args)
                    orig_data = convert_decomp_form(orig_data, d_args{i_decomp}, new_decomp_args);
                    test_data = convert_decomp_form(test_data, d_args{i_decomp}, new_decomp_args);
                end
                
                for i_pts = 2*num_pts:num_pts:4*num_pts
                    %Set up directory to save RF to
                    rf_dir = [exp_dir '/rfs/' num2str(repeat)...
                        '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/' num2str(i_pts) '/'];
                    rf_args.tree_dir = [rf_dir 'trees/'];
                    results_dir = [exp_dir '/results/' num2str(repeat)...
                        '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/' num2str(i_pts) '/'];
                    mkdir(results_dir);

                    %Train
                    rf_args.sampling_args.y = orig_y(1:i_pts);
                    rf_args.sampling_args.X = orig_data(1:i_pts,:);
                    predictor = random_forest_reg_train(rf_args);
                    if isfield(predictor, 'D')             
                        save([rf_dir 'predictor.mat'], 'predictor');  
                        display(['************ FOREST for ' d_args{i_decomp}.decomp_type{1} ' complete!! **************']);

                        %Test
                        predictions = random_forest_reg_predict(predictor, test_data);
                        save([results_dir 'results.mat'], 'predictions');
                    end
                end
            end
        end
    end
    %
    % 1) Original data
    %
    load([exp_dir 'training/' num2str(repeat) '/true_labels.mat'], 'true_*');
    if ismember(1, args.do_tests) && ismember(i_decomp, 1:11);
        for output_type = output_types

            switch output_type{1}                       
                case 'width'
                    
                    orig_y = true_widths;    
                case 'orientation'
                    orig_y = true_oris;
            end

            for win_size = win_sizes

                %Load data for this decomp type
                orig_data = u_load([exp_dir 'training/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
                test_data = u_load([exp_dir 'test/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);

                new_decomp_args = [];
                if win_size == 1
                    new_decomp_args.win_size = win_size;
                end
                if isfield(d_args{i_decomp}, 'num_angles') && (d_args{i_decomp}.num_angles == 18)
                    new_decomp_args.bands = 1:3:16;
                end
                if isfield(d_args{i_decomp}, 'feature_type');
                    new_decomp_args.feature_type = 'conj';
                end 
                if ~isempty(new_decomp_args)
                    orig_data = convert_decomp_form(orig_data, d_args{i_decomp}, new_decomp_args);
                    test_data = convert_decomp_form(test_data, d_args{i_decomp}, new_decomp_args);
                end
                
                %Set up directory to save RF to
                rf_dir = [exp_dir '/rfs/' num2str(repeat)...
                    '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/'];
                rf_args.tree_dir = [rf_dir 'trees/'];
                results_dir = [exp_dir '/results/' num2str(repeat)...
                    '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/'];
                mkdir(results_dir);

                %Train
                rf_args.sampling_args.y = orig_y;
                rf_args.sampling_args.X = orig_data;
                predictor = random_forest_reg_train(rf_args);
                if isfield(predictor, 'D')             
                    save([rf_dir 'predictor.mat'], 'predictor');         
                    display(['************ FOREST for ' d_args{i_decomp}.decomp_type{1} ' complete!! **************']);

                    %Test
                    predictions = random_forest_reg_predict(predictor, test_data);
                    save([results_dir 'results.mat'], 'predictions');
                end
            end
        end
    end
    %
    % 2) Complex representations
    %
    %Dual-tree tests
    if ismember(2, args.do_tests) && ismember(i_decomp, [1 3]);
        for output_type = output_types

            switch output_type{1}                       
                case 'width'
                    rf_args.sampling_args.y = true_widths;    
                case 'orientation'
                    rf_args.sampling_args.y = true_oris;
            end

            %Load data for this decomp type
            orig_data = u_load([exp_dir 'training/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);

            %Load data for this decomp type
            test_data = u_load([exp_dir 'test/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);

            for win_size = win_sizes
                % 1) Complex forms: all (mag/phase), real/imag, conj, mag, phase
                for feature_type = {'real_imag', 'imag', 'real'}%, 'conj', 'mag', 'phase', 'all'}

                    rf_dir = [exp_dir '/rfs/' num2str(repeat) '/' output_type{1} '/'...
                        d_args{i_decomp}.decomp_type{1} '/feature_types/' feature_type{1} '/' num2str(win_size) '/'];
                    results_dir = [exp_dir '/results/' num2str(repeat) '/' output_type{1} '/'...
                        d_args{i_decomp}.decomp_type{1} '/feature_types/' feature_type{1} '/' num2str(win_size) '/'];
                    mkdir(results_dir);

                    %Args to reform data
                    new_decomp_args = [];
                    new_decomp_args.feature_type = feature_type{1};
                    new_decomp_args.win_size = win_size;

                    if isfield(d_args{i_decomp}, 'num_angles') && (d_args{i_decomp}.num_angles == 18)
                        new_decomp_args.bands = 1:3:16;
                    end

                    %Train
                    rf_args.sampling_args.X = convert_decomp_form(orig_data, d_args{i_decomp}, new_decomp_args);
                    rf_args.tree_dir = [rf_dir 'trees/'];
                    predictor = random_forest_reg_train(rf_args);    
                    if isfield(predictor, 'D')             
                        save([rf_dir 'predictor.mat'], 'predictor');         

                        %Test
                        test_data_i = convert_decomp_form(test_data, d_args{i_decomp}, new_decomp_args);
                        predictions = random_forest_reg_predict(predictor, test_data_i);
                        save([results_dir 'results.mat'], 'predictions');
                    end
                end
            end
        end
    end
    %
    % 3) Different number of levels
    %
    if ismember(3, args.do_tests) && ismember(i_decomp, [1 3 5 6]);
        
        for output_type = output_types

            switch output_type{1}                       
                case 'width'
                    rf_args.sampling_args.y = true_widths;    
                case 'orientation'
                    rf_args.sampling_args.y = true_oris;
            end
            
            %Load data for this decomp type
            orig_data = u_load([exp_dir 'training/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);

            %Load data for this decomp type
            test_data = u_load([exp_dir 'test/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);

            for win_size = win_sizes

                for i_level = 4:5

                    rf_dir = [exp_dir '/rfs/' num2str(repeat) '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/levels/' num2str(i_level) '/' num2str(win_size) '/'];
                    results_dir = [exp_dir '/results/' num2str(repeat) '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/levels/' num2str(i_level) '/' num2str(win_size) '/'];
                    mkdir(results_dir);


                    %Args to reform data
                    new_decomp_args = [];
                    new_decomp_args.feature_type = 'conj';
                    new_decomp_args.levels = i_level;
                    new_decomp_args.win_size = win_size;

                    %Train       
                    rf_args.sampling_args.X = convert_decomp_form(orig_data, d_args{i_decomp}, new_decomp_args);
                    rf_args.tree_dir = [rf_dir 'trees/'];
                    predictor = random_forest_reg_train(rf_args);    
                    if isfield(predictor, 'D')             
                        save([rf_dir 'predictor.mat'], 'predictor');         

                        %Test
                        test_data_i = convert_decomp_form(test_data, d_args{i_decomp}, new_decomp_args);
                        predictions = random_forest_reg_predict(predictor, test_data_i);
                        save([results_dir 'results.mat'], 'predictions');
                    end
                end

            end
        end
    end
    % 
    % 5) Different ways of stitching Gaussian filters together
    % G" + G'
    if ismember(5, args.do_tests) && i_decomp == 5
        for output_type = output_types

            switch output_type{1}                       
                case 'width'
                    rf_args.sampling_args.y = true_widths;    
                case 'orientation'
                    rf_args.sampling_args.y = true_oris;
            end

%             g1_data = u_load([exp_dir 'training/' num2str(repeat) '/responses_g1d.mat']);
%             g2_data = u_load([exp_dir 'training/' num2str(repeat) '/responses_g2d.mat']);
%             rf_args.sampling_args.X = [g1_data g2_data]; clear g1_data g2_data;
% 
%             g1_data = u_load([exp_dir 'test/' num2str(repeat) '/responses_g1d.mat']);
%             g2_data = u_load([exp_dir 'test/' num2str(repeat) '/responses_g2d.mat']);
%             test_data = [g1_data g2_data]; clear g1_data g2_data;
% 
%             rf_dir = [exp_dir '/rfs/' num2str(repeat) '/' output_type{1} '/g12d/'];
%             results_dir = [exp_dir '/results/' num2str(repeat) '/' output_type{1} '/g12d/'];
%             mkdir(results_dir);
% 
%             rf_args.tree_dir = [rf_dir 'trees/'];
%             predictor = random_forest_reg_train(rf_args);    
%             if isfield(predictor, 'D')             
%                 save([rf_dir 'predictor.mat'], 'predictor');  
%   
%                 predictions = random_forest_reg_predict(predictor, test_data, 0);
%                 save([results_dir 'results.mat'], 'predictions');
%             end

            %G" + H"
            for win_size = win_sizes
                new_decomp_args = [];
                new_decomp_args.win_size = win_size;
                h2_data = convert_decomp_form(...
                    u_load([exp_dir 'training/' num2str(repeat) '/responses_h2da.mat']),...
                    d_args{2}, new_decomp_args);
                g2_data = convert_decomp_form(...
                    u_load([exp_dir 'training/' num2str(repeat) '/responses_g2da.mat']),...
                    d_args{9}, new_decomp_args);
                rf_args.sampling_args.X = [h2_data g2_data]; clear h2_data g2_data;

                h2_data = convert_decomp_form(...
                    u_load([exp_dir 'test/' num2str(repeat) '/responses_h2da.mat']),...
                    d_args{2}, new_decomp_args);
                g2_data = convert_decomp_form(...
                    u_load([exp_dir 'test/' num2str(repeat) '/responses_g2da.mat']),...
                    d_args{9}, new_decomp_args);
                test_data = [h2_data g2_data]; clear h2_data g2_data;

                rf_dir = [exp_dir '/rfs/' num2str(repeat) '/' output_type{1} '/gh2da/' num2str(win_size)  '/'];
                results_dir = [exp_dir '/results/' num2str(repeat) '/' output_type{1} '/gh2da/' num2str(win_size)  '/'];
                mkdir(results_dir);

                rf_args.tree_dir = [rf_dir 'trees/'];
                predictor = random_forest_reg_train(rf_args);    
                if isfield(predictor, 'D')             
                    save([rf_dir 'predictor.mat'], 'predictor');         

                    predictions = random_forest_reg_predict(predictor, test_data, 0);
                    save([results_dir 'results.mat'], 'predictions');
                end
            end
        end
    end
    
    if ismember(6, args.do_tests)
        for win_size = win_sizes
            %set common parameters for all decomp types and then compute vector sizes
            all_args = d_args;
            all_decomps = [1 3 4 5 6 8];%
            D_all = zeros(max(all_decomps),1);
            for i_d = all_decomps
                all_args{i_d}.win_size = win_size;
                all_args{i_d}.normalise = 0;
                all_args{i_d}.pca = [];
                if isfield(d_args{i_d}, 'feature_type');
                    all_args{i_d}.feature_type = 'conj';
                end

                d_args{i_d}.win_size = 3;
                d_args{i_d}.normalise = 0;
                d_args{i_d}.pca = [];

                D_all(i_d) = get_samples_per_channel(all_args{i_d});
            end

            rf_args.sampling_args.X = zeros(num_pts, sum(D_all));
            rf_args.sampling_args.y = true_oris;
            
            test_data = zeros(num_pts, sum(D_all));
            curr_col = 0;
            for i_d = all_decomps
                cols = curr_col + (1:D_all(i_d));
                curr_col = cols(end);

                %Load data for this decomp type
                training_data_i = u_load([exp_dir 'training/' num2str(repeat) '/responses_' d_args{i_d}.decomp_type{1} '.mat']);
                test_data_i = u_load([exp_dir 'test/' num2str(repeat) '/responses_' d_args{i_d}.decomp_type{1} '.mat']);
                
                new_decomp_args = [];
                if win_size == 1
                    new_decomp_args.win_size = 1; 
                end
                if isfield(d_args{i_d}, 'feature_type');
                    new_decomp_args.feature_type = 'conj';
                end 
                rf_args.sampling_args.X(:, cols) =...
                    convert_decomp_form(training_data_i, d_args{i_d}, new_decomp_args);
                test_data(:, cols) =...
                    convert_decomp_form(test_data_i, d_args{i_d}, new_decomp_args);
            end

            %Set up directory to save RF to
            rf_dir = [exp_dir '/rfs/' num2str(repeat) '/orientation/all/' num2str(win_size) '/'];
            rf_args.tree_dir = [rf_dir 'trees/'];
            results_dir = [exp_dir '/results/' num2str(repeat) '/orientation/all/' num2str(win_size) '/'];
            mkdir(results_dir);

            %Train
            predictor = random_forest_reg_train(rf_args);    
            if isfield(predictor, 'D')             
                save([rf_dir 'predictor.mat'], 'predictor');         

                predictions = random_forest_reg_predict(predictor, test_data, 0);
                save([results_dir 'results.mat'], 'predictions');
            end
        end
    end
    
    if ismember(7, args.do_tests) && i_decomp == 2
        for output_type = output_types

            switch output_type{1}                       
                case 'width'
                    rf_args.sampling_args.y = true_widths;    
                case 'orientation'
                    rf_args.sampling_args.y = true_oris;
            end

            %G" + H"
            h2_data = u_load([exp_dir 'training/' num2str(repeat) '/responses_h2da.mat']);
            g2_data = u_load([exp_dir 'training/' num2str(repeat) '/responses_g2da.mat']);
            training_data = ...
                convert_complex_representation_back([g2_data h2_data], 'real_imag'); 
            clear h2_data g2_data;

            h2_data = u_load([exp_dir 'test/' num2str(repeat) '/responses_h2da.mat']);
            g2_data = u_load([exp_dir 'test/' num2str(repeat) '/responses_g2da.mat']);
            test_data = ...
                convert_complex_representation_back([g2_data h2_data], 'real_imag'); 
            clear h2_data g2_data;
            
            for feature_type = {'all'} %'phase', 'mag', 'conj', 'real_imag'}

                rf_dir = [exp_dir '/rfs/' num2str(repeat) '/' output_type{1} '/gh2d/feature_types/' feature_type{1} '/'];
                results_dir = [exp_dir '/results/' num2str(repeat) '/' output_type{1} '/gh2d/feature_types/' feature_type{1} '/'];
                mkdir(results_dir);

                %Args to reform data
                rf_args.sampling_args.X = ...
                    convert_complex_representation(training_data, feature_type{1}, 3);
                test_data_i = ...
                    convert_complex_representation(test_data, feature_type{1}, 3);

                rf_args.tree_dir = [rf_dir 'trees/'];
                predictor = random_forest_reg_train(rf_args);    
                if isfield(predictor, 'D')             
                    save([rf_dir 'predictor.mat'], 'predictor');         

                    predictions = random_forest_reg_predict(predictor, test_data_i, 0);
                    save([results_dir 'results.mat'], 'predictions');
                end
            end
        end
    end
end
%%
%--------------------------------------------------------------------------
%**************************************************************************
%--------------------------------------------------------------------------
