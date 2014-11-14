function comparing_representations_experiment_csf(i_decomp, repeat, varargin)

args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'noise_level',  unixenv('NOISE_LEVEL', 2), ...
    'num_ims',      unixenv('NUM_IMS', 2), ...
    'num_trees',    unixenv('NUM_TREES', 100), ...
    'pts_per_im',   unixenv('PTS_PER_IM', 10), ...
    'dx',           unixenv('DX',64),...
    'dy',           unixenv('DY',64),...
    'cx',           unixenv('CX',32),...
    'cy',           unixenv('CY',32),...
    'make_data',    unixenv('MAKE_DATA',1), ...
    'add_edge',    unixenv('ADD_EDGE',1), ...
    'add_circles',    unixenv('ADD_CIRCLES',1), ...
    'num_circles', unixenv('NUM_CIRCLES',32), ...
    'do_orientation', unixenv('DO_ORIENTATION',1), ...
    'do_detection', unixenv('DO_DETECTION',1), ...
    'do_width', unixenv('DO_WIDTH',1), ...
    'do_tests', unixenv('DO_TESTS',1), ...
    'win_sizes', unixenv('WIN_SIZES',[1 3]) ...
);

noise_level = args.noise_level;
num_ims = args.num_ims;
pts_per_im = args.pts_per_im;
win_sizes = args.win_sizes;
dx = args.dx;
dy = args.dy;
cx = args.cx;
cy = args.cy;

if args.add_edge
    exp_dir = [asymmetryroot 'experiments/synthetic_lines/line_on_edge/'];
elseif args.add_circles
    exp_dir = [asymmetryroot 'experiments/synthetic_lines/circles/'];
else
    exp_dir = [asymmetryroot 'experiments/synthetic_lines/comparing_representations/'];
end

warning('off', 'ASYM:unexpectedArgument');
%
%Make a load of data for bothing training and testing purposes
% Set up arguments for each decomposition type
d_args{1}.decomp_type = {'dt'};
d_args{1}.levels = 1:4;
d_args{1}.feature_shape = 'rect';
d_args{1}.feature_type = 'complex';
d_args{1}.do_max = 0;
d_args{1}.rotate = 0;
d_args{1}.use_nag = 0;

d_args{2}.decomp_type = {'linop'};
d_args{2}.num_levels = 4;
d_args{2}.num_angles = 6;
d_args{2}.do_max = 0;
d_args{2}.rotate = 0;       

d_args{3}.decomp_type = {'gabor'};
d_args{3}.num_angles = 6;
d_args{3}.sigma_range = [1 2 4 8];	
d_args{3}.do_max = 0;
d_args{3}.rotate = 0;
d_args{3}.feature_type = 'complex';

d_args{4}.decomp_type = {'mono'};
d_args{4}.num_levels = 4;
d_args{4}.min_wavelength = 4;
d_args{4}.onf = 0.65;

d_args{5}.decomp_type = {'g1d'};
d_args{5}.sigma_range = [1 2 4 8];
            
d_args{6}.decomp_type = {'g2d'};
d_args{6}.sigma_range = [1 2 4 8];
            
d_args{7}.decomp_type = {'g2di'};
d_args{7}.sigma_range = [1 4];
            
d_args{8}.decomp_type = {'h2d'};
d_args{8}.sigma_range = [1 2 4 8];

d_args{9}.decomp_type = {'g2da'};
d_args{9}.sigma_range = [1 2 4 8];
d_args{9}.num_angles = 6;
d_args{9}.do_max = 0;
d_args{9}.rotate = 0;

d_args{10}.decomp_type = {'gabori'};
d_args{10}.sigma_range = [1 4];
d_args{10}.num_angles = 6;
d_args{10}.do_max = 0;
d_args{10}.feature_type = 'complex';

num_decomps = 10;

%set common parameters for all decomp types and then compute vector sizes
d_args{i_decomp}.win_size = 3;
d_args{i_decomp}.normalise = 0;
d_args{i_decomp}.pca = [];

D = zeros(num_decomps,1);
D(i_decomp) = get_samples_per_channel(d_args{i_decomp});

%--------------------------------------------------------------------------
% ***************** GENERATE DATA *****************************************
%--------------------------------------------------------------------------

if args.make_data
    for data_type = {'test'}%{'training', 'test'}

        mkdir([exp_dir data_type{1} '/rician_' num2str(noise_level) '/' num2str(repeat)]);

        %Pre-allocate space for the true line parameters
        true_oris = zeros(num_ims*pts_per_im,1);
        true_cons = zeros(num_ims*pts_per_im,1);
        true_widths = zeros(num_ims*pts_per_im,1);
        centre_idx = false(num_ims*pts_per_im,1);
        edge_idx = [];
        circle_idx = [];
        if args.add_edge
            edge_idx = false(num_ims*pts_per_im,1);         
        end
        if args.add_circles
            circle_idx = false(num_ims*pts_per_im,1);         
        end

        %Pre-allocate space for the responses
        eval(['responses_' d_args{i_decomp}.decomp_type{1} ' = zeros(num_ims*pts_per_im, D(i_decomp));']);
        eval(['bg_responses_' d_args{i_decomp}.decomp_type{1} ' = zeros(num_ims*pts_per_im, D(i_decomp));']);

        if strcmpi(data_type{1}, 'training')
            %rng(1000, 'twister');
            rand('twister', 1000 + repeat);
        else
            %rng(2000, 'twister');
            rand('twister', 2000 + repeat);
        end

        for ii = 1:num_ims

            display(['Testing image ' num2str(ii)]);

            %Sample properties of line
            line_width = sample_uniform([1 8]);
            line_contrast = sample_uniform([1 8]);
            line_ori = sample_uniform([0 360]);
            line_rad = pi * line_ori / 180;

            %Generate line
            [line, label, label_centre] =...
                create_ellipse_bar(line_width/2, line_contrast, line_ori, dy, dx, cx, cy);
            
            if args.add_edge
                edge_width = sample_uniform([1 8]);
                edge_contrast = sample_uniform([4 8]);
                edge_ori = sample_from_normal(line_ori + 90, 30^2, 1);
                
                %Generate edge
                [edge, edge_label] = create_sin_step(edge_width/2, edge_contrast, edge_ori, dy, dx, cx, cy);
                line = line + edge;
            end
            if args.add_circles
                circle_label = false(dy,dx);
                shuffle = randperm(dx*dy);
                r_idx = shuffle(1:args.num_circles);
                [r_rows r_cols] = ind2sub([dy dx], r_idx);
                radii = sample_uniform([1 2], [args.num_circles 1]);
                circ_contrast = sample_uniform([1 8], [args.num_circles 1]);
                for i_circ = 1:args.num_circles
                    [blob blob_label] = create_gauss_blob(...
                        radii(i_circ), circ_contrast(i_circ), dy, dx, r_cols(i_circ), r_rows(i_circ));
                    circle_label(blob_label) = 1;
                    line = line + blob;
                end
                
            end
            test_image = ricernd(1 + line, noise_level);
            
            if ispc && ii < 5
                figure; 
                subplot(1,2,1); imgray(test_image);
                subplot(1,2,2); imgray(label);
            end

            %For ecah decomposition compute the responses
            im_responses = compute_filter_responses(test_image, d_args{i_decomp});

            fov_mask = true(size(label));
            fov_mask([1:16 end-15:end],:) = 0;
            fov_mask(:, [1:16 end-15:end]) = 0;

            %Select some random pixels in the signal image
            shuffle = randperm(sum(label(:) & fov_mask(:)));
            idx = find(label & fov_mask);
            r_idx = idx(shuffle(1:pts_per_im));
            [r_rows r_cols] = ind2sub([dy dx], r_idx);
            sample_idx = ((ii-1)*pts_per_im + 1):(ii*pts_per_im);

            %Save the line parameters
            true_oris(sample_idx,:) = line_rad;
            true_cons(sample_idx,:) = line_contrast;
            true_widths(sample_idx,:) = line_width;
            centre_idx(sample_idx,:) = label_centre(r_idx);

            eval(['responses_' d_args{i_decomp}.decomp_type{1} '(sample_idx, :) = '...
                'sample_image_features(im_responses, r_rows(:), r_cols(:), d_args{i_decomp});']);

            %Select some random bg pixels in the signal image
            shuffle = randperm(sum(~label(:) & fov_mask(:)));
            idx = find(~label & fov_mask);
            r_idx = idx(shuffle(1:pts_per_im));
            [r_rows r_cols] = ind2sub([dy dx], r_idx);
            sample_idx = ((ii-1)*pts_per_im + 1):(ii*pts_per_im);
            if args.add_edge
                edge_idx(sample_idx,:) = edge_label(r_idx);
            end
            if args.add_circles
                circle_idx(sample_idx,:) = circle_label(r_idx);
            end

            eval(['bg_responses_' d_args{i_decomp}.decomp_type{1} '(sample_idx, :) = '...
                'sample_image_features(im_responses, r_rows(:), r_cols(:), d_args{i_decomp});']);
        end
        save([exp_dir data_type{1} '/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat'],...
            ['responses_' d_args{i_decomp}.decomp_type{1}]);
        save([exp_dir data_type{1} '/rician_' num2str(noise_level) '/' num2str(repeat) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat'],...
                ['bg_responses_' d_args{i_decomp}.decomp_type{1}]);
        save([exp_dir data_type{1} '/rician_' num2str(noise_level) '/' num2str(repeat) '/true_labels.mat'], 'true_*',...
            'centre_idx', 'edge_idx', 'circle_idx');

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
    rf_args.overwrite = 1;
    rf_args.minimise_size = 0;
    rf_args.split_criterion = 'gdi';
    rf_args.var_criterion = 'mabs';

    rf_args.sampling_args.sampling_method = 'sample_saved_training_data';
    rf_args.decomposition_args = [];
    
    %
    % 1) Original data
    %
    if ismember(0, args.do_tests) && ismember(i_decomp, 1:10);

        for win_size = win_sizes

            fg_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
            bg_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(repeat) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
            test_data = [fg_data; bg_data]; 
            clear fg_data bg_data;
            num_pts = num_ims*pts_per_im;
            orig_data = zeros(20*num_pts, size(test_data,2));
            
            for i_repeat = 0:9 
                %Load data for this decomp type
                fg_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/'...
                    num2str(repeat+i_repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
                bg_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/'...
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
                rf_dir = [exp_dir '/rfs/rician_' num2str(noise_level) '/' num2str(repeat)...
                    '/detection/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/' num2str(i_pts) '/'];
                rf_args.tree_dir = [rf_dir 'trees/'];
                results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(repeat)...
                    '/detection/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/' num2str(i_pts) '/'];
                mkdir(results_dir);

                %Train
                rf_args.sampling_args.y = [true(i_pts,1); false(i_pts,1)];
                rf_args.sampling_args.X = orig_data([1:i_pts 10*num_pts+(1:i_pts)],:);
                predictor = random_forest_class_train(rf_args);
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
    rf_args.sampling_args.y = [true(num_ims*pts_per_im,1); false(num_ims*pts_per_im,1)];
    %
    if ismember(1, args.do_tests) && ismember(i_decomp, 1:10);

        for win_size = win_sizes

            %Load data for this decomp type
            fg_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
            bg_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/' num2str(repeat) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
            orig_data = [fg_data; bg_data]; 
            clear fg_data bg_data;
            fg_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
            bg_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(repeat) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
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
            rf_dir = [exp_dir '/rfs/rician_' num2str(noise_level) '/' num2str(repeat)...
                '/detection/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/'];
            rf_args.tree_dir = [rf_dir 'trees/'];
            results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(repeat)...
                '/detection/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/'];
            mkdir(results_dir);

            %Train
            rf_args.sampling_args.X = orig_data;
            predictor = random_forest_class_train(rf_args);
            save([rf_dir 'predictor.mat'], 'predictor');
            display(['************ FOREST for ' d_args{i_decomp}.decomp_type{1} ' complete!! **************']);

            %Test
            [~, votes, all_votes] = random_forest_class_predict(predictor, test_data);
            predicted_lines = votes(:,2) / rf_args.n_trees;
            save([results_dir 'all_votes.mat'], 'all_votes');
            save([results_dir 'results1.mat'], 'predicted_lines');
        end
    end
    
    %
    % 2) Complex representations
    %
    if ismember(2, args.do_tests) && ismember(i_decomp, [1 3]);

        fg_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
        bg_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(repeat) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
        test_data = [fg_data; bg_data]; 
        clear fg_data bg_data;
        rf_args.sampling_args.y = [true(4*num_ims*pts_per_im,1); false(4*num_ims*pts_per_im,1)];
        
        num_pts = num_ims*pts_per_im;
        orig_data = zeros(8*num_pts, size(test_data,2));

        for i_repeat = 0:3 
            %Load data for this decomp type
            fg_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/'...
                num2str(repeat+i_repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
            bg_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/'...
                num2str(repeat+i_repeat) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
            fg_rows = (1:num_pts) + i_repeat*num_pts;
            bg_rows = fg_rows + 4*num_pts;
            orig_data(fg_rows,:) = fg_data;
            orig_data(bg_rows,:) = bg_data; 
            clear fg_data bg_data;
        end

        % 1) Complex forms: all (mag/phase), real/imag, conj, mag, phase
        for feature_type = {'all' 'real_imag', 'conj', 'mag', 'phase', 'real', 'imag', 'real_abs_imag'}

            rf_dir = [exp_dir '/rfs/rician_' num2str(noise_level) '/' num2str(repeat) '/detection/' d_args{i_decomp}.decomp_type{1} '/feature_types/' feature_type{1} '/'];
            results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(repeat) '/detection/' d_args{i_decomp}.decomp_type{1} '/feature_types/' feature_type{1} '/'];
            mkdir(results_dir);

            %Args to reform data
            new_decomp_args = [];
            new_decomp_args.feature_type = feature_type{1};

            if isfield(d_args{i_decomp}, 'num_angles') && (d_args{i_decomp}.num_angles == 18)
                new_decomp_args.bands = 1:3:16;
            end

            %Train
            rf_args.sampling_args.X = convert_decomp_form(orig_data, d_args{i_decomp}, new_decomp_args);
            rf_args.tree_dir = [rf_dir 'trees/'];
            predictor = random_forest_class_train(rf_args);    
            save([rf_dir 'predictor.mat'], 'predictor');

            %Test
            test_data_i = convert_decomp_form(test_data, d_args{i_decomp}, new_decomp_args);
            [dummy, votes] = random_forest_class_predict(predictor, test_data_i);
            predicted_lines = votes(:,2) / rf_args.n_trees;
            save([results_dir 'results.mat'], 'predicted_lines');

        end
    end

    % 
    % 3) Different number of levels
    %
    if ismember(3, args.do_tests) && ismember(i_decomp, [1 3 5 6]);

        %Load data for this decomp type
        fg_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
        bg_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/' num2str(repeat) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
        orig_data = [fg_data; bg_data]; 
        clear fg_data bg_data;
        fg_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
        bg_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(repeat) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
        test_data = [fg_data; bg_data]; 
        clear fg_data bg_data;
                
        for i_level = 1:4

            for win_size = win_sizes

                rf_dir = [exp_dir '/rfs/rician_' num2str(noise_level) '/' num2str(repeat) '/detection/' d_args{i_decomp}.decomp_type{1} '/levels/' num2str(i_level) '/' num2str(win_size) '/'];
                results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(repeat) '/detection/' d_args{i_decomp}.decomp_type{1} '/levels/' num2str(i_level) '/' num2str(win_size) '/'];
                mkdir(results_dir);

                %Args to reform data
                new_decomp_args = [];
                new_decomp_args.feature_type = 'conj';
                new_decomp_args.levels = i_level;
                new_decomp_args.win_size = win_size;
                if isfield(d_args{i_decomp}, 'num_angles') && (d_args{i_decomp}.num_angles == 18)
                    new_decomp_args.bands = 1:3:16;
                end

                %Train       
                rf_args.sampling_args.X = convert_decomp_form(orig_data, d_args{i_decomp}, new_decomp_args);
                rf_args.tree_dir = [rf_dir 'trees/'];
                predictor = random_forest_class_train(rf_args);    
                save([rf_dir 'predictor.mat'], 'predictor');

                %Test
                test_data_i = convert_decomp_form(test_data, d_args{i_decomp}, new_decomp_args);
                [dummy, votes] = random_forest_class_predict(predictor, test_data_i);
                predicted_lines = votes(:,2) / rf_args.n_trees;
                save([results_dir 'results.mat'], 'predicted_lines');
            end

        end
    end
    %
    % 4) Rotate and do max
    %
    if ismember(4, args.do_tests) && ismember(i_decomp, [1 3 9])

        fg_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
        bg_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(repeat) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
        test_data = [fg_data; bg_data]; 
        clear fg_data bg_data;
        rf_args.sampling_args.y = [true(4*num_ims*pts_per_im,1); false(4*num_ims*pts_per_im,1)];
        
        num_pts = num_ims*pts_per_im;
        orig_data = zeros(8*num_pts, size(test_data,2));

        for i_repeat = 0:3 
            %Load data for this decomp type
            fg_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/'...
                num2str(repeat+i_repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
            bg_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/'...
                num2str(repeat+i_repeat) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
            fg_rows = (1:num_pts) + i_repeat*num_pts;
            bg_rows = fg_rows + 4*num_pts;
            orig_data(fg_rows,:) = fg_data;
            orig_data(bg_rows,:) = bg_data; 
            clear fg_data bg_data;
        end
    
        for i_level = 5

            for win_size = win_sizes

                for reform_type = 1:2

                    if win_size==1 && reform_type==2
                        continue;
                    end

                    new_decomp_args = [];
                    if reform_type == 1
                        rf_dir = [exp_dir '/rfs/rician_' num2str(noise_level) '/' num2str(repeat) '/detection/' d_args{i_decomp}.decomp_type{1} '/rotate/L' num2str(i_level) '_W' num2str(win_size)  '/'];
                        results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(repeat) '/detection/' d_args{i_decomp}.decomp_type{1} '/rotate/L' num2str(i_level) '_W' num2str(win_size) '/'];
                        new_decomp_args.rotate = 1;
                    else
                        rf_dir = [exp_dir '/rfs/rician_' num2str(noise_level) '/' num2str(repeat) '/detection/' d_args{i_decomp}.decomp_type{1} '/do_max/L' num2str(i_level) '_W' num2str(win_size)  '/'];
                        results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(repeat) '/detection/' d_args{i_decomp}.decomp_type{1} '/do_max/L' num2str(i_level) '_W' num2str(win_size) '/'];
                        new_decomp_args.do_max = 1;
                    end
                    mkdir(results_dir);

                    %Args to reform data
                    new_decomp_args.feature_type = 'conj';
                    if i_level < 5
                        new_decomp_args.levels = i_level;
                    else
                        new_decomp_args.levels = 1:4;
                    end

                    new_decomp_args.win_size = win_size;
                    if isfield(d_args{i_decomp}, 'num_angles') && (d_args{i_decomp}.num_angles == 18)
                        new_decomp_args.bands = 1:3:16;
                    end

                    %Train       
                    rf_args.sampling_args.X = convert_decomp_form(orig_data, d_args{i_decomp}, new_decomp_args);
                    rf_args.tree_dir = [rf_dir 'trees/'];
                    predictor = random_forest_class_train(rf_args);    
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

        g1_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_g1d.mat']);
        g2_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_g2d.mat']);
        g1_data_bg = u_load([exp_dir 'training/rician_' num2str(noise_level) '/' num2str(repeat) '/bg_responses_g1d.mat']);
        g2_data_bg = u_load([exp_dir 'training/rician_' num2str(noise_level) '/' num2str(repeat) '/bg_responses_g2d.mat']);
        rf_args.sampling_args.X = [g1_data g2_data; g1_data_bg g2_data_bg]; 
        clear g1_data g2_data;

        g1_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_g1d.mat']);
        g2_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_g2d.mat']);
        g1_data_bg = u_load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(repeat) '/bg_responses_g1d.mat']);
        g2_data_bg = u_load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(repeat) '/bg_responses_g2d.mat']);
        test_data = [g1_data g2_data; g1_data_bg g2_data_bg]; 
        clear g1_data g2_data;

        rf_dir = [exp_dir '/rfs/rician_' num2str(noise_level) '/' num2str(repeat) '/detection/g12d/'];
        results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(repeat) '/detection/g12d/'];
        mkdir(results_dir);

        rf_args.tree_dir = [rf_dir 'trees/'];
        predictor = random_forest_class_train(rf_args);    
        save([rf_dir 'predictor.mat'], 'predictor');

        [dummy, votes] = random_forest_class_predict(predictor, test_data);
        predicted_lines = votes(:,2) / rf_args.n_trees;
        save([results_dir 'results.mat'], 'predicted_lines');

        %G" + H"
        h2_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_h2d.mat']);
        g2_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_g2d.mat']);
        h2_data_bg = u_load([exp_dir 'training/rician_' num2str(noise_level) '/' num2str(repeat) '/bg_responses_h2d.mat']);
        g2_data_bg = u_load([exp_dir 'training/rician_' num2str(noise_level) '/' num2str(repeat) '/bg_responses_g2d.mat']);
        rf_args.sampling_args.X = [h2_data g2_data; h2_data_bg g2_data_bg];
        clear h2_data g2_data;

        h2_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_h2d.mat']);
        g2_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_g2d.mat']);
        h2_data_bg = u_load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(repeat) '/bg_responses_h2d.mat']);
        g2_data_bg = u_load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(repeat) '/bg_responses_g2d.mat']);
        test_data = [h2_data g2_data; h2_data_bg g2_data_bg];
        clear h2_data g2_data h2_data_bg g2_data_bg;

        rf_dir = [exp_dir '/rfs/rician_' num2str(noise_level) '/' num2str(repeat) '/detection/gh2d/'];
        results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(repeat) '/detection/gh2d/'];
        mkdir(results_dir);

        rf_args.tree_dir = [rf_dir 'trees/'];
        predictor = random_forest_class_train(rf_args);    
        save([rf_dir 'predictor.mat'], 'predictor');

        [dummy, votes] = random_forest_class_predict(predictor, test_data);
        predicted_lines = votes(:,2) / rf_args.n_trees;
        save([results_dir 'results.mat'], 'predicted_lines');
    end
% %
% % 6) Different number of angles for Gabor filters
% %
% if i_decomp == 3
%       
%     fg_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
%     bg_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/' num2str(repeat) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
%     orig_data = [fg_data; bg_data]; 
%     clear fg_data bg_data;
%     fg_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
%     bg_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(repeat) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
%     test_data = [fg_data; bg_data]; 
%     clear fg_data bg_data;
%     
%     for num_angles = [3 9 18]
%         
%         rf_dir = [exp_dir '/rfs/rician_' num2str(noise_level) '/' num2str(repeat) '/detection/' d_args{i_decomp}.decomp_type{1} '/num_angles/' num2str(num_angles) '/'];
%         results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(repeat) '/detection/' d_args{i_decomp}.decomp_type{1} '/num_angles/' num2str(num_angles) '/'];
%         mkdir(results_dir);
% 
%         %Args to reform data
%         new_decomp_args = [];
%         new_decomp_args.feature_type = 'conj';
%         band_step = 18 / num_angles;
%         new_decomp_args.bands = 1:band_step:18;
% 
%         %Train       
%         rf_args.sampling_args.X = convert_decomp_form(orig_data, d_args{i_decomp}, new_decomp_args);
%         rf_args.tree_dir = [rf_dir 'trees/'];
%         predictor = random_forest_class_train(rf_args);    
%         save([rf_dir 'predictor.mat'], 'predictor');
% 
%         %Test
%         test_data_i = convert_decomp_form(test_data, d_args{i_decomp}, new_decomp_args);
%         [dummy, votes] = random_forest_class_predict(predictor, test_data_i);
%         predicted_lines = votes(:,2) / rf_args.n_trees;
%         save([results_dir 'results.mat'], 'predicted_lines');
% 
%     end
% end
end

%--------------------------------------------------------------------------
% *********************** ORIENTATION + WIDTH *****************************
%--------------------------------------------------------------------------
%
output_types = cell(1,0);

if args.do_orientation
    output_types{1,end+1} = 'orientation';
end
if args.do_width
    output_types{1,end+1} = 'width';
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
    rf_args.overwrite = 1;
    rf_args.minimise_size = 0;
    rf_args.split_criterion = 'ssq';
    rf_args.var_criterion = 'ssq';

    rf_args.sampling_args.sampling_method = 'sample_saved_training_data';
    rf_args.decomposition_args = [];
    %
    % 0) Original data with increasing num points
    %
    if ismember(0, args.do_tests) && ismember(i_decomp, 1:10);
        for output_type = output_types
            for win_size = win_sizes
                
                %Load data for this decomp type
                test_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);

                num_pts = num_ims*pts_per_im;
                orig_data = zeros(4*num_pts, size(test_data,2));
                orig_y = zeros(4*num_pts, 1);
                for i_repeat = 0:3 
                    %Load data for this decomp type
                    fg_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/'...
                        num2str(repeat+i_repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
                    fg_rows = (1:num_pts) + i_repeat*num_pts;
                    orig_data(fg_rows,:) = fg_data;
                    clear fg_data;
                    
                    switch output_type{1}                       
                        case 'width'
                            load([exp_dir 'training/rician_' num2str(noise_level) '/' num2str(repeat) '/true_labels.mat'], 'true_widths');
                            orig_y(fg_rows,:) = true_widths;    
                        case 'orientation'
                            load([exp_dir 'training/rician_' num2str(noise_level) '/' num2str(repeat) '/true_labels.mat'], 'true_oris');
                            orig_y(fg_rows,:) = complex(cos(2*true_oris), sin(2*true_oris));
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
                    rf_dir = [exp_dir '/rfs/rician_' num2str(noise_level) '/' num2str(repeat)...
                        '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/' num2str(i_pts) '/'];
                    rf_args.tree_dir = [rf_dir 'trees/'];
                    results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(repeat)...
                        '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/' num2str(i_pts) '/'];
                    mkdir(results_dir);

                    %Train
                    rf_args.sampling_args.y = orig_y(1:i_pts);
                    rf_args.sampling_args.X = orig_data(1:i_pts,:);
                    predictor = random_forest_reg_train(rf_args);
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
    % 1) Original data
    %
    if ismember(1, args.do_tests) && ismember(i_decomp, 1:10);
        for output_type = output_types

            switch output_type{1}                       
                case 'width'
                    orig_y = true_widths;    
                case 'orientation'
                    orig_y = complex(cos(2*true_oris), sin(2*true_oris));
            end

            for win_size = win_sizes

                %Load data for this decomp type
                orig_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
                test_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);

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
                rf_dir = [exp_dir '/rfs/rician_' num2str(noise_level) '/' num2str(repeat)...
                    '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/'];
                rf_args.tree_dir = [rf_dir 'trees/'];
                results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(repeat)...
                    '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/'];
                mkdir(results_dir);

                %Train
                rf_args.sampling_args.y = orig_y;
                rf_args.sampling_args.X = orig_data;
                predictor = random_forest_reg_train(rf_args);
                save([rf_dir 'predictor.mat'], 'predictor');
                display(['************ FOREST for ' d_args{i_decomp}.decomp_type{1} ' complete!! **************']);

                %Test
                predictions = random_forest_reg_predict(predictor, test_data);
                save([results_dir 'results.mat'], 'predictions');
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
                    rf_args.sampling_args.y = complex(cos(2*true_oris), sin(2*true_oris));
            end

            %Load data for this decomp type
            orig_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);

            %Load data for this decomp type
            test_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);

            % 1) Complex forms: all (mag/phase), real/imag, conj, mag, phase
            for feature_type = {'all' 'real_imag', 'conj', 'mag', 'phase', 'imag', 'real', 'real_abs_imag'}

                rf_dir = [exp_dir '/rfs/rician_' num2str(noise_level) '/' num2str(repeat) '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/feature_types/' feature_type{1} '/'];
                results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(repeat) '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/feature_types/' feature_type{1} '/'];
                mkdir(results_dir);

                %Args to reform data
                new_decomp_args = [];
                new_decomp_args.feature_type = feature_type{1};

                if isfield(d_args{i_decomp}, 'num_angles') && (d_args{i_decomp}.num_angles == 18)
                    new_decomp_args.bands = 1:3:16;
                end

                %Train
                rf_args.sampling_args.X = convert_decomp_form(orig_data, d_args{i_decomp}, new_decomp_args);
                rf_args.tree_dir = [rf_dir 'trees/'];
                predictor = random_forest_reg_train(rf_args);    
                save([rf_dir 'predictor.mat'], 'predictor');

                %Test
                test_data_i = convert_decomp_form(test_data, d_args{i_decomp}, new_decomp_args);
                predictions = random_forest_reg_predict(predictor, test_data_i);
                save([results_dir 'results.mat'], 'predictions');
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
                    rf_args.sampling_args.y = complex(cos(2*true_oris), sin(2*true_oris));
            end
            
            %Load data for this decomp type
            orig_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);

            %Load data for this decomp type
            test_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);

            for win_size = win_sizes

                for i_level = 1:4

                    rf_dir = [exp_dir '/rfs/rician_' num2str(noise_level) '/' num2str(repeat) '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/levels/' num2str(i_level) '/' num2str(win_size) '/'];
                    results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(repeat) '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/levels/' num2str(i_level) '/' num2str(win_size) '/'];
                    mkdir(results_dir);


                    %Args to reform data
                    new_decomp_args = [];
                    new_decomp_args.feature_type = 'conj';
                    new_decomp_args.levels = i_level;
                    new_decomp_args.win_size = win_size;
                    if isfield(d_args{i_decomp}, 'num_angles') && (d_args{i_decomp}.num_angles == 18)
                        new_decomp_args.bands = 1:3:16;
                    end

                    %Train       
                    rf_args.sampling_args.X = convert_decomp_form(orig_data, d_args{i_decomp}, new_decomp_args);
                    rf_args.tree_dir = [rf_dir 'trees/'];
                    predictor = random_forest_reg_train(rf_args);    
                    save([rf_dir 'predictor.mat'], 'predictor');

                    %Test
                    test_data_i = convert_decomp_form(test_data, d_args{i_decomp}, new_decomp_args);
                    predictions = random_forest_reg_predict(predictor, test_data_i);
                    save([results_dir 'results.mat'], 'predictions');
                end

            end
        end
    end
    % 
    % 4) Different ways of stitching Gaussian filters together
    % G" + G'
    if ismember(4, args.do_tests) && i_decomp == 5
        for output_type = output_types

            switch output_type{1}                       
                case 'width'
                    rf_args.sampling_args.y = true_widths;    
                case 'orientation'
                    rf_args.sampling_args.y = complex(cos(2*true_oris), sin(2*true_oris));
            end

            g1_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_g1d.mat']);
            g2_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_g2d.mat']);
            rf_args.sampling_args.X = [g1_data g2_data]; clear g1_data g2_data;

            g1_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_g1d.mat']);
            g2_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_g2d.mat']);
            test_data = [g1_data g2_data]; clear g1_data g2_data;

            rf_dir = [exp_dir '/rfs/rician_' num2str(noise_level) '/' num2str(repeat) '/' output_type{1} '/g12d/'];
            results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(repeat) '/' output_type{1} '/g12d/'];
            mkdir(results_dir);

            rf_args.tree_dir = [rf_dir 'trees/'];
            predictor = random_forest_reg_train(rf_args);    
            save([rf_dir 'predictor.mat'], 'predictor');

            predictions = random_forest_reg_predict(predictor, test_data, 0);
            save([results_dir 'results.mat'], 'predictions');

            %G" + H"
            h2_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_h2d.mat']);
            g2_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_g2d.mat']);
            rf_args.sampling_args.X = [h2_data g2_data]; clear h2_data g2_data;

            h2_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_h2d.mat']);
            g2_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_g2d.mat']);
            test_data = [h2_data g2_data]; clear h2_data g2_data;

            rf_dir = [exp_dir '/rfs/rician_' num2str(noise_level) '/' num2str(repeat) '/' output_type{1} '/gh2d/'];
            results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(repeat) '/' output_type{1} '/gh2d/'];
            mkdir(results_dir);

            rf_args.tree_dir = [rf_dir 'trees/'];
            predictor = random_forest_reg_train(rf_args);    
            save([rf_dir 'predictor.mat'], 'predictor');

            predictions = random_forest_reg_predict(predictor, test_data, 0);
            save([results_dir 'results.mat'], 'predictions');
        end
    end

% % 
% % 6) Different number of angles for Gabor filters
% %
% for output_type = {'orientation'}%, 'width'}
%     
%     switch output_type{1}                       
%         case 'width'
%             rf_args.sampling_args.y = true_widths;    
%         case 'orientation'
%             rf_args.sampling_args.y = complex(cos(2*true_oris), sin(2*true_oris));
%     end
%     if ismember(i_decomp, 3);
% 
%         Load data for this decomp type
%         orig_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
% 
%         Load data for this decomp type
%         test_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
% 
% 
%         for num_angles = 18%[3 9 18]
%         
%             rf_dir = [exp_dir '/rfs/rician_' num2str(noise_level) '/' num2str(repeat) '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/num_angles/' num2str(num_angles) '/'];
%             results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(repeat) '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/num_angles/' num2str(num_angles) '/'];
%             mkdir(results_dir);
% 
% 
%             Args to reform data
%             new_decomp_args = [];
%             new_decomp_args.feature_type = 'conj';
%             band_step = 18 / num_angles;
%             new_decomp_args.bands = 1:band_step:18;
% 
%             Train       
%             rf_args.sampling_args.X = convert_decomp_form(orig_data, d_args{i_decomp}, new_decomp_args);
%             rf_args.tree_dir = [rf_dir 'trees/'];
%             predictor = random_forest_reg_train(rf_args);    
%             save([rf_dir 'predictor.mat'], 'predictor');
% 
%             Test
%             test_data_i = convert_decomp_form(test_data, d_args{i_decomp}, new_decomp_args);
%             predictions = random_forest_reg_predict(predictor, test_data_i);
%             save([results_dir 'results.mat'], 'predictions');
% 
%         end
%     end
% end
end
%%
%--------------------------------------------------------------------------
%**************************************************************************
%--------------------------------------------------------------------------
