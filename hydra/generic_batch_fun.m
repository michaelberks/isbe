function generic_batch_fun(idx, varargin)

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
d_args{i_decomp}.win_size = 3;
d_args{i_decomp}.normalise = 0;
d_args{i_decomp}.pca = [];
%%
% --------------------------------------------------------------------------
% ********************** DETECTION ****************************************
% --------------------------------------------------------------------------
if args.do_detection
    %
    if ismember(1, args.do_tests) && ismember(i_decomp, 1:11);

        for win_size = win_sizes
            
            %Set up directory to save RF to
            rf_dir = [exp_dir '/rfs/' num2str(repeat)...
                '/detection/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/'];
            results_dir = [exp_dir '/results/' num2str(repeat)...
                '/detection/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/'];
            
            dir(rf_dir);
            dir(results_dir);
            
            if exist([rf_dir 'predictor.mat'], 'file') && ~exist([results_dir 'results1.mat'], 'file')
                display('Re-testing');
                fg_data = u_load([exp_dir 'test/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
                bg_data = u_load([exp_dir 'test/' num2str(repeat) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
                test_data = [fg_data; bg_data]; 
                clear fg_data bg_data;
            
                load([rf_dir 'predictor.mat'], 'predictor');

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

        % 1) Complex forms: all (mag/phase), real/imag, conj, mag, phase
        for feature_type = {'real_imag', 'conj', 'mag', 'phase', 'real', 'imag', 'all'}

            rf_dir = [exp_dir '/rfs/' num2str(repeat) '/detection/' d_args{i_decomp}.decomp_type{1} '/feature_types/' feature_type{1} '/'];
            results_dir = [exp_dir '/results/' num2str(repeat) '/detection/' d_args{i_decomp}.decomp_type{1} '/feature_types/' feature_type{1} '/'];
            if exist([rf_dir 'predictor.mat'], 'file') && ~exist([results_dir 'results.mat'], 'file')
                fg_data = u_load([exp_dir 'test/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
                bg_data = u_load([exp_dir 'test/' num2str(repeat) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
                test_data = [fg_data; bg_data];
                
                load([rf_dir 'predictor.mat'], 'predictor');
                
                %Args to reform data
                new_decomp_args = [];
                new_decomp_args.feature_type = feature_type{1};

                %Test
                test_data_i = convert_decomp_form(test_data, d_args{i_decomp}, new_decomp_args);
                [dummy, votes] = random_forest_class_predict(predictor, test_data_i);
                predicted_lines = votes(:,2) / rf_args.n_trees;
                save([results_dir 'results.mat'], 'predicted_lines');
            end
        end
    end

    % 
    % 3) Different number of levels
    %
    if ismember(3, args.do_tests) && ismember(i_decomp, [1 3 5 6]);                
        for i_level = 1:5

            for win_size = win_sizes

                rf_dir = [exp_dir '/rfs/' num2str(repeat) '/detection/' d_args{i_decomp}.decomp_type{1} '/levels/' num2str(i_level) '/' num2str(win_size) '/'];
                results_dir = [exp_dir '/results/' num2str(repeat) '/detection/' d_args{i_decomp}.decomp_type{1} '/levels/' num2str(i_level) '/' num2str(win_size) '/'];
                
                if exist([rf_dir 'predictor.mat'], 'file') && ~exist([results_dir 'results.mat'], 'file')
                    fg_data = u_load([exp_dir 'test/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
                    bg_data = u_load([exp_dir 'test/' num2str(repeat) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
                    test_data = [fg_data; bg_data]; 
                    load([rf_dir 'predictor.mat'], 'predictor');

                    %Args to reform data
                    new_decomp_args = [];
                    new_decomp_args.feature_type = 'conj';
                    new_decomp_args.levels = i_level;
                    new_decomp_args.win_size = win_size;
                    %Test
                    test_data_i = convert_decomp_form(test_data, d_args{i_decomp}, new_decomp_args);
                    [dummy, votes] = random_forest_class_predict(predictor, test_data_i);
                    predicted_lines = votes(:,2) / rf_args.n_trees;
                    save([results_dir 'results.mat'], 'predicted_lines');
                end
            end
        end
    end
    
    if ismember(6, args.do_tests) && i_decomp == 7
        for win_size = win_sizes
            %Set up directory to save RF to
            rf_dir = [exp_dir '/rfs/' num2str(repeat) '/detection/all/' num2str(win_size) '/'];
            results_dir = [exp_dir '/results/' num2str(repeat) '/detection/all/' num2str(win_size) '/'];
            
            if exist([rf_dir 'predictor.mat'], 'file') && ~exist([results_dir 'results1.mat'], 'file')
                
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

                test_data = zeros(2*num_pts, sum(D_all));
                curr_col = 0;
                for i_d = all_decomps
                    cols = curr_col + (1:D_all(i_d));
                    curr_col = cols(end);

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
                    test_data(:, cols) =...
                        convert_decomp_form(test_data_i, d_args{i_d}, new_decomp_args);
                end

                %Test
                load([rf_dir 'predictor.mat'], 'predictor');
                [~, votes] = random_forest_class_predict(predictor, test_data);
                predicted_lines = votes(:,2) / rf_args.n_trees;
                save([results_dir 'results1.mat'], 'predicted_lines');
            end
        end
    end
    
    if ismember(7, args.do_tests) && i_decomp == 2

        for feature_type = {'phase', 'mag', 'conj', 'real_imag'}

            rf_dir = [exp_dir '/rfs/' num2str(repeat) '/detection/gh2d/feature_types/' feature_type{1} '/'];
            results_dir = [exp_dir '/results/' num2str(repeat) '/detection/gh2d/feature_types/' feature_type{1} '/'];
            
            dir(rf_dir);
            dir(results_dir);
            
            if exist([rf_dir 'predictor.mat'], 'file') && ~exist([results_dir 'results1.mat'], 'file')
                display('Re-testing');
                
                h2_data = u_load([exp_dir 'test/' num2str(repeat) '/responses_h2da.mat']);
                g2_data = u_load([exp_dir 'test/' num2str(repeat) '/responses_g2da.mat']);
                h2_data_bg = u_load([exp_dir 'test/' num2str(repeat) '/bg_responses_h2da.mat']);
                g2_data_bg = u_load([exp_dir 'test/' num2str(repeat) '/bg_responses_g2da.mat']);
                test_data = ...
                    convert_complex_representation_back([g2_data h2_data; g2_data_bg h2_data_bg], 'real_imag'); 
                clear h2_data g2_data g2_data_bg h2_data_bg;

                %Args to reform data
                test_data_i = ...
                    convert_complex_representation(test_data, feature_type{1}, 3);

                load([rf_dir 'predictor.mat'], 'predictor');

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

    %
    % 1) Original data
    %
    if ismember(1, args.do_tests) && ismember(i_decomp, 1:11);
        for output_type = output_types
            for win_size = win_sizes               
                %Set up directory to save RF to
                rf_dir = [exp_dir '/rfs/' num2str(repeat)...
                    '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/'];
                results_dir = [exp_dir '/results/' num2str(repeat)...
                    '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/'];
                
                display(rf_dir);
                dir(rf_dir);
                display(results_dir);
                dir(results_dir);
            
                if exist([rf_dir 'predictor.mat'], 'file') && ~exist([results_dir 'results.mat'], 'file')
                    display('Re-testing');
                    
                    test_data = u_load([exp_dir 'test/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
                    display(['Original data size = ' num2str(size(test_data))]);
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
                        test_data = convert_decomp_form(test_data, d_args{i_decomp}, new_decomp_args);
                    end
                    display(['Converted data size = ' num2str(size(test_data))]);
                    load([rf_dir 'predictor.mat'], 'predictor');
                    display(['Forest size = ' num2str(predictor.D)]);
                    
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

            % 1) Complex forms: all (mag/phase), real/imag, conj, mag, phase
            for feature_type = {'real_imag', 'conj', 'mag', 'phase', 'imag', 'real', 'all'}

                rf_dir = [exp_dir '/rfs/' num2str(repeat) '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/feature_types/' feature_type{1} '/'];
                results_dir = [exp_dir '/results/' num2str(repeat) '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/feature_types/' feature_type{1} '/'];
                
                if exist([rf_dir 'predictor.mat'], 'file') && ~exist([results_dir 'results.mat'], 'file')
                    
                    %Load data for this decomp type
                    test_data = u_load([exp_dir 'test/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);

                    %Args to reform data
                    new_decomp_args = [];
                    new_decomp_args.feature_type = feature_type{1};

                    if isfield(d_args{i_decomp}, 'num_angles') && (d_args{i_decomp}.num_angles == 18)
                        new_decomp_args.bands = 1:3:16;
                    end
                    load([rf_dir 'predictor.mat'], 'predictor');         

                    %Test
                    test_data_i = convert_decomp_form(test_data, d_args{i_decomp}, new_decomp_args);
                    predictions = random_forest_reg_predict(predictor, test_data_i);
                    save([results_dir 'results.mat'], 'predictions');
                end
            end
        end
    end
    %
    % 3) Different number of levels
    %
    if ismember(3, args.do_tests) && ismember(i_decomp, [1 3 5 6]);
        
        for output_type = output_types

            for win_size = win_sizes

                for i_level = 1:5

                    rf_dir = [exp_dir '/rfs/' num2str(repeat) '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/levels/' num2str(i_level) '/' num2str(win_size) '/'];
                    results_dir = [exp_dir '/results/' num2str(repeat) '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/levels/' num2str(i_level) '/' num2str(win_size) '/'];
                    
                    if exist([rf_dir 'predictor.mat'], 'file') && ~exist([results_dir 'results.mat'], 'file')

                        %Load data for this decomp type
                        test_data = u_load([exp_dir 'test/' num2str(repeat) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);

                        %Args to reform data
                        new_decomp_args = [];
                        new_decomp_args.feature_type = 'conj';
                        new_decomp_args.levels = i_level;
                        new_decomp_args.win_size = win_size;
                        
                        load([rf_dir 'predictor.mat'], 'predictor');         

                        %Test
                        test_data_i = convert_decomp_form(test_data, d_args{i_decomp}, new_decomp_args);
                        predictions = random_forest_reg_predict(predictor, test_data_i);
                        save([results_dir 'results.mat'], 'predictions');
                    end
                end

            end
        end
    end
    
    
    if ismember(6, args.do_tests)&& i_decomp == 7
        for win_size = win_sizes
            
            %Set up directory to save RF to
            rf_dir = [exp_dir '/rfs/' num2str(repeat) '/orientation/all/' num2str(win_size) '/'];
            results_dir = [exp_dir '/results/' num2str(repeat) '/orientation/all/' num2str(win_size) '/'];
            
            if exist([rf_dir 'predictor.mat'], 'file') && ~exist([results_dir 'results.mat'], 'file')
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
                    test_data(:, cols) =...
                        convert_decomp_form(test_data_i, d_args{i_d}, new_decomp_args);
                end
                
                load([rf_dir 'predictor.mat'], 'predictor');         

                predictions = random_forest_reg_predict(predictor, test_data, 0);
                save([results_dir 'results.mat'], 'predictions');
            end
        end
    end
    
    if ismember(7, args.do_tests) && i_decomp == 2
        for output_type = output_types
            for feature_type = {'phase', 'mag', 'conj', 'real_imag'}

                rf_dir = [exp_dir '/rfs/' num2str(repeat) '/' output_type{1} '/gh2d/feature_types/' feature_type{1} '/'];
                results_dir = [exp_dir '/results/' num2str(repeat) '/' output_type{1} '/gh2d/feature_types/' feature_type{1} '/'];
                
                display(rf_dir);
                display(results_dir);
            
                if exist([rf_dir 'predictor.mat'], 'file') && ~exist([results_dir 'results.mat'], 'file')
                    display('Re-testing');
                    
                    h2_data = u_load([exp_dir 'test/' num2str(repeat) '/responses_h2da.mat']);
                    g2_data = u_load([exp_dir 'test/' num2str(repeat) '/responses_g2da.mat']);
                    test_data = ...
                        convert_complex_representation_back([g2_data h2_data], 'real_imag'); 
                    clear h2_data g2_data;

                    test_data_i = ...
                        convert_complex_representation(test_data, feature_type{1}, 3);
                    load([rf_dir 'predictor.mat'], 'predictor');         

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
