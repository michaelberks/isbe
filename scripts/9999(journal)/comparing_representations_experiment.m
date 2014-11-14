%--------------------------------------------------------------------------
% Script for producing results for journal submission
% Experiment predicting orientation of synthetic lines in the presence of increasing noise 

%--------------------------------------------------------------------------
%%
%1. Generic arguments to initialise
clear, clc
data_dir = 'C:\isbe\asymmetry_project\data\synthetic_lines\comparing_representations\';
exp_dir = 'C:\isbe\asymmetry_project\experiments\synthetic_lines\comparing_representations\';

dx = 64;
dy = 64;
cx = 32;
cy = 32;
warning('off', 'ASYM:unexpectedArgument');
%%
%Make a load of data for bothing training and testing purposes
% Set up arguments for each decomposition type
d_args{1}.decomp_type = {'dt'};
d_args{1}.levels = 1:5;
d_args{1}.feature_shape = 'rect';
d_args{1}.feature_type = 'complex';
d_args{1}.do_max = 0;
d_args{1}.rotate = 0;
d_args{1}.use_nag = 0;

d_args{2}.decomp_type = {'linop'};
d_args{2}.num_levels = 5;
d_args{2}.num_angles = 6;
d_args{2}.do_max = 0;
d_args{2}.rotate = 0;       

d_args{3}.decomp_type = {'gabor'};
d_args{3}.num_angles = 18;
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

num_decomps = 10;

%set common parameters for all decomp types and then compute vector sizes
D = zeros(num_decomps,1);
for i_decomp = 1:num_decomps
    d_args{i_decomp}.win_size = 3;
    d_args{i_decomp}.normalise = 0;
    d_args{i_decomp}.pca = [];
    
    D(i_decomp) = get_samples_per_channel(d_args{i_decomp});
end

%%
%--------------------------------------------------------------------------
% ***************** GENERATE DATA *****************************************
%--------------------------------------------------------------------------
num_ims = 2000;
pts_per_im = 10;

decomps = 5%1:num_decomps;

for data_type = {'training', 'test'}
    
    for noise_level = 4%[4 2 1];
        mkdir([exp_dir data_type{1} '/rician_' num2str(noise_level)]);

        %Pre-allocate space for the true line parameters
        true_oris = zeros(num_ims*pts_per_im,1);
        true_cons = zeros(num_ims*pts_per_im,1);
        true_widths = zeros(num_ims*pts_per_im,1);
        centre_idx = false(num_ims*pts_per_im,1);

        for i_decomp = decomps
            %Pre-allocate space for the responses
            eval(['responses_' d_args{i_decomp}.decomp_type{1} ' = zeros(num_ims*pts_per_im, D(i_decomp));']);
        end

        if strcmpi(data_type{1}, 'training')
            rng(1000, 'twister');
        else
            rng(2000, 'twister');
        end
        
        for ii = 1:num_ims

            display(['Testing image ' num2str(ii)]);

            %Sample properties of line
            line_width = sample_uniform([1 8]);
            line_contrast = sample_uniform([1 8]);
            line_ori = sample_uniform([0 180]);
            line_rad = pi * line_ori / 180;

            %Generate line
            [line, label, label_centre] =...
                create_ellipse_bar(line_width/2, line_contrast, line_ori, dy, dx, cx, cy);
            test_image = ricernd(1 + line, noise_level);

            label([1 end],:) = 0;
            label(:, [1 end]) = 0;

            %Select some random pixels in the signal image
            shuffle = randperm(sum(label(:)));
            idx = find(label);
            r_idx = idx(shuffle(1:pts_per_im));
            [r_rows r_cols] = ind2sub([dy dx], r_idx);
            sample_idx = ((ii-1)*pts_per_im + 1):(ii*pts_per_im);

            %Save the line parameters
            true_oris(sample_idx,:) = line_rad;
            true_cons(sample_idx,:) = line_contrast;
            true_widths(sample_idx,:) = line_width;
            centre_idx(sample_idx,:) = label_centre(r_idx);

            %For ecah decomposition compute the responses
            for i_decomp = decomps
                im_responses = compute_filter_responses(test_image, d_args{i_decomp});

                eval(['responses_' d_args{i_decomp}.decomp_type{1} '(sample_idx, :) = '...
                    'sample_image_features(im_responses, r_rows(:), r_cols(:), d_args{i_decomp});']);
            end
        end
        for i_decomp = decomps
            save([exp_dir data_type{1} '/rician_' num2str(noise_level) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat'],...
                ['responses_' d_args{i_decomp}.decomp_type{1}]);
        end
        save([exp_dir data_type{1} '/rician_' num2str(noise_level) '/true_labels.mat'], 'true_*', 'centre_idx');

    end
end
%
num_ims = 2000;
pts_per_im = 10;

decomps = 5%1:num_decomps;

for data_type = {'training', 'test'}
    
    for noise_level = 4%[4 2 1];
        mkdir([exp_dir data_type{1} '/rician_' num2str(noise_level)]);

        for i_decomp = decomps
            %Pre-allocate space for the responses
            eval(['responses_' d_args{i_decomp}.decomp_type{1} ' = zeros(num_ims*pts_per_im, D(i_decomp));']);
        end

        if strcmpi(data_type{1}, 'training')
            rng(1000, 'twister');
        else
            rng(2000, 'twister');
        end
        
        for ii = 1:num_ims

            display(['Testing image ' num2str(ii)]);

            %Sample properties of line
            line_width = sample_uniform([1 8]);
            line_contrast = sample_uniform([1 8]);
            line_ori = sample_uniform([0 180]);
            line_rad = pi * line_ori / 180;

            %Generate line
            [line, label, label_centre] =...
                create_ellipse_bar(line_width/2, line_contrast, line_ori, dy, dx, cx, cy);
            test_image = ricernd(1 + line, noise_level);

            label([1 end],:) = 1;
            label(:, [1 end]) = 1;

            %Select some random pixels in the signal image
            shuffle = randperm(sum(~label(:)));
            idx = find(~label);
            r_idx = idx(shuffle(1:pts_per_im));
            [r_rows r_cols] = ind2sub([dy dx], r_idx);
            sample_idx = ((ii-1)*pts_per_im + 1):(ii*pts_per_im);

            %For ecah decomposition compute the responses
            for i_decomp = decomps
                im_responses = compute_filter_responses(test_image, d_args{i_decomp});

                eval(['responses_' d_args{i_decomp}.decomp_type{1} '(sample_idx, :) = '...
                    'sample_image_features(im_responses, r_rows(:), r_cols(:), d_args{i_decomp});']);
            end
        end
        for i_decomp = decomps
            save([exp_dir data_type{1} '/rician_' num2str(noise_level) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat'],...
                ['responses_' d_args{i_decomp}.decomp_type{1}]);
        end

    end
end
%%
%--------------------------------------------------------------------------
%Now make a bunch of different forests for each decomp type
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% *********************** ORIENTATION + WIDTH *****************************
%--------------------------------------------------------------------------

%First build a forest for each original type
load([exp_dir 'training/rician_' num2str(noise_level) '/true_labels.mat'], 'true_*');
rf_args.prediction_type = 'rf_regression';
rf_args.n_trees = 200;
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
% 1) Original data
%
for output_type = {'orientation'}%'width', 'orientation'}
    
    switch output_type{1}                       
        case 'width'
            rf_args.sampling_args.y = true_widths;    
        case 'orientation'
            rf_args.sampling_args.y = complex(cos(2*true_oris), sin(2*true_oris));
    end

    for i_decomp = 5%1:10;

        %Load data for this decomp type
        orig_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
        test_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);

        %Set up directory to save RF to
        rf_dir = [exp_dir '/rf/rician_' num2str(noise_level) '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/orig/'];
        results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/orig/'];
        rf_args.tree_dir = [rf_dir 'trees/'];
        mkdir(results_dir);

        new_decomp_args = [];
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

        %Train
        rf_args.sampling_args.X = orig_data;
        predictor = random_forest_reg_train(rf_args);
        save([rf_dir 'predictor.mat'], 'predictor');
        display(['************ FOREST for ' d_args{i_decomp}.decomp_type{1} ' complete!! **************']);

        %Test
        predictions = random_forest_reg_predict(predictor, test_data);
        save([results_dir 'results.mat'], 'predictions');
    end
end
%%
%
% 2) Complex representations
%
for output_type = {'width', 'orientation'}
    
    switch output_type{1}                       
        case 'width'
            rf_args.sampling_args.y = true_widths;    
        case 'orientation'
            rf_args.sampling_args.y = complex(cos(2*true_oris), sin(2*true_oris));
    end
    %Dual-tree tests
    for i_decomp = [1 3];

        %Load data for this decomp type
        orig_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);

        %Load data for this decomp type
        test_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);

        % 1) Complex forms: all (mag/phase), real/imag, conj, mag, phase
        for feature_type = {'all' 'real_imag', 'conj', 'mag', 'phase'}

            rf_dir = [exp_dir '/rf/rician_' num2str(noise_level) '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/feature_types/' feature_type{1} '/'];
            results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/feature_types/' feature_type{1} '/'];
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
for output_type = {'width', 'orientation'}
    
    switch output_type{1}                       
        case 'width'
            rf_args.sampling_args.y = true_widths;    
        case 'orientation'
            rf_args.sampling_args.y = complex(cos(2*true_oris), sin(2*true_oris));
    end
    for i_decomp = [1 3 6];

        %Load data for this decomp type
        orig_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);

        %Load data for this decomp type
        test_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);


        for i_level = 1:4

            rf_dir = [exp_dir '/rf/rician_' num2str(noise_level) '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/levels/' num2str(i_level) '/'];
            results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/levels/' num2str(i_level) '/'];
            mkdir(results_dir);


            %Args to reform data
            new_decomp_args = [];
            new_decomp_args.feature_type = 'conj';
            new_decomp_args.levels = i_level;
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
% 4) Different ways of stitching Gaussian filters together
% G" + G'
if i_decomp == 5
    for output_type = {'width', 'orientation'}
    
        switch output_type{1}                       
            case 'width'
                rf_args.sampling_args.y = true_widths;    
            case 'orientation'
                rf_args.sampling_args.y = complex(cos(2*true_oris), sin(2*true_oris));
        end
        
        g1_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/responses_g1d.mat']);
        g2_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/responses_g2d.mat']);
        rf_args.sampling_args.X = [g1_data g2_data]; clear g1_data g2_data;

        g1_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/responses_g1d.mat']);
        g2_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/responses_g2d.mat']);
        test_data = [g1_data g2_data]; clear g1_data g2_data;

        rf_dir = [exp_dir '/' output_type{1} '/rf/rician_' num2str(noise_level) '/g12d/'];
        results_dir = [exp_dir '/' output_type{1} '/results/rician_' num2str(noise_level) '/g12d/'];
        mkdir(results_dir);

        rf_args.tree_dir = [rf_dir 'trees/'];
        predictor = random_forest_reg_train(rf_args);    
        save([rf_dir 'predictor.mat'], 'predictor');

        predictions = random_forest_reg_predict(predictor, test_data, 0);
        save([results_dir 'results.mat'], 'predictions');
        
        %G" + H"
        h2_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/responses_h2d.mat']);
        g2_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/responses_g2d.mat']);
        rf_args.sampling_args.X = [h2_data g2_data]; clear h2_data g2_data;

        h2_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/responses_h2d.mat']);
        g2_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/responses_g2d.mat']);
        test_data = [h2_data g2_data]; clear h2_data g2_data;

        rf_dir = [exp_dir '/' output_type{1} '/rf/rician_' num2str(noise_level) '/gh2d/'];
        results_dir = [exp_dir '/' output_type{1} '/results/rician_' num2str(noise_level) '/gh2d/'];
        mkdir(results_dir);

        rf_args.tree_dir = [rf_dir 'trees/'];
        predictor = random_forest_reg_train(rf_args);    
        save([rf_dir 'predictor.mat'], 'predictor');

        predictions = random_forest_reg_predict(predictor, test_data, 0);
        save([results_dir 'results.mat'], 'predictions');
    end
end
%
% 6) Different number of angles for Gabor filters
%
for output_type = {'width', 'orientation'}
    
    switch output_type{1}                       
        case 'width'
            rf_args.sampling_args.y = true_widths;    
        case 'orientation'
            rf_args.sampling_args.y = complex(cos(2*true_oris), sin(2*true_oris));
    end
    if ismember(i_decomp, 3);

        %Load data for this decomp type
        orig_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);

        %Load data for this decomp type
        test_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);


        for num_angles = [3 9 18]
        
            rf_dir = [exp_dir '/rf/rician_' num2str(noise_level) '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/num_angles/' num2str(num_angles) '/'];
            results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/num_angles/' num2str(num_angles) '/'];
            mkdir(results_dir);


            %Args to reform data
            new_decomp_args = [];
            new_decomp_args.feature_type = 'conj';
            band_step = 18 / num_angles;
            new_decomp_args.bands = 1:band_step:18;

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

%%
%--------------------------------------------------------------------------
% ********************** DETECTION ****************************************
%--------------------------------------------------------------------------
rf_args.prediction_type = 'rf_regression';
rf_args.n_trees = 200;
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
rf_args.sampling_args.y = [true(20000,1); false(20000,1)];
%
% 1) Original data
%
for i_decomp = 5%1:10;

    %Load data for this decomp type
    fg_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
    bg_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
    orig_data = [fg_data; bg_data]; 
    clear fg_data bg_data;
    fg_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
    bg_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
    test_data = [fg_data; bg_data]; 
    clear fg_data bg_data;

    %Set up directory to save RF to
    rf_dir = [exp_dir '/rf/rician_' num2str(noise_level) '/detection' d_args{i_decomp}.decomp_type{1} '/orig/'];
    rf_args.tree_dir = [rf_dir 'trees/'];
    mkdir(results_dir);

    new_decomp_args = [];
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

    %Train
    rf_args.sampling_args.X = orig_data;
    predictor = random_forest_class_train(rf_args);
    save([rf_dir 'predictor.mat'], 'predictor');
    display(['************ FOREST for ' d_args{i_decomp}.decomp_type{1} ' complete!! **************']);

    %Test
    [~, ~, all_votes] = random_forest_class_predict(predictor, test_data);
    for i_trees = [50 75 125 150 175 200]
        results_dir = [exp_dir '/results/rician_' num2str(noise_level)...
            '/detection' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(i_trees) '/'];
 
        predicted_lines = sum(all_votes(:,2,1:i_trees),3) / i_trees;
        save([results_dir 'results.mat'], 'predicted_lines');
    end
end
%%
%
% 2) Complex representations
%
for i_decomp = [1 3];

    %Load data for this decomp type
    fg_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
    bg_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
    orig_data = [fg_data; bg_data]; 
    clear fg_data bg_data;
    fg_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
    bg_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
    test_data = [fg_data; bg_data]; 
    clear fg_data bg_data;

    % 1) Complex forms: all (mag/phase), real/imag, conj, mag, phase
    for feature_type = {'all' 'real_imag', 'conj', 'mag', 'phase'}

        rf_dir = [exp_dir '/rf/rician_' num2str(noise_level) '/detection/' d_args{i_decomp}.decomp_type{1} '/feature_types/' feature_type{1} '/'];
        results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/detection/' d_args{i_decomp}.decomp_type{1} '/feature_types/' feature_type{1} '/'];
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
for i_decomp = [1 3 6];
    
    %Load data for this decomp type
    fg_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
    bg_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
    orig_data = [fg_data; bg_data]; 
    clear fg_data bg_data;
    fg_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
    bg_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
    test_data = [fg_data; bg_data]; 
    clear fg_data bg_data;
    
    for i_level = 1:4
        
        rf_dir = [exp_dir '/rf/rician_' num2str(noise_level) '/detection/' d_args{i_decomp}.decomp_type{1} '/levels/' num2str(i_level) '/'];
        results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/detection/' d_args{i_decomp}.decomp_type{1} '/levels/' num2str(i_level) '/'];
        mkdir(results_dir);

        %Args to reform data
        new_decomp_args = [];
        new_decomp_args.feature_type = 'conj';
        new_decomp_args.levels = i_level;
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
% 4) Rotate and do max
%
for i_decomp = [1 3 9]
    
    %Load data for this decomp type
    fg_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
    bg_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
    orig_data = [fg_data; bg_data]; 
    clear fg_data bg_data;
    fg_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
    bg_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
    test_data = [fg_data; bg_data]; 
    clear fg_data bg_data;
    
    for i_level = 1:6
        
        for win_size = [1 3]
            
            for reform_type = 1:2
                
                if win_size==1 && reform_type==2
                    continue;
                end
                
                new_decomp_args = [];
                if reform_type == 1
                    rf_dir = [exp_dir '/rf/rician_' num2str(noise_level) '/detection/' d_args{i_decomp}.decomp_type{1} '/rotate/L' num2str(i_level) '_W' num2str(win_size)  '/'];
                    results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/detection/' d_args{i_decomp}.decomp_type{1} '/rotate/L' num2str(i_level) '_W' num2str(win_size) '/'];
                    new_decomp_args.rotate = 1;
                else
                    rf_dir = [exp_dir '/rf/rician_' num2str(noise_level) '/detection/' d_args{i_decomp}.decomp_type{1} '/do_max/L' num2str(i_level) '_W' num2str(win_size)  '/'];
                    results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/detection/' d_args{i_decomp}.decomp_type{1} '/do_max/L' num2str(i_level) '_W' num2str(win_size) '/'];
                    new_decomp_args.do_max = 1;
                end
                mkdir(results_dir);


                %Args to reform data
                new_decomp_args.feature_type = 'conj';
                if i_level < 6
                    new_decomp_args.levels = i_level;
                else
                    new_decomp_args.levels = 1:5;
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
if i_decomp == 5
        
    g1_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/responses_g1d.mat']);
    g2_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/responses_g2d.mat']);
    rf_args.sampling_args.X = [g1_data g2_data]; clear g1_data g2_data;

    g1_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/responses_g1d.mat']);
    g2_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/responses_g2d.mat']);
    test_data = [g1_data g2_data]; clear g1_data g2_data;

    rf_dir = [exp_dir 'detection/rf/rician_' num2str(noise_level) '/g12d/'];
    results_dir = [exp_dir 'detection/results/rician_' num2str(noise_level) '/g12d/'];
    mkdir(results_dir);

    rf_args.tree_dir = [rf_dir 'trees/'];
    predictor = random_forest_class_train(rf_args);    
    save([rf_dir 'predictor.mat'], 'predictor');

    [dummy, votes] = random_forest_class_predict(predictor, test_data);
    predicted_lines = votes(:,2) / rf_args.n_trees;
    save([results_dir 'results.mat'], 'predicted_lines');

    %G" + H"
    h2_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/responses_h2d.mat']);
    g2_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/responses_g2d.mat']);
    rf_args.sampling_args.X = [h2_data g2_data]; clear h2_data g2_data;

    h2_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/responses_h2d.mat']);
    g2_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/responses_g2d.mat']);
    test_data = [h2_data g2_data]; clear h2_data g2_data;

    rf_dir = [exp_dir 'detection/rf/rician_' num2str(noise_level) '/gh2d/'];
    results_dir = [exp_dir 'detection/results/rician_' num2str(noise_level) '/gh2d/'];
    mkdir(results_dir);

    rf_args.tree_dir = [rf_dir 'trees/'];
    predictor = random_forest_class_train(rf_args);    
    save([rf_dir 'predictor.mat'], 'predictor');

    [dummy, votes] = random_forest_class_predict(predictor, test_data);
    predicted_lines = votes(:,2) / rf_args.n_trees;
    save([results_dir 'results.mat'], 'predicted_lines');
end
%
% 6) Different number of angles for Gabor filters
%
if i_decomp == 3
      
    fg_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
    bg_data = u_load([exp_dir 'training/rician_' num2str(noise_level) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
    orig_data = [fg_data; bg_data]; 
    clear fg_data bg_data;
    fg_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
    bg_data = u_load([exp_dir 'test/rician_' num2str(noise_level) '/bg_responses_' d_args{i_decomp}.decomp_type{1} '.mat']);
    test_data = [fg_data; bg_data]; 
    clear fg_data bg_data;
    
    for num_angles = [3 9 18]
        
        rf_dir = [exp_dir '/rf/rician_' num2str(noise_level) '/detection/' d_args{i_decomp}.decomp_type{1} '/num_angles/' num2str(num_angles) '/'];
        results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/detection/' d_args{i_decomp}.decomp_type{1} '/num_angles/' num2str(num_angles) '/'];
        mkdir(results_dir);

        %Args to reform data
        new_decomp_args = [];
        new_decomp_args.feature_type = 'conj';
        band_step = 18 / num_angles;
        new_decomp_args.bands = 1:band_step:18;

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
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% **************** Analyse all the results ********************************
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% *********************** ORIENTATION + WIDTH *****************************
%--------------------------------------------------------------------------
%%
%
% 1) Original data
%
for output_type = {'orientation'}%{'width'}% 
    
    for i_decomp = 1:10;
        for win_size = [1 3]
            
            mae = zeros(num_repeats,1);
            for i_repeat = 1:num_repeats
                results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(i_repeat) '/'...
                    output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/20000/'];
                load([results_dir 'results.mat'], 'predictions');
                

                switch output_type{1}
                    case 'orientation'
                        load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(i_repeat) '/true_labels.mat'], 'true_oris');
                        complex_orientations = complex(cos(2*true_oris), sin(2*true_oris));
                        
                        [~, error_stats] = ori_error(complex_orientations, predictions);
                        mae(i_repeat) = error_stats.abs_mean;
                    case 'width'
                        load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(i_repeat) '/true_labels.mat'], 'true_widths');
                        width_errors = abs(true_widths - predictions);
                        mae(i_repeat) = mean(width_errors);
                        
                end
            end
            display([results_dir ':']);
            display(['MAE: ' num2str(mean(mae),4) ' +/- ' num2str(std(mae),3)]);
        end
    end
end


%%
%
% 2) Complex representations
%
for output_type = {'orientation'}%{'width'}%

    %Dual-tree tests
    for i_decomp = [1 3];

        % 1) Complex forms: all (mag/phase), real/imag, conj, mag, phase
        for feature_type = {'all', 'real_imag', 'conj', 'mag', 'phase', 'real', 'imag', 'real_abs_imag'}

            mae = zeros(num_repeats,1);
            for i_repeat = 1:num_repeats
                results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(i_repeat) '/'...
                    output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/feature_types/' feature_type{1} '/'];
                load([results_dir 'results.mat'], 'predictions');
                

                switch output_type{1}
                    case 'orientation'
                        load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(i_repeat) '/true_labels.mat'], 'true_oris');
                        complex_orientations = complex(cos(2*true_oris), sin(2*true_oris));
                        
                        [~, error_stats] = ori_error(complex_orientations, predictions);
                        mae(i_repeat) = error_stats.abs_mean;
                    case 'width'
                        load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(i_repeat) '/true_labels.mat'], 'true_widths');
                        width_errors = abs(true_widths - predictions);
                        mae(i_repeat) = mean(width_errors);
                        
                end
            end
            display([results_dir ':']);
            display(['MAE: ' num2str(mean(mae),4) ' +/- ' num2str(std(mae),3)]);
        end
    end
end
%%
%
% 3) Different number of levels
%
for output_type = {'orientation'}%{'width'}% 
    
    for i_decomp = [1 3 6];

       for i_level = 1:4
           
           for win_size = [1 3]

                mae = zeros(num_repeats,1);
                for i_repeat = 1:num_repeats
                    results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(i_repeat) '/'...
                        output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/levels/' num2str(i_level) '/' num2str(win_size) '/'];
                    load([results_dir 'results.mat'], 'predictions');


                    switch output_type{1}
                        case 'orientation'
                            load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(i_repeat) '/true_labels.mat'], 'true_oris');
                            complex_orientations = complex(cos(2*true_oris), sin(2*true_oris));

                            [~, error_stats] = ori_error(complex_orientations, predictions);
                            mae(i_repeat) = error_stats.abs_mean;
                        case 'width'
                            load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(i_repeat) '/true_labels.mat'], 'true_widths');
                            width_errors = abs(true_widths - predictions);
                            mae(i_repeat) = mean(width_errors);

                    end
                end
                display([results_dir ':']);
                display(['MAE: ' num2str(mean(mae),4) ' +/- ' num2str(std(mae),3)]);
           end
        end
    end
end
%%
% 4) Different ways of stitching Gaussian filters together
 

% G" + G'
mae = zeros(10,1);
for i_repeat = 1:11
    results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(i_repeat) '/orientation/g12d/'];
    load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(i_repeat) '/true_labels.mat'], 'true_oris');
    load([results_dir 'results.mat'], 'predictions');
    complex_orientations = complex(cos(2*true_oris), sin(2*true_oris));
    [~, error_stats] = ori_error(complex_orientations, predictions);
    mae(i_repeat) = error_stats.abs_mean;
end
display([results_dir ':']);
display(['MAE: ' num2str(mean(mae),4) ' +/- ' num2str(std(mae),3)]);

%G" + H"
mae = zeros(10,1);
for i_repeat = 1:11
    results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(i_repeat) '/orientation/gh2d/'];
    load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(i_repeat) '/true_labels.mat'], 'true_oris');
    load([results_dir 'results.mat'], 'predictions');
    complex_orientations = complex(cos(2*true_oris), sin(2*true_oris));
    [~, error_stats] = ori_error(complex_orientations, predictions);
    mae(i_repeat) = error_stats.abs_mean;
end
display([results_dir ':']);
display(['MAE: ' num2str(mean(mae),4) ' +/- ' num2str(std(mae),3)]);
    %%
% G" + G'
mae = zeros(10,1);
for i_repeat = 1:11
    results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(i_repeat) '/width/g12d/'];
    load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(i_repeat) '/true_labels.mat'], 'true_widths');
    load([results_dir 'results.mat'], 'predictions');
    width_errors = abs(true_widths - predictions);
    mae(i_repeat) = mean(width_errors);
end
display([results_dir ':']);
display(['MAE: ' num2str(mean(mae),4) ' +/- ' num2str(std(mae),3)]);

%G" + H"
mae = zeros(10,1);
for i_repeat = 1:11
    results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(i_repeat) '/width/gh2d/'];
    load([exp_dir 'test/rician_' num2str(noise_level) '/' num2str(i_repeat) '/true_labels.mat'], 'true_widths');
    load([results_dir 'results.mat'], 'predictions');
    width_errors = abs(true_widths - predictions);
    mae(i_repeat) = mean(width_errors);
end
display([results_dir ':']);
display(['MAE: ' num2str(mean(mae),4) ' +/- ' num2str(std(mae),3)]);

%%
%
% 6) Different number of angles for Gabor filters
%
for output_type = {'width', 'orientation'}   
    for i_decomp = 3
        for num_angles = [3 9 18]        
            results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/num_angles/' num2str(num_angles) '/'];
            load([results_dir 'results.mat'], 'predictions');
        end
    end
end

%%
%--------------------------------------------------------------------------
% ********************** DETECTION ****************************************
%--------------------------------------------------------------------------
class_labels = [true(20000,1); false(20000,1)];
operating_pts = (-1:101)/100;
num_repeats = 1;
%%
%
% 1) Original data
%
noise_level = 2;
for i_decomp = 1:10;
    for win_size = [1 3]
        auc = zeros(num_repeats,1);
        for i_repeat = 1:num_repeats
            results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(i_repeat) '/detection/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/'];
            load([results_dir 'results.mat'], 'predicted_lines');
            [~, auc(i_repeat)] =...
                calculate_roc_curve(predicted_lines,class_labels,operating_pts);
        end
        display([results_dir ':'])
        display(['ROC Az = ' num2str(mean(auc), 3) ' +/- ' num2str(std(auc), 2)]);
    end
    display([]);
end
%%
noise_level = 2;
for i_decomp = 1:10;
    for win_size = [1 3]
        figure; hold all;
        title(['Decomp: ' d_args{i_decomp}.decomp_type{1} ', w = ' num2str(win_size)]);
        leg_text = cell(0,1);
        for i_pts = [10000:2000:20000]% 25000:5000:40000]
            auc = zeros(200,1);
            results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/1'...
                '/detection/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/' num2str(i_pts) '/'];
            
            try
                load([results_dir 'all_votes.mat'], 'all_votes'); 
            catch
                continue;
            end
            for i_trees = 5:5:100

                predicted_lines = sum(all_votes(:,2,1:i_trees),3) / i_trees;
                [~, auc(i_trees)] =...
                    calculate_roc_curve(predicted_lines,class_labels,operating_pts);
            end
            plot(5:5:100, auc(5:5:100));
            leg_text{end+1} = ['N =' num2str(i_pts)];
        end
        legend(leg_text, 'location', 'southeast');
        set(gca, 'ylim', [0.9 0.99]);
    end
end
%%
%
% 2) Complex representations
%
for i_decomp = [1 3];
    %Complex forms: all (mag/phase), real/imag, conj, mag, phase
    for feature_type = {'all' 'real_imag', 'conj', 'mag', 'phase', 'real', 'imag', 'real_abs_imag'}
        auc = zeros(num_repeats,1);
        for i_repeat = 1:num_repeats
            results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(i_repeat) '/detection/' d_args{i_decomp}.decomp_type{1} '/feature_types/' feature_type{1} '/'];
            load([results_dir 'results.mat'], 'predicted_lines');
            [~, auc(i_repeat)] =...
                calculate_roc_curve(predicted_lines,class_labels,operating_pts);
        end
        
        display([results_dir ':'])
        display(['ROC Az = ' num2str(mean(auc), 3) ' +/- ' num2str(std(auc), 2)]);

    end
end
%%
%
% 3) Different number of levels
%
for i_decomp = [1 3 6];
    
    for i_level = 1:4
        for win_size = [1 3]
            
            auc = zeros(num_repeats,1);
            for i_repeat = 1:num_repeats
                results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(i_repeat) '/detection/'...
                    d_args{i_decomp}.decomp_type{1} '/levels/' num2str(i_level) '/' num2str(win_size) '/'];
                load([results_dir 'results.mat'], 'predicted_lines');
                [~, auc(i_repeat)] =...
                    calculate_roc_curve(predicted_lines,class_labels,operating_pts);
            end

            display([results_dir ':'])
            display(['ROC Az = ' num2str(mean(auc), 3) ' +/- ' num2str(std(auc), 2)]);
        end

    end
end
%%
%
% 4) Rotate and do max
%
for i_decomp = [1 3 9]    
    for i_level = 1:5        
        for win_size = [1 3]            
            for reform_type = 1:2                
                if win_size==1 && reform_type==2
                    continue;
                end
                
                auc = zeros(num_repeats,1);
                for i_repeat = 1:num_repeats
                    if reform_type == 1
                        results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(i_repeat) '/detection/' d_args{i_decomp}.decomp_type{1}...
                            '/rotate/L' num2str(i_level) '_W' num2str(win_size) '/'];
                    else
                        results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(i_repeat) '/detection/' d_args{i_decomp}.decomp_type{1}...
                            '/do_max/L' num2str(i_level) '_W' num2str(win_size) '/'];
                    end
                    load([results_dir 'results.mat'], 'predicted_lines');
                    [~, auc(i_repeat)] =...
                        calculate_roc_curve(predicted_lines,class_labels,operating_pts);
                end
                display([results_dir ':'])
                display(['ROC Az = ' num2str(mean(auc), 3) ' +/- ' num2str(std(auc), 2)]);
            end
        end

    end
end
%%
%
% 5) Different ways of stitching Gaussian filters together
%
% G" + G'    
auc = zeros(num_repeats,1);
for i_repeat = 1:num_repeats
    results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(i_repeat) '/detection/g12d/'];
    load([results_dir 'results.mat'], 'predicted_lines');
    [~, auc(i_repeat)] =...
        calculate_roc_curve(predicted_lines,class_labels,operating_pts);
end

display([results_dir ':'])
display(['ROC Az = ' num2str(mean(auc), 3) ' +/- ' num2str(std(auc), 2)]);

%G" + H"
auc = zeros(num_repeats,1);
for i_repeat = 1:num_repeats
    results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/' num2str(i_repeat) '/detection/gh2d/'];
    load([results_dir 'results.mat'], 'predicted_lines');
    [~, auc(i_repeat)] =...
        calculate_roc_curve(predicted_lines,class_labels,operating_pts);
end

display([results_dir ':'])
display(['ROC Az = ' num2str(mean(auc), 3) ' +/- ' num2str(std(auc), 2)]);
%%
%
% 6) Different number of angles for Gabor filters
%
if i_decomp == 3
    for num_angles = [3 9 18]        
        results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/detection/' d_args{i_decomp}.decomp_type{1} '/num_angles/' num2str(num_angles) '/'];
        load([results_dir 'results.mat'], 'predicted_lines');

    end
end
%%
%---------------------------------------------------------------------------
% Analyse analytic predictions for orientation and width
data_dir = [exp_dir 'test\rician_0\9\'];
responses_g2d = u_load([data_dir 'responses_g2d.mat']);

responses_g2d = reshape(responses_g2d, [], 9, 4, 3);
responses_g2d = responses_g2d(:, 5, :, :);

load('true_labels.mat', 'true_oris', 'true_widths');

%--------------------------------------------------------------------------
%%
% 
% %2. Make a set of test images
% mkdir([data_dir 'images']);
% mkdir([data_dir 'orientations']);
% mkdir([data_dir 'widths']);
% mkdir([data_dir 'labels']);
% mkdir([data_dir 'params']);
% 
% for ii = 1:100
%     line_width = sample_uniform([1 8]);
%     line_contrast = sample_uniform([4 8]);
%     line_ori = sample_uniform([0 180]);
%     
%     [line, line_label, ~, orientation_map] =...
%         create_ellipse_bar(line_width/2, line_contrast, line_ori, dy, dx, cx, cy);
%     orientation_map = complex(cosd(2*orientation_map), sind(2*orientation_map));
%     width_map = line_label * line_width;
%     
%     test_image = ricernd(line + 1, snr);
%     
%     params.size = [dy dx];
%     params.line_centre = [cx cy];
%     params.line_type = 'ellipse';
%     params.width = line_width;
%     params.contrast = line_contrast;
%     params.orientation = line_ori;
%     parmas.noise_type = 'rician';
%     params.noise_params = snr;
%     
%     save([data_dir 'images\image' zerostr(ii, 4) '.mat'],...
%         'test_image');
%     save([data_dir 'orientations\image' zerostr(ii, 4) '_ori.mat'],...
%         'orientation_map');
%     save([data_dir 'widths\image' zerostr(ii, 4) '_width.mat'],...
%         'width_map');
%     save([data_dir 'labels\image' zerostr(ii, 4) '_label.mat'],...
%         'line_label');
%     save([data_dir 'params\image' zerostr(ii, 4) '_params.mat'],...
%         'params');
% end
% %%
% %******************* DETECTION ********************************************
% %1. Train a set of random forests on line images
% %
% DECOMP_TYPE="dt" WIN_SIZE=3 OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
% DECOMP_TYPE="g2d" WIN_SIZE=3 OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
% DECOMP_TYPE="g12d" WIN_SIZE=3 OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
% DECOMP_TYPE="mono" WIN_SIZE=3 OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
% DECOMP_TYPE="linop" WIN_SIZE=3 OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
% DECOMP_TYPE="haar" WIN_SIZE=3 OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
% DECOMP_TYPE="gabor" WIN_SIZE=3 OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
% 
% MODEL_ROOT="models/synthetic_lines/detection/rf_classification" MODEL_PATH="'10400'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10400 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
% MODEL_ROOT="models/synthetic_lines/detection/rf_classification" MODEL_PATH="'10401'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10401 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
% MODEL_ROOT="models/synthetic_lines/detection/rf_classification" MODEL_PATH="'10402'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10402 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
% MODEL_ROOT="models/synthetic_lines/detection/rf_classification" MODEL_PATH="'10403'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10403 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
% MODEL_ROOT="models/synthetic_lines/detection/rf_classification" MODEL_PATH="'10821'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10821 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
% MODEL_ROOT="models/synthetic_lines/detection/rf_classification" MODEL_PATH="'10405'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10405 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
% MODEL_ROOT="models/synthetic_lines/detection/rf_classification" MODEL_PATH="'10406'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10406 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
% %%
% %2. Predict line presence on the images using the CSF
% MODEL_ROOT="models/synthetic_lines/detection/rf_classification" MODEL_PATH="'10400'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/comparing_representations" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10407 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
% MODEL_ROOT="models/synthetic_lines/detection/rf_classification" MODEL_PATH="'10401'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/comparing_representations" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10408 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
% MODEL_ROOT="models/synthetic_lines/detection/rf_classification" MODEL_PATH="'10402'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/comparing_representations" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10409 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
% MODEL_ROOT="models/synthetic_lines/detection/rf_classification" MODEL_PATH="'10403'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/comparing_representations" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10410 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
% MODEL_ROOT="models/synthetic_lines/detection/rf_classification" MODEL_PATH="'10821'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/comparing_representations" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10822 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
% MODEL_ROOT="models/synthetic_lines/detection/rf_classification" MODEL_PATH="'10405'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/comparing_representations" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10412 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
% MODEL_ROOT="models/synthetic_lines/detection/rf_classification" MODEL_PATH="'10406'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/comparing_representations" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10413 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
% 
% %%
% %3. Analyse results
% pred_dir = [data_dir 'predictions\detection\rf_classification\'];
% label_dir = [data_dir 'labels\'];
% fov_mask_dir = [data_dir 'fov_masks\'];
% 
% rf_codes = cell(1, 4);
% rf_codes( 1,:) = {'10400', 'dt', '3', 'orig'};
% rf_codes( 2,:) = {'10401', 'g2d', '3', 'orig'};
% rf_codes( 3,:) = {'10402', 'g12d', '3', 'orig'};
% rf_codes( 4,:) = {'10403', 'mono', '3', 'orig'};
% rf_codes( 5,:) = {'10821', 'linop', '3', 'orig'};
% rf_codes( 6,:) = {'10405', 'haar', '3', 'orig'};
% rf_codes( 7,:) = {'10406', 'pixel', '25', 'orig'};
% 
% %
% %3.a compute ROC curves and Az values for the whole images
% roc_pts = [];
% auc = [];
% for ii = [1 2 3 4 6]
%     [roc_pts(:,:,ii), auc(ii,1)] = compute_roc_image_set([pred_dir rf_codes{ii,1} '\'], label_dir, [], 1); %#ok
% end
% %%
% %3.b Visually inspect the prediction images
% for ii = 1:20
%     %test_im = u_load([data_dir 'images\image' zerostr(ii,4) '.mat']);
%     figure;
%     for jj = 1:1
%         pred_im = load_uint8([data_dir 'predictions\detection\rf_classification\' rf_codes{jj,1} '\image' zerostr(ii,4) '_pred.mat']);
%         subplot(2,3,jj); imgray(pred_im); caxis([0 1]);
%         title(['Decomp: ' rf_codes{jj,2} ', w = ' rf_codes{jj,3} ', ' rf_codes{jj,4}]); 
%     end
%     %colormap(jet(256));
%     %exportfig([exp_dir 'predictions ' zerostr(ii,3) '_mp.png']);
% end
% %%
% %4. Inspect tree outputs
% %
% counts = zeros(1,48);
% for ii = 1:100
%     load(['C:\isbe\asymmetry_project\data\models\synthetic_lines\detection\rf_classification\dt\01_trees\rf_tree' zerostr(ii,4) '.mat']);
%     c = hist(tree.var(tree.var > 0), 1:48);
%     counts = counts + c;
% end
% figure; bar(1:48, counts);
% %%
% for ii = 7%[1:4 6]
%     predictor = u_load([data_dir 'rfs\' rf_codes{ii,1} '\predictor.mat']);
%     predictor.tree_root = [data_dir 'rfs\'];
%     save([data_dir 'rfs\' rf_codes{ii,1} '\predictor.mat'], 'predictor');
% end
% %%
% ori_errors = [];
% error_stats = cell(0,1);
% 
% for ii = [1 2 3 4 6]
% 
%     display(['Errors for ' rf_codes{ii,2} ', w = ' rf_codes{ii,3} ', sampling = ' rf_codes{ii,4}]); 
%     [ori_errors(:,ii), ~, error_stats{ii}] =...
%         compute_image_orientation_errors([pred_dir rf_codes{ii,1} '\'], fg_mask_dir,...
%         'label_dir', label_dir); %#ok;
%     display(error_stats{ii});
%     
% end
% %%
% ori_errors_centre = [];
% error_stats_centre = cell(0,1);
% for ii = [1 2 3 4 6]
% 
%     display(['Errors for ' rf_codes{ii,2} ', w = ' rf_codes{ii,3} ', sampling = ' rf_codes{ii,4}]); 
%     [ori_errors_centre(:,ii), ~, error_stats_centre{ii}] =...
%         compute_image_orientation_errors([pred_dir rf_codes{ii,1} '\'], fg_mask_dir,...
%         'label_dir', label_dir, 'centre_only', 1); %#ok;
%     display(error_stats_centre{ii});
%     
% end
% %%
% %************************** WIDTH *****************************************
% 
% %1. Traing models on the CSF
% DECOMP_TYPE="dt" WIN_SIZE=3 OUTPUT_TYPE="width" PREDICTION_TYPE="rf_regression" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
% 
% MODEL_ROOT="models/synthetic_lines/width/rf_regression" MODEL_PATH="'10438'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10438 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
% 
% %2. Predict line width on the images using the CSF
% MODEL_ROOT="models/synthetic_lines/width/rf_regression" MODEL_PATH="'10438'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/comparing_representations" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10439 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
% %%
% %3. Analyse results
% pred_dir = [data_dir 'predictions\width\rf_regression\'];
% label_dir = [data_dir 'widths\'];
% fg_mask_dir = [data_dir 'labels\'];
% 
% rf_codes = cell(1, 4);
% rf_codes( 1,:) = {'10438', 'dt', '3', 'orig'};
% 
% %
% %3.a compute ROC curves and Az values for the whole images
% width_errors = [];
% for ii = 1
%     width_errors(:,ii) = compute_image_width_errors([pred_dir rf_codes{ii,1} '\'], fg_mask_dir, 'label_dir', label_dir); %#ok
% end
