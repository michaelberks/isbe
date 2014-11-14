data_name = 'DRIVE_clean';
args.image_dir = [asymmetryroot,'data/retinograms/' data_name '/test/images/'];
args.vessel_mask_dir = [asymmetryroot,'data/retinograms/' data_name '/test/vessel_masks/'];
args.fov_mask_dir = [asymmetryroot,'data/retinograms/' data_name '/test/fov_masks/'];
args.ori_dir = [asymmetryroot,'data/retinograms/' data_name '/test/orientations/'];
warning('off', 'load_uint8:missing_variables');

i_decomp = 1;
repeat = 1;

num_pts = 5000;
win_sizes = 3;
exp_dir = [asymmetryroot 'experiments/' data_name '/comparing_representations/interp/'];

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
d_args{1}.win_size = 3;
d_args{1}.normalise = 0;
d_args{1}.pca = [];

D = get_samples_per_channel(d_args{1});
%%
%--------------------------------------------------------------------------
% ***************** GENERATE DATA *****************************************
%--------------------------------------------------------------------------
for unwrap = 0:1
    for mag_phase = 0:1
        
        d_args{1}.unwrap_phase = unwrap;
        d_args{1}.interp_mag_phase = mag_phase;

        rand('twister', 1);
        randn('state', 1);

        display(['Generating data for ' d_args{i_decomp}.decomp_type{1}]);

        image_list = dir([args.image_dir '/*.mat']);
        fov_list = dir([args.fov_mask_dir '/*.mat']);
        vessel_list = dir([args.vessel_mask_dir '/*.mat']);
        ori_list = dir([args.ori_dir '/*.mat']);
        
        %Pre-allocate space for the true line parameters
        true_oris_test = zeros(num_pts,1);
        true_centre_test = false(num_pts,1);

        %Pre-allocate space for the responses
        responses_test = zeros(num_pts, D);
        bg_responses_test = zeros(num_pts, D);

        %Check which images are selected - we'll assume if images are selected the
        %user has managed to index images within the corrcet range
        selected_images = 1:length(image_list);
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
            true_centre_test(sample_idx) = centre_mask(v_idx(1:end/2));

            %For each decomposition compute the responses
            im_responses = compute_filter_responses(ret, d_args{i_decomp}); 
            responses_test(sample_idx, :) = ...
                sample_image_features(im_responses, v_rows(1:end/2), v_cols(1:end/2), d_args{i_decomp});
            bg_responses_test(sample_idx, :) = ...
                sample_image_features(im_responses, b_rows(1:end/2), b_cols(1:end/2), d_args{i_decomp});

            %Update the current sample count
            curr_sample = curr_sample + num_samples_image;

        end
        mkdir([exp_dir 'test/' num2str(repeat)]);

        save([exp_dir 'test/' num2str(repeat) '/responses_dt_u' num2str(unwrap) '_i' num2str(mag_phase) '.mat'],...
            'responses_test');
        save([exp_dir 'test/' num2str(repeat) '/bg_responses_dt_u' num2str(unwrap) '_i' num2str(mag_phase) '.mat'],...
            'bg_responses_test');
        save([exp_dir 'test/' num2str(repeat) '/true_labels.mat'],...
            'true*test');

    end
end
%%
% --------------------------------------------------------------------------
% ********************** DETECTION ****************************************
% --------------------------------------------------------------------------

%
% 2) Complex representations
%
for unwrap = 0:1
    for mag_phase = 0:1
        pack;        

        % 1) Complex forms: all (mag/phase), real/imag, conj, mag, phase
        for feature_type = {'conj', 'mag', 'phase', 'real_imag'}

            rf_dir = [exp_dir '/rfs/' num2str(repeat) '/detection/'...
                d_args{i_decomp}.decomp_type{1} '/feature_types/' feature_type{1} '/'];
            
            results_dir = [exp_dir '/results/' num2str(repeat) '/detection/'...
                d_args{i_decomp}.decomp_type{1} '/feature_types/' feature_type{1} '/u' num2str(unwrap) '_i' num2str(mag_phase) '/'];
            mkdir(results_dir);

            %Args to reform data
            new_decomp_args = [];
            new_decomp_args.feature_type = feature_type{1};

            %Load predictor
            predictor = u_load([rf_dir 'predictor.mat']);    
            predictor.tree_dir = [rf_dir 'trees/'];
            
            %Load data
            fg_data = u_load([exp_dir 'test/' num2str(repeat) '/responses_dt_u' num2str(unwrap) '_i' num2str(mag_phase) '.mat']);
            fg_data = convert_decomp_form(fg_data, d_args{i_decomp}, new_decomp_args);
                
            bg_data = u_load([exp_dir 'test/' num2str(repeat) '/bg_responses_dt_u' num2str(unwrap) '_i' num2str(mag_phase) '.mat']);
            bg_data = convert_decomp_form(bg_data, d_args{i_decomp}, new_decomp_args);
            
            test_data = [fg_data; bg_data]; 
            clear fg_data bg_data;

            %Test
            [dummy, votes] = random_forest_class_predict(predictor, test_data);
            predicted_lines = votes(:,2) / length(predictor.trees);
            save([results_dir 'results.mat'], 'predicted_lines');
        end
    end
end
%%
class_labels = [true(num_pts,1); false(num_pts,1)];
operating_pts = (-1:101)/100;
num_repeats = 1;


for feature_type = {'conj', 'mag', 'phase', 'real_imag'}
    for unwrap = 0:1
        for mag_phase = 0:1
            auc = nan(num_repeats,1);
            for i_repeat = 1:num_repeats
                results_dir = [exp_dir '/results/' num2str(repeat) '/detection/'...
                    d_args{i_decomp}.decomp_type{1} '/feature_types/' feature_type{1} '/u' num2str(unwrap) '_i' num2str(mag_phase) '/'];

                try
                    load([results_dir 'results.mat'], 'predicted_lines'); 
                catch
                    continue;
                end            
                [~, auc(i_repeat)] =...
                    calculate_roc_curve(predicted_lines,class_labels,operating_pts);
            end
            display([results_dir ':'])
            display(['ROC Az = ' num2str(naNmean(auc), 3) ' +/- ' num2str(naNstd(auc), '%10.1e')  ' ' num2str(sum(~isnan(auc))) ' tests complete']);
        end
    end
end