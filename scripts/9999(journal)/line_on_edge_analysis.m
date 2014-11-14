%--------------------------------------------------------------------------
% Script for producing results for journal submission
% Experiment predicting orientation of synthetic lines in the presence of increasing noise 

%--------------------------------------------------------------------------
%%
%1. Generic arguments to initialise
warning('off', 'ASYM:unexpectedArgument');
clear, clc
exp_dir = 'C:\isbe\asymmetry_project\experiments\synthetic_lines\line_on_edge\';
exp_dir_no_edge = [asymmetryroot 'experiments/synthetic_lines\comparing_representations/']; 

dx = 64;
dy = 64;
num_ims = 2000;
pts_per_im = 10;
noise_level = 2;
num_repeats = 1;

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
load([exp_dir 'test/rician_' num2str(noise_level) '/1/true_labels.mat'], '*idx');
%%
%
% 1) Original data
%
for i_decomp = 1:10;
    for win_size = [1 3]
        i_pts = 20000;
        results_dir_no_edge = [exp_dir_no_edge '/results/rician_' num2str(noise_level) '/1'...
                '/detection/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/' num2str(i_pts) '/'];
        load([results_dir_no_edge 'all_votes.mat'], 'all_votes');
        predicted_lines_no_edge = sum(all_votes(:,2,:),3) / 200;
        
        results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/1'...
                '/detection/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/' num2str(i_pts) '/'];
        load([results_dir 'all_votes.mat'], 'all_votes');
        predicted_lines = sum(all_votes(:,2,:),3) / 100;
        
        [~, auc_orig] = calculate_roc_curve(...
            predicted_lines_no_edge,...
            class_labels, operating_pts);
        [~, auc_edge] = calculate_roc_curve(...
            predicted_lines,...
            class_labels, operating_pts); 
        [~, auc_roi] = calculate_roc_curve(...
            predicted_lines2([true(20000,1); edge_idx]),...
            class_labels([true(20000,1); edge_idx]), operating_pts); 
        display([results_dir2 ':']);
        display(['Az (orig) = ' num2str(auc_orig,4) ...
            ', Az (edge) = ' num2str(auc_edge,4) ...
            ', Az (roi) = ' num2str(auc_roi,4)]); 
    end
end
%%
%
% 1) Original data - increasing numbers of points
%
noise_level = 2;
for i_decomp = 1:10;
    for win_size = [1 3]
        figure; hold all;
        title(['Decomp: ' d_args{i_decomp}.decomp_type{1} ', w = ' num2str(win_size)]);
        leg_text = cell(0,1);
        for i_pts = 10000:2000:20000
            auc = zeros(200,1);
            results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/1'...
                '/detection/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/' num2str(i_pts) '/'];
            
            try
                load([results_dir 'all_votes.mat'], 'all_votes'); 
            catch
                continue;
            end
            predicted_lines = sum(all_votes(:,2,1:i_trees),3) / i_trees;
            [~, auc] =...
                calculate_roc_curve(predicted_lines,class_labels,operating_pts);
            display(['N =' num2str(i_pts) ', Az =' num2str(auc,4)]);
        end
    end
end
%%
%
% 1) Original data - increasing numbers of trees
%
noise_level = 2;
for i_decomp = 1:10;
    for win_size = [1 3]
        figure; hold all;
        title(['Decomp: ' d_args{i_decomp}.decomp_type{1} ', w = ' num2str(win_size)]);
        auc = zeros(100,1);
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
    end
end
%%
%
% 2) Complex representations
%
for i_decomp = [1];
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