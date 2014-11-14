%--------------------------------------------------------------------------
% Script for producing results for journal submission
% Experiment predicting orientation of synthetic lines in the presence of increasing noise 

%--------------------------------------------------------------------------
%%
%1. Generic arguments to initialise
warning('off', 'ASYM:unexpectedArgument');
warning('off', 'ori_error:nans');
clear, clc
exp_dir = 'C:\isbe\asymmetry_project\experiments\fibre\comparing_representations\';

num_repeats = 10;
num_pts = 100000;

%Make a load of data for bothing training and testing purposes
% Set up arguments for each decomposition type
d_args{1}.decomp_type = {'dt'};

d_args{2}.decomp_type = {'h2da'};
       
d_args{3}.decomp_type = {'gabor'};

d_args{4}.decomp_type = {'mono'};

d_args{5}.decomp_type = {'g1d'};
            
d_args{6}.decomp_type = {'g2d'};
            
d_args{7}.decomp_type = {'g2di'};
            
d_args{8}.decomp_type = {'h2d'};

d_args{9}.decomp_type = {'g2da'};

d_args{10}.decomp_type = {'gabori'};

num_decomps = 10;

class_labels = [true(num_pts,1); false(num_pts,1)];
operating_pts = (-1:101)/100;
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
    
    for i_decomp = 4%1:10;
        for win_size = [1 3]
            
            mae = nan(num_repeats,1);
            for i_repeat = 1:num_repeats
                results_dir = [exp_dir '/results/' num2str(i_repeat) '/'...
                    output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/'];
                try
                    load([results_dir 'results.mat'], 'predictions');
                catch
                    continue;
                end
                

                switch output_type{1}
                    case 'orientation'
                        load([exp_dir 'test/' num2str(i_repeat) '/true_labels.mat'], 'true_oris');
                        
                        [~, error_stats] = ori_error(true_oris, predictions);
                        mae(i_repeat) = error_stats.abs_median;
                    case 'width'
                        load([exp_dir 'test/' num2str(i_repeat) '/true_labels.mat'], 'true_widths');
                        width_errors = abs(true_widths - predictions);
                        mae(i_repeat) = mean(width_errors);
                        
                end
            end
            display([results_dir ':']);
            display(['MAE: ' num2str(naNmean(mae),3) ' +/- ' num2str(naNstd(mae),2) ' ' num2str(sum(~isnan(mae))) ' tests complete']);
        end
    end
end


%%
%
% 2) Complex representations
%
for output_type = {'orientation'}%{'width'}%

    %Dual-tree tests
    for i_decomp = 1%[1 3];

        % 1) Complex forms: all (mag/phase), real/imag, conj, mag, phase
        for feature_type = {'all', 'real', 'imag', 'real_imag', 'conj', 'mag', 'phase', 'real_abs_imag'}

            mae = nan(num_repeats,1);
            for i_repeat = 1:num_repeats
                results_dir = [exp_dir '/results/' num2str(i_repeat) '/'...
                    output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/feature_types/' feature_type{1} '/1/'];
                try
                    load([results_dir 'results.mat'], 'predictions');
                catch
                    continue;
                end
                

                switch output_type{1}
                    case 'orientation'
                        load([exp_dir 'test/' num2str(i_repeat) '/true_labels.mat'], 'true_oris');
                        
                        [~, error_stats] = ori_error(true_oris, predictions);
                        mae(i_repeat) = error_stats.abs_median;
                    case 'width'
                        load([exp_dir 'test/' num2str(i_repeat) '/true_labels.mat'], 'true_widths');
                        width_errors = abs(true_widths - predictions);
                        mae(i_repeat) = mean(width_errors);
                        
                end
            end
            display([results_dir ':']);
            display(['MAE: ' num2str(naNmean(mae),3) ' +/- ' num2str(naNstd(mae),2) ' ' num2str(sum(~isnan(mae))) ' tests complete']);
        end
    end
end
%%
for output_type = {'orientation'}%{'width'}%

    % 1) Complex forms: all (mag/phase), real/imag, conj, mag, phase
    for feature_type = {'all', 'real_imag', 'conj', 'mag', 'phase'}

        mae = nan(num_repeats,1);
        for i_repeat = 1:num_repeats
            results_dir = [exp_dir '/results/' num2str(i_repeat) '/'...
                output_type{1} '/gh2d/feature_types/' feature_type{1} '/'];
            try
                load([results_dir 'results.mat'], 'predictions');
            catch
                continue;
            end


            switch output_type{1}
                case 'orientation'
                    load([exp_dir 'test/' num2str(i_repeat) '/true_labels.mat'], 'true_oris');

                    [~, error_stats] = ori_error(true_oris, predictions);
                    mae(i_repeat) = error_stats.abs_median;
                case 'width'
                    load([exp_dir 'test/' num2str(i_repeat) '/true_labels.mat'], 'true_widths');
                    width_errors = abs(true_widths - predictions);
                    mae(i_repeat) = mean(width_errors);

            end
        end
        display([results_dir ':']);
        display(['MAE: ' num2str(naNmean(mae),3) ' +/- ' num2str(naNstd(mae),2) ' ' num2str(sum(~isnan(mae))) ' tests complete']);
    end
end
%%
%
% 3) Different number of levels
%
for output_type = {'orientation'}%{'width'}% 
    
    for i_decomp = 3%[1 3 6];

       for i_level = 1:5
           
           for win_size = 3

                mae = nan(num_repeats,1);
                for i_repeat = 1:num_repeats
                    results_dir = [exp_dir '/results/' num2str(i_repeat) '/'...
                        output_type{1} '/' d_args{i_decomp}.decomp_type{1} '/levels/' num2str(i_level) '/' num2str(win_size) '/'];
                    try
                        load([results_dir 'results.mat'], 'predictions');
                    catch
                        continue;
                    end


                    switch output_type{1}
                        case 'orientation'
                            load([exp_dir 'test/' num2str(i_repeat) '/true_labels.mat'], 'true_oris');
                            
                            [~, error_stats] = ori_error(true_oris, predictions);
                            mae(i_repeat) = error_stats.abs_median;
                        case 'width'
                            load([exp_dir 'test/' num2str(i_repeat) '/true_labels.mat'], 'true_widths');
                            width_errors = abs(true_widths - predictions);
                            mae(i_repeat) = mean(width_errors);

                    end
                end
                display([results_dir ':']);
                display(['MAE: ' num2str(naNmean(mae),3) ' +/- ' num2str(naNstd(mae),2) ' ' num2str(sum(~isnan(mae))) ' tests complete']);
           end
        end
    end
end
%%
% 4) Different ways of stitching Gaussian filters together
 

% G" + G'
mae = nan(num_repeats,1);
for i_repeat = 1:num_repeats
    results_dir = [exp_dir '/results/' num2str(i_repeat) '/orientation/g12d/'];
    try
        load([exp_dir 'test/' num2str(i_repeat) '/true_labels.mat'], 'true_oris');
        load([results_dir 'results.mat'], 'predictions');
    catch
        continue;
    end 
    [~, error_stats] = ori_error(true_oris, predictions);
    mae(i_repeat) = error_stats.abs_median;
end
display([results_dir ':']);
display(['MAE: ' num2str(naNmean(mae),4) ' +/- ' num2str(naNstd(mae),3) ' ' num2str(sum(~isnan(mae))) ' tests complete']);
%%
%G" + H"
mae = nan(num_repeats,1);
for i_repeat = 1:num_repeats
    results_dir = [exp_dir '/results/' num2str(i_repeat) '/orientation/gh2da/1/'];
    try
        load([exp_dir 'test/' num2str(i_repeat) '/true_labels.mat'], 'true_oris');
        load([results_dir 'results.mat'], 'predictions');
    catch
        continue;
    end
    [~, error_stats] = ori_error(true_oris, predictions);
    mae(i_repeat) = error_stats.abs_median;
end
display([results_dir ':']);
display(['MAE: ' num2str(naNmean(mae),4) ' +/- ' num2str(naNstd(mae),3) ' ' num2str(sum(~isnan(mae))) ' tests complete']);
%%
% G" + G'
mae = nan(num_repeats,1);
for i_repeat = 1:num_repeats
    results_dir = [exp_dir '/results/' num2str(i_repeat) '/width/g12d/'];
    try
        load([exp_dir 'test/' num2str(i_repeat) '/true_labels.mat'], 'true_widths');
        load([results_dir 'results.mat'], 'predictions');
    catch
        continue;
    end 
    width_errors = abs(true_widths - predictions);
    mae(i_repeat) = mean(width_errors);
end
display([results_dir ':']);
display(['MAE: ' num2str(naNmean(mae),4) ' +/- ' num2str(naNstd(mae),3) ' ' num2str(sum(~isnan(mae))) ' tests complete']);

%G" + H"
mae = nan(num_repeats,1);
for i_repeat = 1:num_repeats
    results_dir = [exp_dir '/results/' num2str(i_repeat) '/width/gh2d/'];    
    try
        load([exp_dir 'test/' num2str(i_repeat) '/true_labels.mat'], 'true_widths');
        load([results_dir 'results.mat'], 'predictions');
    catch
        continue;
    end  
    width_errors = abs(true_widths - predictions);
    mae(i_repeat) = mean(width_errors);
end
display([results_dir ':']);
display(['MAE: ' num2str(naNmean(mae),4) ' +/- ' num2str(naNstd(mae),3) ' ' num2str(sum(~isnan(mae))) ' tests complete']);
%%
% Combining all the filters together
%
for win_size = [1 3]
    mae = nan(num_repeats,1);
    for i_repeat = 1:num_repeats
        results_dir = [exp_dir 'results/' num2str(i_repeat) '/orientation/all/' num2str(win_size) '/'];

        try
            load([exp_dir 'test/' num2str(i_repeat) '/true_labels.mat'], 'true_oris');
            load([results_dir 'results.mat'], 'predictions');
        catch
            continue;
        end 
        [~, error_stats] = ori_error(true_oris, predictions);
        mae(i_repeat) = error_stats.abs_median;          
    end
    display([results_dir ':'])
    display(['MAE: ' num2str(naNmean(mae),4) ' +/- ' num2str(naNstd(mae),3) ' ' num2str(sum(~isnan(mae))) ' tests complete']);
end
%%
% Gabor filtering at more scales and angles
%
for exp_name = {'angles', 'scales', 'all', 'orig'}
    for win_size = [1 3]
        mae = nan(num_repeats,1);
        for i_repeat = 1:num_repeats
            results_dir = [exp_dir 'results/' num2str(i_repeat) '/orientation/gabor_all/' exp_name{1} '/' num2str(win_size) '/'];
            try
                load([exp_dir 'test/' num2str(i_repeat) '/true_labels.mat'], 'true_oris');
                load([results_dir 'results.mat'], 'predictions');
            catch
                continue;
            end 
            [~, error_stats] = ori_error(true_oris, predictions);
            mae(i_repeat) = error_stats.abs_median;          
        end
        display([results_dir ':'])
        display(['MAE: ' num2str(naNmean(mae),3) ' +/- ' num2str(naNstd(mae),3) ' ' num2str(sum(~isnan(mae))) ' tests complete']);
    end
end
%%

%%
%--------------------------------------------------------------------------
% ********************** DETECTION ****************************************
%--------------------------------------------------------------------------
%%
%
% 1) Original data - repeats
%

for i_decomp = 1:10%1:10;
    for win_size = [1 3]
        auc = nan(num_repeats,1);
        for i_repeat = 1:num_repeats
            results_dir = [exp_dir 'results/' num2str(i_repeat) '/detection/'...
                d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/'];
            
            try
                load([results_dir 'results1.mat'], 'predicted_lines'); 
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
%%
%
% 2) Complex representations
%
for i_decomp = [1 3];
    %Complex forms: all (mag/phase), real/imag, conj, mag, phase
    for feature_type = {'all', 'real', 'imag' 'real_imag', 'conj', 'mag', 'phase', 'real', 'imag', 'real_abs_imag'}
        auc = nan(num_repeats,1);
        for i_repeat = 1:num_repeats
            results_dir = [exp_dir '/results/' num2str(i_repeat) '/detection/'...
                d_args{i_decomp}.decomp_type{1} '/feature_types/' feature_type{1} '/1/'];
            
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
%%
%Complex forms: all (mag/phase), real/imag, conj, mag, phase
for feature_type = {'all' 'real_imag', 'conj', 'mag', 'phase'}
    auc = nan(num_repeats,1);
    for i_repeat = 1:num_repeats
        results_dir = [exp_dir '/results/' num2str(i_repeat) '/detection/gh2d/feature_types/' feature_type{1} '/'];

        try
            load([results_dir 'results1.mat'], 'predicted_lines'); 
        catch
            continue;
        end            
        [~, auc(i_repeat)] =...
            calculate_roc_curve(predicted_lines,class_labels,operating_pts);
    end

    display([results_dir ':'])
    display(['ROC Az = ' num2str(naNmean(auc), 3) ' +/- ' num2str(naNstd(auc), '%10.1e')  ' ' num2str(sum(~isnan(auc))) ' tests complete']);

end

%%
%
% 3) Different number of levels
%
for i_decomp = 3%[1 3 6];
    
    for i_level = 1:5
        for win_size = 3
            
            auc = nan(num_repeats,1);
            for i_repeat = 1:num_repeats
                results_dir = [exp_dir '/results/' num2str(i_repeat) '/detection/'...
                    d_args{i_decomp}.decomp_type{1} '/levels/' num2str(i_level) '/' num2str(win_size) '/'];
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
%%
%
% 4) Rotate and do max
%
for i_decomp = [1 3 9]%3             
    for win_size = 3            
        for reform_type = 1:2                
            if win_size==1 && reform_type==2
                continue;
            end

            auc = nan(num_repeats,1);
            for i_repeat = 1:num_repeats
                if reform_type == 1
                    results_dir = [exp_dir '/results/' num2str(i_repeat) '/detection/' d_args{i_decomp}.decomp_type{1}...
                        '/rotate/' num2str(win_size) '/'];
                else
                    results_dir = [exp_dir '/results/' num2str(i_repeat) '/detection/' d_args{i_decomp}.decomp_type{1}...
                        '/do_max/' num2str(win_size) '/'];
                end
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
%%
%
% 5) Different ways of stitching Gaussian filters together
%
% G" + G'    
auc = nan(num_repeats,1);
for i_repeat = 1:num_repeats
    results_dir = [exp_dir '/results/' num2str(i_repeat) '/detection/g12d/'];
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
%%
%G" + H"
auc = nan(num_repeats,1);
for i_repeat = 1:num_repeats
    results_dir = [exp_dir '/results/' num2str(i_repeat) '/detection/gh2da/1/'];
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
%%
for win_size = [1 3]
    auc = nan(num_repeats,1);
    for i_repeat = 1:num_repeats
        results_dir = [exp_dir 'results/' num2str(i_repeat) '/detection/all/' num2str(win_size) '/'];

        try
            load([results_dir 'results1.mat'], 'predicted_lines'); 
        catch
            continue;
        end            
        [~, auc(i_repeat)] =...
            calculate_roc_curve(predicted_lines,class_labels,operating_pts);
    end
    display([results_dir ':'])
    display(['ROC Az = ' num2str(naNmean(auc), 3) ' +/- ' num2str(naNstd(auc), '%10.1e')  ' ' num2str(sum(~isnan(auc))) ' tests complete']);
end
%%
for exp_name = {'angles', 'scales', 'all', 'orig'}
    for win_size = [1 3]
        auc = nan(num_repeats,1);
        for i_repeat = 1:num_repeats
            results_dir = [exp_dir '/results/' num2str(i_repeat) '/detection/gabor_all/' exp_name{1} '/' num2str(win_size) '/'];

            try
                load([results_dir 'results1.mat'], 'predicted_lines'); 
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
%%
%%
%
% 1) Original data
%
for i_decomp = 1:10;
    for win_size = [1 3]
        i_pts = 20000;
        
        results_dir = [exp_dir '/results/1'...
                '/detection/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/' num2str(i_pts) '/'];
        load([results_dir 'all_votes.mat'], 'all_votes');
        predicted_lines = sum(all_votes(:,2,:),3) / 100;
        
        [~, auc] = calculate_roc_curve(...
            predicted_lines,...
            class_labels, operating_pts);
        display([results_dir ':']);
        display(['Az (orig) = ' num2str(auc,4)]); 
    end
end
%%
%
% 1) Original data - increasing numbers of points
%
for i_decomp = 1:10;
    for win_size = [1 3]
        figure; hold all;
        title(['Decomp: ' d_args{i_decomp}.decomp_type{1} ', w = ' num2str(win_size)]);
        leg_text = cell(0,1);
        for i_pts = 10000:2000:20000
            results_dir = [exp_dir '/results/1'...
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

for i_decomp = 1:10;
    for win_size = [1 3]
        figure; hold all;
        title(['Decomp: ' d_args{i_decomp}.decomp_type{1} ', w = ' num2str(win_size)]);
        auc = nan(100,1);
        results_dir = [exp_dir '/results/1'...
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
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Some tree analysis
%--------------------------------------------------------------------------
exp_dir = 'C:\isbe\asymmetry_project\experiments\DRIVE\comparing_representations\rfs\1\';
for i_output = {'detection', 'orientation'}
    counts = zeros(1,192);
    for i_tree = 1:100
        load([exp_dir '\' i_output{1} '\gabor_all\martha\1\trees\rf_tree' zerostr(i_tree,4) '.mat']);
        leaves = ~tree.var;
        counts = counts + hist(tree.var(~leaves), 1:192);
    end
    %
    figure; a1 = gca; bar(1:192, counts); title([i_output{1} ': all bands']);
    figure; a2 = gca; hold all; title([i_output{1} ': magnitude']);
    figure; a3 = gca; hold all; title([i_output{1} ': magnitude']);

    colors = lines(16);
    cols = 1:6;
    for i_level = 1:16
        bar(a2, cols, counts(cols), 'facecolor', colors(i_level,:));
        bar(a3, cols, counts(cols+6), 'facecolor', colors(i_level,:));
        cols = cols + 12;
    end
end
%%
for i_output = {'detection', 'orientation'}
    counts = zeros(1,180);
    for i_tree = 1:100
        load([exp_dir '\' i_output{1} '\gabor_all\doris\1\trees\rf_tree' zerostr(i_tree,4) '.mat']);
        leaves = ~tree.var;
        counts = counts + hist(tree.var(~leaves), 1:180);
    end
    %
    figure; a1 = gca; bar(1:180, counts); title(i_output{1});
    figure; a2 = gca; hold all; title([i_output{1} ': magnitude']);
    figure; a3 = gca; hold all; title([i_output{1} ': phase']);

    colors = lines(5);
    cols = 1:18;
    for i_level = 1:5
        bar(a2, cols, counts(cols), 'facecolor', colors(i_level,:));
        bar(a3, cols, counts(cols+18), 'facecolor', colors(i_level,:));
        cols = cols + 36;
    end
end
%%
MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="[3]" LEVELS="[2 4 6 8 10 12 14 16]" NUM_ANGLES=18 EXP_NAME="mike" qsub -V -t 1 -l twoday matlab_code/trunk/hydra/cuc/generic_batch_fun.sh
MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="[1]" LEVELS="[2 4 6 8 10 12 14 16]" NUM_ANGLES=18 EXP_NAME="mike" qsub -V -t 1 -l twoday matlab_code/trunk/hydra/cuc/generic_batch_fun.sh
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="[1]" LEVELS="[2 4 6 8 10 12 14 16]" NUM_ANGLES=18 EXP_NAME="mike" qsub -V -t 1 -l twoday matlab_code/trunk/hydra/cuc/generic_batch_fun.sh
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="[3]" LEVELS="[2 4 6 8 10 12 14 16]" NUM_ANGLES=18 EXP_NAME="mike" qsub -V -t 1 -l twoday matlab_code/trunk/hydra/cuc/generic_batch_fun.sh

MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="[1]" LEVELS="[1 2 4 8 16]" NUM_ANGLES=18 EXP_NAME="doris" qsub -V -t 1 -l twoday -hold_jid 49607 matlab_code/trunk/hydra/cuc/generic_batch_fun.sh

MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="[3]" LEVELS="[1 2 4 8 16]" NUM_ANGLES=18 EXP_NAME="doris" qsub -V -t 1 -l twoday -hold_jid 49607 matlab_code/trunk/hydra/cuc/generic_batch_fun.sh

MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="[3]" LEVELS="[1 2 4 8 16]" NUM_ANGLES=18 EXP_NAME="doris" qsub -V -t 1 -l twoday -hold_jid 49607 matlab_code/trunk/hydra/cuc/generic_batch_fun.sh

MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="[1]" LEVELS="[1 2 4 8 16]" NUM_ANGLES=18 EXP_NAME="doris" qsub -V -t 1 -l twoday -hold_jid 49607 matlab_code/trunk/hydra/cuc/generic_batch_fun.sh

MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="[1]" LEVELS="[1 2 4 8 16]" NUM_ANGLES=9 EXP_NAME="esme" qsub -V -t 1 -l twoday -hold_jid 49607 matlab_code/trunk/hydra/cuc/generic_batch_fun.sh

MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="[3]" LEVELS="[1 2 4 8 16]" NUM_ANGLES=9 EXP_NAME="esme" qsub -V -t 1 -l twoday -hold_jid 49607 matlab_code/trunk/hydra/cuc/generic_batch_fun.sh

MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="[3]" LEVELS="[1 2 4 8 16]" NUM_ANGLES=9 EXP_NAME="esme" qsub -V -t 1 -l twoday -hold_jid 49607 matlab_code/trunk/hydra/cuc/generic_batch_fun.sh

MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="[1]" LEVELS="[1 2 4 8 16]" NUM_ANGLES=9 EXP_NAME="esme" qsub -V -t 1 -l twoday -hold_jid 49607 matlab_code/trunk/hydra/cuc/generic_batch_fun.sh

