% Script to test that all predictors build for all output types
% successfully.

close_all_timebars();

outputs = {'detection', 'width', 'orientation'};

%All the old valid types before we swapped to the new cellstr format
features = {'pixel',... % dumbest
          'linop', ... % template matching
          'g1d', 'g2d', 'g2di', 'g2da', 'haar', 'g2dg', ... % derivs
          'dt', 'mono', 'g12d', 'g2dh', 'gabor', 'gaborg',... % with phase
          'dtg2' % suggested by dumb reviewer
         };

classifiers = {'rf_classification'};
regressors = {'linear_regression','rf_regression','boosted_regression'};
all_predictors = [classifiers regressors];

success = nan(length(outputs),length(features),length(all_predictors));
error_strs = cell(0,2);

for i_output = 1:length(outputs)
    output_type = outputs{i_output};
    
    switch output_type
        case {'detection'}
            predictors = classifiers;
            split_criterion = 'gdi';
        case {'orientation','width'}
            predictors = regressors;
            split_criterion = 'mabs';
        otherwise
            error(['Output type ',output_type,' not recognized']);
    end
    
    feature_type = 'dt';
    i_feature = strcmp(features, 'dt');
    
    for i_prediction = 1:length(predictors)
        prediction_type = predictors{i_prediction};
        i_all_pred = find(strcmp(all_predictors, prediction_type));
        success(i_output, i_feature, i_all_pred) = false;
            build_predictor('output_type', output_type, ...
                            'decomposition_type',feature_type,...
                            'prediction_type', prediction_type,...
                            'split_criterion',split_criterion,...
                            'quiet', 1);
        success(i_output, i_feature, i_all_pred) = true;
    end
end

for i_feature = 1:length(features)
        feature_type = features{i_feature};
        success(i_output, i_feature, i_all_pred) = false;
        build_predictor('output_type', output_type, ...
                                'decomposition_type',feature_type,...
                                'prediction_type', prediction_type,...
                                'split_criterion',split_criterion,...
                                'quiet', 1);
        success(i_output, i_feature, i_all_pred) = true;
end

close_all_timebars();

fprintf('\n\n');
for i = 1:size(error_strs,1)
    fprintf('%s:\n%s\n\n',error_strs{i,1},error_strs{i,2});
end

n_total = sum(~isnan(success(:)));
n_success = sum(success(:)==1);
fprintf('%d out of %d tests successful\n', n_success, n_total);
