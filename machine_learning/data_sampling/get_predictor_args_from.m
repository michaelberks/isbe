function [args_out] = get_predictor_args_from(args_in, args_out)
% Given a set of arguments, return only those arguments that are 
% relevant to the statistical predictor (or append them if args_out is supplied)

% deal with undefined inputs
if ~exist('args_out','var'), args_out = []; end

% Arguments common to all predictors
args_out.prediction_type = args_in.prediction_type;

% Models are stored in a directory structure under
%   model_root/output_type/prediction_type/job_id/
args_out.model_dir = [args_in.model_root '/' args_in.output_type '/' ...
                      args_in.prediction_type '/' args_in.job_id];
args_out.model_dir = prettypath(args_out.model_dir);

% Filename of this task's predictor
predictor_path = [args_in.predictor_name sprintf('%02d.mat', args_in.task_id)];
args_out.predictor_path = ...
    prettypath([args_out.model_dir '/' predictor_path]);

% Data path containing information on sampled data
sampled_data_subdir = sprintf('%02d_sampled_data/', args_in.task_id);
args_out.sampled_data_dir = ...
    prettypath([args_out.model_dir '/' sampled_data_subdir]);

% Arguments common to classifiers vs regressors
switch args_in.prediction_type
    case{'linear_classification', 'logistic_classification', ...
         'boosted_classification', ...
         'tree_classification', 'rf_classification'}
     
    case{'linear_regression', 'logistic_regression', ...
         'boosted_regression', ...
         'tree_regression', 'rf_regression'}
     
    otherwise
        error(['Prediction type ', args_in.prediction_type, ' not recognized']);
end

% Arguments specific to a particular type of classifier/regressor
switch args_in.prediction_type
    case{'linear_classification', 'linear_regression'}

    case{'logistic_classification', 'logistic_regression'}

    case{'boosted_classification', 'boosted_regression'}
        fields_to_copy = {...
            'boost_n_levels'; 'boost_weak_learner'; 'boost_output_type';
            'boost_n_bins'; 'boost_shrinkage';
        };
        args_out = get_substructure(args_in, fields_to_copy, args_out);

    case{'tree_classification', 'tree_regression'}
        args_out.n_trees = 1;
        
    case{'rf_classification', 'rf_regression'}
        % Random forest parameters
        fields_to_copy = { ...
            'n_trees'; 'd'; 'w_prior'; 'impure_thresh'; 
            'split_min'; 'end_cut_min'; 'do_ubound'; 'quiet'; 
            'do_circular'; 'overwrite'; 'minimise_size';
        };
        args_out = get_substructure(args_in, fields_to_copy, args_out);
        
        %Use a different default split/var criterion for regression and
        %classification
        if strcmpi(args_in.prediction_type, 'rf_classification')
            args_out.split_criterion = args_in.split_criterion_c;
            args_out.var_criterion = args_in.var_criterion_c; 
        else
            args_out.split_criterion = args_in.split_criterion_r;
            args_out.var_criterion = args_in.var_criterion_r; 
        end

        % Folder in which to store individual trees
        tree_subdir = sprintf('%02d_trees/', args_in.task_id);
        args_out.tree_dir = ...
            prettypath([args_out.model_dir '/' tree_subdir]);
        
    otherwise
        error(['Prediction type ', args_in.prediction_type, ' not recognized']);
end        
